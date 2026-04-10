#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>
#include <memory>
#include <map>
#include "bsc/alignment.hpp"
#include "bsc/minimizers.hpp"

struct Sequence {
    std::string name;
    std::string data;
};

std::string reverse_complement(const std::string& seq) {
    std::string rc = seq;
    std::reverse(rc.begin(), rc.end());
    for (char& c : rc) {
        if (c == 'A') c = 'T';
        else if (c == 'T') c = 'A';
        else if (c == 'C') c = 'G';
        else if (c == 'G') c = 'C';
        // else leave as is
    }
    return rc;
}

std::vector<std::unique_ptr<Sequence>> parse_fasta(const std::string& path) {
    std::vector<std::unique_ptr<Sequence>> sequences;
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::string line;
    std::string name;
    std::string data;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) {
                sequences.push_back(std::make_unique<Sequence>(Sequence{name, data}));
            }
            // Take until first space
            size_t space = line.find(' ');
            name = (space != std::string::npos) ? line.substr(1, space - 1) : line.substr(1);
            data.clear();
        } else {
            data += line;
        }
    }
    if (!name.empty()) {
        sequences.push_back(std::make_unique<Sequence>(Sequence{name, data}));
    }
    return sequences;
}

int main(int argc, char* argv[]) {
    std::string version = "v0.1.0";
    std::string help = "Usage: bsc_mapper [options] <reference.fasta> <fragments.fasta/q>\n"
                       "Options:\n"
                       "  -h, --help     Show this help message\n"
                       "  --version      Show version\n"
                       "  -a <type>      Alignment type: global, local, semi_global (default: global)\n"
                       "  -m <score>     Match score (default: 1)\n"
                       "  -n <score>     Mismatch score (default: -1)\n"
                       "  -g <score>     Gap score (default: -1)\n"
                       "  -k <len>       K-mer length (default: 15)\n"
                       "  -w <len>       Window length (default: 5)\n"
                       "  -f <frac>      Frequency threshold (default: 0.001)\n"
                       "  -c             Output CIGAR strings\n";

    bsc::AlignmentType atype = bsc::AlignmentType::GLOBAL;
    int match_score = 1;
    int mismatch_score = -1;
    int gap_score = -1;
    unsigned int kmer_len = 15;
    unsigned int window_len = 5;
    double freq_threshold = 0.001;
    bool output_cigar = false;

    int optind = 1;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout << help << std::endl;
            return 0;
        } else if (arg == "--version") {
            std::cout << version << std::endl;
            return 0;
        } else if (arg == "-a") {
            if (i + 1 < argc) {
                std::string type = argv[++i];
                if (type == "global") atype = bsc::AlignmentType::GLOBAL;
                else if (type == "local") atype = bsc::AlignmentType::LOCAL;
                else if (type == "semi_global") atype = bsc::AlignmentType::SEMI_GLOBAL;
                else {
                    std::cerr << "Invalid alignment type" << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Missing argument for -a" << std::endl;
                return 1;
            }
        } else if (arg == "-m") {
            if (i + 1 < argc) match_score = std::stoi(argv[++i]);
            else {
                std::cerr << "Missing argument for -m" << std::endl;
                return 1;
            }
        } else if (arg == "-n") {
            if (i + 1 < argc) mismatch_score = std::stoi(argv[++i]);
            else {
                std::cerr << "Missing argument for -n" << std::endl;
                return 1;
            }
        } else if (arg == "-g") {
            if (i + 1 < argc) gap_score = std::stoi(argv[++i]);
            else {
                std::cerr << "Missing argument for -g" << std::endl;
                return 1;
            }
        } else if (arg == "-k") {
            if (i + 1 < argc) kmer_len = std::stoi(argv[++i]);
            else {
                std::cerr << "Missing argument for -k" << std::endl;
                return 1;
            }
        } else if (arg == "-w") {
            if (i + 1 < argc) window_len = std::stoi(argv[++i]);
            else {
                std::cerr << "Missing argument for -w" << std::endl;
                return 1;
            }
        } else if (arg == "-f") {
            if (i + 1 < argc) freq_threshold = std::stod(argv[++i]);
            else {
                std::cerr << "Missing argument for -f" << std::endl;
                return 1;
            }
        } else if (arg == "-c") {
            output_cigar = true;
        } else {
            // Assume file
            if (optind == 1) {
                optind = i;
                break;
            }
        }
    }

    if (argc - optind != 2) {
        std::cerr << "Error: Exactly two files required.\n" << help << std::endl;
        return 1;
    }

    std::string ref_file = argv[optind];
    std::string frag_file = argv[optind + 1];

    // Parse reference (FASTA)
    auto ref_sequences = parse_fasta(ref_file);

    std::cerr << "Reference sequences:" << std::endl;
    for (const auto& seq : ref_sequences) {
        std::cerr << seq->name << " " << seq->data.size() << std::endl;
    }

    // Parse fragments (FASTA)
    auto frag_sequences = parse_fasta(frag_file);

    size_t num_seqs = frag_sequences.size();
    if (num_seqs == 0) {
        std::cerr << "No sequences in fragments file." << std::endl;
        return 1;
    }

    std::vector<size_t> lengths;
    for (const auto& seq : frag_sequences) {
        lengths.push_back(seq->data.size());
    }

    size_t total_len = std::accumulate(lengths.begin(), lengths.end(), 0ULL);
    double avg_len = static_cast<double>(total_len) / num_seqs;
    size_t min_len = *std::min_element(lengths.begin(), lengths.end());
    size_t max_len = *std::max_element(lengths.begin(), lengths.end());

    // N50
    std::sort(lengths.rbegin(), lengths.rend());
    size_t cumulative = 0;
    size_t n50 = 0;
    for (size_t len : lengths) {
        cumulative += len;
        if (cumulative >= total_len / 2) {
            n50 = len;
            break;
        }
    }

    std::cerr << "Fragments statistics:" << std::endl;
    std::cerr << "Number of sequences: " << num_seqs << std::endl;
    std::cerr << "Average length: " << avg_len << std::endl;
    std::cerr << "N50: " << n50 << std::endl;
    std::cerr << "Min length: " << min_len << std::endl;
    std::cerr << "Max length: " << max_len << std::endl;

    // Compute minimizers
    std::map<std::string, int> ref_minimizer_freq;
    std::string ref_seq = ref_sequences[0]->data;
    auto ref_minimizers = bsc::Minimize(ref_seq.c_str(), ref_seq.size(), kmer_len, window_len);
    for (auto& t : ref_minimizers) {
        unsigned int pos1, pos2;
        bool strand;
        std::tie(pos1, pos2, strand) = t;
        std::string minimizer = strand ? ref_seq.substr(pos1, kmer_len) : reverse_complement(ref_seq).substr(pos2, kmer_len);
        ref_minimizer_freq[minimizer]++;
    }
    int num_distinct_ref = ref_minimizer_freq.size();

    std::map<std::string, int> frag_minimizer_freq;
    for (auto& seq : frag_sequences) {
        auto minims = bsc::Minimize(seq->data.c_str(), seq->data.size(), kmer_len, window_len);
        for (auto& t : minims) {
            unsigned int pos1, pos2;
            bool strand;
            std::tie(pos1, pos2, strand) = t;
            std::string minimizer = strand ? seq->data.substr(pos1, kmer_len) : reverse_complement(seq->data).substr(pos2, kmer_len);
            frag_minimizer_freq[minimizer]++;
        }
    }
    int num_distinct_frag = frag_minimizer_freq.size();
    int num_singletons = 0;
    int max_freq = 0;
    for (auto& p : frag_minimizer_freq) {
        if (p.second == 1) num_singletons++;
        if (p.second > max_freq) max_freq = p.second;
    }
    double fraction_singletons = static_cast<double>(num_singletons) / num_distinct_frag;

    std::cerr << "Reference minimizers:" << std::endl;
    std::cerr << "Number of distinct minimizers: " << num_distinct_ref << std::endl;
    std::cerr << "Fragments minimizers:" << std::endl;
    std::cerr << "Number of distinct minimizers: " << num_distinct_frag << std::endl;
    std::cerr << "Fraction of singletons: " << fraction_singletons << std::endl;
    std::cerr << "Number of occurrences of the most frequent minimizer: " << max_freq << std::endl;

    // Align two random sequences <=5000
    std::vector<size_t> candidates;
    for (size_t idx = 0; idx < num_seqs; ++idx) {
        if (frag_sequences[idx]->data.size() <= 5000) {
            candidates.push_back(idx);
        }
    }
    if (candidates.size() < 2) {
        std::cerr << "Not enough sequences <=5000 for alignment." << std::endl;
        return 1;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, candidates.size() - 1);
    size_t idx1 = candidates[dis(gen)];
    size_t idx2 = candidates[dis(gen)];
    while (idx2 == idx1) idx2 = candidates[dis(gen)];

    const auto& seq1 = frag_sequences[idx1];
    const auto& seq2 = frag_sequences[idx2];

    std::string cigar;
    unsigned int target_begin;
    int score = bsc::Align(
        seq1->data.c_str(), seq1->data.size(),
        seq2->data.c_str(), seq2->data.size(),
        atype, match_score, mismatch_score, gap_score,
        &cigar, &target_begin);

    std::cerr << "Alignment of sequence " << seq1->name << " and " << seq2->name << ":" << std::endl;
    std::cerr << "Score: " << score << std::endl;
    std::cerr << "CIGAR: " << cigar << std::endl;
    std::cerr << "Target begin: " << target_begin << std::endl;

    return 0;
}