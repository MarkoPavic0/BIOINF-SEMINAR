#include "bsc/minimizers.hpp"

#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <algorithm>

namespace bsc {

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

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len) {

    std::string seq(sequence, sequence_len);
    std::string rc = reverse_complement(seq);
    std::set<std::tuple<unsigned int, unsigned int, bool>> unique_minimizers;

    // For original strand
    for (unsigned int start = 0; start + window_len + kmer_len - 1 <= sequence_len; ++start) {
        std::string min_kmer = seq.substr(start, kmer_len);
        unsigned int min_pos = start;
        for (unsigned int i = 1; i < window_len; ++i) {
            unsigned int pos = start + i;
            std::string kmer = seq.substr(pos, kmer_len);
            if (kmer < min_kmer) {
                min_kmer = kmer;
                min_pos = pos;
            }
        }
        unique_minimizers.insert(std::make_tuple(min_pos, min_pos, true));
    }

    // For reverse complement
    for (unsigned int start = 0; start + window_len + kmer_len - 1 <= sequence_len; ++start) {
        std::string min_kmer = rc.substr(start, kmer_len);
        unsigned int min_pos = start;
        for (unsigned int i = 1; i < window_len; ++i) {
            unsigned int pos = start + i;
            std::string kmer = rc.substr(pos, kmer_len);
            if (kmer < min_kmer) {
                min_kmer = kmer;
                min_pos = pos;
            }
        }
        unique_minimizers.insert(std::make_tuple(min_pos, min_pos, false));
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> result(unique_minimizers.begin(), unique_minimizers.end());
    return result;
}

} // namespace bsc