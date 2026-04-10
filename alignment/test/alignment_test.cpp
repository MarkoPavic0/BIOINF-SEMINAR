extern "C" {
  struct gzFile_s;
  typedef struct gzFile_s* gzFile;
  gzFile gzopen(const char* path, const char* mode);
  int gzread(gzFile file, void* buf, unsigned len);
  long gzseek(gzFile file, long offset, int whence);
  int gzclose(gzFile file);
}

#include "alignment/alignment.hpp"
#include "bioparser/fasta_parser.hpp"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

struct Sequence {
  std::string name;
  std::string data;

  Sequence(const char* name_, std::uint32_t name_len,
           const char* data_, std::uint32_t data_len)
      : name(name_, name_len), data(data_, data_len) {}
};

namespace fs = std::filesystem;

static alignment::AlignmentType PathToAlignmentType(const std::string& path) {
  auto lowered = path;
  for (auto& c : lowered) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (lowered.find("lokalno") != std::string::npos) {
    return alignment::AlignmentType::LOCAL;
  }
  if (lowered.find("poluglobalno") != std::string::npos ||
      lowered.find("pologlobalno") != std::string::npos ||
      lowered.find("semi_global") != std::string::npos) {
    return alignment::AlignmentType::SEMI_GLOBAL;
  }
  return alignment::AlignmentType::GLOBAL;
}

static uint32_t AlignmentBlockLength(const std::string& cigar) {
  uint32_t length = 0;
  uint32_t value = 0;
  for (char c : cigar) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      value = value * 10 + (c - '0');
    } else {
      if (c == 'M' || c == 'I' || c == 'D') {
        length += value;
      }
      value = 0;
    }
  }
  return length;
}

static uint32_t TargetConsumed(const std::string& cigar) {
  uint32_t consumed = 0;
  uint32_t value = 0;
  for (char c : cigar) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      value = value * 10 + (c - '0');
    } else {
      if (c == 'M' || c == 'D') {
        consumed += value;
      }
      value = 0;
    }
  }
  return consumed;
}

static uint32_t QueryConsumed(const std::string& cigar) {
  uint32_t consumed = 0;
  uint32_t value = 0;
  for (char c : cigar) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      value = value * 10 + (c - '0');
    } else {
      if (c == 'M' || c == 'I') {
        consumed += value;
      }
      value = 0;
    }
  }
  return consumed;
}

static uint32_t CountMatches(
    const std::string& cigar,
    const std::string& query,
    const std::string& target,
    uint32_t query_begin,
    uint32_t target_begin) {
  uint32_t matches = 0;
  uint32_t qpos = query_begin;
  uint32_t tpos = target_begin;
  uint32_t value = 0;

  for (char c : cigar) {
    if (std::isdigit(static_cast<unsigned char>(c))) {
      value = value * 10 + (c - '0');
      continue;
    }

    if (c == 'M') {
      for (uint32_t i = 0; i < value; ++i) {
        if (qpos + i < query.size() && tpos + i < target.size() &&
            query[qpos + i] == target[tpos + i]) {
          matches += 1;
        }
      }
      qpos += value;
      tpos += value;
    } else if (c == 'I') {
      qpos += value;
    } else if (c == 'D') {
      tpos += value;
    }
    value = 0;
  }

  return matches;
}

int main() {
  const fs::path source_dir = fs::path(__FILE__).parent_path();
  const fs::path results_dir = source_dir / "results";

  if (!fs::exists(results_dir)) {
    fs::create_directories(results_dir);
  }

  const std::vector<std::string> file_names = {
    "1_primjer_globalno_poravnanje.fasta.txt",
    "1_primjer_globalno_poravnanje2.fasta.txt",
    "2_primjer_poluGlobalno_poravnanje.fasta.txt",
    "3_primjer_lokalno_poravnanje.fasta.txt"
  };

  bool any_failed = false;
  for (const auto& file_name : file_names) {
    const fs::path file_path = source_dir / "data" / file_name;
    std::cout << "Processing " << file_name << "\n";

    if (!fs::exists(file_path)) {
      std::cerr << "File not found: " << file_path << '\n';
      any_failed = true;
      continue;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> parser;
    try {
      parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(file_path.string());
    } catch (const std::exception& error) {
      std::cerr << "Unable to open file: " << error.what() << '\n';
      any_failed = true;
      continue;
    }

    std::vector<std::unique_ptr<Sequence>> sequences;
    try {
      sequences = parser->Parse(std::numeric_limits<std::uint64_t>::max());
    } catch (const std::exception& error) {
      std::cerr << "Error parsing file " << file_path << ": " << error.what() << '\n';
      any_failed = true;
      continue;
    }

    if (sequences.size() < 2) {
      std::cerr << "Skipping " << file_name << " because it contains fewer than 2 sequences." << '\n';
      any_failed = true;
      continue;
    }

    const auto& query = *sequences[0];
    const auto& target = *sequences[1];
    const auto type = PathToAlignmentType(file_name);
    const auto& query_data = query.data;
    const auto& target_data = target.data;

    std::string cigar;
    unsigned int target_begin = 0;

    std::cout << "Aligning " << query.name << " (length " << query_data.size() << ") to "
              << target.name << " (length " << target_data.size() << ") using "
              << (type == alignment::AlignmentType::GLOBAL ? "global" :
                  type == alignment::AlignmentType::LOCAL ? "local" : "semi-global")
              << " alignment...\n";

    const int score = alignment::Align(
        query_data.data(), static_cast<unsigned int>(query_data.size()),
        target_data.data(), static_cast<unsigned int>(target_data.size()),
        type,
        2,  // match
        -1, // mismatch
        -1, // gap
        &cigar,
        &target_begin);

    const uint32_t target_consumed = TargetConsumed(cigar);
    const uint32_t query_consumed = QueryConsumed(cigar);
    const uint32_t matches = (type == alignment::AlignmentType::LOCAL)
        ? 0
        : CountMatches(cigar, query_data, target_data, 0u, target_begin);

    const fs::path output_path = results_dir / (file_path.stem().string() + ".paf");
    std::ofstream output(output_path);
    if (!output) {
      std::cerr << "Unable to create output file: " << output_path << '\n';
      any_failed = true;
      continue;
    }

    const uint32_t query_start = 0;
    const uint32_t query_end = static_cast<uint32_t>(query_data.size());
    const uint32_t target_end = target_begin + target_consumed;
    const uint32_t block_length = AlignmentBlockLength(cigar);

    output << query.name << '\t'
           << query_data.size() << '\t'
           << query_start << '\t'
           << query_end << '\t'
           << '+' << '\t'
           << target.name << '\t'
           << target_data.size() << '\t'
           << target_begin << '\t'
           << target_end << '\t'
           << matches << '\t'
           << block_length << '\t'
           << "255" << '\t'
           << "cg:Z:" << cigar << '\n';

    output.close();
    std::cout << "Wrote " << output_path.string() << " score=" << score << " cigar=" << cigar << '\n';
  }

  return any_failed ? 1 : 0;
}
