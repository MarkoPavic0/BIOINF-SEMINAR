#ifndef BSC_MINIMIZERS_HPP
#define BSC_MINIMIZERS_HPP

#include <vector>
#include <tuple>
#include <string>

namespace bsc {

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

} // namespace bsc

#endif // BSC_MINIMIZERS_HPP