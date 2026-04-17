#include "bsc/minimizers.hpp"

#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include <set>
#include <algorithm>
#include <functional>

namespace bsc
{
    std::string reverse_complement(const std::string &seq)
    {
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        for (char &c : rc)
        {
            if (c == 'A')
                c = 'T';
            else if (c == 'T')
                c = 'A';
            else if (c == 'C')
                c = 'G';
            else if (c == 'G')
                c = 'C';
            // else leave as is
        }
        return rc;
    }

    // Funciton returns hashed calue of canonical k-mer (the lexicographically smaller of the k-mer and its reverse complement)
    std::tuple<std::string, bool> get_canonical(const std::string& kmer) {
        std::string rc = reverse_complement(kmer);
        if (kmer <= rc) return {kmer, true}; 
        return {rc, false};
    }

    std::tuple<std::string, unsigned int, bool> find_window_minimizer(const std::string &seq, unsigned int start_pos, unsigned int kmer_len, unsigned int window_len)
    {
        std::string minimizer = seq.substr(start_pos, kmer_len);
        bool is_reversed;
        std::tie(minimizer, is_reversed) = get_canonical(minimizer);

        unsigned int min_pos = start_pos;
        for (unsigned int i = 1; i < window_len; ++i)
        {
            unsigned int pos = start_pos + i;
            std::string kmer = seq.substr(pos, kmer_len);
            bool it_is_reversed;

            std::tie(kmer, it_is_reversed) = get_canonical(kmer);

            if (kmer < minimizer)
            {
                minimizer = kmer;
                min_pos = pos;
                is_reversed = it_is_reversed;
            }
        }
        return std::make_tuple(minimizer, min_pos, is_reversed);
    }

    std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
        const char *sequence, unsigned int sequence_len,
        unsigned int kmer_len,
        unsigned int window_len)
    {

        std::string seq(sequence, sequence_len);
        std::vector<std::tuple<unsigned int, unsigned int, bool>> result;
        std::hash<std::string> hasher;

        unsigned int last_win_minimizer_position;
        std::string last_win_minimizer;

        // First windows
        bool is_reversed;
        std::tie(last_win_minimizer, last_win_minimizer_position, is_reversed) = find_window_minimizer(seq, 0, kmer_len, window_len);
        result.push_back(std::make_tuple(hasher(last_win_minimizer), last_win_minimizer_position, is_reversed));

        // Other Windows iteration
        for (unsigned int start = 1; start + window_len + kmer_len - 1 <= sequence_len; ++start)
        {
            // Check if new kmer is smallest
            bool last_kmer_is_reversed;
            unsigned int last_kmer_in_window_pos = start + window_len - 1;
            std::string last_kmer_in_window = seq.substr(last_kmer_in_window_pos, kmer_len);
            std::tie(last_kmer_in_window, last_kmer_is_reversed) = get_canonical(last_kmer_in_window);

            // If new kmer is smaller or equal to current minimizer, add it to the result
            if (last_kmer_in_window <= last_win_minimizer)
            {
                last_win_minimizer = last_kmer_in_window;
                last_win_minimizer_position = last_kmer_in_window_pos;
                result.push_back(std::make_tuple(hasher(last_win_minimizer), last_win_minimizer_position, last_kmer_is_reversed));
            }

            // If last windows minimizer is 'pushed out' of the windows, find new one
            else if (last_win_minimizer_position < start)
            {
                bool is_reversed;

                std::tie(last_win_minimizer, last_win_minimizer_position, is_reversed) = find_window_minimizer(seq, start, kmer_len, window_len);

                result.push_back(std::make_tuple(hasher(last_win_minimizer), last_win_minimizer_position, is_reversed));
            }
        }

        return result;
    }

} // namespace bsc