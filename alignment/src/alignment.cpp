#include "alignment/alignment.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace alignment {

enum class BacktrackDir : unsigned char {
    NONE = 0,
    DIAG = 1,
    UP = 2,
    LEFT = 3
};

unsigned int getIndex(unsigned int row, unsigned int col, unsigned int cols) {
    return row * cols + col;
}

std::string compressCigar(const std::string& ops) {
    std::string result;
    if (ops.empty()) {
        return result;
    }

    char prev = ops[0];
    unsigned int count = 1;
    for (size_t i = 1; i < ops.size(); ++i) {
        if (ops[i] == prev) {
            ++count;
        } else {
            result += std::to_string(count);
            result.push_back(prev);
            prev = ops[i];
            count = 1;
        }
    }
    result += std::to_string(count);
    result.push_back(prev);
    return result;
}


int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar,
    unsigned int* target_begin) {

    unsigned int m = query_len;
    unsigned int n = target_len;
    unsigned int cols = n + 1;
    unsigned int cells = (m + 1) * cols;

    std::vector<int> dp(cells, 0);
    std::vector<unsigned char> backtrack(cells, static_cast<unsigned char>(BacktrackDir::NONE));

    // Initialization
    dp[getIndex(0, 0, cols)] = 0;
    for (unsigned int i = 1; i <= m; ++i) {
        if (type == AlignmentType::GLOBAL) {
            dp[getIndex(i, 0, cols)] = dp[getIndex(i - 1, 0, cols)] + gap;
        } else {
            dp[getIndex(i, 0, cols)] = 0;
        }
        backtrack[getIndex(i, 0, cols)] = static_cast<unsigned char>(BacktrackDir::UP);
    }
    for (unsigned int j = 1; j <= n; ++j) {
        if (type == AlignmentType::GLOBAL) {
            dp[getIndex(0, j, cols)] = dp[getIndex(0, j - 1, cols)] + gap;
        } else {
            dp[getIndex(0, j, cols)] = 0;
        }
        backtrack[getIndex(0, j, cols)] = static_cast<unsigned char>(BacktrackDir::LEFT);
    }

    int best_score = std::numeric_limits<int>::min();
    unsigned int best_i = 0;
    unsigned int best_j = 0;

    if (type == AlignmentType::GLOBAL) {
        best_score = dp[getIndex(m, n, cols)];
        best_i = m;
        best_j = n;
    } else if (type == AlignmentType::LOCAL) {
        best_score = 0;
    }

    // Fill DP matrix
    for (unsigned int i = 1; i <= m; ++i) {
        for (unsigned int j = 1; j <= n; ++j) {
            int score_diag = dp[getIndex(i - 1, j - 1, cols)] + ((query[i - 1] == target[j - 1]) ? match : mismatch);
            int score_up = dp[getIndex(i - 1, j, cols)] + gap;
            int score_left = dp[getIndex(i, j - 1, cols)] + gap;

            int score = score_diag;
            BacktrackDir dir = BacktrackDir::DIAG;
            if (score_up > score) {
                score = score_up;
                dir = BacktrackDir::UP;
            }
            if (score_left > score) {
                score = score_left;
                dir = BacktrackDir::LEFT;
            }
            if (type == AlignmentType::LOCAL && score < 0) {
                score = 0;
                dir = BacktrackDir::NONE;
            }

            dp[getIndex(i, j, cols)] = score;
            backtrack[getIndex(i, j, cols)] = static_cast<unsigned char>(dir);

            if (type == AlignmentType::LOCAL && score > best_score) {
                best_score = score;
                best_i = i;
                best_j = j;
            }
            //std::cout << " " << score << " (" << (dir == BacktrackDir::DIAG ? "M" : dir == BacktrackDir::UP ? "I" : dir == BacktrackDir::LEFT ? "D" : "0") << ")";
        }
        //std::cout << "\n";
    }

    if (type == AlignmentType::SEMI_GLOBAL) {
        best_score = std::numeric_limits<int>::min();
        best_i = 0;
        best_j = 0;
        for (unsigned int j = 0; j <= n; ++j) {
            int score = dp[getIndex(m, j, cols)];
            if (score > best_score) {
                best_score = score;
                best_i = m;
                best_j = j;
            }
        }
        for (unsigned int i = 0; i <= m; ++i) {
            int score = dp[getIndex(i, n, cols)];
            if (score > best_score) {
                best_score = score;
                best_i = i;
                best_j = n;
            }
        }
    }

    if (type == AlignmentType::GLOBAL) {
        best_score = dp[getIndex(m, n, cols)];
        best_i = m;
        best_j = n;
    }

    if (type == AlignmentType::LOCAL && best_score == 0) {
        if (cigar) {
            *cigar = std::string();
        }
        if (target_begin) {
            *target_begin = 0;
        }
        return 0;
    }

    unsigned int i = best_i;
    unsigned int j = best_j;
    std::string operations;
    operations.reserve(i + j);

    while (i > 0 || j > 0) {
        if (type == AlignmentType::LOCAL && dp[getIndex(i, j, cols)] == 0) {
            break;
        }

        BacktrackDir dir = static_cast<BacktrackDir>(backtrack[getIndex(i, j, cols)]);
        if (dir == BacktrackDir::DIAG) {
            operations.push_back('M');
            --i;
            --j;
        } else if (dir == BacktrackDir::UP) {
            operations.push_back('I');
            --i;
        } else if (dir == BacktrackDir::LEFT) {
            operations.push_back('D');
            --j;
        } else {
            break;
        }
    }

    if (type == AlignmentType::GLOBAL || type == AlignmentType::SEMI_GLOBAL) {
        while (i > 0) {
            operations.push_back('I');
            --i;
        }
        while (j > 0) {
            operations.push_back('D');
            --j;
        }
    }

    if (target_begin) {
        *target_begin = j;
    }

    if (cigar) {
        std::reverse(operations.begin(), operations.end());
        *cigar = compressCigar(operations);
    }

    return best_score;
}

}
