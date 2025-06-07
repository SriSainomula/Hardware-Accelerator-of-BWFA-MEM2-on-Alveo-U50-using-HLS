#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;


// Function to perform the Needleman-Wunsch alignment algorithm
// for the calculation of the partial alignment score (i.e., scoring for one row of the matrix)
vector<int> needleman_wunsch_scoring(const string& seq1, const string& seq2, int gap_penalty, int match_score, int mismatch_penalty) {
    int m = seq1.size();
    int n = seq2.size();
    vector<int> score(n + 1);

    // Initialize the first row of the dynamic programming table
    score[0] = 0;
    for (int j = 1; j <= n; ++j) {
        score[j] = score[j - 1] + gap_penalty;
    }

    vector<int> new_score(n + 1);

    // Calculate the score row by row
    for (int i = 1; i <= m; ++i) {
        new_score[0] = score[0] + gap_penalty;
        for (int j = 1; j <= n; ++j) {
            int match_or_mismatch = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;
            new_score[j] = max({new_score[j - 1] + gap_penalty, score[j] + gap_penalty, score[j - 1] + match_or_mismatch});
        }
        swap(score, new_score);
    }
    
    return score;
}

// BWFA's Algorithm for Sequence Alignment in Linear Space
//pair is used to store two values of possibly different types as a single unit.
//These values are stored as members of the pair and can be accessed using the first and second members.
pair<string, string> bwfa_alignment(const string& seq1, const string& seq2, int gap_penalty, int match_score, int mismatch_penalty) {
    if (seq1.empty()) return {string(seq2.size(), '-'), seq2};
    if (seq2.empty()) return {seq1, string(seq1.size(), '-')};

    if (seq1.size() == 1 || seq2.size() == 1) {
        //This line creates a 2D vector dp with dimensions(seq1.size() + 1) x (seq2.size() + 1), 
        //where each element is initialized to 0. 
        //It uses the constructor of the outer vector vector<vector<int>> to initialize the rows,
        //each containing a vector of size (seq2.size() + 1).
        vector<vector<int>> dp(seq1.size() + 1, vector<int>(seq2.size() + 1));
        
        // Fill the dynamic programming table
        for (int i = 0; i <= seq1.size(); ++i) dp[i][0] = i * gap_penalty;
        for (int j = 0; j <= seq2.size(); ++j) dp[0][j] = j * gap_penalty;

        for (int i = 1; i <= seq1.size(); ++i) {
            for (int j = 1; j <= seq2.size(); ++j) {
                int match_or_mismatch = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;
                dp[i][j] = max({dp[i - 1][j - 1] + match_or_mismatch, dp[i - 1][j] + gap_penalty, dp[i][j - 1] + gap_penalty});
            }
        }

        // Backtrack to get the optimal alignment
        string aligned_seq1, aligned_seq2;
        int i = seq1.size(), j = seq2.size();
        while (i > 0 && j > 0) {
            if (dp[i][j] == dp[i - 1][j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty)) {
                aligned_seq1.push_back(seq1[i - 1]);
                aligned_seq2.push_back(seq2[j - 1]);
                --i;
                --j;
            } else if (dp[i][j] == dp[i - 1][j] + gap_penalty) {
                aligned_seq1.push_back(seq1[i - 1]);
                aligned_seq2.push_back('-');
                --i;
            } else {
                aligned_seq1.push_back('-');
                aligned_seq2.push_back(seq2[j - 1]);
                --j;
            }
        }

        // Fill in the remaining parts
        while (i > 0) {
            aligned_seq1.push_back(seq1[i - 1]);
            aligned_seq2.push_back('-');
            --i;
        }
        while (j > 0) {
            aligned_seq1.push_back('-');
            aligned_seq2.push_back(seq2[j - 1]);
            --j;
        }

        reverse(aligned_seq1.begin(), aligned_seq1.end());
        reverse(aligned_seq2.begin(), aligned_seq2.end());

        return {aligned_seq1, aligned_seq2};
    }

    int mid = seq1.size() / 2;
    vector<int> score1 = needleman_wunsch_scoring(seq1.substr(0, mid), seq2, gap_penalty, match_score, mismatch_penalty);
    vector<int> score2 = needleman_wunsch_scoring(string(seq1.rbegin(), seq1.rend()).substr(0, mid), string(seq2.rbegin(), seq2.rend()), gap_penalty, match_score, mismatch_penalty);

    int best_split = 0;
    int best_score = score1[0] + score2[seq2.size()];
    for (int i = 1; i <= seq2.size(); ++i) {
        int split_score = score1[i] + score2[seq2.size() - i];
        if (split_score >= best_score) {
            best_score = split_score;
            best_split = i;
        }
    }

    string seq1_left = seq1.substr(0, mid);
    string seq1_right = seq1.substr(mid);
    string seq2_left = seq2.substr(0, best_split);
    string seq2_right = seq2.substr(best_split);

    auto left_align = bwfa_alignment(seq1_left, seq2_left, gap_penalty, match_score, mismatch_penalty);
    auto right_align = bwfa_alignment(seq1_right, seq2_right, gap_penalty, match_score, mismatch_penalty);

    return {left_align.first + right_align.first, left_align.second + right_align.second};
}

int main() {
    //string seq1 = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
    //string seq2 = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
    
    string seq1 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGA";
    string seq2 = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACG";
    
    //string seq1 = "ACTTAATT";
    //string seq2 = "GAGCAATT";
    
    int gap_penalty = 0;
    int match_score = 1;
    int mismatch_penalty = 0;
    
    //The auto keyword in C++ is used for type inference. 
    //It allows the compiler to automatically deduce the type of a variable from 
    //its initializer expression. 
    auto result = bwfa_alignment(seq1, seq2, gap_penalty, match_score, mismatch_penalty);

    cout << "Aligned Sequences:" << endl;
    cout << result.first << endl;
    cout << result.second << endl;

    return 0;
}
