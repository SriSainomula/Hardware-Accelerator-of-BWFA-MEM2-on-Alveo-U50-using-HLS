#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <ctime>
#include <chrono>

using namespace std;



// Function to perform the Needleman-Wunsch alignment algorithm
// for the calculation of the partial alignment score (i.e., scoring for one row of the matrix)
vector<int> needleman_wunsch_scoring_affine_linear_space(const string& seq1, const string& seq2, int gap_open_penalty, int gap_extend_penalty, int match_score, int mismatch_penalty) {
    int n = seq1.size();
    int m = seq2.size();

    // Previous and current rows for match, insertion, and deletion matrices
    vector<int> prev_D(n + 1, 0), prev_P(n + 1, 0), prev_Q(n + 1, 0);
    vector<int> curr_D(n + 1, 0), curr_P(n + 1, 0), curr_Q(n + 1, 0);

    // Initialize the first column (when aligning seq1 to gaps in seq2)
    //prev_insertion[0] = gap_open_penalty;  // First gap in insertion row (gap opening)
    //prev_deletion[0] = gap_open_penalty;   // First gap in deletion row (gap opening)
    //prev_match[0] = gap_open_penalty;      // First match row score (gap opening)
    
    //cout << "M[0]:" << prev_match[0] << ",";
    for (int j = 1; j <= n; ++j) {
        // For insertion, start with gap opening penalty, followed by gap extension penalties
        //prev_insertion[j] = gap_open_penalty + j * gap_extend_penalty;
        prev_P[j] = -1000;
        // For deletion, start with gap opening penalty, followed by gap extension penalties
        //prev_deletion[j] = gap_open_penalty + j * gap_extend_penalty;
        prev_Q[j] = -1000;
        // For match, this is based on aligning seq1[0] to seq2[j-1], but for first row, it's just gap penalties
        prev_D[j] = gap_open_penalty + j * gap_extend_penalty;
        //cout << "M[" << j << "]:" << prev_match[j] << ",";
    }
    //cout << endl;
    
    //cout << "P      Q      D" << endl;

    // Fill the matrices column by column
    for (int i = 1; i <= m; ++i) {
        // First row for each column (aligning gaps in seq1 or seq2)
        curr_P[0] = -1000;
        curr_Q[0] = -1000;
        curr_D[0] = gap_open_penalty + i * gap_extend_penalty;
          
        
        for (int j = 1; j <= n; ++j) {
            int match_or_mismatch = (seq1[j - 1] == seq2[i - 1]) ? match_score : mismatch_penalty;

            curr_P[j] = max(
                curr_P[j - 1] + gap_extend_penalty, // gap extension in seq2
                curr_D[j - 1] + gap_open_penalty + gap_extend_penalty           // gap opening in seq1
            );

            curr_Q[j] = max(
                prev_Q[j] + gap_extend_penalty,     // gap extension in seq2
                prev_D[j] + gap_open_penalty + gap_extend_penalty           // gap opening in seq2
            );
            
            // Calculate the current row's values for match, insertion, and deletion
            curr_D[j] = max({prev_D[j - 1] + match_or_mismatch, curr_P[j], curr_Q[j]});
            //cout << curr_P[j] << "    " << curr_Q[j] << "   " << curr_D[j] << endl;

        }
        //cout << endl;

        // Swap the current and previous rows
        swap(prev_D, curr_D);
        swap(prev_P, curr_P);
        swap(prev_Q, curr_Q);
    }

    // Return the final row of scores from the match matrix
    return prev_D;
}

// BWFA's Algorithm for Sequence Alignment in Linear Space
// Returns Two Srings
pair<string, string> bwfa_alignment(const string& seq1, const string& seq2, int gap_penalty, int gap_extend, int match_score, int mismatch_penalty) {
    if (seq1.empty()){
      cout << "Aligned Seq-1:" << string(seq2.size(), '-') << endl;
      cout << "Aligned Seq-2:" << seq2 << endl;
      return {string(seq2.size(), '-'), seq2}; //When seq1 is empty, we return a string of gaps of length seq2.size() and seq2 itself.
    }
    
    if (seq2.empty()) return {seq1, string(seq1.size(), '-')};

    if (seq1.size() == 1 || seq2.size() == 1) {
        cout << "Reached Size=1" << endl;
        // initializes a 2D DP table with dimensions (seq1.size() + 1) x (seq2.size() + 1) to store the alignment scores.
        vector<vector<int>> dp(seq1.size() + 1, vector<int>(seq2.size() + 1));
        
        // Fill the dynamic programming table
        dp[0][0] = 0;
        for (int i = 1; i <= seq1.size(); ++i) dp[i][0] = gap_penalty + i * gap_extend;
        for (int j = 1; j <= seq2.size(); ++j) dp[0][j] = gap_penalty + j * gap_extend;

        for (int i = 1; i <= seq1.size(); ++i) {
            for (int j = 1; j <= seq2.size(); ++j) {
                int match_or_mismatch = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_penalty;
                dp[i][j] = max({dp[i - 1][j - 1] + match_or_mismatch, dp[i - 1][j] + gap_extend, dp[i][j - 1] + gap_extend});
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
            } else if (dp[i][j] == dp[i - 1][j] + gap_extend) {
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
        
        cout << "Aligned Seq-1:" << aligned_seq1 << endl;
        cout << "Aligned Seq-2:" << aligned_seq2 << endl;

        return {aligned_seq1, aligned_seq2};
    }

    int mid= seq2.size() / 2;
    
    cout << "Seq-1:" << seq1 << endl;
    cout << "Seq-2:" << seq2 << endl;
    cout << "Mid = " << mid << endl;
    cout << "Scoring-1:" << endl;
    vector<int> score1 = needleman_wunsch_scoring_affine_linear_space(seq1, seq2.substr(0, mid), gap_penalty, gap_extend, match_score, mismatch_penalty);
    cout << "Scoring-2:" << endl;
    vector<int> score2 = needleman_wunsch_scoring_affine_linear_space(string(seq1.rbegin(), seq1.rend()), string(seq2.rbegin(), seq2.rend()).substr(0, mid), gap_penalty, gap_extend, match_score, mismatch_penalty);
    
    
    int best_split = 0;
    int best_score = score1[0] + score2[seq1.size()];
    for (int i = 1; i <= seq1.size(); ++i) {
        int split_score = score1[i] + score2[seq1.size() - i];
        if (split_score > best_score) {
            best_score = split_score;
            best_split = i;
        }
    }
    cout << "Best Score = " << best_score << ", Split = " << best_split << endl;
    
    string seq1_left = seq1.substr(0, best_split);
    string seq1_right = seq1.substr(best_split);
    string seq2_left = seq2.substr(0, mid);
    string seq2_right = seq2.substr(mid);
    
    cout << "Subproblem-1:" << seq1_left << "," << seq2_left << endl;
    cout << "Subproblem-2:" << seq1_right << "," << seq2_right << endl;
      
    
    auto left_align = bwfa_alignment(seq1_left, seq2_left, gap_penalty, gap_extend, match_score, mismatch_penalty);
    auto right_align = bwfa_alignment(seq1_right, seq2_right, gap_penalty, gap_extend, match_score, mismatch_penalty);

    return {left_align.first + right_align.first, left_align.second + right_align.second};
}

int main() 
{
    //string seq1 = "CG";
    //string seq2 = "CCGA";
    
    //string seq1 = "acgtcact";
    //string seq2 = "ccatcctt";
    
    string seq1 = "AGTGATGG";
    string seq2 = "AGTGATGT";
    
    //string seq1 = "TCTATACTGCGCGTTTGGAGAAATAAAATAGTCCGTGGCCCCGGGCGCCTTCCTGGGCCTGA";
    //string seq2 = "TCTTTACTCGCGCGTTGGAGAAATACAATAGTATGGTGCTGTCTCCTGCCGACAAGACCAAC";
    
    //string seq1 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGA";
    //string seq2 = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACG";
    
    //string seq1 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGAGT";
    //string seq2 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGAGT";
    
    
    //string seq1 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGACCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGACCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGATCTATACTGCGCGTTTGGAGAAATAAAATAGTCCGTGGCCCCGGGCGCCTTCCTGGGCCTGAATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGTCTTTACTCGCGCGTTGGAGAAATACAATAGTATGGTGCTGTCTCCTGCCGACAAGACCAAC";
    //string seq2 = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGTCTTTACTCGCGCGTTGGAGAAATACAATAGTATGGTGCTGTCTCCTGCCGACAAGACCAACCCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGACCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGACCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGATCTATACTGCGCGTTTGGAGAAATAAAATAGTCCGTGGCCCCGGGCGCCTTCCTGGGCCTGA";
    
    //string seq1 = "CCGTGGCCCCGGGCGCCTTCCTGGGCCTGAAGGCGCTGCGATGGCTGGACCTGTCCCACATGGATCTATACTGCGCGTTTGGAGAAATAAAATAGTCCGTGGCCCCGGGCGCCTTCCTGGGCCTGA";
    //string seq2 = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGTCTTTACTCGCGCGTTGGAGAAATACAATAGTATGGTGCTGTCTCCTGCCGACAAGACCAAC";
  
    int gap_penalty = -2;
    int gap_extend = -1;
    int match_score = 1;
    int mismatch_penalty = -2;
    
    //clock_t start = clock();
    auto start = std::chrono::high_resolution_clock::now();
    
    auto result = bwfa_alignment(seq1, seq2, gap_penalty,gap_extend, match_score, mismatch_penalty);
    
    //clock_t end = clock();
    //double cpu_time_used = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start).count();
    
    cout << "Aligned Sequences:" << endl;
    cout << result.first << endl;
    cout << result.second << endl;
    
    std::cout << "Length of Sequences: " << seq1.size() << endl;
    
    //std::cout << "CPU time used: " << cpu_time_used << " seconds\n";
    std::cout << "Time taken: " << duration << " ms\n";
    
    return 0;
}