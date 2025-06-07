// 1. Finding length of Longest Common Subsequence (LCS)
/*
#include <bits/stdc++.h>
using namespace std;

// Function to find the length of the Longest Common Subsequence (LCS)
int lcsUtil(string& s1, string& s2, int ind1, int ind2, vector<vector<int>>& dp) {
    // Base case: If either string reaches the end, return 0
    if (ind1 < 0 || ind2 < 0)
        return 0;

    // If the result for this pair of indices is already calculated, return it
    if (dp[ind1][ind2] != -1)
        return dp[ind1][ind2];

    // If the characters at the current indices match, increment the LCS length
    if (s1[ind1] == s2[ind2])
        return dp[ind1][ind2] = 1 + lcsUtil(s1, s2, ind1 - 1, ind2 - 1, dp);
    else
        // If the characters don't match, consider two options: moving either left or up in the strings
        return dp[ind1][ind2] = max(lcsUtil(s1, s2, ind1, ind2 - 1, dp), lcsUtil(s1, s2, ind1 - 1, ind2, dp));
}

// Function to calculate the length of the Longest Common Subsequence
int lcs(string s1, string s2) {
    int n = s1.size();
    int m = s2.size();

    vector<vector<int>> dp(n, vector<int>(m, -1)); // Create a DP table
    return lcsUtil(s1, s2, n - 1, m - 1, dp);
}

int main() {
    string s1 = "acd";
    string s2 = "ced";

    // Call the function to find and output the length of the Longest Common Subsequence
    cout << "The Length of Longest Common Subsequence is " << lcs(s1, s2) << endl;

    return 0; // Return 0 to indicate successful program execution
}
*/

// 2. Printing lcs

#include <bits/stdc++.h>
using namespace std;

void lcs(string s1, string s2) {

  int n = s1.size();
  int m = s2.size();

  vector < vector < int >> dp(n + 1, vector < int > (m + 1, 0));
  for (int i = 0; i <= n; i++) {
    dp[i][0] = 0;
  }
  for (int i = 0; i <= m; i++) {
    dp[0][i] = 0;
  }

  for (int ind1 = 1; ind1 <= n; ind1++) {
    for (int ind2 = 1; ind2 <= m; ind2++) {
      if (s1[ind1 - 1] == s2[ind2 - 1])
        dp[ind1][ind2] = 1 + dp[ind1 - 1][ind2 - 1];
      else
        dp[ind1][ind2] = 0 + max(dp[ind1 - 1][ind2], dp[ind1][ind2 - 1]);
    }
  }

  int len = dp[n][m];
  int i = n;
  int j = m;

  int index = len - 1;
  string str = "";
  for (int k = 1; k <= len; k++) {
    str += "$"; // dummy string
  }

  while (i > 0 && j > 0) {
    if (s1[i - 1] == s2[j - 1]) {
      str[index] = s1[i - 1];
      index--;
      i--;
      j--;
    } else if (s1[i - 1] > s2[j - 1]) {
      i--;
    } else j--;
  }
  cout << str;
}

int main() {

  string s1 = "ABCABBA";
  string s2 = "CBABAC";

  cout << "The Longest Common Subsequence is ";
  lcs(s1, s2);
}
