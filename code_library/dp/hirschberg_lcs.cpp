/***
 *
 * Hirschberg LCS
 * One longest common subsequence of two sequences, in linear memory
 *
 * Complexity: O(n m) time, O(n + m) memory
 *
 * lcs_sequence(a, b): one longest common subsequence of a and b, in their Sequence type (a string for strings)
 *
 * Notes:
 *   - Elements need only ==, a and b any random-access containers with push_back
 *   - Only the length: bit_string_lcs.cpp in O(n m / 64), or hunt_szymanski.cpp when few positions match
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
struct Hirschberg{
    const Sequence &a, &b;
    Sequence res;

    Hirschberg(const Sequence& a, const Sequence& b) : a(a), b(b){
        solve(0, a.size(), 0, b.size());
    }

    void solve(int a0, int a1, int b0, int b1){
        if (a0 >= a1 || b0 >= b1) return;
        if (a1 - a0 == 1){
            for (int j = b0; j < b1; j++){
                if (a[a0] == b[j]){
                    res.push_back(a[a0]);
                    return;
                }
            }

            return;
        }

        /// An optimal subsequence splits b where the first half of a ends; the rows are freed before recursing
        int mid = (a0 + a1) / 2, split = b0;
        {
            vector<int> left = forward_row(a0, mid, b0, b1), right = backward_row(mid, a1, b0, b1);
            for (int j = 0, best = -1; j <= b1 - b0; j++){
                if (left[j] + right[j] > best) best = left[j] + right[j], split = b0 + j;
            }
        }

        solve(a0, mid, b0, split);
        solve(mid, a1, split, b1);
    }

    /// row[j] = LCS length of a[a0..a1) and b[b0..b0 + j)
    vector<int> forward_row(int a0, int a1, int b0, int b1) const{
        vector<int> row(b1 - b0 + 1, 0);
        for (int i = a0; i < a1; i++){
            int diag = 0;
            for (int j = b0; j < b1; j++){
                int up = row[j - b0 + 1];
                row[j - b0 + 1] = a[i] == b[j] ? diag + 1 : max(up, row[j - b0]);
                diag = up;
            }
        }

        return row;
    }

    /// row[j] = LCS length of a[a0..a1) and b[b0 + j..b1)
    vector<int> backward_row(int a0, int a1, int b0, int b1) const{
        vector<int> row(b1 - b0 + 1, 0);
        for (int i = a1 - 1; i >= a0; i--){
            int diag = 0;
            for (int j = b1 - 1; j >= b0; j--){
                int down = row[j - b0];
                row[j - b0] = a[i] == b[j] ? diag + 1 : max(down, row[j - b0 + 1]);
                diag = down;
            }
        }

        return row;
    }
};

template <typename Sequence>
Sequence lcs_sequence(const Sequence& a, const Sequence& b){
    /// The rows run over the second sequence, so make it the shorter one
    if (a.size() < b.size()) return Hirschberg<Sequence>(b, a).res;
    return Hirschberg<Sequence>(a, b).res;
}

int main(){
    /// "GTAB" is the only longest common subsequence of these two strings
    assert(lcs_sequence(string("AGGTAB"), string("GXTXAYB")) == "GTAB");
    assert((lcs_sequence(vector<int>{5, 1, 2, 3}, vector<int>{1, 9, 2, 3, 7}) == vector<int>{1, 2, 3}));
    assert(lcs_sequence(string("same"), string("same")) == "same");
    assert(lcs_sequence(string("abc"), string("xyz")).empty());
    assert(lcs_sequence(string("abc"), string("")).empty());

    return 0;
}
