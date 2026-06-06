/***
 *
 * Suffix Array using DC3 (Difference Cover 3) Algorithm
 * Linear time construction with LCP array
 * Constructs lexicographically sorted array of all suffixes
 *
 * Complexity: O(n) for suffix array construction, O(n) for LCP array
 *
 * Highly optimized DC3 algorithm uses recursive radix sort with difference cover modulo 3
 *
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Working memory for DC3 algorithm
namespace SuffixArrayDC3 {
    vector<int> s0, sa0, bucket, mem;

    void radixsort(int* source, int* dest, int* val, int n, int lim){
        bucket.assign(lim, 0);
        for (int i = 0; i < n; i++) bucket[val[source[i]]]++;

        int s = 0;
        for (int i = 0; i < lim; i++){
            int x = bucket[i];
            bucket[i] = s;
            s += x;
        }
        for (int i = 0; i < n; i++) dest[bucket[val[source[i]]]++] = source[i];
    }

    void DC3(int* ar, int* sa, int n, int lim, int ptr){
        int *s12, *sa12;
        int allc = (n / 3) << 1, n0 = (n + 2) / 3;
        int i, j, k, l, c, d, p, t, m, r, counter;

        s12 = &mem[ptr], ptr += (allc + 5);
        sa12 = &mem[ptr], ptr += (allc + 5);

        c = 0, m = 0, r = n + ((n % 3) == 1);
        for (i = 0; i < r; i++, m++){
            if (m == 3) m = 0;
            if (m) s12[c++] = i;
        }
        s12[c] = sa12[c] = s12[c + 1] = sa12[c + 1] = s12[c + 2] = sa12[c + 2] = 0;

        radixsort(s12, sa12, ar + 2, c, lim + 1);
        radixsort(sa12, s12, ar + 1, c, lim + 1);
        radixsort(s12, sa12, ar, c, lim + 1);

        counter = 0, j = -1;
        for (i = 0; i < c; i++){
            if ((ar[sa12[i]] != j) || (ar[sa12[i] + 1] != k) || (ar[sa12[i] + 2] != l)){
                counter++;
                j = ar[sa12[i]], k = ar[sa12[i] + 1], l = ar[sa12[i] + 2];
            }
            if ((sa12[i] % 3) == 1) s12[sa12[i] / 3] = counter;
            else s12[(sa12[i] / 3) + n0] = counter;
        }

        if (counter == c){
            for (i = 0; i < c; i++) sa12[s12[i] - 1] = i;
        }
        else{
            DC3(s12, sa12, c, counter, ptr);
            for (i = 0; i < c; i++) s12[sa12[i]] = i + 1;
        }

        for (i = 0, d = 0; i < c; i++){
            if (sa12[i] < n0) s0[d++] = (sa12[i] * 3);
        }
        radixsort(&s0[0], &sa0[0], ar, d, lim + 1);

        for (k = 0, l = ((n % 3) == 1), r = 0; r < n; r++){
            j = sa0[k];
            i = ((sa12[l] < n0) ? (sa12[l] * 3) + 1 : ((sa12[l] - n0) * 3) + 2);
            if (l == c) sa[r] = sa0[k++];
            else if (k == d) sa[r] = i, l++;
            else{
                if (sa12[l] < n0){
                    if ((ar[i] < ar[j]) || (ar[i] == ar[j] && s12[sa12[l] + n0] <= s12[j / 3])) sa[r] = i, l++;
                    else sa[r] = j, k++;
                }
                else{
                    if ((ar[i] < ar[j]) || (ar[i] == ar[j] && ar[i + 1] < ar[j + 1]) ||
                        (ar[i] == ar[j] && ar[i + 1] == ar[j + 1] && s12[sa12[l] - n0 + 1] <= s12[(j / 3) + n0]))
                        sa[r] = i, l++;
                    else sa[r] = j, k++;
                }
            }
        }
    }

    template<typename Container>
    void build_lcp(const Container& str, int n, const vector<int>& sa, vector<int>& lcp){
        vector<int> rank(n);
        for (int i = 0; i < n; i++) rank[sa[i]] = i;

        int k = 0;
        for (int i = 0; i < n; i++, k ? k-- : 0){
            if (rank[i] == n - 1){
                k = 0;
            }
            else{
                int j = sa[rank[i] + 1];
                while ((i + k) < n && (j + k) < n && str[i + k] == str[j + k]) k++;
            }
            lcp[rank[i]] = k;
        }
    }
}

struct SuffixArray {
	int n;
    vector<int> sa;   // sa[i] = starting position of i-th smallest suffix
    vector<int> lcp;  // lcp[i] = longest common prefix of sa[i] and sa[i+1]

    // Count distinct substrings = n*(n+1)/2 - sum(lcp)
    long long distinct_substrings() const {
        long long total = ((long long)n * (n + 1)) / 2;
        for (int i = 0; i + 1 < n; i++) total -= lcp[i];
        return total;
    }
};

/// Build suffix array for any container (string, vector<int>, etc)
template<typename Container>
SuffixArray suffix_array(const Container& c){
    int n = c.size();
    SuffixArray result;
    result.n = n;

    if (n == 0){
        return result;
    }

    vector<int> ar(n + 3);
    int lim = 0;
    for (int i = 0; i < n; i++){
        ar[i] = c[i];
        if (ar[i] > lim) lim = ar[i];
    }
    ar[n] = ar[n + 1] = ar[n + 2] = 0;

    int alloc_size = max(n + 3, ((n / 3) + 10) * 2);
    SuffixArrayDC3::s0.resize(alloc_size);
    SuffixArrayDC3::sa0.resize(alloc_size);
    SuffixArrayDC3::bucket.resize(lim + 10);
    SuffixArrayDC3::mem.resize(n * 4 + 100);

    result.sa.resize(n);
    SuffixArrayDC3::DC3(&ar[0], &result.sa[0], n, lim, 0);

    result.lcp.resize(n);
    SuffixArrayDC3::build_lcp(c, n, result.sa, result.lcp);

    return result;
}

SuffixArray suffix_array(const char* s){
    return suffix_array(string(s));
}

int main(){
    auto sa1 = suffix_array(vector<int>{2, 1, 4, 1, 4, 1});
    assert((sa1.sa == vector<int>{5, 3, 1, 0, 4, 2}));
    assert((sa1.lcp == vector<int>{1, 3, 0, 0, 2, 0}));

    auto sa2 = suffix_array("mississippi");
    assert((sa2.sa == vector<int>{10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2}));
    assert(sa2.distinct_substrings() == 53);

    return 0;
}
