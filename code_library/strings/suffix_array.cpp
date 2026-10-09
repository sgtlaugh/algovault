/***
 *
 * Suffix Array using DC3 (Difference Cover 3) Algorithm
 * Builds the lexicographically sorted array of all suffixes with its LCP array, plus O(1) LCP, O(log n) substring and O(|pattern| log n) pattern occurrence queries
 *
 * Complexity: O(n) construction when the values span at most n + 256, else O(n log n) for the value compression sort
 * LCP array in O(n), SuffixArrayQueries adds an O(n log n) sparse table on top
 * Memory: the sparse table holds ~n log n ints (~76 MB at n = 1e6) and SuffixArrayQueries keeps a copy of the text
 *
 * suffix_array(c) works on any container of integers (string, vector<int>, vector<unsigned long long>), any values including 0 and negatives
 *   sa[i]  = start of the i-th smallest suffix
 *   lcp[i] = longest common prefix of the suffixes at sa[i] and sa[i + 1], lcp[n - 1] = 0
 *
 * SuffixArrayQueries<Container> q(c) has everything above and also:
 *   rank[p]                 = position of the suffix starting at p in sa
 *   common_prefix(i, j)     = longest common prefix of the suffixes starting at i and j, O(1)
 *   occurrences(p, len)     = [lo, hi) such that sa[lo .. hi) are exactly the starts of s[p .. p + len), O(log n)
 *   occurrences(pattern)    = the same range for any pattern, O(|pattern| log n)
 * The occurrence count is hi - lo, an absent pattern gives lo == hi (its insertion point), an empty one gives [0, n)
 *
 * Example:
 *   SuffixArrayQueries q("banana");
 *   q.occurrences("ana");     // {1, 3}: sa[1] = 3, sa[2] = 1
 *   q.common_prefix(1, 3);    // 3, "anana" and "ana"
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Scratch buffers shared by every call, so a call is not reentrant or thread safe
/// Per call locals measured 4-35% slower (noisy) on repeated n = 1e6 builds (each call page faults ~24 MB fresh) and equal on small n
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
        int i, j, k, l, c, d, m, r, counter;

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

        counter = 0, j = k = l = -1;
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
    vector<int> sa;  // sa[i] = starting position of i-th smallest suffix
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

    /// DC3 wants values in [1, lim] with 0 as sentinel and lim buckets: shift small ranges, compress large ones
    /// Values stay in their own type: a cast to long long would wrap unsigned 64-bit values past LLONG_MAX below small ones
    using T = typename Container::value_type;
    using U = make_unsigned_t<T>;
    T lo = *min_element(c.begin(), c.end()), hi = *max_element(c.begin(), c.end());
    auto offset = [&](T x){ return (unsigned long long)U(U(x) - U(lo)); };  /// modular unsigned difference, exact as x >= lo

    vector<int> ar(n + 3, 0);

    int lim;
    if (offset(hi) <= (unsigned long long)n + 256){
        for (int i = 0; i < n; i++) ar[i] = offset(c[i]) + 1;
        lim = offset(hi) + 1;
    }
    else{
        vector<T> values(c.begin(), c.end());
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());
        for (int i = 0; i < n; i++) ar[i] = lower_bound(values.begin(), values.end(), c[i]) - values.begin() + 1;
        lim = values.size();
    }

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

template<typename Container>
struct SuffixArrayQueries : SuffixArray {
    Container text;
    vector<int> rank;
    vector<vector<int>> table;  /// table[k][i] = min(lcp[i .. i + 2^k))

    SuffixArrayQueries(const Container& c) : SuffixArray(suffix_array(c)), text(c), rank(n){
        for (int i = 0; i < n; i++) rank[sa[i]] = i;

        table.push_back(lcp);
        for (int k = 1; (1 << k) <= n; k++){
            vector<int> level(n - (1 << k) + 1);
            for (int i = 0; i < (int)level.size(); i++) level[i] = min(table[k - 1][i], table[k - 1][i + (1 << (k - 1))]);
            table.push_back(move(level));
        }
    }

    int common_prefix(int i, int j) const {
        if (i == j) return n - i;

        int a = rank[i], b = rank[j];
        if (a > b) swap(a, b);
        return lcp_min(a, b);
    }

    /// Requires 0 <= p, 0 <= len, p + len <= n
    pair<int, int> occurrences(int p, int len) const {
        if (len == 0) return {0, n};

        int r = rank[p], lo = 0, hi = r;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (lcp_min(mid, r) >= len) hi = mid;
            else lo = mid + 1;
        }

        int first = lo;
        lo = r, hi = n - 1;
        while (lo < hi){
            int mid = (lo + hi + 1) / 2;
            if (lcp_min(r, mid) >= len) lo = mid;
            else hi = mid - 1;
        }

        return {first, lo + 1};
    }

    pair<int, int> occurrences(const Container& pattern) const {
        return {bound(pattern, false), bound(pattern, true)};
    }

    /// First rank whose suffix cut to |pattern| values is >= pattern, or > pattern when strict
    int bound(const Container& pattern, bool strict) const {
        int lo = 0, hi = n;
        while (lo < hi){
            int mid = (lo + hi) / 2, cmp = compare(sa[mid], pattern);
            if (cmp < 0 || (strict && cmp == 0)) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    /// Sign of the suffix at p cut to |pattern| values against pattern, a suffix that ends first is smaller
    int compare(int p, const Container& pattern) const {
        int m = pattern.size();
        for (int k = 0; k < m; k++){
            if (p + k == n) return -1;
            if (text[p + k] != pattern[k]) return text[p + k] < pattern[k] ? -1 : 1;
        }
        return 0;
    }

    /// min(lcp[l .. r)) for l < r, the lcp of the suffixes at sa[l] and sa[r]
    int lcp_min(int l, int r) const {
        int k = __lg(r - l);
        return min(table[k][l], table[k][r - (1 << k)]);
    }
};

SuffixArrayQueries(const char*) -> SuffixArrayQueries<string>;

int main(){
    auto sa1 = suffix_array(vector<int>{2, 1, 4, 1, 4, 1});
    assert((sa1.sa == vector<int>{5, 3, 1, 0, 4, 2}));
    assert((sa1.lcp == vector<int>{1, 3, 0, 0, 2, 0}));

    auto sa2 = suffix_array("mississippi");
    assert((sa2.sa == vector<int>{10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2}));
    assert(sa2.distinct_substrings() == 53);

    SuffixArrayQueries q1("banana");
    assert((q1.sa == vector<int>{5, 3, 1, 0, 4, 2}));
    assert((q1.rank == vector<int>{3, 2, 5, 1, 4, 0}));
    assert(q1.common_prefix(1, 3) == 3 && q1.common_prefix(2, 4) == 2 && q1.common_prefix(0, 5) == 0);
    assert(q1.common_prefix(0, 0) == 6 && q1.common_prefix(5, 5) == 1);
    assert((q1.occurrences("ana") == pair<int, int>{1, 3}));
    assert((q1.occurrences("na") == pair<int, int>{4, 6}));
    assert((q1.occurrences("banana") == pair<int, int>{3, 4}));
    assert((q1.occurrences("bananas") == pair<int, int>{4, 4}));
    assert((q1.occurrences("x") == pair<int, int>{6, 6}));
    assert((q1.occurrences("") == pair<int, int>{0, 6}));
    assert((q1.occurrences(1, 3) == pair<int, int>{1, 3}));
    assert((q1.occurrences(4, 2) == pair<int, int>{4, 6}));
    assert((q1.occurrences(0, 6) == pair<int, int>{3, 4}));
    assert((q1.occurrences(5, 1) == pair<int, int>{0, 3}));
    assert((q1.occurrences(3, 0) == pair<int, int>{0, 6}));

    SuffixArrayQueries q2("mississippi");
    assert(q2.common_prefix(1, 4) == 4 && q2.common_prefix(2, 5) == 3 && q2.common_prefix(0, 10) == 0);
    assert((q2.occurrences("issi") == pair<int, int>{2, 4}));
    assert((q2.occurrences("ss") == pair<int, int>{9, 11}));
    assert((q2.occurrences(5, 3) == pair<int, int>{9, 11}));
    assert((q2.occurrences(0, 1) == pair<int, int>{4, 5}));

    SuffixArrayQueries<vector<long long>> q3({-5, 0, -5, 0});
    assert((q3.sa == vector<int>{2, 0, 3, 1}));
    assert(q3.common_prefix(0, 2) == 2 && q3.common_prefix(1, 3) == 1);
    assert((q3.occurrences({-5, 0}) == pair<int, int>{0, 2}));
    assert((q3.occurrences({0}) == pair<int, int>{2, 4}));
    assert((q3.occurrences({0, 0}) == pair<int, int>{4, 4}));

    SuffixArrayQueries q4("");
    assert((q4.occurrences("") == pair<int, int>{0, 0}));
    assert((q4.occurrences("a") == pair<int, int>{0, 0}));
    assert((q4.occurrences(0, 0) == pair<int, int>{0, 0}));

    SuffixArrayQueries<vector<unsigned long long>> q5({1ULL << 63, 1, 1ULL << 63, 5});
    assert((q5.sa == vector<int>{1, 3, 0, 2}));
    assert((q5.occurrences({1}) == pair<int, int>{0, 1}));
    assert((q5.occurrences({1ULL << 63}) == pair<int, int>{2, 4}));

    return 0;
}
