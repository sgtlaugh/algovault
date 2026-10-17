/***
 *
 * Suffix Array (SA-IS)
 * Sorted suffixes and their LCP array in O(n), with O(1) common prefix and pattern occurrence queries
 *
 * Complexity:
 *   - suffix_array: O(n) when the values span at most n + 256, else O(n log n) to compress them
 *   - SuffixArrayQueries: O(n log n) more for a sparse table of n log n ints, about 76 MB at n = 1e6, plus a copy of the text
 *   - common_prefix: O(1)
 *   - occurrences(p, len): O(log n)
 *   - occurrences(pattern): O(|pattern| log n)
 *
 * SuffixArray res = suffix_array(c): c is any container of integers (string, vector<int>, vector<unsigned long long>), any values
 * res.sa[i]: start of the i-th smallest suffix
 * res.lcp[i]: longest common prefix of the suffixes at sa[i] and sa[i + 1], lcp[n - 1] = 0
 * res.distinct_substrings(): number of distinct non-empty substrings, O(n)
 *
 * SuffixArrayQueries<Container> q(c): everything in SuffixArray plus the calls below
 * q.rank[p]: position of the suffix starting at p in sa
 * q.common_prefix(i, j): longest common prefix of the suffixes starting at i and j
 * q.occurrences(p, len): [lo, hi) such that sa[lo .. hi) are exactly the starts of s[p .. p + len), 0 <= p, 0 <= len, p + len <= n
 * q.occurrences(pattern): the same range for any pattern
 *
 * The occurrence count is hi - lo, an absent pattern gives lo == hi at its insertion point, an empty one gives [0, n)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Replaced DC3, at n = 5e5 to 1e6: 2-3.5x faster on repetitive text, 10-30% on random letters, 3-8% slower on all distinct values
namespace SuffixArraySAIS{
    /// Sorts the suffixes of text[0 .. n), values in [0, alphabet), into sa[0 .. n) by induced sorting (Nong, Zhang and Chan)
    /// Past the end sits a virtual sentinel smaller than every value, and the recursion runs inside sa, so only buckets are allocated
    template<typename Char>
    void sais(const Char* text, int* sa, int n, int alphabet){
        if (n == 1){
            sa[0] = 0;
            return;
        }

        vector<int> start(alphabet + 1, 0), bucket(alphabet);  /// bucket c spans sa[start[c] .. start[c + 1])
        for (int i = 0; i < n; i++) start[text[i] + 1]++;
        for (int c = 0; c < alphabet; c++) start[c + 1] += start[c];

        /// Suffix i is S-type when smaller than suffix i + 1, else L-type, and LMS when S-type right after an L-type
        /// Visits the LMS positions right to left, n - 1 is L-type as the sentinel is smaller
        auto for_each_lms = [&](auto visit){
            bool next_s = false;
            for (int i = n - 2; i >= 0; i--){
                bool s = text[i] < text[i + 1] || (text[i] == text[i + 1] && next_s);
                if (!s && next_s) visit(i + 1);
                next_s = s;
            }
        };

        /// Seeds at the ends of their buckets induce every L-type suffix left to right, then every S-type one right to left
        /// p is stored as ~p when the other pass must place p - 1, which is S-type iff text[p - 1] < text[p], or equal and p is S-type
        auto induce = [&](){
            copy(start.begin(), start.end() - 1, bucket.begin());
            auto place_l = [&](int p){ sa[bucket[text[p]]++] = p > 0 && text[p - 1] < text[p] ? ~p : p; };
            place_l(n - 1);
            for (int i = 0; i < n; i++){
                if (sa[i] > 0) place_l(sa[i] - 1);
            }

            copy(start.begin() + 1, start.end(), bucket.begin());
            for (int i = n - 1; i >= 0; i--){
                if (sa[i] >= 0) continue;

                sa[i] = ~sa[i];
                int p = sa[i] - 1;
                sa[--bucket[text[p]]] = p > 0 && text[p - 1] <= text[p] ? ~p : p;
            }
        };

        /// Sort the LMS substrings, each running from an LMS position to the next one, by inducing from unsorted seeds
        fill(sa, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        int m = 0;
        for_each_lms([&](int p){
            sa[--bucket[text[p]]] = p;
            m++;
        });
        induce();
        if (m == 0) return;  /// the text never rises, so the L pass alone sorted it

        /// p is LMS iff text[p - 1] > text[p] and the run of text[p] ends below a larger value, the runs are disjoint so this is O(n)
        int sorted = 0;
        for (int i = 0; i < n; i++){
            int p = sa[i];
            if (p == 0 || text[p - 1] <= text[p]) continue;

            int end = p + 1;
            while (end < n && text[end] == text[p]) end++;
            if (end < n && text[end] > text[p]) sa[sorted++] = p;
        }

        /// LMS positions are at least 2 apart, so sa[m + p / 2] is a free slot of its own for each, first holding its substring length
        fill(sa + m, sa + n, 0);
        int next = 0;
        for_each_lms([&](int p){
            sa[m + p / 2] = next ? next - p : 0;  /// 0 for the last one, which ends at the sentinel and equals no other
            next = p;
        });

        /// Name each LMS substring by its rank among the distinct ones, from 1 so 0 still marks a free slot
        int names = 0, prev = -1, prev_len = 0;
        for (int i = 0; i < m; i++){
            int p = sa[i], len = sa[m + p / 2];
            bool same = prev >= 0 && len == prev_len;
            for (int d = 0; same && d <= len; d++) same = text[p + d] == text[prev + d];
            if (!same) names++;

            sa[m + p / 2] = names;
            prev = p, prev_len = len;
        }

        /// The names in text order are a reduced text of length m <= n / 2 that orders its suffixes as the LMS suffixes,
        /// kept in sa[n - m .. n) while its suffix array fills sa[0 .. m)
        int* reduced = sa + n - m;
        int filled = m;
        for (int i = n - 1; i >= m; i--){
            if (sa[i] > 0) reduced[--filled] = sa[i] - 1;
        }

        if (names == m){
            for (int i = 0; i < m; i++) sa[reduced[i]] = i;  /// distinct names are already the ranks
        }
        else{
            sais(reduced, sa, m, names);
        }

        /// Map the sorted reduced suffixes back to LMS positions and seed them in order at their bucket ends
        filled = m;
        for_each_lms([&](int p){ reduced[--filled] = p; });
        for (int i = 0; i < m; i++) sa[i] = reduced[sa[i]];

        fill(sa + m, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        for (int i = m - 1; i >= 0; i--){  /// the i-th smallest lands at or after i, so no unread seed is overwritten
            int p = sa[i];
            sa[i] = 0;
            sa[--bucket[text[p]]] = p;
        }
        induce();
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

struct SuffixArray{
    int n;
    vector<int> sa;
    vector<int> lcp;

    long long distinct_substrings() const {
        long long total = ((long long)n * (n + 1)) / 2;
        for (int i = 0; i + 1 < n; i++) total -= lcp[i];
        return total;
    }
};

template<typename Container>
SuffixArray suffix_array(const Container& c){
    int n = c.size();
    SuffixArray result;
    result.n = n;

    if (n == 0){
        return result;
    }

    /// SA-IS wants values in [0, alphabet) with one bucket each: shift small ranges, compress large ones
    /// Values stay in their own type: a cast to long long would wrap unsigned 64-bit values past LLONG_MAX below small ones
    using T = typename Container::value_type;
    using U = make_unsigned_t<T>;
    T lo = *min_element(c.begin(), c.end()), hi = *max_element(c.begin(), c.end());
    auto offset = [&](T x){ return (unsigned long long)U(U(x) - U(lo)); };  /// modular unsigned difference, exact as x >= lo

    bool compress = offset(hi) > (unsigned long long)n + 256;
    vector<T> values;
    if (compress){
        values.assign(c.begin(), c.end());
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());
    }
    auto value_rank = [&](T x){ return compress ? int(lower_bound(values.begin(), values.end(), x) - values.begin()) : int(offset(x)); };
    int alphabet = compress ? values.size() : offset(hi) + 1;

    result.sa.resize(n);
    auto build = [&](auto text){
        for (int i = 0; i < n; i++) text[i] = value_rank(c[i]);
        SuffixArraySAIS::sais(text.data(), result.sa.data(), n, alphabet);
    };

    /// The induce passes read the text at random, bytes instead of ints measured 30% faster at n = 1e7 with 26 letters
    if (alphabet <= 256) build(vector<unsigned char>(n));
    else build(vector<int>(n));

    result.lcp.resize(n);
    SuffixArraySAIS::build_lcp(c, n, result.sa, result.lcp);

    return result;
}

SuffixArray suffix_array(const char* s){
    return suffix_array(string(s));
}

template<typename Container>
struct SuffixArrayQueries : SuffixArray{
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
    /// Suffixes of banana in order: a, ana, anana, banana, na, nana
    SuffixArray res = suffix_array("banana");
    assert((res.sa == vector<int>{5, 3, 1, 0, 4, 2}));
    assert((res.lcp == vector<int>{1, 3, 0, 0, 2, 0}));
    assert(res.distinct_substrings() == 15);

    /// Any integer container and any values: -5 0 < -5 0 -5 0 < 0 < 0 -5 0
    SuffixArray nums = suffix_array(vector<long long>{-5, 0, -5, 0});
    assert((nums.sa == vector<int>{2, 0, 3, 1}));

    /// Occurrence ranges index sa: "ana" starts at sa[1] = 3 and sa[2] = 1
    SuffixArrayQueries q("banana");
    assert(q.rank[0] == 3);
    assert(q.common_prefix(1, 3) == 3);
    assert((q.occurrences("ana") == pair<int, int>{1, 3}));
    assert((q.occurrences("na") == pair<int, int>{4, 6}));
    assert((q.occurrences("x") == pair<int, int>{6, 6}));
    assert((q.occurrences(1, 3) == pair<int, int>{1, 3}));

    return 0;
}
