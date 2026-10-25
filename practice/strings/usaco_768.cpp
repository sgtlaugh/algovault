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
 * res.rank[p]: position of the suffix starting at p in sa
 * res.lcp[i]: longest common prefix of the suffixes at sa[i] and sa[i + 1], lcp[n - 1] = 0
 * res.distinct_substrings(): number of distinct non-empty substrings, O(n)
 *
 * SuffixArrayQueries q(c): everything in SuffixArray plus the calls below
 * q.common_prefix(i, j): longest common prefix of the suffixes starting at i and j
 * q.occurrences(p, len): [lo, hi) such that sa[lo .. hi) are exactly the starts of s[p .. p + len), 0 <= p, 0 <= len, p + len <= n
 * q.occurrences(pattern): the same range for any pattern
 *
 * The occurrence count is hi - lo, an absent pattern gives lo == hi at its insertion point, an empty one gives [0, n)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Measured against DC3 at n = 5e5 to 1e6: 2-3.5x faster on repetitive text, 10-30% on random letters, 3-8% slower on all distinct values
namespace SuffixArraySAIS{
    /// Sorts the suffixes of text, values in [0, alphabet), as if a sentinel below every value ended it
    template<typename Char>
    void sais(const Char* text, int* sa, int n, int alphabet){
        if (n == 1){
            sa[0] = 0;
            return;
        }

        vector<int> start(alphabet + 1, 0), bucket(alphabet);  /// bucket c spans sa[start[c] .. start[c + 1])
        for (int i = 0; i < n; i++) start[text[i] + 1]++;
        for (int c = 0; c < alphabet; c++) start[c + 1] += start[c];

        /// LMS: an S-type suffix (smaller than the next one) right after an L-type one, visited right to left
        auto for_each_lms = [&](auto visit){
            bool next_s = false;
            for (int i = n - 2; i >= 0; i--){
                bool s = text[i] < text[i + 1] || (text[i] == text[i + 1] && next_s);
                if (!s && next_s) visit(i + 1);
                next_s = s;
            }
        };

        /// ~p marks p - 1 for the other pass, which replaces a type array: p - 1 has the type of p when text[p - 1] == text[p]
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

        /// Inducing from unsorted LMS seeds sorts the LMS substrings
        fill(sa, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        int m = 0;
        for_each_lms([&](int p){
            sa[--bucket[text[p]]] = p;
            m++;
        });
        induce();
        if (m == 0) return;  /// non-increasing text, the L pass sorted it

        int sorted = 0;
        for (int i = 0; i < n; i++){
            int p = sa[i];
            if (p == 0 || text[p - 1] <= text[p]) continue;

            int end = p + 1;
            while (end < n && text[end] == text[p]) end++;
            if (end < n && text[end] > text[p]) sa[sorted++] = p;
        }

        /// LMS positions are at least 2 apart, so sa[m + p / 2] is a free slot for each, holding its length, then its name
        fill(sa + m, sa + n, 0);
        int next = 0;
        for_each_lms([&](int p){
            sa[m + p / 2] = next ? next - p : 0;  /// the last one ends at the sentinel and equals no other
            next = p;
        });

        /// Names start at 1, as 0 marks a free slot
        int names = 0, prev = -1, prev_len = 0;
        for (int i = 0; i < m; i++){
            int p = sa[i], len = sa[m + p / 2];
            bool same = prev >= 0 && len == prev_len;
            for (int d = 0; same && d <= len; d++) same = text[p + d] == text[prev + d];
            if (!same) names++;

            sa[m + p / 2] = names;
            prev = p, prev_len = len;
        }

        /// The names in text order sort like the LMS suffixes, m <= n / 2 so they and their suffix array both fit in sa
        int* reduced = sa + n - m;
        int filled = m;
        for (int i = n - 1; i >= m; i--){
            if (sa[i] > 0) reduced[--filled] = sa[i] - 1;
        }

        if (names == m){
            for (int i = 0; i < m; i++) sa[reduced[i]] = i;
        }
        else{
            sais(reduced, sa, m, names);
        }

        filled = m;
        for_each_lms([&](int p){ reduced[--filled] = p; });
        for (int i = 0; i < m; i++) sa[i] = reduced[sa[i]];

        fill(sa + m, sa + n, 0);
        copy(start.begin() + 1, start.end(), bucket.begin());
        for (int i = m - 1; i >= 0; i--){  /// the i-th smallest lands at or after i, past every unread seed
            int p = sa[i];
            sa[i] = 0;
            sa[--bucket[text[p]]] = p;
        }
        induce();
    }

    template<typename Container>
    vector<int> build_lcp(const Container& c, const vector<int>& sa, const vector<int>& rank){
        int n = sa.size(), k = 0;
        vector<int> lcp(n);
        for (int i = 0; i < n; i++){
            if (rank[i] == n - 1){
                k = 0;
                continue;
            }

            int j = sa[rank[i] + 1];
            while (i + k < n && j + k < n && c[i + k] == c[j + k]) k++;
            lcp[rank[i]] = k;
            if (k > 0) k--;
        }
        return lcp;
    }
}

struct SuffixArray{
    int n;
    vector<int> sa, rank, lcp;

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

    /// Values are ranked in T, a cast to long long would wrap unsigned values past LLONG_MAX
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

    /// Byte text measured 30% faster than int at n = 1e7, the induce passes read it at random
    if (alphabet <= 256) build(vector<unsigned char>(n));
    else build(vector<int>(n));

    result.rank.resize(n);
    for (int i = 0; i < n; i++) result.rank[result.sa[i]] = i;
    result.lcp = SuffixArraySAIS::build_lcp(c, result.sa, result.rank);

    return result;
}

SuffixArray suffix_array(const char* s){
    return suffix_array(string(s));
}

template<typename Container>
struct SuffixArrayQueries : SuffixArray{
    Container text;
    vector<vector<int>> table;  /// table[k][i] = min(lcp[i .. i + 2^k))

    SuffixArrayQueries(const Container& c) : SuffixArray(suffix_array(c)), text(c){
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

    /// First rank whose suffix cut to |pattern| is >= pattern, or > when strict
    int bound(const Container& pattern, bool strict) const {
        int lo = 0, hi = n;
        while (lo < hi){
            int mid = (lo + hi) / 2, cmp = compare(sa[mid], pattern);
            if (cmp < 0 || (strict && cmp == 0)) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    /// Sign of the suffix at p cut to |pattern| minus pattern, a suffix that ends first is smaller
    int compare(int p, const Container& pattern) const {
        int m = pattern.size();
        for (int k = 0; k < m; k++){
            if (p + k == n) return -1;
            if (text[p + k] != pattern[k]) return text[p + k] < pattern[k] ? -1 : 1;
        }
        return 0;
    }

    int lcp_min(int l, int r) const {
        int k = __lg(r - l);
        return min(table[k][l], table[k][r - (1 << k)]);
    }
};

SuffixArrayQueries(const char*) -> SuffixArrayQueries<string>;

int main(){
    freopen("standingout.in", "r", stdin);
    freopen("standingout.out", "w", stdout);

    int cows;
    cin >> cows;
    vector<int> c, owner, end;
    for (int i = 0; i < cows; i++){
        string name;
        cin >> name;
        int stop = c.size() + name.size();
        for (char ch : name){
            c.push_back(ch);
            owner.push_back(i);
            end.push_back(stop);
        }
        c.push_back(-1 - i);
        owner.push_back(-1);
        end.push_back(stop);
    }

    SuffixArray res = suffix_array(c);
    int n = c.size();

    vector<int> below(n, 0);
    for (int i = n - 2; i >= 0; i--){
        int a = owner[res.sa[i]], b = owner[res.sa[i + 1]];
        below[i] = a != b ? res.lcp[i] : min(res.lcp[i], below[i + 1]);
    }

    vector<long long> unique_count(cows, 0);
    for (int i = 0; i < n; i++){
        int p = res.sa[i];
        if (owner[p] < 0) continue;

        int taken = max(below[i], i ? res.lcp[i - 1] : 0);
        unique_count[owner[p]] += max(0, end[p] - p - taken);
    }

    for (long long x : unique_count) printf("%lld\n", x);
    return 0;
}
