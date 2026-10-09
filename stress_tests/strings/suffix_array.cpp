#include "../common.h"

#define main library_main
#include "../../code_library/strings/suffix_array.cpp"
#undef main

/// A permutation of the suffixes is sorted iff each adjacent pair is in order, and the pair's lcp k decides that order,
/// so one naive lcp pass checks both arrays without the quadratic compares of sorting repetitive input
template <typename Container>
void check(const Container& c){
    int n = c.size();
    auto res = suffix_array(c);

    assert((int)res.sa.size() == n);
    vector<int> sorted_sa = res.sa;
    sort(sorted_sa.begin(), sorted_sa.end());
    for (int i = 0; i < n; i++) assert(sorted_sa[i] == i);

    long long distinct = (long long)n * (n + 1) / 2;
    for (int i = 0; i + 1 < n; i++){
        int a = res.sa[i], b = res.sa[i + 1], k = 0;
        while (a + k < n && b + k < n && c[a + k] == c[b + k]) k++;
        assert(res.lcp[i] == k && (a + k == n || (b + k < n && c[a + k] < c[b + k])));  /// suffix a ends first or is smaller at k
        distinct -= k;
    }
    assert(res.distinct_substrings() == distinct);
}

template <typename Container>
int naive_prefix(const Container& c, int i, const Container& pattern, int j, int limit){
    int k = 0;
    while (k < limit && i + k < (int)c.size() && j + k < (int)pattern.size() && c[i + k] == pattern[j + k]) k++;
    return k;
}

/// The range must hold exactly the starts where a direct scan finds the pattern, and an absent pattern must sit at its
/// insertion point: every suffix before lo is smaller than the pattern, every suffix from hi on is larger
template <typename Container>
void check_range(const SuffixArrayQueries<Container>& q, const Container& c, const Container& pattern, pair<int, int> range){
    int n = c.size(), m = pattern.size(), lo = range.first, hi = range.second;
    assert(0 <= lo && lo <= hi && hi <= n);

    vector<int> expected, found(q.sa.begin() + lo, q.sa.begin() + hi);
    for (int i = 0; i < n; i++){
        if (naive_prefix(c, i, pattern, 0, m) == m) expected.push_back(i);
    }
    sort(found.begin(), found.end());
    assert(found == expected);

    auto below = [&](int p){
        int k = naive_prefix(c, p, pattern, 0, m);
        return k < m && (p + k == n || c[p + k] < pattern[k]);
    };
    if (lo > 0) assert(below(q.sa[lo - 1]));
    if (hi < n) assert(!below(q.sa[hi]) && naive_prefix(c, q.sa[hi], pattern, 0, m) < m);
}

/// min - 1 or max + 1 of the text, so the insertion point check also runs on values the text never holds
template <typename Container, typename T = typename Container::value_type>
bool absent_value(const Container& c, T& value){
    auto [lo, hi] = minmax_element(c.begin(), c.end());
    bool below = *lo > numeric_limits<T>::min(), above = *hi < numeric_limits<T>::max();
    if (!below && !above) return false;

    value = below && (!above || stress::rand_int(0, 1)) ? T(*lo - 1) : T(*hi + 1);
    return true;
}

template <typename Container>
void check_queries(const Container& c){
    using T = typename Container::value_type;
    int n = c.size();
    SuffixArrayQueries<Container> q(c);

    for (int i = 0; i < n; i++) assert(q.sa[q.rank[i]] == i);

    int pairs = n <= 12 ? n * n : 60;
    for (int t = 0; t < pairs; t++){
        int i = n <= 12 ? t / n : stress::rand_int(0, n - 1), j = n <= 12 ? t % n : stress::rand_int(0, n - 1);
        assert(q.common_prefix(i, j) == naive_prefix(c, i, c, j, n));
    }

    for (int t = 0; t < (n <= 40 ? 12 : 4); t++){
        int p = stress::rand_int(0, n), len = stress::rand_int(0, min(n - p, n <= 40 ? n : 64));  /// long patterns on repetitive input make the naive scan quadratic
        Container sub(c.begin() + p, c.begin() + p + len);
        check_range(q, c, sub, q.occurrences(p, len));
        check_range(q, c, sub, q.occurrences(sub));

        if (n == 0) continue;
        if (stress::rand_int(0, 1)) sub.push_back(c[stress::rand_int(0, n - 1)]);
        if (!sub.empty() && stress::rand_int(0, 2) == 0) sub[stress::rand_int(0, (int)sub.size() - 1)] = c[stress::rand_int(0, n - 1)];
        check_range(q, c, sub, q.occurrences(sub));

        T outside;
        if (!absent_value(c, outside)) continue;
        if (sub.empty() || stress::rand_int(0, 1)) sub.insert(sub.begin() + stress::rand_int(0, (int)sub.size()), outside);
        else sub[stress::rand_int(0, (int)sub.size() - 1)] = outside;
        check_range(q, c, sub, q.occurrences(sub));
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(0, it % 10 ? 40 : 2000), alphabet = stress::rand_int(1, 4), mode = stress::rand_int(0, 3);
        vector<long long> v(n);
        for (auto& x : v){
            long long s = stress::rand_int(0, alphabet - 1);
            if (mode == 0) x = s;                                  /// zeros, DC3's sentinel value
            else if (mode == 1) x = s - 2;                         /// negatives
            else if (mode == 2) x = s * 1000000000000LL;           /// a range far too large for DC3 buckets
            else x = stress::rand_int(-1000000, 1000000);
        }
        check(v);
        check_queries(v);

        string s(n, 'a');
        for (auto& ch : s) ch = (char)stress::rand_int(0, 2) ? 'a' + stress::rand_int(0, alphabet - 1) : (char)stress::rand_int(-128, 127);
        check(s);
        check_queries(s);

        vector<unsigned long long> u(n);
        for (auto& x : u){                                         /// values past LLONG_MAX next to small ones, 64-bit hashes
            unsigned long long s = stress::rand_int(0, alphabet - 1);
            x = mode < 2 ? s << 63 | s : mode == 2 ? ULLONG_MAX - s : (unsigned long long)stress::rand_int(LLONG_MIN, LLONG_MAX);
        }
        check(u);
        check_queries(u);
    }

    return 0;
}
