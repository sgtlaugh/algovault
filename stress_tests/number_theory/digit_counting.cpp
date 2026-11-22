#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/digit_counting.cpp"
#undef main

struct Tally{
    __int128 numbers = 0, cnt[10] = {};
};

/// Digit DP over the decimal string of n, an independent reference for values far beyond brute force
Tally reference(const string& s, int pos, bool tight, bool started, map<pair<int, bool>, Tally>& memo){
    if (pos == (int)s.size()){
        Tally t;
        t.numbers = 1;
        if (!started) t.cnt[0] = 1;
        return t;
    }
    if (!tight && memo.count({pos, started})) return memo[{pos, started}];

    Tally res;
    int top = tight ? s[pos] - '0' : 9;
    for (int d = 0; d <= top; d++){
        bool now = started || d > 0;
        Tally sub = reference(s, pos + 1, tight && d == top, now, memo);
        res.numbers += sub.numbers;
        for (int k = 0; k < 10; k++) res.cnt[k] += sub.cnt[k];
        if (now) res.cnt[d] += sub.numbers;
    }

    if (!tight) memo[{pos, started}] = res;
    return res;
}

void check_reference(long long n){
    map<pair<int, bool>, Tally> memo;
    Tally t = reference(to_string(n), 0, true, false, memo);
    array<long long, 10> got = digit_counts(n);
    for (int d = 0; d < 10; d++){
        assert(t.cnt[d] <= LLONG_MAX);
        assert(got[d] == (long long)t.cnt[d]);
    }
}

int main(){
    /// Every [l, r] inside 0..1200 against writing each number out, covering 0, l = r and every power of 10 crossing
    const int limit = 1200;
    vector<array<long long, 10>> prefix(limit + 1);
    for (int x = 0; x <= limit; x++){
        if (x > 0) prefix[x] = prefix[x - 1];
        for (char ch : to_string(x)) prefix[x][ch - '0']++;
    }
    for (int r = 0; r <= limit; r++){
        assert(digit_counts(r) == prefix[r]);
        for (int l = 0; l <= r; l++){
            array<long long, 10> got = digit_counts(l, r);
            for (int d = 0; d < 10; d++) assert(got[d] == prefix[r][d] - (l ? prefix[l - 1][d] : 0));
        }
    }
    assert(digit_counts(-1) == (array<long long, 10>{}));

    for (int it = 0; it < stress::scaled(300); it++){
        long long l = stress::rand_int(0, 1000000), r = l + stress::rand_int(0, 3000);
        array<long long, 10> brute{};
        for (long long x = l; x <= r; x++){
            for (char ch : to_string(x)) brute[ch - '0']++;
        }
        assert(digit_counts(l, r) == brute);
    }

    /// Up to the documented 4e18 against the digit DP, around every power of 10 where a bad place value would show
    for (int k = 0; k <= 18; k++){
        long long p = stoll("1" + string(k, '0'));
        for (long long n : {p - 1, p, p + 1, 2 * p - 1, 4 * p}) check_reference(n);
    }
    for (int it = 0; it < stress::scaled(3000); it++){
        long long hi = it % 2 ? (long long)4e18 : 1LL << stress::rand_int(1, 61);
        check_reference(stress::rand_int(0, hi));
    }

    stress::expect_assert("number_theory/digit_counting.cpp", []{ digit_counts(5, 3); });
    stress::expect_assert("number_theory/digit_counting.cpp", []{ digit_counts(-5, 3); });

    return 0;
}
