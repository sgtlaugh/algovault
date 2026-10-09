#include "../common.h"

#define main library_main
#include "../../code_library/strings/tandem_repeats.cpp"
#undef main

/// Expands the triples, checking each one's bounds and the n * ceil(log2 n) triple bound on the way
template <typename Container>
vector<pair<int, int>> squares(const Container& s){
    int n = s.size(), levels = 0;
    while ((1 << levels) < n) levels++;

    vector<pair<int, int>> res;
    long long triples = 0;
    tandem_repeats(s, [&](int first, int last, int l){
        assert(l >= 1 && 0 <= first && first <= last && last + 2 * l <= n);
        triples++;
        for (int i = first; i <= last; i++) res.push_back({i, l});
    });
    assert(triples <= (long long)n * levels);

    sort(res.begin(), res.end());
    return res;
}

/// O(n^3): compares both halves of every substring of even length
template <typename Container>
vector<pair<int, int>> brute_cubic(const Container& s){
    int n = s.size();
    vector<pair<int, int>> res;
    for (int i = 0; i < n; i++){
        for (int l = 1; i + 2 * l <= n; l++){
            if (equal(s.begin() + i, s.begin() + i + l, s.begin() + i + l)) res.push_back({i, l});
        }
    }
    return res;
}

/// O(n^2): for each period l, a run of l consecutive j with s[j] == s[j + l] starting at i is a square at i
template <typename Container>
vector<pair<int, int>> brute_quadratic(const Container& s){
    int n = s.size();
    vector<pair<int, int>> res;
    for (int l = 1; 2 * l <= n; l++){
        for (int j = n - l - 1, run = 0; j >= 0; j--){
            run = s[j] == s[j + l] ? run + 1 : 0;
            if (run >= l) res.push_back({j, l});
        }
    }
    sort(res.begin(), res.end());
    return res;
}

string random_string(int n, int sigma){
    string s(n, 'a');
    for (auto& c : s) c = 'a' + stress::rand_int(0, sigma - 1);
    return s;
}

int main(){
    for (int n = 0; n <= 13; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s;
            for (int i = 0; i < n; i++) s += mask >> i & 1 ? 'b' : 'a';
            assert(squares(s) == brute_cubic(s));
        }
    }

    for (int n = 0; n <= 9; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 3;
        for (int code = 0; code < total; code++){
            string s;
            for (int i = 0, c = code; i < n; i++, c /= 3) s += char("#\0a"[c % 3]);
            assert(squares(s) == brute_cubic(s));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 100 ? 100 : 1000), sigma = stress::rand_int(1, it % 3 ? 2 : 5);
        string s = random_string(n, sigma);
        if (it % 4 == 0 && n){
            string block = random_string(stress::rand_int(1, 6), sigma);
            for (int i = 0; i < n; i++) s[i] = stress::rand_int(0, 30) ? block[i % block.size()] : 'z';
        }
        assert(squares(s) == brute_quadratic(s));

        vector<int> v(n);
        for (int i = 0; i < n; i++) v[i] = s[i] == 'a' ? INT_MIN : s[i] - 'a' - 1;
        assert(squares(v) == brute_quadratic(v));
    }

    string fib_prev = "a", fib = "ab";
    while (fib.size() < 4000){
        string next = fib + fib_prev;
        fib_prev = fib;
        fib = next;
    }
    assert(squares(fib) == brute_quadratic(fib));

    string thue(4096, 'a');
    for (int i = 0; i < 4096; i++) thue[i] = 'a' + __builtin_popcount(i) % 2;
    assert(squares(thue) == brute_quadratic(thue));

    int n = 200000, levels = 18;
    long long total = 0, triples = 0;
    tandem_repeats(string(n, 'a'), [&](int first, int last, int){
        total += last - first + 1;
        triples++;
    });
    assert(total == (long long)n / 2 * (n / 2));
    assert(triples <= (long long)n * levels);

    return 0;
}
