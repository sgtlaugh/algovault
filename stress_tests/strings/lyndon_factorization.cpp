#include "../common.h"

#define main library_main
#include "../../code_library/strings/lyndon_factorization.cpp"
#undef main

long long comparisons = 0;

struct Counted{
    int v;

    bool operator<(const Counted& o) const{
        comparisons++;
        return v < o.v;
    }
};

template <typename Container>
bool is_lyndon(const Container& w){
    int m = w.size();
    for (int i = 1; i < m; i++){
        if (!(w < Container(w.begin() + i, w.end()))) return false;
    }

    return m > 0;
}

template <typename Container>
int brute_min_rotation(const Container& s){
    int n = s.size(), best = 0;
    auto rotation = [&](int i){
        Container r(s.begin() + i, s.end());
        r.insert(r.end(), s.begin(), s.begin() + i);
        return r;
    };

    for (int i = 1; i < n; i++){
        if (rotation(i) < rotation(best)) best = i;
    }

    return best;
}

/// The Lyndon factorization is unique, so these three properties pin down the only correct answer
template <typename Container>
void check(const Container& s){
    int n = s.size();
    auto starts = lyndon_factorization(s);
    assert(starts.empty() == (n == 0));
    if (n) assert(starts[0] == 0);

    Container prev;
    for (int t = 0; t < (int)starts.size(); t++){
        int end = t + 1 < (int)starts.size() ? starts[t + 1] : n;
        assert(starts[t] < end && end <= n);
        Container w(s.begin() + starts[t], s.begin() + end);
        assert(is_lyndon(w));
        if (t) assert(!(prev < w));
        prev = w;
    }

    assert(lyndon_min_rotation(s) == brute_min_rotation(s));
}

/// Each inner step costs at most two < calls, Duval takes at most 4n steps on n elements and min rotation runs on 2n
void check_linear(const vector<int>& v){
    int n = v.size();
    vector<Counted> s(n);
    for (int i = 0; i < n; i++) s[i].v = v[i];

    comparisons = 0;
    lyndon_factorization(s);
    assert(comparisons <= 8LL * n);

    comparisons = 0;
    lyndon_min_rotation(s);
    assert(comparisons <= 16LL * n);
}

int main(){
    for (int n = 0; n <= 12; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s;
            for (int i = 0; i < n; i++) s += mask >> i & 1 ? 'b' : 'a';
            check(s);
        }
    }

    for (int n = 0; n <= 7; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 3;
        for (int code = 0; code < total; code++){
            vector<int> v(n);
            for (int i = 0, c = code; i < n; i++, c /= 3) v[i] = c % 3 - 1;
            check(v);
        }
    }

    for (long long it = 0; it < stress::scaled(10000); it++){
        int n = stress::rand_int(0, it % 50 ? 30 : 150), sigma = stress::rand_int(1, it % 4 ? 3 : 26);
        string s;
        for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, sigma - 1));
        if (it % 5 == 0 && n){
            string block = s.substr(0, stress::rand_int(1, n));
            s.clear();
            while ((int)s.size() < n) s += block;
            if (it % 10 == 0) s.pop_back();
        }
        check(s);
        check(vector<int>(s.begin(), s.end()));
        assert(lyndon_factorization(s.c_str()) == lyndon_factorization(s));  /// the const char* overloads
        assert(lyndon_min_rotation(s.c_str()) == lyndon_min_rotation(s));
    }

    for (long long it = 0; it < stress::scaled(40); it++){
        int n = stress::rand_int(1, 200000), sigma = stress::rand_int(1, 3);
        vector<int> v(n);
        for (auto& x : v) x = stress::rand_int(0, sigma - 1);
        if (it % 2){
            int period = stress::rand_int(1, 50);
            for (int i = period; i < n; i++) v[i] = v[i - period];
            if (it % 4 == 1) v[n - 1] = -1;
        }
        check_linear(v);
    }

    const int big = 500000;
    vector<int> ones(big, 1);
    check_linear(ones);
    assert((int)lyndon_factorization(ones).size() == big);

    vector<int> almost(big, 0);
    almost[big - 1] = 1;
    check_linear(almost);
    assert((lyndon_factorization(almost) == vector<int>{0}));
    assert(lyndon_min_rotation(almost) == 0);

    string ab;
    for (int i = 0; i < big; i++) ab += "ab";
    ab += 'a';
    auto starts = lyndon_factorization(ab);
    assert((int)starts.size() == big + 1);
    for (int i = 0; i <= big; i++) assert(starts[i] == 2 * i);
    assert(lyndon_min_rotation(ab) == 2 * big);

    string late(big, 'a');
    late[123456] = 'b';
    assert(lyndon_min_rotation(late) == 123457);

    return 0;
}
