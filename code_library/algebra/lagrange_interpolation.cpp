/***
 *
 * Lagrange Interpolation
 * Value of a polynomial at any x from its values at x = 0, 1, ..., n - 1, modulo a prime
 *
 * Complexity:
 *   - Lagrange(terms, mod): O(n + log mod)
 *   - interpolate: O(n) time and memory per call
 *   - get_terms(terms, mod, l, r): O(n (r - l + 1) + log mod)
 *   - find_degree: O(n (d + 2)) for the degree d it finds, O(n^2) when it returns -1
 *
 * Lagrange lagrange(terms, mod): the polynomial f of degree below n with f(x) = terms[x] for x = 0..n - 1
 * lagrange.interpolate(x): f(x) in [0, mod), for any long long x, negatives included
 * get_terms(terms, mod, l, r): f(l), f(l + 1), ..., f(r), for any ints l <= r
 * find_degree(terms, mod): the lowest degree d of a polynomial through all the terms, or -1 when fewer than d + 2 terms leave it unconfirmed
 *
 * mod must be a prime below 2^31 with n <= mod, which is asserted
 * Terms may be any ints, negatives included, and are reduced modulo mod
 * The zero polynomial has degree 0, so find_degree needs at least 2 terms of it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Lagrange{
    vector<int> terms, weight;
    int mod, n;

    Lagrange(const vector<int>& terms, int mod) : terms(terms), mod(mod){
        n = terms.size();
        assert(n <= mod);
        for (auto& t : this->terms){
            if ((t %= mod) < 0) t += mod;
        }

        int fact = 1;
        for (int i = 1; i < n; i++) fact = (long long)fact * i % mod;

        vector<int> inv_fact(n, expo(fact, mod - 2));
        for (int i = n - 1; i > 0; i--) inv_fact[i - 1] = (long long)inv_fact[i] * i % mod;

        /// The basis denominator, the product of (i - j) over j != i, is (-1)^(n - 1 - i) i! (n - 1 - i)!
        weight.resize(n);
        for (int i = 0; i < n; i++){
            weight[i] = (long long)inv_fact[i] * inv_fact[n - i - 1] % mod * this->terms[i] % mod;
        }
    }

    int interpolate(long long x){
        if (0 <= x && x < n) return terms[x];
        x %= mod;

        /// prefix[i] * suffix[i] is the product of (x - j) over j != i
        vector<int> prefix(n, 1), suffix(n, 1);
        for (int i = 1; i < n; i++){
            prefix[i] = (long long)prefix[i - 1] * (x - i + 1) % mod;
            if (prefix[i] < 0) prefix[i] += mod;
        }

        for (int i = n - 2; i >= 0; i--){
            suffix[i] = (long long)suffix[i + 1] * (x - i - 1) % mod;
            if (suffix[i] < 0) suffix[i] += mod;
        }

        long long res = 0;
        for (int i = 0; i < n; i++){
            int w = (long long)prefix[i] * suffix[i] % mod * weight[i] % mod;
            if ((n - 1 - i) & 1) w = mod - w;
            res += w;
        }

        return res % mod;
    }

    int expo(int a, int b){
        int res = 1;

        while (b){
            if (b & 1) res = (long long)res * a % mod;
            a = (long long)a * a % mod;
            b >>= 1;
        }

        return res;
    }
};

vector<int> get_terms(const vector<int>& terms, int mod, int l, int r){
    Lagrange lagrange(terms, mod);

    vector<int> res;
    for (int i = l; i <= r; i++){
        res.push_back(lagrange.interpolate(i));
    }

    return res;
}

int find_degree(const vector<int>& terms, int mod){
    int n = terms.size();
    assert(n <= mod);

    vector<int> diff(terms);
    for (auto& x : diff){
        if ((x %= mod) < 0) x += mod;
    }

    /// The k-th differences all vanish exactly when k > d, as the d-th ones equal d! times the leading coefficient and d < n <= mod
    /// Exact but slower than matching prefix fits at a few points, which roots at those points fool: 509 ms against 107 ms at n = 30000, d = n - 1
    for (int k = 0; k < n; k++){
        if (all_of(diff.begin(), diff.end() - k, [](int x){ return x == 0; })){
            int d = max(k - 1, 0);
            return n >= d + 2 ? d : -1;
        }

        for (int i = 0; i + k + 1 < n; i++){
            diff[i] = diff[i + 1] - diff[i];
            if (diff[i] < 0) diff[i] += mod;
        }
    }

    return -1;
}

int main(){
    const int mod = 1000000007;

    /// f(x) = 0^2 + 1^2 + ... + x^2 = x (x + 1) (2x + 1) / 6, a cubic, from its first 5 values
    vector<int> terms = {0, 1, 5, 14, 30};
    Lagrange lagrange(terms, mod);
    assert(lagrange.interpolate(5) == 55);
    assert(lagrange.interpolate(6) == 91);
    assert(get_terms(terms, mod, 4, 6) == (vector<int>{30, 55, 91}));
    assert(find_degree(terms, mod) == 3);

    terms.pop_back();
    assert(find_degree(terms, mod) == -1);  /// 4 terms cannot tell a cubic from a quartic

    /// f(x) = -5 - 2x from negative terms, f(2^40) = -5 - 2^41 = 976759843 modulo 10^9 + 7
    Lagrange line({-5, -7, -9}, mod);
    assert(line.interpolate(-1) == mod - 3);
    assert(line.interpolate(1LL << 40) == 976759843);
    assert(find_degree({-5, -7, -9}, mod) == 1);

    return 0;
}
