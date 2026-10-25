/***
 *
 * Lagrange Interpolation
 * Evaluates a polynomial anywhere from its values at x = 0, 1, ..., n - 1, modulo a prime
 *
 * Complexity:
 *   - constructor: O(n + log mod)
 *   - f(x): O(n)
 *   - find_degree: O(n (d + 1)) for the degree d it returns, O(n^2) when it returns -1
 *
 * Lagrange f(terms, mod): the polynomial of degree below n with f(x) = terms[x]
 * f(x): the value at any long long x
 * find_degree(terms, mod): the lowest degree d that fits, -1 if fewer than d + 2 terms
 *
 * Notes:
 *   - mod is a prime below 2^31 and n <= mod
 *   - Terms and x may be negative, every result is in [0, mod)
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

        weight.resize(n);
        for (int i = 0; i < n; i++){
            weight[i] = (long long)inv_fact[i] * inv_fact[n - i - 1] % mod * this->terms[i] % mod;
        }
    }

    int operator()(long long x){
        if (0 <= x && x < n) return terms[x];
        x %= mod;

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

int find_degree(const vector<int>& terms, int mod){
    int n = terms.size();

    assert(n <= mod);

    vector<int> diff(terms);
    for (auto& x : diff){
        if ((x %= mod) < 0) x += mod;
    }

    /// k-th differences all vanish exactly when k > d: the d-th are d! times the leading coefficient and d < n <= mod
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
    Lagrange f(terms, mod);
    assert(f(5) == 55);
    assert(f(6) == 91);
    assert(find_degree(terms, mod) == 3);

    terms.pop_back();
    assert(find_degree(terms, mod) == -1);  /// 4 terms cannot tell a cubic from a quartic

    /// f(x) = -5 - 2x from negative terms, f(2^40) = -5 - 2^41 = 976759843 modulo 10^9 + 7
    Lagrange line({-5, -7, -9}, mod);
    assert(line(-1) == mod - 3);
    assert(line(1LL << 40) == 976759843);
    assert(find_degree({-5, -7, -9}, mod) == 1);

    return 0;
}
