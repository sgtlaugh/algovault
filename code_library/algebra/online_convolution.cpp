/***
 *
 * Online Convolution
 * A sequence f whose every term depends on the convolution of the earlier terms with a sequence g known in advance
 *
 * Complexity: O(n log^2 n) with an O(n log n) multiply, O(n) calls to it
 *
 * online_convolution(n, g, mod, multiply, next): f_i = next(i, s_i) for i < n, s_i = sum over j < i of f_j g_(i-j)
 *
 * Notes:
 *   - multiply(a, b) returns the product of two vector<long long> modulo mod in [0, mod), like ntt::multiply
 *   - multiply may drop trailing zeros, like Polynomial's operator*, so ntt::mod_multiply or an FFT plug in as well
 *   - next is called once per i in increasing order with s_i in [0, mod), so it may keep state between calls
 *   - 1 <= mod < 2^31, g_0 is never read, missing g terms are 0, and values of g and next are reduced into [0, mod)
 *   - g = f itself, as in f_i = next(i, (f * f)_i), is not supported: g must be known before f is built
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED online_convolution
template<class Multiply, class Next>
vector<long long> online_convolution(int n, vector<long long> g, long long mod, Multiply multiply, Next next){
    assert(n >= 0 && 1 <= mod && mod < (1LL << 31));
    const int LEAF = 32;

    g.resize(n, 0);
    for (auto& x : g) x = (x % mod + mod) % mod;
    vector<long long> f(n, 0), s(n, 0);

    auto solve = [&](auto& self, int l, int r) -> void{
        if (r - l <= LEAF){
            for (int i = l; i < r; i++){
                for (int j = l; j < i; j++) s[i] = (s[i] + f[j] * g[i - j]) % mod;
                f[i] = (next(i, s[i]) % mod + mod) % mod;
            }

            return;
        }

        int mid = (l + r) / 2;
        self(self, l, mid);

        /// f_j for j in [l, mid) times g_(i-j) lands on index i - l of the product, for every i in [mid, r)
        vector<long long> prod = multiply(vector<long long>(f.begin() + l, f.begin() + mid),
                                          vector<long long>(g.begin(), g.begin() + (r - l)));
        /// A multiply that trims trailing zeros, like Polynomial's, returns fewer terms
        for (int i = mid; i < r && i - l < (int)prod.size(); i++) s[i] = (s[i] + prod[i - l]) % mod;

        self(self, mid, r);
    };

    solve(solve, 0, n);
    return f;
}
/// END SHARED online_convolution

int main(){
    /// Fibonacci modulo 7 from f = 1 + x f + x^2 f, a schoolbook multiply standing in for ntt::multiply
    auto multiply = [](const vector<long long>& a, const vector<long long>& b){
        vector<long long> res(a.size() + b.size() - 1, 0);
        for (size_t i = 0; i < a.size(); i++){
            for (size_t j = 0; j < b.size(); j++) res[i + j] = (res[i + j] + a[i] * b[j]) % 7;
        }

        return res;
    };

    auto fib_next = [](int i, long long s){ return i ? s : 1; };
    vector<long long> fib = online_convolution(10, {0, 1, 1}, 7, multiply, fib_next);
    assert((fib == vector<long long>{1, 1, 2, 3, 5, 1, 6, 0, 6, 6}));

    /// 1 / (2 - x) modulo 7: f_0 = 1 / 2 = 4 and f_i = -s_i / 2 = 4 f_(i-1), negative g and next values reduced
    auto inv_next = [](int i, long long s){ return i ? -s * 4 : 4; };
    vector<long long> inv = online_convolution(6, {2, -1}, 7, multiply, inv_next);
    assert((inv == vector<long long>{4, 2, 1, 4, 2, 1}));

    return 0;
}
