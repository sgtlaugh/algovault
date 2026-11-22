/***
 *
 * Sparse Power Series
 * Inverse, log, exp, power and square root of a power series with few nonzero terms, modulo any prime
 *
 * Complexity: O(n k) for n coefficients of a series with k nonzero terms, plus O(log mod), O(log^2 mod) for sqrt
 *
 * SparseSeries f(terms, mod = 998244353): the sum of c x^d over the pairs {d, c} of terms, repeated degrees add up
 * f.exp(n): e^f mod x^n, requires a zero constant term
 * f.inverse(n): 1 / f mod x^n, requires a nonzero constant term
 * f.log(n): ln f mod x^n, requires the constant term 1
 * f.pow(k, n): f^k mod x^n for any k >= 0, f may start with zeros, and 0^0 = 1
 * f.sqrt(n): a square root of f mod x^n, lowest nonzero coefficient the smaller modular root, nullopt when none exists
 * f.terms: the nonzero terms {d, c} by increasing degree, c in [1, mod)
 *
 * Notes:
 *   - mod is a prime below 2^31 with n <= mod, so 1..n - 1 have inverses; sqrt also needs mod > 2
 *   - Degrees are >= 0 and coefficients any long long; results have exactly n coefficients in [0, mod)
 *   - Every call reads only the terms of degree below n
 *   - Modulo an NTT prime, polynomial.cpp does the same in O(n log n)
 *   - At n = 10^6 polynomial.cpp is faster once k passes about 40 for inverse, 60 for log and sqrt, 110 for exp and pow
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SparseSeries{
    vector<pair<int, int>> terms;
    int mod;

    SparseSeries(const vector<pair<int, long long>>& raw, int mod = 998244353) : mod(mod){
        assert(mod >= 2);
        vector<pair<int, long long>> sorted = raw;
        sort(sorted.begin(), sorted.end());

        for (auto [d, c] : sorted){
            assert(d >= 0);
            int r = (c % mod + mod) % mod;
            if (!terms.empty() && terms.back().first == d) terms.back().second = (terms.back().second + (long long)r) % mod;
            else terms.push_back({d, r});
        }
        terms.erase(remove_if(terms.begin(), terms.end(), [](auto& t){ return t.second == 0; }), terms.end());
    }

    vector<int> exp(int n) const{
        assert(n >= 0 && (terms.empty() || terms[0].first > 0));
        vector<int> res(n), inv = inverses(n);
        if (n == 0) return res;

        /// From g' = f' g: i g[i] = sum of d c g[i - d]
        res[0] = 1;
        for (int i = 1; i < n; i++){
            long long sum = 0;
            for (auto [d, c] : terms){
                if (d > i) break;
                sum = (sum + (long long)d * c % mod * res[i - d]) % mod;
            }
            res[i] = sum * inv[i] % mod;
        }

        return res;
    }

    vector<int> inverse(int n) const{
        assert(n >= 0 && !terms.empty() && terms[0].first == 0);
        vector<int> res(n);
        if (n == 0) return res;

        long long c0_inv = power(terms[0].second, mod - 2);
        res[0] = c0_inv;
        for (int i = 1; i < n; i++){
            long long sum = 0;
            for (int j = 1; j < (int)terms.size() && terms[j].first <= i; j++){
                sum = (sum + (long long)terms[j].second * res[i - terms[j].first]) % mod;
            }
            res[i] = (mod - sum) * c0_inv % mod;
        }

        return res;
    }

    vector<int> log(int n) const{
        assert(n >= 0 && !terms.empty() && terms[0] == make_pair(0, 1));
        vector<int> res(n), inv = inverses(n);

        /// q = (ln f)' solves f q = f', so q[m] = (m + 1) f[m + 1] - sum of c q[m - d] over d > 0
        vector<long long> q(max(n - 1, 0), 0);
        for (auto [d, c] : terms){
            if (d > 0 && d < n) q[d - 1] = (long long)d * c % mod;
        }
        for (int m = 0; m < n - 1; m++){
            long long sum = q[m];
            for (int j = 1; j < (int)terms.size() && terms[j].first <= m; j++){
                sum = (sum + (long long)(mod - terms[j].second) * q[m - terms[j].first]) % mod;
            }
            q[m] = sum;
        }

        for (int i = 1; i < n; i++) res[i] = q[i - 1] * inv[i] % mod;
        return res;
    }

    vector<int> pow(long long k, int n) const{
        assert(k >= 0 && n >= 0);
        vector<int> res(n);
        if (n == 0) return res;
        if (k == 0){
            res[0] = 1;
            return res;
        }
        if (terms.empty() || terms[0].first >= n) return res;

        int s = terms[0].first, c0 = terms[0].second;
        if (s > 0 && k >= (n + s - 1) / s) return res;

        int shift = s * k;
        vector<int> root = unit_power(k % mod, n - shift, n);
        long long scale = power(c0, k);
        for (int i = shift; i < n; i++) res[i] = root[i - shift] * scale % mod;
        return res;
    }

    optional<vector<int>> sqrt(int n) const{
        assert(n >= 0 && mod > 2);
        vector<int> res(n);
        if (terms.empty() || terms[0].first >= n) return res;

        int s = terms[0].first, r = mod_sqrt(terms[0].second);
        if (s % 2 || r < 0) return nullopt;

        int shift = s / 2;
        vector<int> root = unit_power(mod / 2 + 1, n - shift, n);
        for (int i = shift; i < n; i++) res[i] = (long long)root[i - shift] * r % mod;
        return res;
    }

    /// (f / (c x^s))^e mod x^len for the lowest term c x^s, e a residue mod p, from T g' = e T' g with T(0) = 1
    vector<int> unit_power(long long e, int len, int n) const{
        int s = terms[0].first, top = min(s + len, n);
        long long c0_inv = power(terms[0].second, mod - 2);
        vector<pair<int, long long>> t;
        for (auto [d, c] : terms){
            if (d > s && d < top) t.push_back({d - s, c * c0_inv % mod});
        }

        vector<int> res(len), inv = inverses(len);
        res[0] = 1;
        for (int i = 1; i < len; i++){
            long long sum = 0;
            for (auto [d, c] : t){
                if (d > i) break;
                long long w = (e * d - i + d) % mod;
                if (w < 0) w += mod;
                sum = (sum + c * w % mod * res[i - d]) % mod;
            }
            res[i] = sum * inv[i] % mod;
        }

        return res;
    }

    /// inv[i] = 1 / i for 1 <= i < n, through inv[i] = -(mod / i) inv[mod % i], which needs a prime mod > i
    vector<int> inverses(int n) const{
        assert(n <= mod);
        vector<int> inv(max(n, 2), 1);
        for (int i = 2; i < n; i++) inv[i] = (long long)(mod - mod / i) * inv[mod % i] % mod;
        return inv;
    }

    long long power(long long b, long long e) const{
        long long res = 1;
        b %= mod;

        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }

        return res;
    }

    /// Tonelli-Shanks, -1 when x is not a quadratic residue, otherwise the smaller of the two roots
    int mod_sqrt(int x) const{
        if (power(x, (mod - 1) / 2) != 1) return -1;

        int s = 0, q = mod - 1, z = 2;
        while (q % 2 == 0){
            q /= 2;
            s++;
        }
        while (power(z, (mod - 1) / 2) != mod - 1) z++;

        long long c = power(z, q), t = power(x, q), r = power(x, (q + 1) / 2);
        while (t != 1){
            int i = 0;
            for (long long u = t; u != 1; u = u * u % mod) i++;

            long long b = c;
            for (int j = 0; j < s - i - 1; j++) b = b * b % mod;
            s = i;
            c = b * b % mod;
            t = t * c % mod;
            r = r * b % mod;
        }

        return min<long long>(r, mod - r);
    }
};

int main(){
    const int mod = 998244353, inv2 = 499122177, inv6 = 166374059, minus_inv8 = 124780544;

    /// Partitions into parts 1 and 2: 1 / ((1 - x)(1 - x^2)) = 1 / (1 - x - x^2 + x^3) = 1 + x + 2x^2 + 2x^3 + 3x^4
    SparseSeries f({{0, 1}, {1, -1}, {2, -1}, {3, 1}});
    assert((f.inverse(5) == vector<int>{1, 1, 2, 2, 3}));

    /// ln(1 - x) = -x - x^2 / 2 - ..., e^x = 1 + x + x^2 / 2 + x^3 / 6
    SparseSeries one_minus_x({{0, 1}, {1, -1}});
    assert((one_minus_x.log(3) == vector<int>{0, mod - 1, mod - inv2}));

    SparseSeries x({{1, 1}});
    assert((x.exp(4) == vector<int>{1, 1, inv2, inv6}));

    /// (x + x^2)^3 = x^3 + 3x^4 + 3x^5 + x^6, sqrt(1 + x) = 1 + x / 2 - x^2 / 8, and x has no square root
    SparseSeries g({{1, 1}, {2, 1}});
    assert((g.pow(3, 6) == vector<int>{0, 0, 0, 1, 3, 3}));

    SparseSeries one_plus_x({{0, 1}, {1, 1}});
    assert((one_plus_x.sqrt(3) == vector<int>{1, inv2, minus_inv8}));
    assert((!x.sqrt(3).has_value()));

    /// Any prime works: 1 / (1 - 2x) = 1 + 2x + 4x^2 modulo 7 wraps 8 to 1
    SparseSeries mod7({{0, 1}, {1, -2}}, 7);
    assert((mod7.inverse(4) == vector<int>{1, 2, 4, 1}));

    return 0;
}
