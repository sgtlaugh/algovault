#include <bits/stdtr1c++.h>
using namespace std;
namespace nt{
    int extended_gcd(int a, int b, int& x, int& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }
        int g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }
    int mod_inverse(int a, int mod){
        int x, y;
        extended_gcd(a, mod, x, y);
        return (x + mod) % mod;
    }
    int chinese_remainder_theorem(vector<int> ar, vector<int> mods){
        int i, res = 0, mod = 1;
        for (auto m: mods) mod *= m;
        for (i = 0; i < (int)ar.size(); i++){
            long long x = mod / mods[i];
            long long y = mod_inverse(x, mods[i]);
            res = (res + (x * ar[i] % mod) * y) % mod;
        }

        return res;
    }
}
namespace fft{
    struct ComplexNum{
        double real, img;
        inline ComplexNum(double real=0.0, double img=0.0) : real(real), img(img) {}
        inline ComplexNum conjugate(){
            return ComplexNum(real, -img);
        }
        inline ComplexNum operator + (ComplexNum other){
            return ComplexNum(real + other.real, img + other.img);
        }
        inline ComplexNum operator - (ComplexNum other){
            return ComplexNum(real - other.real, img - other.img);
        }
        inline ComplexNum operator * (ComplexNum other){
            return ComplexNum((real * other.real) - (img * other.img), (real * other.img) + (img * other.real));
        }
    };
    vector<ComplexNum> roots = {ComplexNum(), ComplexNum(1)};
    void extend_roots(int len){
        for (int k = roots.size(); k < len; k <<= 1){
            roots.resize(2 * k);
            for (int j = 0; j < k; j++){
                long double theta = acosl(-1.0L) * j / k;
                roots[k + j] = ComplexNum(cosl(theta), sinl(theta));
            }
        }
    }
    void transform(vector<ComplexNum>& a){
        int i, j, k, len = a.size();
        for (i = 1, j = 0; i < len; i++){
            for (k = len >> 1; j & k; k >>= 1) j ^= k;
            j ^= k;
            if (i < j) swap(a[i], a[j]);
        }
        for (k = 1; k < len; k <<= 1){
            for (i = 0; i < len; i += (k << 1)){
                for (j = 0; j < k; j++){
                    auto z = a[i + j + k] * roots[j + k];
                    a[i + j + k] = a[i + j] - z;
                    a[i + j] = a[i + j] + z;
                }
            }
        }
    }
    inline void combine(ComplexNum fi, ComplexNum fj, ComplexNum gi, ComplexNum gj, int len, ComplexNum& u, ComplexNum& v){
        auto a1 = (fi + fj.conjugate()) * ComplexNum(0.5, 0);
        auto a2 = (fi - fj.conjugate()) * ComplexNum(0, -0.5);
        auto b1 = (gi + gj.conjugate()) * ComplexNum(0.5 / len, 0);
        auto b2 = (gi - gj.conjugate()) * ComplexNum(0, -0.5 / len);
        v = a1 * b2 + a2 * b1;
        u = a1 * b1 + a2 * b2 * ComplexNum(0, 1);
    }
    vector<int> mod_multiply(const vector<int>& p1, const vector<int>& p2, int mod, int parity = -1){
        int i, j, a = p1.size(), b = p2.size(), p_len = a + b - 1, step = parity < 0 ? 1 : 2;
        while (a > 1 && p1[a - 1] == 0) a--;
        while (b > 1 && p2[b - 1] == 0) b--;
        bool is_equal = (a == b && equal(p1.begin(), p1.begin() + a, p2.begin()));
        int len = 2;
        while (len < a + b) len <<= 1;
        extend_roots(len);
        vector<ComplexNum> f(len), g(len);
        for (i = 0; i < a; i++) f[i] = ComplexNum((p1[i] % mod) & 32767, (p1[i] % mod) >> 15);
        for (i = 0; i < b; i++) g[i] = ComplexNum((p2[i] % mod) & 32767, (p2[i] % mod) >> 15);
        transform(f);
        if (is_equal) g = f;
        else transform(g);
        for (i = 0; i <= len / 2; i++){
            j = (len - i) & (len - 1);
            auto fi = f[i], fj = f[j], gi = g[i], gj = g[j];
            combine(fi, fj, gi, gj, len, f[j], g[j]);
            if (i != j) combine(fj, fi, gj, gi, len, f[i], g[i]);
        }
        if (parity >= 0){
            len >>= 1;
            for (i = 0; i < len; i++){
                if (parity) f[i] = (f[i] - f[i + len]) * roots[len + i], g[i] = (g[i] - g[i + len]) * roots[len + i];
                else f[i] = f[i] + f[i + len], g[i] = g[i] + g[i + len];
            }
            f.resize(len), g.resize(len);
        }
        transform(f);
        transform(g);
        vector<int> res((p_len - max(parity, 0) + step - 1) / step, 0);
        for (i = 0; i < min(len, (int)res.size()); i++){
            long long x = f[i].real + 0.5, y = g[i].real + 0.5, z = f[i].img + 0.5;
            res[i] = (x + ((y % mod) << 15) + ((z % mod) << 30)) % mod;
        }

        return res;
    }
}
struct LinearRecurrence{
    long long mod;
    vector<int> sequence, recurrence;
    LinearRecurrence() {}
    LinearRecurrence(vector<int> sequence, int mod, vector<int> recurrence={}) : mod(mod), sequence(sequence), recurrence(recurrence){
        int n = sequence.size(), k = recurrence.size();
        assert(n >= k && mod >= 1 && 2LL * mod < INT_MAX);
        for (auto& x : this->sequence) x = ((long long)x % mod + mod) % mod;
        for (auto& c : this->recurrence) c = ((long long)c % mod + mod) % mod;
    }
    int nth_term(long long n){
        return nth_terms(n, 1)[0];
    }
    vector<int> nth_terms(long long n, int k){
        int len = recurrence.size();
        if (n + k <= (long long)sequence.size()) return vector<int>(sequence.begin() + n, sequence.begin() + n + k);
        if (mod == 1 || len == 0) return vector<int>(k, 0);
        assert(1 <= k && k <= len);
        if (k == 1) return {bostan_mori(n)};
        auto p = power_of_x(n);
        vector<int> terms = sequence;
        while ((int)terms.size() < len + k - 1){
            long long val = 0;
            int s = terms.size();
            for (int i = 0; i < len; i++) val = (val + (long long)recurrence[i] * terms[s - 1 - i]) % mod;
            terms.push_back(val);
        }
        vector<int> res(k);
        for (int j = 0; j < k; j++){
            long long val = 0;
            for (int i = 0; i < (int)p.size(); i++) val = (val + (long long)p[i] * terms[i + j]) % mod;
            res[j] = val;
        }

        return res;
    }
    int bostan_mori(long long n){
        int k = recurrence.size();
        vector<int> q(k + 1), v(k + 1);
        q[0] = 1;
        for (int i = 1; i <= k; i++) q[i] = (mod - recurrence[i - 1]) % mod;
        auto p = fft::mod_multiply(vector<int>(sequence.begin(), sequence.begin() + k), q, mod);
        p.resize(k);
        for (; n; n >>= 1){
            for (int i = 0; i <= k; i++) v[i] = (i & 1) && q[i] ? mod - q[i] : q[i];
            p = fft::mod_multiply(p, v, mod, n & 1);
            q = fft::mod_multiply(q, v, mod, 0);
        }

        return p[0];
    }
    vector<int> power_of_x(long long n){
        int k = recurrence.size();
        vector<int> charpoly(k + 1, 1);
        for (int j = 0; j < k; j++) charpoly[j] = (mod - recurrence[k - 1 - j]) % mod;
        vector<int> reversed(charpoly.rbegin(), charpoly.rend());
        vector<int> inv = {1};
        while ((int)inv.size() < k - 1){
            int len = min(2 * (int)inv.size(), k - 1);
            auto t = fft::mod_multiply(inv, vector<int>(reversed.begin(), reversed.begin() + len), mod);
            t.resize(len);
            for (auto& x : t) x = x ? mod - x : 0;
            t[0] = (t[0] + 2) % mod;
            inv = fft::mod_multiply(inv, t, mod);
            inv.resize(len);
        }
        vector<int> res = {1};
        for (int bit = __lg(n); bit >= 0; bit--){
            res = fft::mod_multiply(res, res, mod);
            if ((int)res.size() > k){
                int m = res.size() - k;
                auto quot = fft::mod_multiply(vector<int>(res.rbegin(), res.rbegin() + m), vector<int>(inv.begin(), inv.begin() + m), mod);
                quot.resize(m);
                reverse(quot.begin(), quot.end());
                auto sub = fft::mod_multiply(vector<int>(charpoly.begin(), charpoly.begin() + k), quot, mod);
                res.resize(k);
                for (int j = 0; j < k; j++) res[j] = (res[j] - sub[j] + mod) % mod;
            }
            if (n >> bit & 1){
                res.insert(res.begin(), 0);
                if ((int)res.size() > k){
                    long long top = res.back();
                    res.pop_back();
                    for (int j = 0; j < k; j++) res[j] = (res[j] + (mod - charpoly[j]) * top) % mod;
                }
            }
        }

        return res;
    }
};
#define MOD 1111111111
int main(){
    int n, d, i, j, k;
    scanf("%d %d", &n, &d);
    vector<long long> c(d + 1);
    for (i = 0; i <= d; i++) scanf("%lld", &c[i]);
    vector<long long> poly = {1};
    for (i = 0; i <= d; i++){
        poly.insert(poly.begin(), 0);
        for (j = 0; (j + 1) < (int)poly.size(); j++) poly[j] = (poly[j] - poly[j + 1] + MOD) % MOD;
    }
    k = d + 3;
    vector<long long> charpoly(k + 1, 0);
    for (j = 0; j <= d + 1; j++){
        charpoly[j + 2] = (charpoly[j + 2] + poly[j]) % MOD;
        charpoly[j + 1] = (charpoly[j + 1] + MOD - poly[j]) % MOD;
        charpoly[j] = (charpoly[j] + MOD - poly[j]) % MOD;
    }
    vector<int> rec(k), seq(k, 1);
    for (i = 0; i < k; i++) rec[i] = (MOD - charpoly[k - 1 - i]) % MOD;
    for (i = 2; i < k; i++){
        long long val = 0;
        for (j = d; j >= 0; j--) val = (val * i + c[j]) % MOD;
        seq[i] = (val + seq[i - 1] + seq[i - 2]) % MOD;
    }
    vector<int> mods = {11111, 100001}, res;
    for (auto p: mods){
        LinearRecurrence lr(seq, p, rec);
        res.push_back(lr.nth_term(n));
    }
    printf("%d\n", nt::chinese_remainder_theorem(res, mods));
    return 0;
}
