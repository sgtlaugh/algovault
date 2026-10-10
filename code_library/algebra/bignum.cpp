/***
 *
 * Arbitrary precision signed integer
 *
 * Stored as base 10^9 limbs in little endian order, zero is always non-negative
 * Supports + - * / % with Bignum or long long on either side, comparisons, stream I/O and to_string()
 * Division truncates toward zero and the remainder takes the sign of the dividend, like C++ integers
 *
 * Complexity, with n and m the number of digits:
 *     Addition, subtraction, comparison, I/O: O(n)
 *     Multiplication: O(n * m) for small inputs, Karatsuba O(n^1.58) for mid sizes,
 *                     three prime NTT O((n + m) log(n + m)) once both have 700+ limbs (~6300 digits)
 *     Division: O(n) when the divisor fits in 9 digits, O((n - m) * m) while the divisor or the quotient has under
 *               400 limbs (~3600 digits), O(n log n) beyond with a Newton reciprocal
 *
 * At -O2, multiplying two 2 * 10^6 digit numbers takes ~0.15 s and dividing 2 * 10^6 digits by 10^6 digits ~0.5 s
 * / and % each run the full division
 * Products of more than 2^23 limbs (~7.5 * 10^7 digits) split with Karatsuba until the pieces fit the NTT
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Bignum{
    static const uint32_t BASE = 1000000000;
    static const int KARATSUBA_CUTOFF = 96;
    static const int NTT_CUTOFF = 700;
    static const size_t NTT_MAX = 1 << 23;
    static const int NEWTON_CUTOFF = 400;

    bool neg = false;
    vector<uint32_t> a;

    Bignum(long long v = 0){
        neg = v < 0;
        for (unsigned long long x = neg ? 0ULL - v : v; x; x /= BASE) a.push_back(x % BASE);
    }

    explicit Bignum(const string& s){
        size_t i = 0;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
        assert(i < s.size());

        for (size_t j = s.size(); j > i; j -= min<size_t>(9, j - i)){
            uint32_t limb = 0;
            for (size_t k = j - min<size_t>(9, j - i); k < j; k++){
                assert(isdigit(s[k]));
                limb = limb * 10 + (s[k] - '0');
            }
            a.push_back(limb);
        }

        trim();
    }

    explicit Bignum(const char* s) : Bignum(string(s)) {}

    string to_string() const{
        if (a.empty()) return "0";
        string s = (neg ? "-" : "") + std::to_string(a.back());
        size_t pos = s.size();
        s.resize(pos + 9 * (a.size() - 1));

        for (int i = (int)a.size() - 2; i >= 0; i--, pos += 9){
            for (int k = 8, v = a[i]; k >= 0; k--, v /= 10) s[pos + k] = '0' + v % 10;
        }

        return s;
    }

    friend ostream& operator << (ostream& os, const Bignum& x){
        return os << x.to_string();
    }

    friend istream& operator >> (istream& is, Bignum& x){
        string s;
        if (is >> s) x = Bignum(s);
        return is;
    }

    Bignum operator - () const{
        Bignum res = *this;
        if (!res.a.empty()) res.neg = !res.neg;
        return res;
    }

    friend Bignum operator + (const Bignum& x, const Bignum& y){
        if (x.neg == y.neg) return Bignum(add(x.a, y.a), x.neg);
        if (cmp(x.a, y.a) >= 0) return Bignum(sub(x.a, y.a), x.neg);
        return Bignum(sub(y.a, x.a), y.neg);
    }

    friend Bignum operator - (const Bignum& x, const Bignum& y){
        return x + (-y);
    }

    friend Bignum operator * (const Bignum& x, const Bignum& y){
        return Bignum(mul(x.a, y.a), x.neg != y.neg);
    }

    friend Bignum operator / (const Bignum& x, const Bignum& y){
        return Bignum(divmod(x.a, y.a).first, x.neg != y.neg);
    }

    friend Bignum operator % (const Bignum& x, const Bignum& y){
        return Bignum(divmod(x.a, y.a).second, x.neg);
    }

    Bignum& operator += (const Bignum& y){ return *this = *this + y; }
    Bignum& operator -= (const Bignum& y){ return *this = *this - y; }
    Bignum& operator *= (const Bignum& y){ return *this = *this * y; }
    Bignum& operator /= (const Bignum& y){ return *this = *this / y; }
    Bignum& operator %= (const Bignum& y){ return *this = *this % y; }

    friend bool operator == (const Bignum& x, const Bignum& y){ return x.neg == y.neg && x.a == y.a; }
    friend bool operator != (const Bignum& x, const Bignum& y){ return !(x == y); }
    friend bool operator < (const Bignum& x, const Bignum& y){
        if (x.neg != y.neg) return x.neg;
        return x.neg ? cmp(y.a, x.a) < 0 : cmp(x.a, y.a) < 0;
    }
    friend bool operator > (const Bignum& x, const Bignum& y){ return y < x; }
    friend bool operator <= (const Bignum& x, const Bignum& y){ return !(y < x); }
    friend bool operator >= (const Bignum& x, const Bignum& y){ return !(x < y); }

private:
    typedef vector<uint32_t> Limbs;

    Bignum(Limbs&& limbs, bool negative) : neg(negative), a(move(limbs)){
        trim();
    }

    void trim(){
        while (!a.empty() && !a.back()) a.pop_back();
        if (a.empty()) neg = false;
    }

    static void trim(Limbs& x){
        while (!x.empty() && !x.back()) x.pop_back();
    }

    static int cmp(const Limbs& x, const Limbs& y){
        if (x.size() != y.size()) return x.size() < y.size() ? -1 : 1;
        for (int i = (int)x.size() - 1; i >= 0; i--){
            if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
        }
        return 0;
    }

    static Limbs add(const Limbs& x, const Limbs& y){
        Limbs res(max(x.size(), y.size()) + 1);
        add_to(res, x, 0), add_to(res, y, 0);
        trim(res);
        return res;
    }

    /// Requires x >= y
    static Limbs sub(Limbs x, const Limbs& y){
        sub_from(x, y);
        trim(x);
        return x;
    }

    /// res += y * BASE^offset, res must be large enough to hold the result
    static void add_to(Limbs& res, const Limbs& y, size_t offset){
        uint32_t carry = 0;
        for (size_t i = 0; i < y.size() || carry; i++){
            uint32_t s = res[i + offset] + carry + (i < y.size() ? y[i] : 0);
            carry = s >= BASE;
            res[i + offset] = carry ? s - BASE : s;
        }
    }

    /// x -= y, requires x >= y
    static void sub_from(Limbs& x, const Limbs& y){
        int borrow = 0;
        for (size_t i = 0; i < y.size() || borrow; i++){
            long long d = (long long)x[i] - borrow - (i < y.size() ? y[i] : 0);
            borrow = d < 0;
            x[i] = borrow ? d + BASE : d;
        }
    }

    /// Rows accumulate in 64 bits and fold carries every 18 rows, since 18 * (BASE - 1)^2 < 2^64
    /// Only the window touched since the last fold is normalized, ~1.5x faster than summing whole columns in 128 bits
    static Limbs mul_schoolbook(const Limbs& x, const Limbs& y){
        if (x.empty() || y.empty()) return {};
        size_t n = x.size(), m = y.size();
        vector<unsigned long long> acc(n + m + 1);

        for (size_t i = 0, last = 0; i < n; i++){
            unsigned long long xi = x[i], *r = acc.data() + i;
            for (size_t j = 0; j < m; j++) r[j] += xi * y[j];
            if (i - last == 17 || i == n - 1){
                unsigned long long carry = 0;
                for (size_t k = last; k < i + m || carry; k++){
                    unsigned long long t = acc[k] + carry;
                    acc[k] = t % BASE, carry = t / BASE;
                }
                last = i + 1;
            }
        }

        Limbs res(acc.begin(), acc.end() - 1);
        trim(res);
        return res;
    }

    template<uint32_t MOD>
    static uint32_t pow_mod(unsigned long long b, unsigned long long e){
        unsigned long long res = 1;
        for (b %= MOD; e; e >>= 1, b = b * b % MOD){
            if (e & 1) res = res * b % MOD;
        }
        return res;
    }

    /// rt[len + j] = w^j for the primitive 2 * len-th root of unity w (inverted if inverse), len = 1, 2, 4, ... < n
    template<uint32_t MOD>
    static vector<uint32_t> ntt_roots(size_t n, bool inverse){
        vector<uint32_t> rt(max<size_t>(n, 2));
        for (size_t len = 1; len < n; len <<= 1){
            unsigned long long w = pow_mod<MOD>(3, (MOD - 1) / (2 * len));
            if (inverse) w = pow_mod<MOD>(w, MOD - 2);

            rt[len] = 1;
            for (size_t j = 1; j < len; j++) rt[len + j] = rt[len + j - 1] * w % MOD;
        }

        return rt;
    }

    /// The forward pass (decimation in frequency) leaves f in bit reversed order and the inverse pass (decimation in time)
    /// takes it in that order, so pointwise products between them need no bit reversal permutation
    template<uint32_t MOD>
    static void ntt(vector<uint32_t>& f, const vector<uint32_t>& rt, bool inverse){
        size_t n = f.size();
        if (!inverse){
            for (size_t len = n >> 1; len; len >>= 1){
                for (size_t i = 0; i < n; i += 2 * len){
                    for (size_t j = 0; j < len; j++){
                        uint32_t u = f[i + j], v = f[i + j + len];
                        f[i + j] = u + v >= MOD ? u + v - MOD : u + v;
                        f[i + j + len] = (unsigned long long)(u + MOD - v) * rt[len + j] % MOD;
                    }
                }
            }
            return;
        }

        for (size_t len = 1; len < n; len <<= 1){
            for (size_t i = 0; i < n; i += 2 * len){
                for (size_t j = 0; j < len; j++){
                    uint32_t u = f[i + j], v = (unsigned long long)f[i + j + len] * rt[len + j] % MOD;
                    f[i + j] = u + v >= MOD ? u + v - MOD : u + v;
                    f[i + j + len] = u >= v ? u - v : u + MOD - v;
                }
            }
        }
    }

    /// Cyclic convolution of x and y modulo MOD with n points, n a power of 2 dividing MOD - 1
    template<uint32_t MOD>
    static vector<uint32_t> convolve_mod(const Limbs& x, const Limbs& y, size_t n){
        vector<uint32_t> fx(n), fy(n);
        for (size_t i = 0; i < x.size(); i++) fx[i] = x[i] % MOD;
        for (size_t i = 0; i < y.size(); i++) fy[i] = y[i] % MOD;

        vector<uint32_t> rt = ntt_roots<MOD>(n, false);
        ntt<MOD>(fx, rt, false), ntt<MOD>(fy, rt, false);

        unsigned long long inv_n = pow_mod<MOD>(n, MOD - 2);
        for (size_t i = 0; i < n; i++) fx[i] = (unsigned long long)fx[i] * fy[i] % MOD * inv_n % MOD;

        ntt<MOD>(fx, ntt_roots<MOD>(n, true), true);
        return fx;
    }

    /// Three prime NTT, each with primitive root 3 and 2^23 | p - 1
    /// Every coefficient is below min(n, m) * BASE^2 <= 2^22 * 10^18 < P1 * P2 * P3 ~ 7.9 * 10^25, so the CRT is exact
    static Limbs mul_ntt(const Limbs& x, const Limbs& y){
        const uint32_t P1 = 998244353, P2 = 167772161, P3 = 469762049;
        size_t len = x.size() + y.size() - 1, n = 1;
        while (n < len) n <<= 1;

        vector<uint32_t> r1 = convolve_mod<P1>(x, y, n), r2 = convolve_mod<P2>(x, y, n), r3 = convolve_mod<P3>(x, y, n);

        const unsigned long long inv1 = pow_mod<P2>(P1, P2 - 2), p12 = (unsigned long long)P1 * P2;
        const unsigned long long inv12 = pow_mod<P3>(p12, P3 - 2);
        Limbs res(len + 1);
        unsigned __int128 carry = 0;
        for (size_t i = 0; i < len; i++){
            unsigned long long a1 = r1[i];
            unsigned long long a2 = (r2[i] + P2 - a1 % P2) * inv1 % P2;
            unsigned long long a3 = (r3[i] + P3 - (a1 + a2 * P1) % P3) * inv12 % P3;
            carry += a1 + a2 * P1 + (unsigned __int128)a3 * p12;
            res[i] = carry % BASE, carry /= BASE;
        }

        res[len] = carry;
        trim(res);
        return res;
    }

    static Limbs mul(const Limbs& x, const Limbs& y){
        if (x.size() < y.size()) return mul(y, x);
        size_t n = x.size(), m = y.size();
        if (m <= 1){
            if (!m) return {};
            Limbs res = mul_small(x, y[0]);
            trim(res);
            return res;
        }
        if (m < KARATSUBA_CUTOFF) return mul_schoolbook(x, y);
        if (m >= NTT_CUTOFF && n + m - 1 <= NTT_MAX) return mul_ntt(x, y);

        /// Beyond NTT_MAX the split below keeps recursing until the pieces fit the transform
        Limbs res(n + m + 1);
        if (n >= 2 * m){
            /// Unbalanced, multiply y by each m limb chunk of x
            for (size_t i = 0; i < n; i += m){
                Limbs chunk(x.begin() + i, x.begin() + min(n, i + m));
                trim(chunk);
                add_to(res, mul(chunk, y), i);
            }
            trim(res);
            return res;
        }

        /// Karatsuba: x * y = z2 * B^2h + z1 * B^h + z0 with z1 = (x0 + x1)(y0 + y1) - z0 - z2
        size_t h = n / 2;
        Limbs x0(x.begin(), x.begin() + h), x1(x.begin() + h, x.end());
        Limbs y0(y.begin(), y.begin() + h), y1(y.begin() + h, y.end());
        trim(x0), trim(y0);

        Limbs z0 = mul(x0, y0), z2 = mul(x1, y1), z1 = mul(add(x0, x1), add(y0, y1));
        sub_from(z1, z0), sub_from(z1, z2);
        trim(z1);

        add_to(res, z0, 0), add_to(res, z1, h), add_to(res, z2, 2 * h);
        trim(res);
        return res;
    }

    static uint32_t divmod_small(Limbs& x, uint32_t d){
        unsigned long long rem = 0;
        for (int i = (int)x.size() - 1; i >= 0; i--){
            unsigned long long cur = rem * BASE + x[i];
            x[i] = cur / d, rem = cur % d;
        }

        trim(x);
        return rem;
    }

    static Limbs mul_small(const Limbs& x, uint32_t f){
        Limbs res(x.size() + 1);
        unsigned long long carry = 0;
        for (size_t i = 0; i < x.size(); i++){
            unsigned long long cur = (unsigned long long)x[i] * f + carry;
            res[i] = cur % BASE, carry = cur / BASE;
        }

        res[x.size()] = carry;
        return res;
    }

    /// Returns {|x| / |y|, |x| % |y|}
    static pair<Limbs, Limbs> divmod(const Limbs& x, const Limbs& y){
        assert(!y.empty());
        if (cmp(x, y) < 0) return {{}, x};
        if (y.size() == 1){
            Limbs q = x;
            uint32_t r = divmod_small(q, y[0]);
            return {q, r ? Limbs{r} : Limbs{}};
        }

        size_t quotient_limbs = x.size() - y.size() + 1;
        if (y.size() >= NEWTON_CUTOFF && quotient_limbs >= NEWTON_CUTOFF) return divmod_newton(x, y);
        return divmod_knuth(x, y);
    }

    static Limbs shift_left(const Limbs& x, size_t k){
        Limbs res(k + x.size());
        copy(x.begin(), x.end(), res.begin() + k);
        return res;
    }

    static Limbs shift_right(const Limbs& x, size_t k){
        return k >= x.size() ? Limbs{} : Limbs(x.begin() + k, x.end());
    }

    /// floor(BASE^(2L) / y) for y with L limbs, give or take a few units
    /// R0 = r * BASE^(L - h), with r the reciprocal of the top h limbs, has relative error ~BASE^(1 - h) when the top limb
    /// is 1. The Newton step R = R0 + R0 * (BASE^(2L) - y * R0) / BASE^(2L) squares it, and since R <= BASE^(L + 1),
    /// 2h >= L + 3 leaves about one unit of error, plus one each from truncating e and the floor below
    static Limbs reciprocal(const Limbs& y){
        size_t L = y.size();
        if (L <= NEWTON_CUTOFF){
            Limbs power(2 * L + 1);
            power[2 * L] = 1;
            return divmod_knuth(power, y).first;
        }

        size_t h = L / 2 + 2;
        Limbs r = reciprocal(shift_right(y, L - h));

        /// With R0 = r * BASE^(L - h): R = R0 + r * e / BASE^(2h) where e = BASE^(L + h) - y * r, |e| < ~BASE^(L + 1)
        /// Dropping the low h - 1 limbs of e costs at most 1 unit since r < BASE^(h + 1)
        Limbs power(L + h + 1), t = mul(y, r);
        power[L + h] = 1;
        bool below = cmp(t, power) > 0;
        Limbs e = shift_right(below ? sub(t, power) : sub(power, t), h - 1);
        Limbs correction = shift_right(mul(r, e), h + 1);

        Limbs res = shift_left(r, L - h);
        return below ? sub(res, correction) : add(res, correction);
    }

    static pair<Limbs, Limbs> divmod_newton(const Limbs& x, const Limbs& y){
        size_t n = x.size(), m = y.size(), L = n - m + 2;
        if (m >= L){
            /// Dropping the low m - L limbs of both moves the quotient by at most 1, since it is below BASE^(L - 1)
            /// and the shortened divisor is at least BASE^(L - 1)
            Limbs ys = shift_right(y, m - L);
            return fix_quotient(x, y, shift_right(mul(shift_right(x, m - L), reciprocal(ys)), 2 * L));
        }

        /// Long division in base BASE^m: each window is below y * BASE^m, so its quotient fits in one m limb digit
        Limbs inv = reciprocal(y), q(n), r;
        for (size_t b = (n - 1) / m + 1; b-- > 0;){
            Limbs cur = shift_left(r, m);
            copy(x.begin() + b * m, x.begin() + min(n, (b + 1) * m), cur.begin());
            trim(cur);

            pair<Limbs, Limbs> step = fix_quotient(cur, y, shift_right(mul(cur, inv), 2 * m));
            copy(step.first.begin(), step.first.end(), q.begin() + b * m);
            r = move(step.second);
        }

        trim(q);
        return {q, r};
    }

    /// Steps an estimate q, off by a few units, to the exact quotient using the remainder x - q * y
    static pair<Limbs, Limbs> fix_quotient(const Limbs& x, const Limbs& y, Limbs q){
        Limbs p = mul(q, y);
        while (cmp(p, x) > 0){
            sub_from(q, {1}), trim(q);
            sub_from(p, y), trim(p);
        }

        Limbs r = sub(x, p);
        while (cmp(r, y) >= 0){
            q = add(q, {1});
            sub_from(r, y), trim(r);
        }

        return {q, r};
    }

    /// Knuth's algorithm D, requires |x| >= |y| and y with at least 2 limbs
    static pair<Limbs, Limbs> divmod_knuth(const Limbs& x, const Limbs& y){
        /// Scaling so the top divisor limb is at least BASE / 2 keeps each estimated quotient limb off by at most 2
        uint32_t norm = BASE / (y.back() + 1ULL);
        Limbs u = mul_small(x, norm), v = mul_small(y, norm);
        v.pop_back();

        size_t n = v.size(), m = u.size() - n;
        Limbs q(m);
        for (size_t j = m; j-- > 0;){
            unsigned long long num = (unsigned long long)u[j + n] * BASE + u[j + n - 1];
            unsigned long long qhat = num / v[n - 1], rhat = num % v[n - 1];
            while (qhat >= BASE || qhat * v[n - 2] > rhat * BASE + u[j + n - 2]){
                qhat--, rhat += v[n - 1];
                if (rhat >= BASE) break;
            }

            long long borrow = 0;
            unsigned long long carry = 0;
            for (size_t i = 0; i < n; i++){
                unsigned long long p = qhat * v[i] + carry;
                carry = p / BASE;
                long long d = (long long)u[i + j] - (long long)(p % BASE) - borrow;
                borrow = d < 0;
                u[i + j] = borrow ? d + BASE : d;
            }
            long long top = (long long)u[j + n] - (long long)carry - borrow;

            if (top < 0){
                /// qhat was one too large, add the divisor back
                qhat--;
                uint32_t c = 0;
                for (size_t i = 0; i < n; i++){
                    uint32_t s = u[i + j] + v[i] + c;
                    c = s >= BASE;
                    u[i + j] = c ? s - BASE : s;
                }
                top += c;
            }
            u[j + n] = top;
            q[j] = qhat;
        }

        u.resize(n);
        trim(u), trim(q);
        divmod_small(u, norm);
        return {q, u};
    }
};

int main(){
    Bignum x = 1;
    for (int i = 0; i < 100; i++) x *= 2;
    assert(x.to_string() == "1267650600228229401496703205376");
    assert((x + x).to_string() == "2535301200456458802993406410752");
    assert((x % 1000000007).to_string() == "976371285");
    assert((x / 1000000007).to_string() == "1267650591354675262013");
    assert(x / 1000000007 * 1000000007 + x % 1000000007 == x);
    assert((x - x - x) == -x && (x - x).to_string() == "0");

    assert(Bignum("-0") == 0 && Bignum("+000123") == 123 && Bignum("-000").to_string() == "0");
    assert(Bignum(LLONG_MIN).to_string() == "-9223372036854775808");
    assert(Bignum(LLONG_MIN) - 1 < LLONG_MIN && -Bignum(LLONG_MIN) > LLONG_MAX);

    for (int a = -100; a <= 100; a++){
        for (int b = -100; b <= 100; b++){
            Bignum p = a, q = b;
            assert(p + q == a + b && p - q == a - b && p * q == a * b);
            if (b) assert(p / q == a / b && p % q == a % b);
            assert((p < q) == (a < b) && (p > q) == (a > b) && (p == q) == (a == b));
            assert((p <= q) == (a <= b) && (p >= q) == (a >= b) && (p != q) == (a != b));
        }
    }

    Bignum fact = 1;
    for (int i = 1; i <= 1000; i++) fact *= i;
    assert(fact.to_string().size() == 2568);
    for (int i = 1000; i >= 1; i--) fact /= i;
    assert(fact == 1);

    auto nines = [](int k){ return Bignum(string(k, '9')); };
    auto pow10 = [](int k){ return Bignum("1" + string(k, '0')); };
    assert((nines(20000) * nines(20000)).to_string() == string(19999, '9') + "8" + string(19999, '0') + "1");
    assert(nines(27000) / nines(18000) == pow10(9000) && nines(27000) % nines(18000) == nines(9000));
    assert((-pow10(27000) / nines(9000)).to_string() == "-1" + string(8999, '0') + "1" + string(8999, '0') + "1");
    assert(-pow10(27000) % nines(9000) == -1);

    return 0;
}
