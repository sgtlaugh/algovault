/***
 *
 * Big Integer
 * Arbitrary precision signed integers that behave like C++ integers
 *
 * Complexity, for numbers of n and m digits:
 *   - +, -, comparisons, to_string and stream I/O: O(n + m)
 *   - *: O(n m) below about 860 digits, O(n^1.58) above, O((n + m) log(n + m)) once both have about 6300 digits
 *   - / and %: O(n) for a divisor below 10^9, O((n - m) m) unless divisor and quotient both pass 3600 digits, O(n log n) then
 *
 * Bignum x(v): the long long v, 0 by default
 * Bignum x(s): the decimal string s, optional sign and leading zeros allowed, asserts at least one digit and nothing else
 * x + y: the sum, either side may be a long long
 * x - y: the difference
 * x * y: the product
 * x / y: the quotient truncated toward zero, asserts y != 0
 * x % y: the remainder with the sign of x, asserts y != 0
 * -x: the negation
 * x += y: x + y stored in x, likewise -=, *=, /= and %=
 * x < y: the comparison, likewise <=, >, >=, == and !=
 * x.to_string(): the decimal digits, with a leading - when x is negative
 * cin >> x: reads one whitespace separated decimal token
 * cout << x: writes x.to_string()
 *
 * In loops prefer the compound operators: x += y reuses x's storage, and *=, /= and %= by a long long below 10^9 never allocate
 * x / y and x % y each run the full division, so needing both costs two
 * At -O2, multiplying two 2 * 10^6 digit numbers takes ~0.15 s and dividing 2 * 10^6 digits by 10^6 digits ~0.5 s
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
    typedef vector<uint32_t> Limbs;

    bool neg = false;      /// false for zero, so 0 has one representation
    Limbs limbs;           /// 9 decimal digits each, value = sum limbs[i] * 10^(9i), no zero limbs at the top

    Bignum(long long v = 0){
        neg = v < 0;
        for (unsigned long long x = neg ? 0ULL - v : v; x; x /= BASE) limbs.push_back(x % BASE);
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
            limbs.push_back(limb);
        }

        trim();
    }

    Bignum(Limbs&& x, bool negative) : neg(negative), limbs(move(x)){
        trim();
    }

    friend Bignum operator + (const Bignum& x, const Bignum& y){
        if (x.neg == y.neg) return Bignum(add(x.limbs, y.limbs), x.neg);
        if (cmp(x.limbs, y.limbs) >= 0) return Bignum(sub(x.limbs, y.limbs), x.neg);
        return Bignum(sub(y.limbs, x.limbs), y.neg);
    }

    friend Bignum operator - (const Bignum& x, const Bignum& y){
        return x + (-y);
    }

    friend Bignum operator * (const Bignum& x, const Bignum& y){
        return Bignum(mul(x.limbs, y.limbs), x.neg != y.neg);
    }

    friend Bignum operator / (const Bignum& x, const Bignum& y){
        return Bignum(divmod(x.limbs, y.limbs).first, x.neg != y.neg);
    }

    friend Bignum operator % (const Bignum& x, const Bignum& y){
        return Bignum(divmod(x.limbs, y.limbs).second, x.neg);
    }

    Bignum operator - () const{
        Bignum res = *this;
        if (!res.limbs.empty()) res.neg = !res.neg;
        return res;
    }

    Bignum& operator += (const Bignum& y){ return add_signed(y.limbs, y.neg); }
    Bignum& operator -= (const Bignum& y){ return add_signed(y.limbs, !y.neg); }
    Bignum& operator *= (const Bignum& y){ return *this = *this * y; }
    Bignum& operator /= (const Bignum& y){ return *this = *this / y; }
    Bignum& operator %= (const Bignum& y){ return *this = *this % y; }

    /// A factor or divisor below BASE works on these limbs without allocating, larger ones take the general path
    Bignum& operator *= (long long y){
        if (y <= -(long long)BASE || y >= (long long)BASE) return *this *= Bignum(y);
        if (y < 0) neg = !neg, y = -y;
        mul_small(limbs, y);
        trim();
        return *this;
    }

    Bignum& operator /= (long long y){
        assert(y != 0);
        if (y <= -(long long)BASE || y >= (long long)BASE) return *this /= Bignum(y);
        if (y < 0) neg = !neg, y = -y;
        divmod_small(limbs, y);
        trim();
        return *this;
    }

    Bignum& operator %= (long long y){
        assert(y != 0);
        if (y <= -(long long)BASE || y >= (long long)BASE) return *this %= Bignum(y);
        uint32_t r = divmod_small(limbs, y < 0 ? -y : y);
        limbs.assign(r ? 1 : 0, r);
        trim();
        return *this;
    }

    friend bool operator == (const Bignum& x, const Bignum& y){ return x.neg == y.neg && x.limbs == y.limbs; }
    friend bool operator != (const Bignum& x, const Bignum& y){ return !(x == y); }
    friend bool operator < (const Bignum& x, const Bignum& y){
        if (x.neg != y.neg) return x.neg;
        return x.neg ? cmp(y.limbs, x.limbs) < 0 : cmp(x.limbs, y.limbs) < 0;
    }
    friend bool operator > (const Bignum& x, const Bignum& y){ return y < x; }
    friend bool operator <= (const Bignum& x, const Bignum& y){ return !(y < x); }
    friend bool operator >= (const Bignum& x, const Bignum& y){ return !(x < y); }

    string to_string() const{
        if (limbs.empty()) return "0";
        string s = (neg ? "-" : "") + std::to_string(limbs.back());
        size_t pos = s.size();
        s.resize(pos + 9 * (limbs.size() - 1));

        for (int i = (int)limbs.size() - 2; i >= 0; i--){
            uint32_t v = limbs[i];
            for (int k = 8; k >= 0; k--){
                s[pos + k] = '0' + v % 10;
                v /= 10;
            }
            pos += 9;
        }

        return s;
    }

    friend istream& operator >> (istream& is, Bignum& x){
        string s;
        if (is >> s) x = Bignum(s);
        return is;
    }

    friend ostream& operator << (ostream& os, const Bignum& x){
        return os << x.to_string();
    }

    Bignum& add_signed(const Limbs& y, bool y_neg){
        if (neg == y_neg){
            limbs.resize(max(limbs.size(), y.size()) + 1);
            add_to(limbs, y, 0);
        }
        else if (cmp(limbs, y) >= 0) sub_from(limbs, y);
        else{
            Limbs d = y;
            sub_from(d, limbs);
            limbs = move(d), neg = y_neg;
        }

        trim();
        return *this;
    }

    void trim(){
        trim(limbs);
        if (limbs.empty()) neg = false;
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
        size_t i = offset;
        for (uint32_t v : y){
            uint32_t s = res[i] + v + carry;
            carry = s >= BASE;
            res[i++] = carry ? s - BASE : s;
        }

        for (; carry; i++){
            carry = res[i] == BASE - 1;
            res[i] = carry ? 0 : res[i] + 1;
        }
    }

    /// Requires x >= y
    static void sub_from(Limbs& x, const Limbs& y){
        int borrow = 0;
        size_t i = 0;
        for (uint32_t v : y){
            long long d = (long long)x[i] - v - borrow;
            borrow = d < 0;
            x[i++] = borrow ? d + BASE : d;
        }

        for (; borrow; i++){
            borrow = x[i] == 0;
            x[i] = borrow ? BASE - 1 : x[i] - 1;
        }
    }

    static Limbs mul(const Limbs& x, const Limbs& y){
        if (x.size() < y.size()) return mul(y, x);
        size_t n = x.size(), m = y.size();
        if (m <= 1){
            if (!m) return {};
            Limbs res = x;
            mul_small(res, y[0]);
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

    /// Always appends the carry limb, even 0, since Knuth D's normalized dividend needs that extra top limb
    static void mul_small(Limbs& x, uint32_t f){
        unsigned long long carry = 0;
        for (auto& limb : x){
            unsigned long long cur = (unsigned long long)limb * f + carry;
            limb = cur % BASE, carry = cur / BASE;
        }

        x.push_back(carry);
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

    template<uint32_t MOD>
    static uint32_t pow_mod(unsigned long long b, unsigned long long e){
        unsigned long long res = 1;
        for (b %= MOD; e; e >>= 1, b = b * b % MOD){
            if (e & 1) res = res * b % MOD;
        }

        return res;
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

    static uint32_t divmod_small(Limbs& x, uint32_t d){
        /// One reciprocal turns each 64-bit division into a multiply: inv underestimates 2^64 / d, so q is low by at most 1
        unsigned long long inv = ~0ULL / d, rem = 0;
        for (int i = (int)x.size() - 1; i >= 0; i--){
            unsigned long long cur = rem * BASE + x[i], q = (unsigned __int128)cur * inv >> 64;
            rem = cur - q * d;
            if (rem >= d) q++, rem -= d;
            x[i] = q;
        }

        trim(x);
        return rem;
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

    static Limbs shift_right(const Limbs& x, size_t k){
        return k >= x.size() ? Limbs{} : Limbs(x.begin() + k, x.end());
    }

    static Limbs shift_left(const Limbs& x, size_t k){
        Limbs res(k + x.size());
        copy(x.begin(), x.end(), res.begin() + k);
        return res;
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
        Limbs u = x, v = y;
        mul_small(u, norm), mul_small(v, norm);
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

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
 *
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ... and returns how many it read before the input ended
 * read accepts integers of any width including __int128, float, double, long double, strings, chars and vectors of these
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... separated by single spaces and ends the line, a vector argument on its own line
 * write accepts integers, strings, chars, vectors of these and fio::fixed values
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 * fio::flush(): writes the pending output now, which also happens at normal exit
 *
 * Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 * A char reads one non-space character
 * A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 * Writing a float without fixed() does not compile, so its precision is always chosen
 * Floats print at printf's speed, since the formatting itself dominates; a floating token keeps its first 127 characters
 * Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: stdin and stdout are one stream each per process, so a second reader instance
/// would split the input between two buffers and a second writer would reorder the output
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
            if (buf_len == 0) return EOF;
        }

        return (unsigned char)inbuf[inptr++];
    }

    /// Integers: takes every type the plain overloads below miss, so __int128 needs no is_integral specialization, which breaks C++20
    template <typename T>
    bool read_one(T &x){
        int c = ' ', neg = 0;
        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

        return true;
    }

    bool read_token(){
        int c = ' ', len = 0;
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        if (!read_token()) return false;
        x = strtof(tmpbuf, nullptr);
        return true;
    }

    bool read_one(double &x){
        if (!read_token()) return false;
        x = strtod(tmpbuf, nullptr);
        return true;
    }

    bool read_one(long double &x){
        if (!read_token()) return false;
        x = strtold(tmpbuf, nullptr);
        return true;
    }

    bool read_one(string &s){
        int c = ' ';
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !isspace(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }

        return true;
    }

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T> &v){
        for (auto &&x: v){
            if (!read_one(x)) return false;
        }

        return true;
    }

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
        return true;
    }

    template <typename T>
    int read(T &x){
        return read_one(x);
    }

    template <typename T, typename ...Args>
    int read(T &x, Args &...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(char c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    /// Integers, like read_one: the plain overloads below win for every other type
    template <typename T>
    void write_one(T x){
        int l = 0;
        if (x < 0){
            /// Peel one digit while still negative so that -x cannot overflow for the minimum value
            write_char('-');
            tmpbuf[l++] = '0' - x % 10;
            x = -(x / 10);
        }

        while (x || !l){
            tmpbuf[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmpbuf[--l]);
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char *s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

    void write_one(const Fixed &f){
        int len = snprintf(tmpbuf, sizeof(tmpbuf), "%.*Lf", f.digits, f.value);
        if (len < (int)sizeof(tmpbuf)){
            for (int j = 0; j < len; j++) write_char(tmpbuf[j]);
            return;
        }

        string s(len + 1, '\0');
        snprintf(&s[0], len + 1, "%.*Lf", f.digits, f.value);
        s.resize(len);
        write_one(s);
    }

    template <typename T>
    void write_one(const vector<T> &v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    template <typename T>
    char separator_after(const T &){
        return ' ';
    }

    template <typename T>
    char separator_after(const vector<T> &){
        return '\n';
    }

    template <typename T>
    void write(const T &x){
        write_one(x);
        write_char('\n');
    }

    template <typename T, typename ...Args>
    void write(const T &x, const Args &...args){
        write_one(x);
        write_char(separator_after(x));
        write(args...);
    }
}

using namespace fio;

int cnt[367][10];

int main(){
    int n, d;
    char buf[64];
    string line;
    Bignum fact = 1;

    for (n = 1; n <= 366; n++){
        fact *= n;
        for (char c: fact.to_string()) cnt[n][c - '0']++;
    }

    while (read(n) && n){
        snprintf(buf, sizeof(buf), "%d! --", n);
        write(string(buf));

        for (d = 0; d < 10; d++){
            snprintf(buf, sizeof(buf), "%s(%d)%5d", d % 5 ? "    " : "   ", d, cnt[n][d]);
            line += buf;
            if (d % 5 == 4) write(line), line.clear();
        }
    }

    return 0;
}
