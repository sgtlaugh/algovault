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
 *     Multiplication: O(n * m) for small inputs, Karatsuba O(n^1.58) for large ones
 *     Division: O((n - m) * m), O(n) when the divisor fits in 9 digits
 *
 * Multiplying two 10^5 digit numbers takes ~0.03 s, two 10^6 digit numbers ~1 s
 * Dividing 2 * 10^5 digits by 10^5 digits takes ~0.3 s, quadratic beyond that (Newton iteration with FFT would be needed)
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Bignum{
    static const uint32_t BASE = 1000000000;
    static const int KARATSUBA_CUTOFF = 96;

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

    /// Knuth's algorithm D, returns {|x| / |y|, |x| % |y|}
    static pair<Limbs, Limbs> divmod(const Limbs& x, const Limbs& y){
        assert(!y.empty());
        if (cmp(x, y) < 0) return {{}, x};
        if (y.size() == 1){
            Limbs q = x;
            uint32_t r = divmod_small(q, y[0]);
            return {q, r ? Limbs{r} : Limbs{}};
        }

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

    return 0;
}
