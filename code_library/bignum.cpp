#include <stdio.h>
#include <bits/stdtr1c++.h>

#define SIZE 100

using namespace std;

const int base_len = 8;
const long long base = 0x5F5E100;

struct UnsignedBignum{
    int len;
    long long data[SIZE];

    long long &operator [](int x){
        return(data[x]);
    }
    const long long &operator [](int x) const{
        return(data[x]);
    }

    UnsignedBignum(){
        len = 0;
        memset(data, 0, sizeof(data));
    }

    void clear(){
        for (int i = 1; i <= len; i++) data[i] = 0;
        len = 0;
    }

    inline int compare(const UnsignedBignum &a, const UnsignedBignum &b){
        if (a.len > b.len) return 0;
        if (a.len < b.len) return 1;
        for (int i = a.len; i >= 1; i--){
            if (a.data[i] > b.data[i]) return 0;
            if (a.data[i] < b.data[i]) return 1;
        }
        return 2;
    }

    inline bool operator < (const UnsignedBignum &x){
        return (compare(*this, x) == 1);
    }
    inline bool operator > (const UnsignedBignum &x){
        return (compare(*this, x) == 0);
    }
    inline bool operator <= (const UnsignedBignum &x){
        return (compare(*this, x) >= 1);
    }
    inline bool operator >= (const UnsignedBignum &x){
        return ((compare(*this, x) & 1) == 0);
    }
    inline bool operator != (const UnsignedBignum &x){
        return (compare(*this, x) != 2);
    }
    inline bool operator == (const UnsignedBignum &x){
        return (compare(*this, x) == 2);
    }
    inline UnsignedBignum operator = (const UnsignedBignum& x){
        for (int i = x.len + 1; i <= len; i++) data[i] = 0;
        for (int i = 1; i <= x.len; i++) data[i] = x.data[i];
        len = x.len;
        return *this;
    }

    inline UnsignedBignum operator = (long long x){
        for (int i = 0; i <= len; i++) data[i] = 0;
        len = 0;

        do{
            data[++len] = x % base;
            x /= base;
        }
        while (x);
        return *this;
    }

    inline UnsignedBignum operator += (const UnsignedBignum &b){
        return *this = (*this + b);
    }
    inline UnsignedBignum operator *= (const UnsignedBignum &b){
        return *this = (*this * b);
    }
    inline UnsignedBignum operator -= (const UnsignedBignum &b){
        return *this = (*this - b);
    }
    inline UnsignedBignum operator /= (const UnsignedBignum &b){
        return *this = (*this / b);
    }
    inline UnsignedBignum operator %= (const UnsignedBignum &b){
        return *this = (*this % b);
    }
    inline UnsignedBignum operator *= (long long x){
        return (*this = (*this * x));
    }
    inline UnsignedBignum operator += (long long x){
        return (*this = (*this + x));
    }
    inline UnsignedBignum operator -= (long long x){
        return (*this = (*this - x));
    }
    inline UnsignedBignum operator /= (long long x){
        return (*this = (*this / x));
    }

    UnsignedBignum(long long x){
        len = 0;
        memset(data, 0, sizeof(data));
        (*this) = x;
    }

    UnsignedBignum operator * (const UnsignedBignum& x){
        int i, j;
        UnsignedBignum res = UnsignedBignum(0);
        for (i = 1; i <= len; i++){
            if (data[i]){
                for (j = 1; j <= x.len; j++){
                    if (x.data[j]){
                        res.data[i + j - 1] += data[i] * x.data[j];
                        res.data[i + j] += res.data[i + j - 1] / base;
                        res.data[i + j - 1] %= base;
                    }
                }
            }
        }

        res.len = len + x.len - 1;
        while (res.data[res.len + 1]) res.len++;
        while (res[res.len] == 0 && res.len > 1) res.len--;
        return res;
    }

    UnsignedBignum operator / (long long x){
        assert(x != 0);

        UnsignedBignum res;
        long long y = 0;
        for (int i = len; i >= 1; i--){
            y = y * base + data[i];
            res[i] = y / x;
            y %= x;
        }

        res.len = len;
        while (res[res.len] == 0 && res.len > 1) res.len--;
        return res;
    }

    void divide(const UnsignedBignum& b, UnsignedBignum& rem, UnsignedBignum& quot){
        assert(!(b.len == 1 && b[1] == 0));

        long long x = data[len], y = b[b.len];
        int i, l1 = base_len * (len - 1), l2 = base_len * (b.len - 1);

        while (x) x /= 10, l1++;
        while (y) y /= 10, l2++;

        UnsignedBignum temp = b;
        rem = *this;

        for (i = 1; i * base_len <= (l1 - l2); i++) temp *= base;
        for (i = 1; i <= (l1 - l2) % base_len; i++) temp *= 10;
        for (i = l1 - l2; i >= 0; i--){
            x = 0;
            while (rem >= temp) rem -= temp, x++;
            quot[i / base_len + 1] = quot[i / base_len + 1] * 10 + x;
            temp /= 10;
        }

        quot.len = (l1 - l2) / base_len + 1;
        while(rem.len > 1 && !rem[rem.len]) rem.len--;
        while(quot.len > 1 && !quot[quot.len]) quot.len--;
    }

    UnsignedBignum operator / (const UnsignedBignum& x){
        UnsignedBignum rem, quot;
        divide(x, rem, quot);
        assert(quot.len > 0);
        return quot;
    }

    UnsignedBignum operator % (const UnsignedBignum& x){
        UnsignedBignum rem, quot;
        divide(x, rem, quot);
        assert(rem.len > 0);
        return rem;
    }

    long long operator % (long long x){
        long long res = 0;
        for (int i = len; i >= 1; i--) res = (res * base + data[i]) % x;
        return res;
    }

    UnsignedBignum operator + (const UnsignedBignum& x){
        UnsignedBignum res;
        int i, l = max(len, x.len);
        for (i = 1; i <= l; i++) res[i] = data[i] + x[i];
        for (i = 1; i <= l; i++) res[i + 1] += res[i] / base, res[i] %= base;

        res.len = l;
        if (res[res.len + 1]) res.len++;
        while (res[res.len] == 0 && res.len > 1) res.len--;
        return res;
    }

    UnsignedBignum operator - (const UnsignedBignum& x){
        UnsignedBignum res;
        for (int i = 1; i <= len; i++) res.data[i] = data[i] - x.data[i];
        for (int i = 1; i <= len; i++){
            if (res[i] < 0) res.data[i] += base, res.data[i + 1]--;
        }

        res.len = len;
        while (res[res.len] == 0 && res.len > 1) res.len--;
        return res;
    }

    UnsignedBignum(const char* str, int dummy){
        int l = strlen(str + 1);
        long long val = 0;
        (*this).clear();

        for (int i = 1; i <= (l - 1) / base_len + 1; i++){
            val = 0;
            for (int j = l - base_len * i + 1; j <= l - base_len * i + base_len; j++){
                if (j >= 1) val = val * 10 + str[j] - 48;
            }
            data[++len] = val;
        }
    }

    void nextUnsignedBignum(){
        char str[SIZE * base_len + 10];
        scanf("%s", str + 1);
        UnsignedBignum(str, 666);
    }

    void print(char endChar = '\n'){
        printf("%lld", data[len]);
        for (int i = len - 1; i >= 1; i--) printf("%0*lld", base_len, data[i]);
        putchar(endChar);
    }
};

struct Bignum{
    int sign;
    UnsignedBignum val;

    Bignum(){}
    Bignum(long long x){
        sign = 0;
        val = llabs(x);
        if (x < 0) sign = 1;
    }

    Bignum(UnsignedBignum x){
        val = x;
        sign = 0;
    }

    inline Bignum operator = (long long x){
        sign = 0;
        val = llabs(x);
        if (x < 0) sign = 1;
        return *this;
    }

    Bignum operator * (const Bignum& x){
        UnsignedBignum y = x.val;
        Bignum res;
        res.val = y * this->val;
        res.sign = this->sign ^ x.sign;
        if (res.val == UnsignedBignum(0)) res.sign = 0;
        return res;
    }

    Bignum operator / (const Bignum& x){
        UnsignedBignum y = this->val;
        Bignum res;
        res.val = y / x.val;
        res.sign = this->sign ^ x.sign;
        if (res.val == UnsignedBignum(0)) res.sign = 0;
        return res;
    }

    Bignum operator % (const Bignum& x){
        UnsignedBignum y = this->val;
        Bignum res;
        res.val = y % x.val;
        res.sign = this->sign;
        if (res.val == UnsignedBignum(0)) res.sign = 0;
        return res;
    }

    Bignum operator + (const Bignum& x){
        Bignum res;

        if (this->sign == x.sign){
            res.sign = this->sign;
            res.val = this->val + x.val;
        }
        else{
            UnsignedBignum v = this->val, w = x.val;

            int cmp = v.compare(v, w);
            if (cmp != 1){
                res.sign = this->sign;
                res.val = v - w;
            }
            else{
                res.sign = x.sign;
                res.val = w - v;
            }
        }

        if (res.val == UnsignedBignum(0)) res.sign = 0;
        return res;
    }

    Bignum operator - (const Bignum& x){
        Bignum res = x;
        res.sign ^= 1;
        return *this + res;
    }

    inline int compare(const Bignum &a, const Bignum &b){
        if (a.sign != b.sign){
            if (a.sign == 0) return 0;
            if (b.sign == 0) return 1;
        }

        UnsignedBignum x = a.val, y = b.val;
        if (a.sign == 0) return x.compare(x, y);
        return x.compare(y, x);
    }

    inline bool operator < (const Bignum &x){
        return (compare(*this, x) == 1);
    }
    inline bool operator > (const Bignum &x){
        return (compare(*this, x) == 0);
    }
    inline bool operator <= (const Bignum &x){
        return (compare(*this, x) >= 1);
    }
    inline bool operator >= (const Bignum &x){
        return ((compare(*this, x) & 1) == 0);
    }
    inline bool operator != (const Bignum &x){
        return (compare(*this, x) != 2);
    }
    inline bool operator == (const Bignum &x){
        return (compare(*this, x) == 2);
    }

    inline Bignum operator += (const Bignum &b){
        return *this = (*this + b);
    }
    inline Bignum operator *= (const Bignum &b){
        return *this = (*this * b);
    }
    inline Bignum operator -= (const Bignum &b){
        return *this = (*this - b);
    }
    inline Bignum operator /= (const Bignum &b){
        return *this = (*this / b);
    }
    inline Bignum operator %= (const Bignum &b){
        return *this = (*this % b);
    }
    inline Bignum operator *= (long long x){
        return (*this = (*this * x));
    }
    inline Bignum operator += (long long x){
        return (*this = (*this + x));
    }
    inline Bignum operator -= (long long x){
        return (*this = (*this - x));
    }
    inline Bignum operator /= (long long x){
        return (*this = (*this / x));
    }

    inline void nextBignum(){
        char str[SIZE * base_len + 10];
        scanf("%s", str);
        sign = (str[0] == '-');
        if (!sign){
            int i, l = strlen(str);
            for (i = l; i >= 1; i--) str[i] = str[i - 1];
            str[l + 1] = 0;
        }
        val = UnsignedBignum(str, 666);
    }

    void print(char endChar = 10){
        if (sign) putchar('-');
        val.print(endChar);
    }
};


int main(){
    Bignum x = Bignum(1);
    for (int i = 0; i < 100; i++) x = x * Bignum(2);
    x.print(); /// 1267650600228229401496703205376 (2^100)

    Bignum y = x + x;
    y.print(); /// 2535301200456458802993406410752

    Bignum r = x % 1000000007;
    r.print(); /// 976371285

    Bignum q = x / 1000000007;
    q.print(); /// 1267650591354675262013

    Bignum v = q * 1000000007 + r;
    v.print(); /// 1267650600228229401496703205376

    Bignum c = v - y;
    c.print(); /// -1267650600228229401496703205376

    Bignum chksum = x + c;
    chksum.print(); /// 0
    return 0;
}
