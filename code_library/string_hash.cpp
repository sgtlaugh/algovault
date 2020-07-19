#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

const int PRIMES[] = {2146936003, 2146946617, 2146958423, 2146960229, 2146960807, 2146970347, 2146975387, 2146977451, 2146983689, 2146990819, 2146995293, 2147005727, 2147007413, 2147009827, 2147011109, 2147015821, 2147023189, 2147023981, 2147024147, 2147030671, 2147033533, 2147037509, 2147042083, 2147048129, 2147050289, 2147051521, 2147055829, 2147060693, 2147062429, 2147063411, 2147063497, 2147066567, 2147070421, 2147073067, 2147076623, 2147084963, 2147088469, 2147089541, 2147092657, 2147095961, 2147098991, 2147102401, 2147103983, 2147115797, 2147116787, 2147120981, 2147124677, 2147127991, 2147139413, 2147139721, 2147144119, 2147147797, 2147156203, 2147172509, 2147180779, 2147181097, 2147182469, 2147186579, 2147188543, 2147189749, 2147194691, 2147194759, 2147204693, 2147207407, 2147211307, 2147212129, 2147212853, 2147215639, 2147223263, 2147234723, 2147235547, 2147236447, 2147237669, 2147238547, 2147240659, 2147247163, 2147261713, 2147261723, 2147273507, 2147276561, 2147286569, 2147286721, 2147291453, 2147301719, 2147303281, 2147308561, 2147309809, 2147317297, 2147317427, 2147325461, 2147327191, 2147332529, 2147333003, 2147336467, 2147340733, 2147345821, 2147348179, 2147348257, 2147353051, 2147353097};

/// Get the hash of any substring in O(1) after O(n) pre-processing

struct SimpleHash{
    long long base, mod;
    vector <int> basepow, f_hash, r_hash;

    SimpleHash() {}
    SimpleHash(const char* str, long long base, long long mod): base(base), mod(mod){
        int len = strlen(str);
        basepow.resize(len + 3, 1), f_hash.resize(len + 3, 0), r_hash.resize(len + 3, 0);

        for (int i = 1; i <= len; i++){
            basepow[i] = basepow[i - 1] * base % mod;
            f_hash[i] = (f_hash[i - 1] * base + str[i - 1] + 1007) % mod;
        }

        for (int i = len; i >= 1; i--){
            r_hash[i] = (r_hash[i + 1] * base + str[i - 1] + 1007) % mod;
        }
    }

    /// returns the hash of the substring from indices l to r (inclusive)
    inline int forward_hash(int l, int r){
        int h = f_hash[r + 1] - ((long long)basepow[r - l + 1] * f_hash[l] % mod);
        return h < 0 ? h + mod : h;
    }

    /// returns the hash of the reversed substring from indices l to r (inclusive)
    inline int reverse_hash(int l, int r){;
        int h = r_hash[l + 1] - ((long long)basepow[r - l + 1] * r_hash[r + 2] % mod);
        return h < 0 ? h + mod : h;
    }
};

/// Get the hash of any substring in O(1) after O(n) pre-processing with double hashing

struct StringHash{
    SimpleHash sh1, sh2;

    StringHash(const char* str){
        unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
        srand(seed);  /// never get hacked in Codeforces :D, comment this line for easier debugging

        int q_len = (sizeof(PRIMES) / sizeof(PRIMES[0])) / 4;
        int base_1 = PRIMES[rand() % q_len];
        int mod_1  = PRIMES[rand() % q_len + q_len];
        int base_2 = PRIMES[rand() % q_len + 2 * q_len];
        int mod_2  = PRIMES[rand() % q_len + 3 * q_len];

        sh1 = SimpleHash(str, base_1, mod_1);
        sh2 = SimpleHash(str, base_2, mod_2);
    }

    /// returns the hash of the substring from indices l to r (inclusive)
    long long forward_hash(int l, int r){
        return ((long long)sh1.forward_hash(l, r) << 32) ^ sh2.forward_hash(l, r);
    }

    /// returns the hash of the reversed substring from indices l to r (inclusive)
    long long reverse_hash(int l, int r){
        return ((long long)sh1.reverse_hash(l, r) << 32) ^ sh2.reverse_hash(l, r);
    }
};

int main(){
    auto hasher = StringHash("racecar");
    assert(hasher.forward_hash(0, 6) == hasher.reverse_hash(0, 6));
    assert(hasher.forward_hash(1, 5) != hasher.reverse_hash(0, 4));
    assert(hasher.forward_hash(1, 1) == hasher.reverse_hash(5, 5));
    assert(hasher.forward_hash(1, 1) != hasher.reverse_hash(5, 6));
    assert(hasher.forward_hash(2, 4) == hasher.reverse_hash(2, 4));
    return 0;
}
