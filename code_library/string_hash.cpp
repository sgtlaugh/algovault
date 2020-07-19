/***
 * 
 * String hashing with double hash
 * Get the forward and reverse hash of any substring
 * Complexity - O(n) to build and O(1) for substring hash query
 * 
***/

#include <bits/stdc++.h>

using namespace std;

const int PRIMES[] = {2147462393, 2147462419, 2147462587, 2147462633, 2147462747, 2147463167, 2147463203, 2147463569, 2147463727, 2147463863, 2147464211, 2147464549, 2147464751, 2147465153, 2147465563, 2147465599, 2147465743, 2147465953, 2147466457, 2147466463, 2147466521, 2147466721, 2147467009, 2147467057, 2147467067, 2147467261, 2147467379, 2147467463, 2147467669, 2147467747, 2147468003, 2147468317, 2147468591, 2147468651, 2147468779, 2147468801, 2147469017, 2147469041, 2147469173, 2147469229, 2147469593, 2147469881, 2147469983, 2147470027, 2147470081, 2147470177, 2147470673, 2147470823, 2147471057, 2147471327, 2147471581, 2147472137, 2147472161, 2147472689, 2147472697, 2147472863, 2147473151, 2147473369, 2147473733, 2147473891, 2147473963, 2147474279, 2147474921, 2147474929, 2147475107, 2147475221, 2147475347, 2147475397, 2147475971, 2147476739, 2147476769, 2147476789, 2147476927, 2147477063, 2147477107, 2147477249, 2147477807, 2147477933, 2147478017, 2147478521, 2147478563, 2147478649, 2147479447, 2147479589, 2147480707, 2147480837, 2147480927, 2147480971, 2147481263, 2147481311, 2147481337, 2147481367, 2147481997, 2147482021, 2147482063, 2147482081, 2147482343, 2147482591, 2147483069, 2147483123};

struct SingleHash{
    long long base, mod;
    vector <int> basepow, f_hash, r_hash;

    SingleHash() {}
    SingleHash(const char* str, long long base, long long mod): base(base), mod(mod){
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

    inline int forward_hash(int l, int r){
        int h = f_hash[r + 1] - ((long long)basepow[r - l + 1] * f_hash[l] % mod);
        return h < 0 ? h + mod : h;
    }

    inline int reverse_hash(int l, int r){;
        int h = r_hash[l + 1] - ((long long)basepow[r - l + 1] * r_hash[r + 2] % mod);
        return h < 0 ? h + mod : h;
    }
};

struct StringHash{
    SingleHash sh1, sh2;

    StringHash(const char* str){
        unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
        srand(seed);  /// to avoid getting hacked in CF, comment this line for easier debugging

        int q_len = (sizeof(PRIMES) / sizeof(PRIMES[0])) / 4;
        int base_1 = PRIMES[rand() % q_len];
        int mod_1  = PRIMES[rand() % q_len + q_len];
        int base_2 = PRIMES[rand() % q_len + 2 * q_len];
        int mod_2  = PRIMES[rand() % q_len + 3 * q_len];

        sh1 = SingleHash(str, base_1, mod_1);
        sh2 = SingleHash(str, base_2, mod_2);
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
