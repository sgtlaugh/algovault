/***
 *
 * Generates two different strings with the same hash
 * Hashing is standard polynomial hashing with base and modulo
 * Base and modulo can be any 64 bit positive integer, provided base + modulo does not exceed signed 64 bit int limit

 * Extremely useful for hacking solutions in Codeforces :-D
 * Or if you're the setter and you're in an ominous mood, construct test cases against commonly used base/mod pairs
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct AntiHash{
    /// Magic constants! Feel free to tweak and experiment

    const static int n = 2000;
    const static int lim = 6;
    const static int chunk = 15;

    long long base, mod;

    int len, found, sample[1 << chunk];
    string res[2] = {string(n, '0'), string(n, '0')};
    pair<long long, int> ar[1 << chunk], dp[lim][n + 10];

    AntiHash(long long base, long long mod) : base(base), mod(mod) {
        generate();
    }

    inline long long fast_modmul(long long a, long long b, long long m){
        if (a >= m) a %= m;
        if (b >= m) b %= m;
        if (m < (long long)UINT_MAX) return (uint64_t)a * b % m;

        #ifdef __SIZEOF_INT128__
            return __int128(a) * b % m;
        #endif

        long double x = (long double)a * b;
        long long c = x / m;

        a = (uint64_t)a * b - (uint64_t)c * m;
        if (a >= m) a -= m;
        if (a < 0) a += m;

        return a;
    }

    inline long long get_hash(const string& s){
        long long h = 0;
        for (int i = 0; i < (int)s.length(); i++){
            h = fast_modmul(h, base, mod);
            if (s[i] == '0') h = (h + 17) % mod;
            if (s[i] == '1') h = (h + 1000000009) % mod;
        }
        return h;
    }

    inline void generate(){
        len = 0, found = 0;
        dp[0][len++] = make_pair(1, n - 1);
        for (int i = 1; i < n; i++){
            dp[0][len++] = make_pair(fast_modmul(dp[0][i - 1].first, base, mod), n - i - 1);
        }
        sort(dp[0], dp[0] + len);

        for (int i = 1; i < lim; i++){
            for (int j = 0, k = 0; k < len; j++, k += 2){
                long long x = (dp[i - 1][k + 1].first - dp[i - 1][k].first + mod) % mod;
                dp[i][j] = make_pair(x, j);
            }
            len >>= 1;
            sort(dp[i], dp[i] + len);
        }
    }

    inline void dfs(unsigned int l, int i, int pos){
        if (l == 0) res[pos][dp[0][i].second] = '1';
        else{
            dfs(l - 1, dp[l][i].second * 2 + 1, pos);
            dfs(l - 1, dp[l][i].second * 2, pos ^ 1);
        }
    }

    inline void insert_hash(int i, long long sum, long long mask){
        if (i == chunk) ar[mask] = make_pair(sum, mask);
        else{
            insert_hash(i + 1, sum, mask);
            insert_hash(i + 1, (sum + dp[lim - 1][sample[i]].first) % mod, mask | (1LL << i));
        }
    }

    inline void find_hash(int i, long long sum, long long mask){
        if (found) return;
        if (i == chunk){
            int idx = lower_bound(ar + 1, ar + (1 << chunk), pair<long long, int>(sum, 0)) - ar;

            if (idx < (1 << chunk) && ar[idx].first == sum){
                mask = (mask << (long long)chunk) | ar[idx].second;
                for (i = 0; i < 2 * chunk; i++){
                    if (mask & (1LL << i)){
                        dfs(lim - 1, sample[i], (i < chunk));
                    }
                }

                if (get_hash(res[0]) == get_hash(res[1]) && res[0] != res[1]) found = 1;
                else {
                    for (i = 0; i < n; i++) res[0][i] = res[1][i] = '0';
                }
                return;
            }
        }
        else{
            find_hash(i + 1, sum, mask);
            find_hash(i + 1, (sum + dp[lim - 1][sample[i + chunk]].first) % mod, mask | (1LL << i));
        }
    }

    pair<string, string> solve(){
        while (!found){
            vector <int> v;
            for (int i = 0; i < len; i++){
                v.push_back(i < (2 * chunk));
            }
            random_shuffle(v.begin(), v.end());
            for (int i = 0, k = 0; i < (int)v.size(); i++){
                if (v[i] == 1) sample[k++] = i;
            }
            random_shuffle(sample, sample + chunk * 2);

            insert_hash(0, 0, 0);
            sort(ar, ar + (1 << chunk));
            find_hash(0, 0, 0);
        }

        return {res[0], res[1]};
    }
};

int main(){
    srand(666);
    clock_t start = clock();

    long long base = 1000000007;
    long long mod = (1ULL << 61) - 1; /// A Mersenne prime

    string s1, s2;
    auto ah = AntiHash(base, mod);
    tie(s1, s2) = ah.solve();

    long long h1 = 0, h2 = 0;
    for (int i = 0; i < (int)s1.length(); i++){
        h1 = (ah.fast_modmul(h1, base, mod) + s1[i]) % mod;
    }
    for (int i = 0; i < (int)s2.length(); i++){
        h2 = (ah.fast_modmul(h2, base, mod) + s2[i]) % mod;
    }

    assert(s1.length() == s2.length() && s1 != s2 && h1 == h2);

    /// freopen("out.txt", "w", stdout);
    puts(s1.c_str());
    puts(s2.c_str());

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));   /// 0.002585
    return 0;
}
