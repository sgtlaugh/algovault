/***
 *
 * Anti Polynomial Hash
 * Builds two different strings with the same polynomial hash, to hack solutions or break weak tests
 *
 * Complexity, with chunk = 15 and lim = 1 + log2(n / 60) rounded down, at least 1:
 *   - AntiHash: O(n log n)
 *   - solve: O(2^chunk chunk + n + chunk 2^lim) per randomized trial, trials repeat until a collision, with no bound
 *   - 512 KB in a fixed array inside the struct, plus 16 n lim bytes, about 600 KB at n = 1000
 *
 * AntiHash hack(base, mod, n): the attack on h = (h * base + s[i]) % mod, n >= 30, base + mod < 2^63, both positive
 * hack.solve(): {s, t}, two different n-character strings of '0' and '1' with equal hashes
 *
 * Notes:
 *   - solve averaged under 10 ms at n = 1000 for mod <= 1e14 and at n = 2000 for mod < 2^61
 *   - Too small an n for mod slows solve down fast: n = 1000 with mod near 2^61 took 2 s on average and up to 9 s
 *   - The strings differ by a +-1 combination of base powers, so they collide whatever values '0' and '1' hash to
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct AntiHash{
    /// Magic constant! Each trial meets in the middle over 2^chunk subset sums

    const static int chunk = 15;

    int n, lim, len, found, sample[2 * chunk];
    long long base, mod;
    string res[2];
    pair<long long, int> ar[1 << chunk];
    vector<vector<pair<long long, int>>> dp;

    AntiHash(long long base, long long mod, int n) : n(n), base(base), mod(mod){
        assert(n >= 2 * chunk);
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
        res[0] = res[1] = string(n, '0');

        /// Halving stops before the last level drops under 4 chunk values, as the tuned n = 2000, lim = 6 left it
        for (lim = 1, len = n; len / 2 >= 4 * chunk; lim++) len /= 2;
        dp.assign(lim, vector<pair<long long, int>>(n));

        len = 0, found = 0;
        dp[0][len++] = make_pair(1, n - 1);
        for (int i = 1; i < n; i++){
            dp[0][len++] = make_pair(fast_modmul(dp[0][i - 1].first, base, mod), n - i - 1);
        }
        sort(dp[0].begin(), dp[0].begin() + len);

        for (int i = 1; i < lim; i++){
            for (int j = 0, k = 0; k + 1 < len; j++, k += 2){
                long long x = (dp[i - 1][k + 1].first - dp[i - 1][k].first + mod) % mod;
                dp[i][j] = make_pair(x, j);
            }
            len >>= 1;
            sort(dp[i].begin(), dp[i].begin() + len);
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
            vector<int> v;
            for (int i = 0; i < len; i++){
                v.push_back(i < (2 * chunk));
            }
            shuffle(v.begin(), v.end(), mt19937(rand()));
            for (int i = 0, k = 0; i < (int)v.size(); i++){
                if (v[i] == 1) sample[k++] = i;
            }
            shuffle(sample, sample + chunk * 2, mt19937(rand()));

            insert_hash(0, 0, 0);
            sort(ar, ar + (1 << chunk));
            find_hash(0, 0, 0);
        }

        return {res[0], res[1]};
    }
};

int main(){
    /// Two 1000-character strings colliding under base 3 and mod 1e14, scored with '0' -> 1 and '1' -> 2
    long long base = 3, mod = 100000000000000LL;
    AntiHash hack(base, mod, 1000);

    string s, t;
    tie(s, t) = hack.solve();
    assert(s != t);
    assert(s.size() == 1000);
    assert(t.size() == 1000);

    long long hs = 0, ht = 0;
    for (int i = 0; i < 1000; i++){
        hs = (hack.fast_modmul(hs, base, mod) + (s[i] - '0' + 1)) % mod;
        ht = (hack.fast_modmul(ht, base, mod) + (t[i] - '0' + 1)) % mod;
    }
    assert(hs == ht);

    return 0;
}
