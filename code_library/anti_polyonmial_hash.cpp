#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/// Magic constants, feel free to tweak and experiment
const int n = 2000;
const int lim = 6;
const int chunk = 15;

long long b, m;
int len, found, sample[1 << chunk];
string res[2] = {string(n, '0'), string(n, '0')};
pair<long long, int> ar[1 << chunk], dp[lim][n + 10];

inline long long multiply(long long a, long long b, long long m){ /// a * b % m
    a %= m, b %= m;
    long long res = 0;
    while (b){
        if (b & 1) res = (res + a) % m;
        a = (a << 1) % m, b >>= 1;
    }
    return res;
}

inline long long get_hash(string s, long long b, long long m){
    long long h = 0;
    for (int i = 0; i < (int)s.length(); i++){
        h = multiply(h, b, m);
        if (s[i] == '0') h = (h + 17) % m;
        if (s[i] == '1') h = (h + 1000000009) % m;
    }
    return h;
}

inline void generate(){
    int i, j, k;
    len = 0, found = 0;

    dp[0][len++] = make_pair(1, n - 1);
    for (i = 1; i < n; i++){
        dp[0][len++] = make_pair(multiply(dp[0][i - 1].first, b, m), n - i - 1);
    }
    sort(dp[0], dp[0] + len);

    for (i = 1; i < lim; i++){
        for (j = 0, k = 0; k < len; j++, k += 2){
            long long x = (dp[i - 1][k + 1].first - dp[i - 1][k].first + m) % m;
            dp[i][j] = make_pair(x, j);
        }
        len >>= 1;
        sort(dp[i], dp[i] + len);
    }
}

inline void dfs(int l, int i, int pos){
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
        insert_hash(i + 1, (sum + dp[lim - 1][sample[i]].first) % m, mask | (1LL << i));
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

            long long hash1 = get_hash(res[0], b, m);
            long long hash2 = get_hash(res[1], b, m);
            if (hash1 == hash2 && res[0] != res[1]){
                found = 1;
                puts(res[0].c_str());
                puts(res[1].c_str());

                fprintf(stderr, "\n%lld\n", hash1);
                fprintf(stderr, "%lld\n", hash2);
            }
            else{
                for (i = 0; i < n; i++) res[0][i] = res[1][i] = '0';
            }
			return;
		}
    }
    else{
        find_hash(i + 1, sum, mask);
        find_hash(i + 1, (sum + dp[lim - 1][sample[i + chunk]].first) % m, mask | (1LL << i));
    }
}

void solve(){
	/***
	 * Prints two different strings with the same hash
	 * Hashing is standard polynomial hashing with base and modulo
	 * Base and modulo can be any 64 bit positive integer, provided base + modulo does not exceed signed 64 bit int limit

	 * Extremely useful for hacking solutions in Codeforces :-D
	 * Or if you're the setter and you're in an ominous mood, construct test cases against commonly used base/mod pairs
	***/

    generate();
    res[0] = string(n, '0'), res[1] = string(n, '0');

    while (found ^ 1){
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
}

int main(){
    srand(666);
    /// freopen("out.txt", "w", stdout);

    b = 1000000007;
    m = 1000000000000000003LL; /// 10^18 + 3

    solve();
    return 0;
}
