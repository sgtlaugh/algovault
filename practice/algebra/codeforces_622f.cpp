/***
 * Lagrange polynomial interpolation
 * A polynomial of degree d can be uniquely identified given its values on d + 1 unique points
 *
 * Pre-calculates given the first n points (x=0 to n-1) of the polynomial
 * Then answers each query to interpolate the x'th term
 * All values are done modulo mod, which needs to be a prime as we need its inverse modulo
 *
 * Also includes an additional helper function called find_degree(terms, mod)
 * Given at least the first d+2 points of a polynomial of degree d, it finds d
 *
 * Usage:
 *   Lagrange lagrange(terms, mod);                  // terms[x] = f(x) for x = 0 .. n - 1, any int values (reduced mod mod)
 *   int y = lagrange.interpolate(x);                // f(x) mod mod in [0, mod) for any long long x, negative too
 *   vector<int> w = get_terms(terms, mod, l, r);    // f(l), f(l + 1), ..., f(r) mod mod
 *   int d = find_degree(terms, mod);                // degree d, or -1 if fewer than d + 2 terms are given
 *
 * Complexity:
 *   - O(n + log mod) to pre-calculate, O(n) memory
 *   - O(n) time and memory per interpolate
 *   - O((k n + log mod) log n) for find_degree, k = smallest power with mod^k >= INT_MAX, 2 for mod near 10^9
 *
 * Note, n should not exceed mod due to the way modular inverse is used
 * In such cases, we can use interpolation without modulo in big integers and take the remainder later
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Lagrange{
    vector<int> terms, dp;
    int mod, n;

    Lagrange(const vector<int>& terms, int mod) : terms(terms), mod(mod){
        n = terms.size();
        assert(n <= mod);

        int i, v, f;
        for (auto& t : this->terms){
            if ((t %= mod) < 0) t += mod;
        }

        for (f = 1, i = 1; i < n; i++) f = (long long)f * i % mod;
        v = expo(f, mod - 2);

        vector<int> inv(n, v);
        for (i = n - 1; i > 0; i--){
            inv[i - 1] = (long long)inv[i] * i % mod;
        }

        dp.resize(n, 1);
        for (i = 0; i < n; i++){
            dp[i] = (long long)inv[i] * inv[n - i - 1] % mod;
            dp[i] = (long long)dp[i] * this->terms[i] % mod;
        }
    }

    int expo(int a, int b){
        int res = 1;

        while (b){
            if (b & 1) res = (long long)res * a % mod;
            a = (long long)a * a % mod;
            b >>= 1;
        }

        return res;
    }

    int interpolate(long long x){
        if (0 <= x && x < n) return terms[x];
        x %= mod;

        int i, w;
        vector<int> X(n, 1), Y(n, 1);
        for (i = 1; i < n; i++){
            X[i] = (long long)X[i - 1] * (x - i + 1) % mod;
            if (X[i] < 0) X[i] += mod;
        }

        for (i = n - 2; i >= 0; i--){
            Y[i] = (long long)Y[i + 1] * (x - i - 1) % mod;
            if (Y[i] < 0) Y[i] += mod;
        }

        long long res = 0;
        for (i = 0; i < n; i++){
            w = ((long long)X[i] * Y[i] % mod) * dp[i] % mod;
            if ((n - i + 1) & 1) w = mod - w;
            res += w;
        }

        return res % mod;
    }
};

vector<int> get_terms(const vector<int>& terms, int mod, int l, int r){
    auto lagrange = Lagrange(terms, mod);

    vector<int> res;
    for (int i = l; i <= r; i++){
        res.push_back(lagrange.interpolate(i));
    }

    return res;
}

int find_degree(const vector<int>& terms, int mod){
    long long v = mod;
    int k = 1, n = terms.size();

    while (v < INT_MAX){
        v *= mod;
        k++;
    }

    /// Top residues, a point congruent to a prefix's nodes 0 .. mid - 1 would match it by construction
    int l = mod - k, r = mod - 1;
    auto expected = get_terms(terms, mod, l, r);

    int low = 1, high = n - 1;
    while ((low + 1) < high){
        int mid = (low + high) >> 1;
        vector<int> v(terms.begin(), terms.begin() + mid);
        if (get_terms(v, mod, l, r) == expected) high = mid;
        else low = mid;
    }

    for (int d = low; d <= high; d++){
        vector<int> v(terms.begin(), terms.begin() + d);
        if (get_terms(v, mod, l, r) == expected) return d - 1;
    }

    return -1;
}

int power(int a, int b, int mod){
    int res = 1;
    for (; b; b >>= 1, a = (long long)a * a % mod){
        if (b & 1) res = (long long)res * a % mod;
    }
    return res;
}

int main(){
    const int mod = 1000000007;
    int n, k;

    scanf("%d %d", &n, &k);

    vector<int> terms(k + 2, 0);
    for (int x = 1; x <= k + 1; x++) terms[x] = (terms[x - 1] + power(x, k, mod)) % mod;

    printf("%d\n", Lagrange(terms, mod).interpolate(n));
    return 0;
}
