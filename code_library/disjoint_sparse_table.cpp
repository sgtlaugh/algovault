/***
 *
 * Disjoint Sparse Table For RMQ on a static array
 * Idempotency is necessary for the regular sparse table
 * Because we answer queries in O(1) allowing partial overlaps which doesn't change the result
 * Some examples are finding minimum or maximum
 *
 * Disjoint Sparse Table can answer RMQ on a static array in O(1), even when they are not idempotent
 * The operation needs to be associative
 * Examples are sum, xor, product, matrix multiplication, etc
 * For sum or xor, we can just keep a cumulative array but for many operations like the product, disjoint sparse tables can be useful
 * Further reading - https://codeforces.com/blog/entry/87940
 *
 * Time and space complexity is the same as in regular sparse table
 * O(n log n) to build
 * O(1) to query
 *
 * Note, we can still do this in O(log n) using regular sparse tables
 * And then there are sqrt trees...
 * Disjoint sparse tables should be rarely useful in CP, nonetheless they are quite cool
 *
 * To use, update the combine function and use an appropriate identity element
 * An identity element is any element which leaves any element unchanged when combined with it
 * Or in other words, i is an identity element if combine(i, x) = x for all x
 * It is 0 for xor, 1 for product, -infinity for maximum and so on
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename T>
struct DisjointST{
    vector <T> dp[32];

    /// defined for sum by default, change as required
    T combine(const T& x, const T& y){
        return x + y;
    }

    DisjointST(const vector<T> &ar, const T& identity_elem){
        int i, c, h, l, n = (int)ar.size();

        for (h = 0; h <= __lg(n); h++){
            dp[h].resize(n + 1, identity_elem);

            for (c = l = 1 << h; c < n + l; c += (l << 1)){
                for (i = c + 1; i <= min(n, c + l); i++){
                    dp[h][i] = combine(dp[h][i - 1], ar[i - 1]);
                }

                for (i = min(n, c) - 1; i >= c - l; i--){
                    dp[h][i] = combine(ar[i], dp[h][i + 1]);
                }
            }
        }
    }

    T query(int l, int r){
        int h = __lg(l ^ (r + 1));
        return combine(dp[h][l], dp[h][r + 1]);
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto rmq = DisjointST<int>(v, 0);

    assert(rmq.query(0, 0) == 5);
    assert(rmq.query(0, 2) == 12);
    assert(rmq.query(1, 4) == 27);
    assert(rmq.query(3, 5) == 24);
    assert(rmq.query(4, 7) == 86);
    assert(rmq.query(2, 6) == 34);
    assert(rmq.query(1, 8) == 129);

    v.clear();
    mt19937 rng(0);

    const int n = 2000000;
    for (int i = 0; i < n; i++){
        v.push_back(rng() % 1000000000);
    }

    clock_t start = clock();
    rmq = DisjointST<int>(v, 0);
    fprintf(stderr, "\nTime taken to build = %0.6f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Took 0.105 s locally

    return 0;
}
