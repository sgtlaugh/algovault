/***
 *
 * Square Root Decomposition
 * Range add and "how many values in [l, r] are below x", which a plain segment tree cannot answer
 *
 * Complexity: O(sqrt(n) log n) per operation, O(n log n) to build
 *
 * SqrtDecomposition<T> s(values), 0-based positions, ranges [l, r] inclusive
 * s.add(l, r, v): adds v to every value in [l, r]
 * s.count_less(l, r, x): how many values in [l, r] are < x
 * s.get(i): current value at position i
 * Every current value and every x must stay within half the range of T (|v| <= 1e9 for int, ~4.6e18 for long long),
 * because a block's pending addition is the difference of two such values
 *
 * Each block keeps a sorted copy of its raw values plus one pending addition:
 * whole blocks only change the pending addition, partial blocks are fixed up and re-sorted
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct SqrtDecomposition{
    int n, block;
    vector<T> raw, pending;
    vector<vector<T>> sorted_block;

    SqrtDecomposition(const vector<T>& values) : n(values.size()), raw(values){
        block = max(1, (int)sqrt((double)n));
        int blocks = (n + block - 1) / block;
        pending.assign(blocks, 0);
        sorted_block.resize(blocks);
        for (int b = 0; b < blocks; b++) rebuild(b);
    }

    void rebuild(int b){
        int lo = b * block, hi = min(n, lo + block);
        sorted_block[b].assign(raw.begin() + lo, raw.begin() + hi);
        sort(sorted_block[b].begin(), sorted_block[b].end());
    }

    /// Adds v to [l, r] inside block b, folding its pending addition in first so raw values equal real values
    void add_partial(int b, int l, int r, T v){
        int lo = b * block, hi = min(n, lo + block);
        for (int i = lo; i < hi && pending[b] != 0; i++) raw[i] += pending[b];
        pending[b] = 0;
        for (int i = l; i <= r; i++) raw[i] += v;
        rebuild(b);
    }

    void add(int l, int r, T v){
        assert(0 <= l && l <= r && r < n);
        int bl = l / block, br = r / block;
        if (bl == br) return add_partial(bl, l, r, v);
        add_partial(bl, l, (bl + 1) * block - 1, v);
        for (int b = bl + 1; b < br; b++) pending[b] += v;
        add_partial(br, br * block, r, v);
    }

    int count_less(int l, int r, T x) const{
        assert(0 <= l && l <= r && r < n);
        int bl = l / block, br = r / block, res = 0;
        if (bl == br){
            for (int i = l; i <= r; i++) res += raw[i] + pending[bl] < x;
            return res;
        }
        for (int i = l; i < (bl + 1) * block; i++) res += raw[i] + pending[bl] < x;
        for (int b = bl + 1; b < br; b++){
            T shift = pending[b];  /// compare raw + shift with x, x - shift could overflow
            res += lower_bound(sorted_block[b].begin(), sorted_block[b].end(), x, [&](const T& v, const T& key){ return v + shift < key; }) - sorted_block[b].begin();
        }
        for (int i = br * block; i <= r; i++) res += raw[i] + pending[br] < x;
        return res;
    }

    T get(int i) const{
        assert(0 <= i && i < n);
        return raw[i] + pending[i / block];
    }
};

int main(){
    SqrtDecomposition<long long> s({5, 1, 4, 1, 3, 9, 2, 6, 5});
    assert(s.count_less(0, 8, 4) == 4);
    assert(s.count_less(2, 6, 3) == 2);
    s.add(1, 7, 10);
    assert(s.get(0) == 5 && s.get(1) == 11 && s.get(7) == 16 && s.get(8) == 5);
    assert(s.count_less(0, 8, 10) == 2);
    assert(s.count_less(1, 7, 13) == 3);
    s.add(4, 4, -100);
    assert(s.count_less(0, 8, 0) == 1 && s.get(4) == -87);

    SqrtDecomposition<int> one({7});
    one.add(0, 0, 3);
    assert(one.get(0) == 10 && one.count_less(0, 0, 11) == 1 && one.count_less(0, 0, 10) == 0);
    return 0;
}
