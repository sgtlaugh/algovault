/***
 *
 * Square Root Decomposition
 * Range add, and the count of values below x in a range
 *
 * Complexity: O(sqrt(n) log n) per operation, O(n log n) to build
 *
 * SqrtDecomposition<T> s(values): positions 0..n - 1
 * s.add(l, r, v): adds v to every value in [l, r], 0 <= l <= r < n
 * s.count_less(l, r, x): how many values in [l, r] are < x, 0 <= l <= r < n
 * s.get(i): the current value at position i, 0 <= i < n
 *
 * Every current value and every x must stay within half the range of T, |v| <= 1e9 for int and about 4.6e18 for long long
 * The bound holds because a block's pending addition is the difference of two such values
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct SqrtDecomposition{
    int n, block_size;
    vector<T> raw, pending;            /// a value is raw[i] + pending[i / block_size], whole-block adds only touch pending
    vector<vector<T>> sorted_block;    /// raw values of each block, sorted

    SqrtDecomposition(const vector<T>& values) : n(values.size()), raw(values){
        block_size = max(1, (int)sqrt((double)n));
        int blocks = (n + block_size - 1) / block_size;
        pending.assign(blocks, 0);
        sorted_block.resize(blocks);
        for (int b = 0; b < blocks; b++) rebuild(b);
    }

    void add(int l, int r, T v){
        assert(0 <= l && l <= r && r < n);
        int bl = l / block_size, br = r / block_size;
        if (bl == br) return add_partial(bl, l, r, v);
        add_partial(bl, l, (bl + 1) * block_size - 1, v);
        for (int b = bl + 1; b < br; b++) pending[b] += v;
        add_partial(br, br * block_size, r, v);
    }

    int count_less(int l, int r, T x) const{
        assert(0 <= l && l <= r && r < n);
        int bl = l / block_size, br = r / block_size, res = 0;
        if (bl == br){
            for (int i = l; i <= r; i++) res += raw[i] + pending[bl] < x;
            return res;
        }

        for (int i = l; i < (bl + 1) * block_size; i++) res += raw[i] + pending[bl] < x;
        for (int b = bl + 1; b < br; b++){
            T shift = pending[b];  /// compare raw + shift with x, x - shift could overflow
            res += lower_bound(sorted_block[b].begin(), sorted_block[b].end(), x, [&](const T& v, const T& key){ return v + shift < key; }) - sorted_block[b].begin();
        }
        for (int i = br * block_size; i <= r; i++) res += raw[i] + pending[br] < x;
        return res;
    }

    T get(int i) const{
        assert(0 <= i && i < n);
        return raw[i] + pending[i / block_size];
    }

    void rebuild(int b){
        int lo = b * block_size, hi = min(n, lo + block_size);
        sorted_block[b].assign(raw.begin() + lo, raw.begin() + hi);
        sort(sorted_block[b].begin(), sorted_block[b].end());
    }

    /// Adds v to [l, r] inside block b, folding its pending addition in first so raw values equal real values
    void add_partial(int b, int l, int r, T v){
        int lo = b * block_size, hi = min(n, lo + block_size);
        for (int i = lo; i < hi && pending[b] != 0; i++) raw[i] += pending[b];
        pending[b] = 0;
        for (int i = l; i <= r; i++) raw[i] += v;
        rebuild(b);
    }
};

int main(){
    /// Counts below x, before and after range adds
    SqrtDecomposition<long long> s({5, 1, 4, 1, 3, 9, 2, 6, 5});
    assert(s.count_less(0, 8, 4) == 4);   /// 1, 1, 3 and 2
    assert(s.count_less(2, 6, 3) == 2);
    s.add(1, 7, 10);                      /// 5 11 14 11 13 19 12 16 5
    assert(s.get(1) == 11);
    assert(s.get(8) == 5);
    assert(s.count_less(0, 8, 10) == 2);
    assert(s.count_less(1, 7, 13) == 3);
    s.add(4, 4, -100);                    /// 5 11 14 11 -87 19 12 16 5
    assert(s.get(4) == -87);
    assert(s.count_less(0, 8, 0) == 1);

    return 0;
}
