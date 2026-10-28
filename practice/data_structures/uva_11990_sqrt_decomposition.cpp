/***
 *
 * Square Root Decomposition
 * Range add, and the count of values below x in a range
 *
 * Complexity:
 *   - build: O(n log n)
 *   - add: O(sqrt(n)) amortized, the first add covering part of a block sorts that block once
 *   - count_less: O(sqrt(n) log n)
 *   - get: O(1)
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
    vector<vector<int>> sorted_index;  /// sorted_index[b][j] is the position holding sorted_block[b][j], empty until block b's first partial add
    vector<T> moved_value;             /// add_partial's scratch, kept to skip an allocation per add
    vector<int> moved_index;

    /// Sorting values alone builds about 25% faster than sorting positions by value, so the positions wait for the first partial add
    SqrtDecomposition(const vector<T>& values) : n(values.size()), raw(values){
        block_size = max(1, (int)sqrt((double)n));
        int blocks = (n + block_size - 1) / block_size;
        pending.assign(blocks, 0);
        sorted_block.resize(blocks);
        sorted_index.resize(blocks);
        for (int b = 0; b < blocks; b++){
            sorted_block[b].assign(raw.begin() + b * block_size, raw.begin() + min(n, (b + 1) * block_size));
            sort(sorted_block[b].begin(), sorted_block[b].end());
        }
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

    /// Adds v to [l, r] inside block b in O(block_size): the touched and untouched entries each stay sorted, so one merge restores the order
    void add_partial(int b, int l, int r, T v){
        int lo = b * block_size, hi = min(n, lo + block_size);
        vector<T>& value = sorted_block[b];
        vector<int>& index = sorted_index[b];
        if (pending[b] != 0){
            for (int i = lo; i < hi; i++) raw[i] += pending[b];
            for (auto& x : value) x += pending[b];
            pending[b] = 0;
        }
        for (int i = l; i <= r; i++) raw[i] += v;
        if (index.empty()) return rebuild(b);

        int m = value.size(), kept = 0;
        moved_value.clear(), moved_index.clear();
        for (int j = 0; j < m; j++){
            if (index[j] < l || index[j] > r) value[kept] = value[j], index[kept] = index[j], kept++;
            else moved_value.push_back(value[j] + v), moved_index.push_back(index[j]);
        }

        int i = kept - 1, k = (int)moved_value.size() - 1;
        for (int w = m - 1; k >= 0; w--){
            if (i >= 0 && value[i] > moved_value[k]) value[w] = value[i], index[w] = index[i], i--;
            else value[w] = moved_value[k], index[w] = moved_index[k], k--;
        }
    }

    void rebuild(int b){
        int lo = b * block_size, hi = min(n, lo + block_size);
        vector<int>& index = sorted_index[b];
        index.resize(hi - lo);
        iota(index.begin(), index.end(), lo);
        sort(index.begin(), index.end(), [&](int i, int j){ return raw[i] < raw[j]; });
        for (int j = 0; j < hi - lo; j++) sorted_block[b][j] = raw[index[j]];
    }
};

struct Fenwick{
    vector<int> tree;

    Fenwick(int n) : tree(n + 1, 0) {}

    void add(int i, int v){
        for (i++; i < (int)tree.size(); i += i & -i) tree[i] += v;
    }

    int prefix(int i) const{
        int res = 0;
        for (; i > 0; i -= i & -i) res += tree[i];
        return res;
    }
};

int main(){
    int n, m;

    while (scanf("%d %d", &n, &m) == 2){
        vector<int> a(n), pos(n + 1);
        for (int i = 0; i < n; i++) scanf("%d", &a[i]), pos[a[i]] = i;

        long long res = 0;
        Fenwick alive_pos(n), alive_val(n + 1);
        for (int i = 0; i < n; i++){
            res += i - alive_val.prefix(a[i]);
            alive_val.add(a[i], 1), alive_pos.add(i, 1);
        }

        SqrtDecomposition<int> s(a);
        while (m--){
            printf("%lld\n", res);
            int v, p;
            scanf("%d", &v);
            p = pos[v];
            int less_before = p > 0 ? s.count_less(0, p - 1, v) : 0;
            int greater_before = alive_pos.prefix(p) - less_before;
            int less_after = alive_val.prefix(v) - less_before;
            res -= greater_before + less_after;
            s.add(p, p, n + 1 - v);
            alive_pos.add(p, -1), alive_val.add(v, -1);
        }
    }

    return 0;
}
