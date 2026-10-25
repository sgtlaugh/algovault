/***
 *
 * Run-Length Implicit Treap
 * A sequence of up to 4.6e18 values kept as runs of equal values, edited and queried by position
 *
 * Complexity: expected O(log r) per operation, r = number of runs, each operation adds at most 3 runs
 *
 * RunLengthImplicitTreap<Runs = SumRuns> seq(n, x): n copies of x, 0 <= n <= 4.6e18
 * RunLengthImplicitTreap<Runs> seq(runs): built from {value, count} pairs in O(r), every count >= 1
 * seq.assign(l, r, x): sets positions l..r to x
 * seq.erase(l, r): removes positions l..r
 * seq.insert(pos, x, count): count copies of x before position pos, pos = seq.size() appends
 * seq.reverse(l, r): reverses the order of positions l..r
 * seq.update(l, r, lz): applies the lazy update lz to positions l..r
 * seq.for_each_run(l, r, f): f(x, count) for each run inside l..r, left to right
 * seq.get(pos): the value at position pos
 * seq.max_right(l, ok): last r with ok(seq.query(l, r)) true, l - 1 if none, 0 <= l <= size()
 * seq.min_left(r, ok): first l with ok(seq.query(l, r)) true, r + 1 if none, -1 <= r < size()
 * seq.query(l, r): the Agg of positions l..r
 * seq.size(): the length of the sequence
 *
 * Notes:
 *   - Positions are 0-based and ranges need 0 <= l <= r < size()
 *   - ok must stay false once false as the range grows, it is never called on an empty range
 *   - For fixed positions known in advance, a segment tree over compressed coordinates is simpler
 *   - Runs says what a run holds, SumRuns below (long long values, sums, add as the lazy) is the model to copy:
 *       - Value, Agg: what a run stores and what query returns
 *       - run_agg(x, count): the Agg of count copies of x, x * count for sums
 *       - combine(a, b): the Agg of two neighbouring pieces, a on the left
 *       - Lazy, apply(x, lz), apply(a, lz, count), compose(old, new): only for seq.update, compose is old then new
 *       - apply(x, lz) updates one value, apply(a, lz, count) an Agg covering count positions
 *       - reversed(x, count), reversed(a): only for seq.reverse, the run and the Agg read backwards
 *       - cut(x, k): the first k values and the rest of a run, {x, x} when every run is all equal
 *       - Calling update without Lazy or reverse without reversed fails to compile
 *   - Every Agg, the whole sequence's included, must fit its type: SumRuns needs every sum within long long
 *   - 0/1 flags with flip (Codeforces 817F): Value int, Agg {ones, count}, Lazy int, apply flips, compose(a, b) = a ^ b
 *   - Arithmetic progressions: Value {a, d}, cut(x, k) = {x, {a + k * d, d}}, reversed(x, c) = {a + (c - 1) * d, -d}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SumRuns{
    using Value = long long;
    using Agg = long long;
    using Lazy = long long;

    static Agg run_agg(Value x, long long count){ return x * count; }
    static Agg combine(Agg a, Agg b){ return a + b; }
    static Value apply(Value x, Lazy lz){ return x + lz; }
    static Agg apply(Agg a, Lazy lz, long long count){ return a + lz * count; }
    static Lazy compose(Lazy old_lz, Lazy new_lz){ return old_lz + new_lz; }
    static Value reversed(Value x, long long){ return x; }
    static Agg reversed(Agg a){ return a; }
    static pair<Value, Value> cut(Value x, long long){ return {x, x}; }
};

/// Optional Runs parts are detected, not defaulted, so a missing Lazy or reversed fails update or reverse to compile
template<typename R, typename = void>
struct LazyOf{
    using type = char;
    static constexpr bool exists = false;
};

template<typename R>
struct LazyOf<R, void_t<typename R::Lazy>>{
    using type = typename R::Lazy;
    static constexpr bool exists = true;
};

template<typename R, typename = void>
struct HasReversed : false_type {};

template<typename R>
struct HasReversed<R, void_t<decltype(R::reversed(declval<typename R::Agg>()))>> : true_type {};

template<typename Runs = SumRuns>
struct RunLengthImplicitTreap{
    using Value = typename Runs::Value;
    using Agg = typename Runs::Agg;
    using Lazy = typename LazyOf<Runs>::type;

    struct Node{
        Value x;
        Agg agg;
        Lazy lz;
        long long count, size;
        unsigned int priority;
        int l, r;
        bool has_lz, rev;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    RunLengthImplicitTreap(long long n, const Value& x){
        assert(0 <= n && n <= LLONG_MAX / 2);
        if (n > 0) root = new_node(x, n);
    }

    RunLengthImplicitTreap(const vector<pair<Value, long long>>& runs){
        /// Cartesian tree on random priorities with a stack, O(n)
        vector<int> stack;
        for (auto& [x, count] : runs){
            assert(count >= 1);
            int t = new_node(x, count), last = -1;
            while (!stack.empty() && nodes[stack.back()].priority < nodes[t].priority){
                last = stack.back();
                stack.pop_back();
                pull(last);
            }
            nodes[t].l = last;
            if (!stack.empty()) nodes[stack.back()].r = t;
            stack.push_back(t);
        }

        if (!stack.empty()) root = stack[0];
        while (!stack.empty()){
            pull(stack.back());
            stack.pop_back();
        }
        assert(size() <= LLONG_MAX / 2);
    }

    void assign(long long l, long long r, const Value& x){
        on_range(l, r, [&](int& t){ t = new_node(x, r - l + 1); });
    }

    void erase(long long l, long long r){
        on_range(l, r, [&](int& t){ t = -1; });
    }

    void insert(long long pos, const Value& x, long long count){
        assert(0 <= pos && pos <= size() && count >= 1 && count <= LLONG_MAX / 2 - size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, new_node(x, count)), b);
    }

    void reverse(long long l, long long r){
        static_assert(HasReversed<Runs>::value, "define reversed(value, count) and reversed(agg) to use reverse");
        on_range(l, r, [&](int& t){ apply_reverse(t); });
    }

    void update(long long l, long long r, const Lazy& lz){
        static_assert(LazyOf<Runs>::exists, "define Lazy, apply(value, lz), apply(agg, lz, count) and compose to use update");
        on_range(l, r, [&](int& t){ apply_lazy(t, lz); });
    }

    template<typename F>
    void for_each_run(long long l, long long r, F f){
        on_range(l, r, [&](int& t){
            vector<int> stack;
            for (int u = t; u != -1 || !stack.empty(); ){
                for (; u != -1; u = nodes[u].l) push(u), stack.push_back(u);
                u = stack.back();
                stack.pop_back();
                f(nodes[u].x, nodes[u].count);
                u = nodes[u].r;
            }
        });
    }

    Value get(long long pos){
        Value res{};
        on_range(pos, pos, [&](int& t){ res = nodes[t].x; });
        return res;
    }

    template<typename Ok>
    long long max_right(long long l, Ok ok){
        assert(0 <= l && l <= size());
        int a, b;
        split(root, l, a, b);
        long long taken = descend(b, ok, false);
        root = merge(a, b);
        return l + taken - 1;
    }

    template<typename Ok>
    long long min_left(long long r, Ok ok){
        assert(-1 <= r && r < size());
        int a, b;
        split(root, r + 1, a, b);
        long long taken = descend(a, ok, true);
        root = merge(a, b);
        return r - taken + 1;
    }

    Agg query(long long l, long long r){
        Agg res{};
        on_range(l, r, [&](int& t){ res = nodes[t].agg; });
        return res;
    }

    long long size() const{
        return size(root);
    }

    long long size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int new_node(const Value& x, long long count){
        nodes.push_back({x, Runs::run_agg(x, count), Lazy(), count, count, (unsigned int)rng(), -1, -1, false, false});
        return nodes.size() - 1;
    }

    void pull(int t){
        Node& u = nodes[t];
        u.size = u.count + size(u.l) + size(u.r);
        Agg a = Runs::run_agg(u.x, u.count);
        if (u.l != -1) a = Runs::combine(nodes[u.l].agg, a);
        if (u.r != -1) a = Runs::combine(a, nodes[u.r].agg);
        u.agg = a;
    }

    /// f may replace the treap of [l, r], with -1 or a new node
    template<typename F>
    void on_range(long long l, long long r, F f){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        f(b);
        root = merge(merge(a, b), c);
    }

    /// First k positions into a, the rest into b; new_node may reallocate nodes, so children go through locals
    void split(int t, long long k, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        push(t);
        long long left = size(nodes[t].l), count = nodes[t].count;
        if (k <= left){
            int p, q;
            split(nodes[t].l, k, p, q);
            nodes[t].l = q;
            a = p, b = t;
        }
        else if (k >= left + count){
            int p, q;
            split(nodes[t].r, k - left - count, p, q);
            nodes[t].r = p;
            a = t, b = q;
        }
        else{
            /// The new node draws a fresh priority: pieces of one run sharing a priority would chain into a path
            auto [first, rest] = Runs::cut(nodes[t].x, k - left);
            int u = new_node(rest, count - (k - left)), right = nodes[t].r;
            nodes[t].r = -1;
            nodes[t].x = first, nodes[t].count = k - left;
            a = t, b = merge(u, right);
        }
        pull(t);
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            push(a);
            int r = merge(nodes[a].r, b);
            nodes[a].r = r;
            pull(a);
            return a;
        }
        push(b);
        int l = merge(a, nodes[b].l);
        nodes[b].l = l;
        pull(b);
        return b;
    }

    void apply_reverse(int t){
        if constexpr (HasReversed<Runs>::value){
            if (t == -1) return;
            Node& u = nodes[t];
            swap(u.l, u.r);
            u.x = Runs::reversed(u.x, u.count);
            u.agg = Runs::reversed(u.agg);
            u.rev ^= 1;
        }
    }

    void apply_lazy(int t, const Lazy& lz){
        if constexpr (LazyOf<Runs>::exists){
            if (t == -1) return;
            Node& u = nodes[t];
            u.x = Runs::apply(u.x, lz);
            u.agg = Runs::apply(u.agg, lz, u.size);
            u.lz = u.has_lz ? Runs::compose(u.lz, lz) : lz;
            u.has_lz = true;
        }
    }

    void push(int t){
        Node& u = nodes[t];

        if (u.has_lz){
            apply_lazy(u.l, u.lz);
            apply_lazy(u.r, u.lz);
            u.has_lz = false;
        }
        if (u.rev){
            apply_reverse(u.l);
            apply_reverse(u.r);
            u.rev = false;
        }
    }

    /// Count of positions of t from its left end (right end if backwards) whose Agg passes ok, never on an empty range
    template<typename Ok>
    long long descend(int t, Ok& ok, bool backwards){
        long long taken = 0;
        bool empty = true;
        Agg acc{};
        auto extend = [&](const Agg& a){ return empty ? a : backwards ? Runs::combine(a, acc) : Runs::combine(acc, a); };

        while (t != -1){
            push(t);
            int near = backwards ? nodes[t].r : nodes[t].l, far = backwards ? nodes[t].l : nodes[t].r;
            if (near != -1){
                Agg a = extend(nodes[near].agg);
                if (!ok(a)){
                    t = near;
                    continue;
                }
                acc = a, empty = false, taken += nodes[near].size;
            }

            Value x = nodes[t].x;
            long long count = nodes[t].count;
            Agg whole = extend(Runs::run_agg(x, count));
            if (ok(whole)){
                acc = whole, empty = false, taken += count;
                t = far;
                continue;
            }

            /// The answer ends inside this run: binary search the longest part of it that still passes
            long long lo = 0, hi = count - 1;
            while (lo < hi){
                long long mid = lo + (hi - lo + 1) / 2;
                auto parts = Runs::cut(x, backwards ? count - mid : mid);
                if (ok(extend(Runs::run_agg(backwards ? parts.second : parts.first, mid)))) lo = mid;
                else hi = mid - 1;
            }
            return taken + lo;
        }

        return taken;
    }
};

int main(){
    /// Edits by position on a short sequence
    RunLengthImplicitTreap<SumRuns> seq(10, 0);           /// 0 0 0 0 0 0 0 0 0 0
    seq.assign(2, 5, 3);                                  /// 0 0 3 3 3 3 0 0 0 0
    seq.insert(4, 7, 2);                                  /// 0 0 3 3 7 7 3 3 0 0 0 0
    seq.update(5, 9, 1);                                   /// 0 0 3 3 7 8 4 4 1 1 0 0
    assert(seq.query(0, 11) == 31);
    seq.reverse(2, 6);                                    /// 0 0 4 8 7 3 3 4 1 1 0 0
    assert(seq.get(3) == 8);
    assert(seq.query(2, 4) == 19);

    /// A sequence far too long for an array
    RunLengthImplicitTreap<SumRuns> big(1000000000000000000LL, 1);
    big.assign(10, 19, 5);                                /// ten 1s, ten 5s, then 1s
    big.erase(0, 4);                                      /// five 1s, ten 5s, then 1s
    assert(big.size() == 999999999999999995LL);
    assert(big.query(0, 24) == 5 + 50 + 10);
    assert(big.max_right(0, [](long long s){ return s <= 20; }) == 7);   /// 1 1 1 1 1 5 5 5 sums to 20

    return 0;
}
