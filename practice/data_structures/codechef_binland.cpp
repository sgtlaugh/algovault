/***
 *
 * Sliding Window Aggregation
 * A queue that folds its elements in queue order under any associative operation, commutative or not
 *
 * Complexity: O(1) amortized per operation, a pop that refills the front stack takes O(n)
 *
 * SlidingWindowAggregation q(identity, op): an empty queue
 * q.push(x): appends x
 * q.pop(): removes the oldest element, the queue must not be empty
 * q.fold(): op(...op(op(a_0, a_1), a_2)..., a_k) over the elements oldest first, identity when empty
 * q.size(): the number of elements
 * q.empty(): true if the queue has no elements
 * MonotonicQueue<T, Compare = less<T>> q: an empty queue with the same push, pop, size and empty
 * q.best(): the minimum under Compare of the elements, the queue must not be empty
 *
 * op only needs associativity, and identity only needs op(identity, x) = x since it is never combined on the right
 * For min or max use MonotonicQueue: no identity and fewer comparisons, greater<T> gives the maximum
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op>
struct SlidingWindowAggregation{
    T identity;
    Op combine;
    vector<T> front, back;  /// front.back() is the fold of the whole front stack, oldest element first
    T back_fold;            /// op over the back stack, oldest element first

    SlidingWindowAggregation(T identity, Op combine) : identity(identity), combine(combine), back_fold(identity) {}

    bool empty() const{
        return front.empty() && back.empty();
    }

    T fold() const{
        if (front.empty()) return back_fold;
        return back.empty() ? front.back() : combine(front.back(), back_fold);
    }

    void pop(){
        if (front.empty()){
            for (int i = (int)back.size() - 1; i >= 0; i--){
                front.push_back(front.empty() ? back[i] : combine(back[i], front.back()));
            }
            back.clear();
            back_fold = identity;
        }
        front.pop_back();
    }

    void push(const T& x){
        back.push_back(x);
        back_fold = combine(back_fold, x);
    }

    int size() const{
        return front.size() + back.size();
    }
};

template <typename T, typename Compare = less<T>>
struct MonotonicQueue{
    Compare cmp;
    deque<pair<T, long long>> candidates;  /// increasing push index, best value at the front
    long long pushed = 0, popped = 0;

    T best() const{
        return candidates.front().first;
    }

    bool empty() const{
        return pushed == popped;
    }

    void pop(){
        if (candidates.front().second == popped) candidates.pop_front();
        popped++;
    }

    void push(const T& x){
        while (!candidates.empty() && !cmp(candidates.back().first, x)) candidates.pop_back();
        candidates.emplace_back(x, pushed++);
    }

    int size() const{
        return pushed - popped;
    }
};

const unsigned MOD = 1000000007;

int n;

/// Path counts from each column of the top row to each column of the bottom row; m stays empty for a single row
struct Strip{
    int rows, top, bottom;
    vector<unsigned> m;

    unsigned at(int i, int j) const{
        return m.empty() ? i == j : m[i * n + j];
    }
};

bool linked(int upper, int lower, int j, int k){
    return abs(j - k) <= 1 && ((upper >> j) & 1) == ((lower >> k) & 1);
}

Strip join(const Strip& a, const Strip& b){
    if (!a.rows) return b;
    if (!b.rows) return a;

    Strip c = {a.rows + b.rows, a.top, b.bottom, vector<unsigned>(n * n, 0)};
    vector<pair<int, int>> links;
    for (int j = 0; j < n; j++){
        for (int k = max(0, j - 1); k <= min(n - 1, j + 1); k++){
            if (linked(a.bottom, b.top, j, k)) links.emplace_back(j, k);
        }
    }

    vector<unsigned> through(n * n, 0);  /// through[i][k]: paths from a's top column i to b's top column k
    if (a.m.empty()){
        for (auto& l : links) through[l.first * n + l.second] = 1;
    }
    else{
        for (int i = 0; i < n; i++){
            const unsigned* from = &a.m[i * n];
            unsigned* to = &through[i * n];
            for (auto& l : links) to[l.second] += from[l.first];  /// at most 3 terms below MOD, fits in 32 bits
            for (int k = 0; k < n; k++) to[k] %= MOD;
        }
    }

    if (b.m.empty()){
        c.m = through;
        return c;
    }
    for (int i = 0; i < n; i++){
        unsigned long long s[20] = {};
        for (int k = 0; k < n; k++){
            unsigned long long x = through[i * n + k];
            if (!x) continue;
            for (int d = 0; d < n; d++){
                s[d] += x * b.m[k * n + d];
                if (s[d] >> 63) s[d] %= MOD;  /// each term is below 2^60, so this never overflows
            }
        }
        for (int d = 0; d < n; d++) c.m[i * n + d] = s[d] % MOD;
    }
    return c;
}

int main(){
    int q;
    scanf("%d %d", &n, &q);

    SlidingWindowAggregation<Strip, Strip(*)(const Strip&, const Strip&)> grid(Strip{0, 0, 0, {}}, join);
    char cmd[16], row[32];
    while (q--){
        scanf("%s", cmd);
        if (cmd[0] == 'a'){
            scanf("%s", row);
            int mask = 0;
            for (int i = 0; i < n; i++) mask |= (row[i] - '0') << i;
            grid.push(Strip{1, mask, mask, {}});
        }
        else if (cmd[0] == 'r'){
            grid.pop();
        }
        else{
            int c, d;
            scanf("%d %d", &c, &d);
            c--, d--;

            /// fold() would multiply two full matrices, O(n^3); one entry of that product only needs O(n)
            unsigned long long res = 0;
            if (grid.front.empty()) res = grid.back_fold.at(c, d);
            else if (grid.back.empty()) res = grid.front.back().at(c, d);
            else{
                const Strip& up = grid.front.back();
                const Strip& down = grid.back_fold;
                for (int j = 0; j < n; j++){
                    for (int k = max(0, j - 1); k <= min(n - 1, j + 1); k++){
                        if (linked(up.bottom, down.top, j, k)) res = (res + (unsigned long long)up.at(c, j) * down.at(k, d)) % MOD;
                    }
                }
            }
            printf("%llu\n", res);
        }
    }
    return 0;
}
