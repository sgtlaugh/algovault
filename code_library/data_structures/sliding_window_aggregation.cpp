/***
 *
 * Sliding Window Aggregation
 * A queue that folds its elements in queue order under any associative operation, commutative or not
 *
 * Complexity: O(1) amortized per push, pop and fold, each element is combined at most twice
 * A single pop can take O(n) when it refills the front stack, so the bound is amortized, not worst case
 *
 * SlidingWindowAggregation(identity, op): two stacks, pushes go to the back, pops come from the front
 *   push(x) appends x, pop() removes the oldest element (the queue must not be empty)
 *   fold() returns op(...op(op(a_0, a_1), a_2)..., a_k) over the elements oldest first, identity when empty
 *   op only needs associativity and identity only needs op(identity, x) = x (identity is never combined on the right),
 *   so composition of functions, matrix products and string concatenation work as well as sum, min, gcd
 *   Composition f_k(...f_0(x)) of linear maps f(x) = a x + b: op(f, g) = (f.a * g.a, f.b * g.a + g.b)
 *
 * MonotonicQueue<T, Compare = less<T>>: the same queue API restricted to min (or max with greater<T>)
 *   best() returns the minimum under Compare of the current elements, best() and pop() need a non-empty queue
 *   It keeps only the elements that can still become the answer, fewer comparisons and no identity needed
 *
 * Sliding minimum of every window of size k:
 *   MonotonicQueue<int> q;
 *   for (int i = 0; i < n; i++){
 *       q.push(a[i]);
 *       if (i >= k) q.pop();
 *       if (i >= k - 1) res.push_back(q.best());
 *   }
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op>
struct SlidingWindowAggregation{
    T identity;
    Op op;
    vector<T> front, back;  /// front.back() is the fold of the whole front stack, oldest element first
    T back_fold;

    SlidingWindowAggregation(T identity, Op op) : identity(identity), op(op), back_fold(identity) {}

    bool empty() const{
        return front.empty() && back.empty();
    }

    T fold() const{
        if (front.empty()) return back_fold;
        return back.empty() ? front.back() : op(front.back(), back_fold);
    }

    void pop(){
        if (front.empty()){
            for (int i = (int)back.size() - 1; i >= 0; i--){
                front.push_back(front.empty() ? back[i] : op(back[i], front.back()));
            }
            back.clear();
            back_fold = identity;
        }
        front.pop_back();
    }

    void push(const T& x){
        back.push_back(x);
        back_fold = op(back_fold, x);
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

int main(){
    using Linear = pair<long long, long long>;
    auto compose = [](const Linear& f, const Linear& g){
        return Linear(f.first * g.first, f.second * g.first + g.second);
    };
    auto eval = [](const Linear& f, long long x){
        return f.first * x + f.second;
    };

    SlidingWindowAggregation maps(Linear(1, 0), compose);
    assert(maps.empty() && eval(maps.fold(), 4) == 4);
    maps.push({2, 3});
    maps.push({5, 1});
    assert(eval(maps.fold(), 4) == 56);
    maps.pop();
    maps.push({3, 7});
    assert(maps.size() == 2 && eval(maps.fold(), 4) == 70);
    maps.pop();
    assert(eval(maps.fold(), 4) == 19);

    SlidingWindowAggregation words(string(), [](const string& a, const string& b){ return a + b; });
    for (string w : {"slow", "ly", "but"}) words.push(w);
    assert(words.fold() == "slowlybut");
    words.pop();
    words.push("ter");
    assert(words.fold() == "lybutter");
    words.pop();
    words.pop();
    assert(words.fold() == "ter");
    words.pop();
    assert(words.empty() && words.fold() == "");

    SlidingWindowAggregation newest(-1, [](int, int b){ return b; });
    newest.push(1);
    newest.push(2);
    newest.pop();
    assert(newest.fold() == 2);
    newest.push(3);
    assert(newest.fold() == 3);
    newest.pop();
    newest.pop();
    assert(newest.fold() == -1);

    vector<int> a = {4, 2, 12, 3, 8, 1, 7}, lows, highs;
    MonotonicQueue<int> low;
    MonotonicQueue<int, greater<int>> high;
    for (int i = 0; i < (int)a.size(); i++){
        low.push(a[i]), high.push(a[i]);
        if (i >= 3) low.pop(), high.pop();
        if (i >= 2) lows.push_back(low.best()), highs.push_back(high.best());
    }
    assert((lows == vector<int>{2, 2, 3, 1, 1}));
    assert((highs == vector<int>{12, 12, 12, 8, 8}));

    MonotonicQueue<int> ties;
    for (int x : {5, 5, 9}) ties.push(x);
    ties.pop();
    assert(ties.size() == 2 && ties.best() == 5);
    ties.pop();
    assert(ties.best() == 9);
    ties.pop();
    assert(ties.empty());

    return 0;
}
