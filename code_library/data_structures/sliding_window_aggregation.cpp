/***
 *
 * Sliding Window Aggregation
 * A queue that folds its elements in queue order under any associative operation, commutative or not
 *
 * Complexity: O(1) amortized per operation, a pop that refills the front stack takes O(n)
 *
 * SlidingWindowAggregation q(identity, op): an empty queue
 * q.pop(): removes the oldest element, the queue must not be empty
 * q.push(x): appends x
 * q.empty(): true if the queue has no elements
 * q.fold(): op(...op(op(a_0, a_1), a_2)..., a_k) over the elements oldest first, identity when empty
 * q.size(): the number of elements
 * MonotonicQueue<T, Compare = less<T>> q: an empty queue with the same push, pop, size and empty
 * q.best(): the minimum under Compare of the elements, the queue must not be empty
 *
 * Notes:
 *   - op only needs associativity, and identity only needs op(identity, x) = x
 *   - For min or max use MonotonicQueue: no identity and fewer comparisons, greater<T> gives the maximum
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

    void pop(){
        assert(!empty());
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

    bool empty() const{
        return front.empty() && back.empty();
    }

    T fold() const{
        if (front.empty()) return back_fold;
        return back.empty() ? front.back() : combine(front.back(), back_fold);
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

    void pop(){
        assert(!empty());
        if (candidates.front().second == popped) candidates.pop_front();
        popped++;
    }

    void push(const T& x){
        while (!candidates.empty() && !cmp(candidates.back().first, x)) candidates.pop_back();
        candidates.emplace_back(x, pushed++);
    }

    T best() const{
        assert(!empty());
        return candidates.front().first;
    }

    bool empty() const{
        return pushed == popped;
    }

    int size() const{
        return pushed - popped;
    }
};

int main(){
    /// f(x) = a x + b as (a, b), composition does not commute so the fold order shows
    using Linear = pair<long long, long long>;
    auto compose = [](const Linear& f, const Linear& g){
        return Linear(f.first * g.first, f.second * g.first + g.second);
    };

    SlidingWindowAggregation maps(Linear(1, 0), compose);
    maps.push({2, 3});
    maps.push({5, 1});
    assert(maps.fold() == Linear(10, 16));  /// 5 (2x + 3) + 1, the oldest map applies first
    maps.pop();
    maps.push({3, 7});
    assert(maps.fold() == Linear(15, 10));  /// 3 (5x + 1) + 7
    maps.pop();
    maps.pop();
    assert(maps.empty());
    assert(maps.fold() == Linear(1, 0));

    /// Minimum and maximum of every window of 3
    vector<int> a = {4, 2, 12, 3, 8, 1, 7}, lows, highs;
    MonotonicQueue<int> low;

    MonotonicQueue<int, greater<int>> high;
    for (int i = 0; i < (int)a.size(); i++){
        low.push(a[i]);
        high.push(a[i]);
        if (i >= 3){
            low.pop();
            high.pop();
        }
        if (i >= 2){
            lows.push_back(low.best());
            highs.push_back(high.best());
        }
    }
    assert((lows == vector<int>{2, 2, 3, 1, 1}));
    assert((highs == vector<int>{12, 12, 12, 8, 8}));

    return 0;
}
