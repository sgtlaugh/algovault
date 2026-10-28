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

int main(){
    int n, l;

    scanf("%d %d", &n, &l);

    MonotonicQueue<int> window;
    for (int i = 0; i < n; i++){
        int x;
        scanf("%d", &x);
        window.push(x);
        if (i >= l) window.pop();
        if (i >= l - 1) printf("%d%c", window.best(), i == n - 1 ? '\n' : ' ');
    }

    return 0;
}
