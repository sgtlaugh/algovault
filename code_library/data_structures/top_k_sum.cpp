/***
 *
 * Top-K Sum
 * A multiset that keeps the sum of its k largest values while values come and go and k changes
 *
 * Complexity: O(log n) amortized per insert and erase, O(|k' - k| log n) for set_k(k'), n the number of calls
 *
 * TopKSum s(k): an empty multiset summing its k largest values, as long long
 * TopKSum<T, Compare> s(k): sums the first k values under Compare instead, less<T> gives the k smallest
 * s.erase(x): removes one copy of x, which must be present
 * s.insert(x): adds one copy of x
 * s.set_k(k): sums the k largest from now on
 * s.size(): the number of values
 * s.sum: the sum of the min(k, size()) largest values, 0 when k = 0 or the multiset is empty
 *
 * Notes:
 *   - erase of a value never inserted is not detected and corrupts every later answer
 *   - Compare may tie only equal values: an order by |x| ties 3 and -3, and erase then corrupts the sum
 *   - T is both the value and the sum type, so int values with large sums need TopKSum<long long>
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T = long long, typename Compare = greater<T>>
struct TopKSum{
    struct Reversed{
        bool operator()(const T& a, const T& b) const{
            return Compare()(b, a);
        }
    };

    /// Erased values wait in a second heap until they reach the top: 5-8x faster than two multisets
    template <typename Cmp>
    struct LazyHeap{
        priority_queue<T, vector<T>, Cmp> values, erased;
        int size = 0;

        void erase(const T& x){
            erased.push(x);
            size--;
        }

        void pop(){
            top();
            values.pop();
            size--;
        }

        void push(const T& x){
            values.push(x);
            size++;
        }

        T top(){
            while (!erased.empty() && values.top() == erased.top()){
                values.pop();
                erased.pop();
            }

            return values.top();
        }
    };

    int k;
    T sum = 0;
    LazyHeap<Compare> top;    /// the counted values, the last of them under Compare on top
    LazyHeap<Reversed> rest;  /// every other value, the first of them under Compare on top

    TopKSum(int k) : k(k){
        assert(k >= 0);
    }

    /// An x equal to the last counted value is taken from top even if a copy sits in rest: the sums agree
    void erase(const T& x){
        assert(size() > 0);
        if (top.size > 0 && !Compare()(top.top(), x)){
            top.erase(x);
            sum -= x;
        }
        else rest.erase(x);

        rebalance();
    }

    /// top.size < k implies rest is empty, so x may join top without checking rest
    void insert(const T& x){
        if (top.size < k || (top.size > 0 && Compare()(x, top.top()))){
            top.push(x);
            sum += x;
        }
        else rest.push(x);

        rebalance();
    }

    void set_k(int new_k){
        assert(new_k >= 0);
        k = new_k;
        rebalance();
    }

    int size() const{
        return top.size + rest.size;
    }

    void rebalance(){
        while (top.size > k){
            T x = top.top();
            top.pop();
            sum -= x;
            rest.push(x);
        }

        while (top.size < k && rest.size > 0){
            T x = rest.top();
            rest.pop();
            sum += x;
            top.push(x);
        }
    }
};

int main(){
    /// Values {5, 1, 8, 3}: the 2 largest are 8 and 5
    TopKSum s(2);
    s.insert(5);
    s.insert(1);
    s.insert(8);
    s.insert(3);
    assert(s.sum == 13);

    s.erase(8);
    assert(s.sum == 8);  /// 5 + 3

    s.set_k(3);
    assert(s.sum == 9);  /// 5 + 3 + 1

    s.set_k(10);
    s.insert(5);
    assert(s.sum == 14);  /// every value counts: 5 + 5 + 3 + 1
    assert(s.size() == 4);

    s.set_k(0);
    assert(s.sum == 0);

    /// The 2 smallest of {4, 4, 9, 2}: 2 and 4
    TopKSum<int, less<int>> low(2);
    for (int x : {4, 4, 9, 2}) low.insert(x);
    assert(low.sum == 6);

    low.erase(2);
    assert(low.sum == 8);  /// 4 + 4

    return 0;
}
