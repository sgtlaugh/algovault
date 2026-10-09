/***
 *
 * Persistent Leftist Heap
 * Meldable priority queue where every operation returns a new version and old versions stay valid
 *
 * Complexity: O(log n) time and at most 2 * log2(n + 1) + 2 new nodes per push, pop and meld,
 *   where n is the size of the largest heap involved
 *
 * A version is an int root into the object's node pool, EMPTY (-1) is the empty heap
 * Versions from one LeftistHeap object can only be passed back to that object
 *   LeftistHeap<T> heap: min-heap, LeftistHeap<T, greater<T>> is a max-heap
 *   push(h, x), pop(h), meld(a, b): return the new version, h, a and b are unchanged
 *   top(h): first value under Compare (the smallest for the default less<T>), pop and top need h != EMPTY
 *   meld(h, h) is valid and holds every value twice
 *
 * Nodes are never freed, so a run of q operations holds O(q log n) nodes
 * Usual use: k shortest walks (Eppstein), where each vertex's heap extends its parent's
 *
 *   LeftistHeap<int> heap;
 *   int a = heap.push(heap.push(LeftistHeap<int>::EMPTY, 20), 5);
 *   int b = heap.meld(a, heap.push(a, 1));
 *   // heap.top(a) = 5, heap.top(b) = 1, b holds {1, 5, 5, 20, 20}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
struct LeftistHeap{
    static constexpr int EMPTY = -1;

    struct Node{
        T val;
        int left, right, rank;
    };

    vector<Node> nodes;
    Compare cmp;

    int push(int h, const T& x){
        nodes.push_back({x, EMPTY, EMPTY, 1});
        return meld(h, (int)nodes.size() - 1);
    }

    int pop(int h){
        assert(h != EMPTY);
        return meld(nodes[h].left, nodes[h].right);
    }

    T top(int h) const{
        assert(h != EMPTY);
        return nodes[h].val;
    }

    /// Copies only the nodes on the merged right spines, whose length is bounded by the ranks
    int meld(int a, int b){
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        if (cmp(nodes[b].val, nodes[a].val)) swap(a, b);

        Node copy = nodes[a];
        nodes.push_back(copy);
        int c = (int)nodes.size() - 1;

        int right = meld(copy.right, b);  /// the call can reallocate nodes, assign after it
        nodes[c].right = right;
        if (rank(nodes[c].left) < rank(right)) swap(nodes[c].left, nodes[c].right);
        nodes[c].rank = rank(nodes[c].right) + 1;
        return c;
    }

    int rank(int h) const{
        return h == EMPTY ? 0 : nodes[h].rank;
    }
};

int main(){
    using Heap = LeftistHeap<int>;
    Heap heap;
    auto drain = [](auto& any_heap, int h){
        vector<decay_t<decltype(any_heap.top(h))>> res;
        for (; h != Heap::EMPTY; h = any_heap.pop(h)) res.push_back(any_heap.top(h));
        return res;
    };

    int h0 = Heap::EMPTY;
    int h1 = heap.push(h0, 20);
    int h2 = heap.push(h1, 10);
    int h3 = heap.push(h2, 30);
    int h4 = heap.push(heap.push(h1, 5), 15);
    int h5 = heap.meld(h4, h3);

    assert(heap.top(h1) == 20);
    assert(heap.top(h2) == 10);
    assert(heap.top(h3) == 10);
    assert(heap.top(h4) == 5);
    assert(heap.top(h5) == 5);
    assert((drain(heap, h5) == vector<int>{5, 10, 15, 20, 20, 30}));
    assert((drain(heap, h3) == vector<int>{10, 20, 30}));
    assert((drain(heap, h4) == vector<int>{5, 15, 20}));
    assert((drain(heap, h1) == vector<int>{20}));
    assert(heap.pop(h1) == Heap::EMPTY);
    assert(heap.meld(h0, h0) == Heap::EMPTY);
    assert(heap.meld(h0, h2) == h2);

    int a = heap.push(heap.push(Heap::EMPTY, 20), 5);
    int b = heap.meld(a, heap.push(a, 1));
    assert(heap.top(a) == 5 && heap.top(b) == 1);
    assert((drain(heap, b) == vector<int>{1, 5, 5, 20, 20}));
    assert((drain(heap, heap.meld(h2, h2)) == vector<int>{10, 10, 20, 20}));
    assert((drain(heap, heap.pop(heap.pop(h5))) == vector<int>{15, 20, 20, 30}));

    LeftistHeap<long long, greater<long long>> max_heap;
    int m = max_heap.EMPTY;
    for (long long x : {0LL, LLONG_MIN, LLONG_MAX, -7LL, LLONG_MAX, 0LL}) m = max_heap.push(m, x);
    assert(max_heap.top(m) == LLONG_MAX);
    assert((drain(max_heap, m) == vector<long long>{LLONG_MAX, LLONG_MAX, 0, 0, -7, LLONG_MIN}));

    LeftistHeap<int> chain;
    int c = Heap::EMPTY;
    for (int i = 100000; i >= 1; i--) c = chain.push(c, i);
    assert(chain.top(c) == 1);
    vector<int> sorted = drain(chain, c);
    for (int i = 0; i < 100000; i++) assert(sorted[i] == i + 1);

    return 0;
}
