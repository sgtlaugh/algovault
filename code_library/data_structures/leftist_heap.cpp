/***
 *
 * Persistent Leftist Heap
 * Meldable priority queue where every operation returns a new version and old versions stay valid
 *
 * Complexity: O(log n) time and at most 2 * log2(n + 1) + 2 new nodes per push, pop and meld, n = size of the largest heap involved
 *
 * LeftistHeap<T> heap: a pool of min-heaps, each version is an int and EMPTY (-1) is the empty heap
 * LeftistHeap<T, greater<T>> heap: the same pool of max-heaps
 * heap.push(h, x): a new version holding the values of h plus x
 * heap.pop(h): a new version holding the values of h without its top, h != EMPTY
 * heap.top(h): first value of h under Compare, the smallest for less<T>, h != EMPTY
 * heap.meld(a, b): a new version holding the values of a and b, meld(h, h) holds every value twice
 *
 * No call changes an existing version, and a version only works with the object that made it
 * Nodes are never freed, so q operations hold O(q log n) nodes
 * Usual use: k shortest walks (Eppstein), where each vertex's heap extends its parent's
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED LeftistHeap
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
/// END SHARED LeftistHeap

int main(){
    /// Versions branch off each other and stay valid
    using Heap = LeftistHeap<int>;
    Heap heap;
    int h1 = heap.push(Heap::EMPTY, 20);
    int h2 = heap.push(h1, 10);
    int h3 = heap.push(h1, 30);
    assert(heap.top(h2) == 10);
    assert(heap.top(h3) == 20);
    assert(heap.top(h1) == 20);

    int both = heap.meld(h2, h3);  /// {10, 20, 20, 30}
    assert(heap.top(both) == 10);
    int rest = heap.pop(both);
    assert(heap.top(rest) == 20);
    assert(heap.pop(h1) == Heap::EMPTY);

    /// Max-heap
    using MaxHeap = LeftistHeap<int, greater<int>>;
    MaxHeap max_heap;
    int m1 = max_heap.push(MaxHeap::EMPTY, 7);
    int m2 = max_heap.push(m1, 42);
    assert(max_heap.top(m2) == 42);

    return 0;
}
