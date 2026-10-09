#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/leftist_heap.cpp"
#undef main

int log2_floor(long long n){
    return 63 - __builtin_clzll(n);
}

/// Pops the whole version, which must not disturb it or any other version
template <typename Heap, typename T>
void check_drain(Heap& heap, int h, const vector<T>& expected){
    auto it = expected.begin();
    for (; h != Heap::EMPTY; h = heap.pop(h), it++){
        assert(it != expected.end());
        assert(heap.top(h) == *it);
    }
    assert(it == expected.end());
}

/// Random push, pop and meld on random old versions, each version checked against its own sorted multiset
template <typename T, typename Compare>
void run(int ops, long long lo, long long hi, size_t max_size, int drain_every){
    LeftistHeap<T, Compare> heap;
    vector<int> versions = {LeftistHeap<T, Compare>::EMPTY};
    vector<vector<T>> expected(1);
    Compare cmp;

    for (int op = 0; op < ops; op++){
        int i = stress::rand_int(0, (int)versions.size() - 1);
        int j = stress::rand_int(0, 3) ? stress::rand_int(0, (int)versions.size() - 1) : i;
        int kind = stress::rand_int(0, 5);
        size_t before = heap.nodes.size();
        long long budget;
        vector<T> next = expected[i];
        int h;

        if (kind <= 2 || expected[i].empty()){
            T x = (T)stress::rand_int(lo, hi);
            h = heap.push(versions[i], x);
            next.insert(upper_bound(next.begin(), next.end(), x, cmp), x);
            budget = log2_floor(expected[i].size() + 1) + 2;
        }
        else if (kind <= 4 || expected[i].size() + expected[j].size() > max_size){
            h = heap.pop(versions[i]);
            next.erase(next.begin());
            budget = 2 * log2_floor(expected[i].size() + 1);
        }
        else{
            h = heap.meld(versions[i], versions[j]);
            next.clear();
            merge(expected[i].begin(), expected[i].end(), expected[j].begin(), expected[j].end(), back_inserter(next), cmp);
            budget = 2 * log2_floor(max(expected[i].size(), expected[j].size()) + 1);
        }

        assert((long long)(heap.nodes.size() - before) <= budget);
        assert((h == LeftistHeap<T, Compare>::EMPTY) == next.empty());
        if (!next.empty()) assert(heap.top(h) == *next.begin());
        versions.push_back(h);
        expected.push_back(next);

        int k = stress::rand_int(0, (int)versions.size() - 1);
        if (!expected[k].empty()) assert(heap.top(versions[k]) == *expected[k].begin());
        if (op % drain_every == 0) check_drain(heap, versions[k], expected[k]);
    }

    for (size_t k = 0; k < versions.size(); k++) check_drain(heap, versions[k], expected[k]);
}

int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int ops = it < 300 ? it + 1 : stress::rand_int(1, 120);
        if (it % 3 == 0) run<int, less<int>>(ops, -5, 5, 40, 1);
        else if (it % 3 == 1) run<int, greater<int>>(ops, 0, 3, 40, 1);
        else run<long long, less<long long>>(ops, LLONG_MIN, LLONG_MAX, 60, 3);
    }

    for (long long it = 0; it < stress::scaled(10); it++){
        run<long long, less<long long>>(3000, -1000000000, 1000000000, 1000, 50);
        run<long long, greater<long long>>(3000, LLONG_MIN, LLONG_MAX, 1000, 50);
    }

    /// One long chain of versions with large heaps: sorted input is the worst case for the right spine
    for (long long it = 0; it < stress::scaled(2); it++){
        const int n = 100000;
        LeftistHeap<int> heap;
        vector<int> roots = {LeftistHeap<int>::EMPTY};
        for (int i = 0; i < n; i++){
            int x = it ? n - i : i;
            size_t before = heap.nodes.size();
            roots.push_back(heap.push(roots.back(), x));
            assert((long long)(heap.nodes.size() - before) <= log2_floor(i + 1) + 2);
        }

        int whole = heap.meld(roots[n], roots[n / 2]);
        vector<int> got;
        for (int h = whole; h != LeftistHeap<int>::EMPTY; h = heap.pop(h)) got.push_back(heap.top(h));
        vector<int> want;
        for (int i = 0; i < n; i++) want.push_back(it ? n - i : i);
        for (int i = 0; i < n / 2; i++) want.push_back(it ? n - i : i);
        sort(want.begin(), want.end());
        assert(got == want);

        for (int i : {1, n / 3, n - 1}) assert(heap.top(roots[i]) == (it ? n - i + 1 : 0));
    }

    return 0;
}
