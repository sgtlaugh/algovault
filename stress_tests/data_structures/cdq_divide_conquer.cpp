#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/cdq_divide_conquer.cpp"
#undef main

struct Item{
    long long x, y, z, w;
    bool query;
};

vector<long long> brute(const vector<Item>& items){
    vector<long long> res;
    for (const Item& q : items){
        if (!q.query) continue;
        long long sum = 0;
        for (const Item& p : items){
            if (!p.query && p.x <= q.x && p.y <= q.y && p.z <= q.z) sum += p.w;
        }
        res.push_back(sum);
    }

    return res;
}

vector<long long> run(const vector<Item>& items){
    DominanceSum<long long> ds;
    int next = 0;
    for (const Item& it : items){
        if (it.query) assert(ds.add_query(it.x, it.y, it.z) == next++);
        else ds.add_point(it.x, it.y, it.z, it.w);
    }

    vector<long long> res = ds.solve();
    assert(ds.solve() == res);
    return res;
}

/// Points and queries interleaved in random order, coordinates from [lo, hi] so ties on every axis are common
vector<Item> random_items(int n, long long lo, long long hi, long long max_w){
    vector<Item> items(n);
    for (Item& it : items){
        it = {stress::rand_int(lo, hi), stress::rand_int(lo, hi), stress::rand_int(lo, hi), stress::rand_int(-max_w, max_w),
              stress::rand_int(0, 1) == 1};
    }

    return items;
}

int main(){
    assert(run({}).empty());
    assert((run({{5, 5, 5, 0, true}}) == vector<long long>{0}));
    assert((run({{5, 5, 5, 7, false}, {5, 5, 5, 0, true}, {4, 5, 5, 0, true}}) == vector<long long>{7, 0}));

    for (int it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 40);
        long long range = stress::rand_int(0, 4);
        vector<Item> items = random_items(n, -range, range, 9);
        assert(run(items) == brute(items));
    }

    /// Full long long range, where a coordinate difference would overflow
    for (int it = 0; it < stress::scaled(200); it++){
        vector<Item> items = random_items(stress::rand_int(1, 60), LLONG_MIN, LLONG_MAX, 1000000000);
        for (int k = 0; k < 6; k++) items.push_back({LLONG_MAX, LLONG_MAX, LLONG_MAX, 0, true});
        items.push_back({LLONG_MIN, LLONG_MIN, LLONG_MIN, 3, false});
        assert(run(items) == brute(items));
    }

    /// The header's recipe: a query at every point counts the others it dominates once its own weight is subtracted
    for (int it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(1, 50);
        vector<Item> points = random_items(n, 0, 3, 1);
        DominanceSum<long long> ds;
        for (Item& p : points) p.w = 1, ds.add_point(p.x, p.y, p.z);
        for (const Item& p : points) ds.add_query(p.x, p.y, p.z);
        vector<long long> res = ds.solve();
        for (int i = 0; i < n; i++){
            int others = 0;
            for (int j = 0; j < n; j++){
                if (j != i && points[j].x <= points[i].x && points[j].y <= points[i].y && points[j].z <= points[i].z) others++;
            }
            assert(res[i] - 1 == others);
        }
    }

    /// Larger inputs cross many recursion levels
    for (int it = 0; it < stress::scaled(4); it++){
        vector<Item> items = random_items(stress::rand_int(1000, 3000), 0, stress::rand_int(10, 1000000), 1000000);
        assert(run(items) == brute(items));
    }

    return 0;
}
