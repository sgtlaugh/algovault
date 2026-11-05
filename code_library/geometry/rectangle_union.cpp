/***
 *
 * Rectangle Union
 * Area and perimeter of the union of axis-aligned rectangles
 *
 * Complexity: O(n log n)
 *
 * rectangle_union_area(rects): area of the union
 * rectangle_union_perimeter(rects): perimeter of the union, hole boundaries included, edges between touching ones not
 * rectangle_union_sweep(rects): {area, length of the vertical boundary} of the union
 *
 * Notes:
 *   - A rectangle is (x1, y1, x2, y2), empty ones (x1 >= x2 or y1 >= y2) are ignored
 *   - Coordinates within [-1e9, 1e9] keep every result inside long long (the area is at most 4e18)
 *   - rectangle_union_perimeter sweeps twice, rectangle_union_sweep gets the area and vertical boundary in one pass
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Rectangle{
    long long x1, y1, x2, y2;
};

/// Every removal undoes an earlier add of the same range, so a count above 0 means the node's whole range is covered
struct CoverTree{
    struct Node{
        int count;
        long long length, covered;
    };

    int n;
    vector<Node> tree;

    /// A power of two keeps every leaf on one level, so tree[1] spans all of ys
    CoverTree(const vector<long long>& ys) : n(1){
        while (n < (int)ys.size() - 1) n *= 2;
        tree.assign(2 * n, Node{});
        for (int i = 0; i + 1 < (int)ys.size(); i++) tree[n + i].length = ys[i + 1] - ys[i];
        for (int i = n - 1; i >= 1; i--) tree[i].length = tree[2 * i].length + tree[2 * i + 1].length;
    }

    long long covered() const{
        return tree[1].covered;
    }

    void update(int l, int r, int delta){
        int lo = l + n, hi = r + n;
        for (l += n, r += n + 1; l < r; l /= 2, r /= 2){
            if (l & 1) add(l++, delta);
            if (r & 1) add(--r, delta);
        }

        /// Every changed node's parent is an ancestor of leaf lo or hi; above where the two paths merge, one walk does
        for (lo /= 2, hi /= 2; lo != hi; lo /= 2, hi /= 2){
            pull(lo);
            pull(hi);
        }
        for (; lo >= 1; lo /= 2) pull(lo);
    }

    void add(int node, int delta){
        tree[node].count += delta;
        pull(node);
    }

    void pull(int node){
        Node& v = tree[node];
        v.covered = v.count > 0 ? v.length : node >= n ? 0 : tree[2 * node].covered + tree[2 * node + 1].covered;
    }
};

/***
 *
 * Returns {union area, total length of the union's vertical boundary}
 * Events at the same x open before they close: a seam where one rectangle ends and another
 * starts then adds and removes the same length, so it never counts as boundary
 *
***/

pair<long long, long long> rectangle_union_sweep(const vector<Rectangle>& rects){
    vector<long long> ys;
    vector<tuple<long long, int, long long, long long>> events;

    for (const auto& r: rects){
        if (r.x1 >= r.x2 || r.y1 >= r.y2) continue;
        ys.push_back(r.y1), ys.push_back(r.y2);
        events.emplace_back(r.x1, -1, r.y1, r.y2);
        events.emplace_back(r.x2, 1, r.y1, r.y2);
    }
    if (events.empty()) return {0, 0};

    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());
    sort(events.begin(), events.end());
    CoverTree tree(ys);

    long long area = 0, boundary = 0, last_x = get<0>(events[0]);
    for (auto [x, type, y1, y2]: events){
        long long before = tree.covered();
        area += before * (x - last_x);
        last_x = x;

        int l = lower_bound(ys.begin(), ys.end(), y1) - ys.begin();
        int r = lower_bound(ys.begin(), ys.end(), y2) - ys.begin() - 1;
        tree.update(l, r, -type);
        boundary += abs(tree.covered() - before);
    }

    return {area, boundary};
}

long long rectangle_union_area(const vector<Rectangle>& rects){
    return rectangle_union_sweep(rects).first;
}

long long rectangle_union_perimeter(const vector<Rectangle>& rects){
    vector<Rectangle> transposed;
    for (const auto& r: rects) transposed.push_back({r.y1, r.x1, r.y2, r.x2});
    return rectangle_union_sweep(rects).second + rectangle_union_sweep(transposed).second;
}

int main(){
    /// Two 2 x 2 squares overlapping in a unit square: 4 + 4 - 1
    vector<Rectangle> rects = {{0, 0, 2, 2}, {1, 1, 3, 3}};
    assert(rectangle_union_area(rects) == 7);
    assert(rectangle_union_perimeter(rects) == 12);

    /// A shared edge is not boundary: two unit squares side by side form a 2 x 1 rectangle
    assert(rectangle_union_perimeter({{0, 0, 1, 1}, {1, 0, 2, 1}}) == 6);

    /// A 3 x 3 frame around a unit hole: the hole's boundary counts
    vector<Rectangle> frame = {{0, 0, 3, 1}, {0, 2, 3, 3}, {0, 1, 1, 2}, {2, 1, 3, 2}};
    assert(rectangle_union_area(frame) == 8);
    assert(rectangle_union_perimeter(frame) == 16);

    /// One sweep gives the area and the vertical boundary only: 1 + 1 + 3 + 3
    assert(rectangle_union_sweep({{0, 0, 1, 1}, {5, 5, 7, 8}}) == make_pair(7LL, 8LL));
    return 0;
}
