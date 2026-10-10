/***
 *
 * Rectangle Union
 * Area and perimeter of the union of axis-aligned rectangles
 *
 * Complexity: O(n log n), O(n) memory
 *
 * Sweeps a vertical line over x, a segment tree over the compressed y coordinates keeps
 * the minimum cover count and the length covered by that minimum, so the covered length is
 * the total minus that length whenever the minimum is 0
 *
 * A rectangle is (x1, y1, x2, y2) with x1 < x2 and y1 < y2, empty ones (x1 >= x2 or y1 >= y2) are ignored
 * Coordinates within [-1e9, 1e9] keep every result inside long long (the area is at most 4e18)
 * The perimeter counts hole boundaries too and never counts edges shared by touching rectangles
 * rectangle_union_perimeter sweeps twice (original and transposed), call rectangle_union_sweep
 * directly to get {area, vertical boundary length} in one pass
 *
 * Example:
 *   vector<Rectangle> rects = {{0, 0, 2, 2}, {1, 1, 3, 3}};
 *   rectangle_union_area(rects);       // 7
 *   rectangle_union_perimeter(rects);  // 12
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Rectangle{
    long long x1, y1, x2, y2;
};

/// Counts never go negative, so the root's minimum is 0 exactly when some cell is uncovered
struct CoverTree{
    int m;
    vector<long long> ys, len;
    vector<int> low, add;

    CoverTree(const vector<long long>& ys) : m(ys.size() - 1), ys(ys), len(4 * m), low(4 * m), add(4 * m){
        build(1, 0, m - 1);
    }

    long long covered() const{
        return ys.back() - ys.front() - (low[1] == 0 ? len[1] : 0);
    }

    void update(int l, int r, int delta){
        update(1, 0, m - 1, l, r, delta);
    }

    void build(int node, int l, int r){
        if (l == r){
            len[node] = ys[l + 1] - ys[l];
            return;
        }

        int mid = (l + r) / 2;
        build(2 * node, l, mid);
        build(2 * node + 1, mid + 1, r);
        pull(node);
    }

    /// No push down: a node's add applies to its whole range, so it is folded into low on the way up
    void pull(int node){
        int a = 2 * node, b = 2 * node + 1;
        int child_low = min(low[a], low[b]);
        low[node] = child_low + add[node];
        len[node] = (low[a] == child_low ? len[a] : 0) + (low[b] == child_low ? len[b] : 0);
    }

    void update(int node, int l, int r, int ql, int qr, int delta){
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr){
            low[node] += delta;
            add[node] += delta;
            return;
        }

        int mid = (l + r) / 2;
        update(2 * node, l, mid, ql, qr, delta);
        update(2 * node + 1, mid + 1, r, ql, qr, delta);
        pull(node);
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

    /// A rectangle inside another adds nothing
    assert(rectangle_union_area({{0, 0, 10, 10}, {2, 2, 3, 3}}) == 100);

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
