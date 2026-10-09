#include "../common.h"

#define main library_main
#include "../../code_library/geometry/rectangle_union.cpp"
#undef main

/// Paints every rectangle onto the grid of compressed cells, then sums painted cell areas and painted/unpainted cell borders
pair<__int128, __int128> paint(const vector<Rectangle>& rects){
    vector<long long> xs, ys;
    for (const auto& r: rects){
        if (r.x1 >= r.x2 || r.y1 >= r.y2) continue;
        xs.push_back(r.x1), xs.push_back(r.x2), ys.push_back(r.y1), ys.push_back(r.y2);
    }
    sort(xs.begin(), xs.end()), xs.erase(unique(xs.begin(), xs.end()), xs.end());
    sort(ys.begin(), ys.end()), ys.erase(unique(ys.begin(), ys.end()), ys.end());
    int nx = xs.size(), ny = ys.size();

    vector<vector<int>> diff(nx + 1, vector<int>(ny + 1, 0));
    for (const auto& r: rects){
        if (r.x1 >= r.x2 || r.y1 >= r.y2) continue;
        int a = lower_bound(xs.begin(), xs.end(), r.x1) - xs.begin(), b = lower_bound(xs.begin(), xs.end(), r.x2) - xs.begin();
        int c = lower_bound(ys.begin(), ys.end(), r.y1) - ys.begin(), d = lower_bound(ys.begin(), ys.end(), r.y2) - ys.begin();
        diff[a][c]++, diff[b][c]--, diff[a][d]--, diff[b][d]++;
    }
    for (int i = 0; i <= nx; i++){
        for (int j = 0; j <= ny; j++){
            if (i) diff[i][j] += diff[i - 1][j];
            if (j) diff[i][j] += diff[i][j - 1];
            if (i && j) diff[i][j] -= diff[i - 1][j - 1];
        }
    }

    auto painted = [&](int i, int j){ return i >= 0 && j >= 0 && i + 1 < nx && j + 1 < ny && diff[i][j] > 0; };
    __int128 area = 0, perimeter = 0;
    for (int i = -1; i + 1 < nx; i++){
        for (int j = -1; j + 1 < ny; j++){
            if (painted(i, j)) area += (__int128)(xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]);
            if (j >= 0 && painted(i, j) != painted(i + 1, j)) perimeter += ys[j + 1] - ys[j];
            if (i >= 0 && painted(i, j) != painted(i, j + 1)) perimeter += xs[i + 1] - xs[i];
        }
    }

    return {area, perimeter};
}

vector<Rectangle> random_rects(int n, long long lo, long long hi){
    vector<Rectangle> rects(n);
    for (auto& r: rects){
        r = {stress::rand_int(lo, hi), stress::rand_int(lo, hi), stress::rand_int(lo, hi), stress::rand_int(lo, hi)};
        if (stress::rand_int(0, 9)){
            if (r.x1 > r.x2) swap(r.x1, r.x2);
            if (r.y1 > r.y2) swap(r.y1, r.y2);
        }
    }
    return rects;
}

void check(const vector<Rectangle>& rects){
    auto [area, perimeter] = paint(rects);
    assert(rectangle_union_area(rects) == area);
    assert(rectangle_union_perimeter(rects) == perimeter);
}

int main(){
    const long long B = 1000000000;

    /// Tiny coordinates force shared edges, seams, nesting, duplicates and holes
    for (long long it = 0; it < stress::scaled(20000); it++){
        int range = stress::rand_int(1, 6);
        check(random_rects(stress::rand_int(0, 8), 0, range));
    }

    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = stress::rand_int(1, 60);
        long long range = it % 3 ? 30 : 1000;
        check(random_rects(n, -range, range));
    }

    /// Full coordinate range, where the area approaches 4e18
    for (long long it = 0; it < stress::scaled(500); it++){
        auto rects = random_rects(stress::rand_int(1, 40), -B, B);
        if (it % 5 == 0) rects.push_back({-B, -B, B, B});
        if (it % 5 == 1) rects.push_back({-B, stress::rand_int(-B, B - 1), B, B});
        check(rects);
    }

    /// Large n: random rectangles inside the cells of a 200 x 200 grid of disjoint, non touching boxes,
    /// so the union area and perimeter are the sums over the distinct boxes
    for (long long it = 0; it < stress::scaled(1); it++){
        const int k = 200;
        const long long cell = 2 * B / k;
        vector<Rectangle> rects;
        long long area = 0, perimeter = 0;
        for (int i = 0; i < k; i++){
            for (int j = 0; j < k; j++){
                long long x0 = -B + i * cell, y0 = -B + j * cell;
                long long x1 = x0 + stress::rand_int(1, cell / 2 - 1), x2 = x0 + stress::rand_int(cell / 2, cell - 1);
                long long y1 = y0 + stress::rand_int(1, cell / 2 - 1), y2 = y0 + stress::rand_int(cell / 2, cell - 1);
                area += (x2 - x1) * (y2 - y1), perimeter += 2 * (x2 - x1 + y2 - y1);
                rects.push_back({x1, y1, x2, y2});
                if (stress::rand_int(0, 3) == 0) rects.push_back({x1, y1, x2, y2});
                if (stress::rand_int(0, 3) == 0) rects.push_back({x1 + 1, y1, x2, y2 - 1});
            }
        }

        shuffle(rects.begin(), rects.end(), stress::rng());
        assert(rectangle_union_area(rects) == area);
        assert(rectangle_union_perimeter(rects) == perimeter);
    }

    return 0;
}
