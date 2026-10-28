#include "../common.h"

#define main library_main
#include "../../code_library/graphs/weighted_graph_matching.cpp"
#undef main

#define main hungarian_main
#include "../../code_library/graphs/hungarian_algorithm.cpp"
#undef main

using Matrix = vector<vector<long long>>;

/// Best (size, weight) by always deciding the lowest unmatched vertex, lexicographic when by_size is set
pair<int, long long> brute(const Matrix& w, bool by_size){
    int n = w.size();
    vector<pair<int, long long>> memo(1 << n, {-1, -1});
    function<pair<int, long long>(int)> best = [&](int used){
        if (used == (1 << n) - 1) return pair<int, long long>{0, 0};
        auto& res = memo[used];
        if (res.first != -1) return res;

        int x = __builtin_ctz(~used);
        auto cur = best(used | 1 << x);
        for (int y = x + 1; y < n; y++){
            if (w[x][y] <= 0 || (used >> y & 1)) continue;
            auto sub = best(used | 1 << x | 1 << y);
            pair<int, long long> cand = {sub.first + 1, sub.second + w[x][y]};
            if (by_size ? cand > cur : cand.second > cur.second) cur = cand;
        }

        return res = cur;
    };

    return best(0);
}

/// mate must pair distinct vertices along positive edges and sum to the claimed total
void check_mates(const WeightedBlossom& wb, const Matrix& w, long long total){
    long long sum = 0;
    for (int v = 0; v < wb.n; v++){
        int u = wb.mate[v];
        if (u == -1) continue;
        assert(u >= 0 && u < wb.n && u != v && wb.mate[u] == v && w[v][u] > 0);
        if (u > v) sum += w[v][u];
    }
    assert(sum == total);
}

/// The final duals must certify optimality by LP duality, so a wrong answer on graphs too large to brute force fails
void check_duals(const WeightedBlossom& wb, const Matrix& w){
    int n = wb.n;
    vector<vector<int>> blossoms_of(n);
    for (int b = n; b < 2 * n; b++){
        if (wb.base[b] == -1) continue;
        assert(wb.dual[b] >= 0);

        vector<int> stack = {b}, inside;
        while (!stack.empty()){
            int c = stack.back();
            stack.pop_back();
            if (c < n) inside.push_back(c);
            else stack.insert(stack.end(), wb.children[c].begin(), wb.children[c].end());
        }
        int matched_inside = 0;
        for (int v : inside){
            blossoms_of[v].push_back(b);
            int u = wb.mate[v];
            if (u != -1 && count(inside.begin(), inside.end(), u)) matched_inside++;
        }
        /// A blossom with a positive dual must be full, (size - 1) / 2 of its edges matched
        if (wb.dual[b] > 0) assert(matched_inside == (int)inside.size() - 1);
    }

    for (int v = 0; v < n; v++){
        assert(wb.dual[v] >= 0);
        if (wb.mate[v] == -1) assert(wb.dual[v] == 0);
        for (int u = v + 1; u < n; u++){
            if (w[v][u] <= 0) continue;
            long long s = wb.dual[v] + wb.dual[u] - 2 * w[v][u];
            for (int b : blossoms_of[v]){
                if (count(blossoms_of[u].begin(), blossoms_of[u].end(), b)) s += 2 * wb.dual[b];
            }
            assert(s >= 0);
            if (wb.mate[v] == u) assert(s == 0);
        }
    }
}

Matrix random_graph(int n, int density, long long max_w, bool allow_bad){
    Matrix w(n, vector<long long>(n, 0));
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (stress::rand_int(0, 99) < density) w[u][v] = w[v][u] = stress::rand_int(allow_bad ? -2 : 1, max_w);
        }
    }

    return w;
}

/// Feeds w to the struct in random order, some edges twice with a smaller copy, nonpositive edges included
WeightedBlossom build(const Matrix& w){
    int n = w.size();
    vector<pair<int, int>> order;
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++) if (w[u][v] != 0) order.push_back({u, v});
    }
    shuffle(order.begin(), order.end(), stress::rng());

    WeightedBlossom wb(n);
    for (auto [u, v] : order){
        if (stress::rand_int(0, 1)) swap(u, v);
        if (stress::rand_int(0, 7) == 0) wb.add_edge(u, v, w[u][v] - stress::rand_int(1, 5));
        wb.add_edge(u, v, w[u][v]);
    }

    return wb;
}

void check_against_brute(const Matrix& w){
    WeightedBlossom wb = build(w);
    long long got = wb.maximum_matching();
    assert(got == brute(w, false).second);
    check_mates(wb, w, got);
    check_duals(wb, w);
}

/// Every graph on up to 4 vertices with weights 0 to 3 (0 is no edge) and on 5 with 0 or 1: ties force zero duals
void exhaustive_small(){
    for (int n = 0; n <= 5; n++){
        int slots = n * (n - 1) / 2, levels = n <= 4 ? 4 : 2, total = 1;
        for (int i = 0; i < slots; i++) total *= levels;

        for (int code = 0; code < total; code++){
            Matrix w(n, vector<long long>(n, 0));
            int rest = code;
            for (int u = 0; u < n; u++){
                for (int v = u + 1; v < n; v++) w[u][v] = w[v][u] = rest % levels, rest /= levels;
            }
            check_against_brute(w);
        }
    }
}

/// Random graphs up to 14 vertices, small weights for ties and huge ones at the 1e18 / n bound for overflow
void random_small(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 14), density = stress::rand_int(10, 100);
        long long max_w = it % 4 == 0 ? (long long)1e18 / n : it % 4 == 1 ? 3 : stress::rand_int(1, 1000);
        check_against_brute(random_graph(n, density, max_w, it % 2 == 1));
    }
}

/// The largest-size recipe from the header: weights shifted by B, then the brute force ranks by size first
void largest_size_recipe(){
    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(2, 12);
        long long max_w = stress::rand_int(1, 50);
        Matrix w = random_graph(n, stress::rand_int(20, 80), max_w, false);

        long long shift = n * max_w / 2 + 1;
        Matrix shifted = w;
        for (auto& row : shifted){
            for (auto& x : row) if (x > 0) x += shift;
        }
        WeightedBlossom wb = build(shifted);
        long long got = wb.maximum_matching();
        auto [size, best] = brute(w, true);
        assert(got == size * shift + best);
    }
}

/// Bipartite graphs too large for the brute force, against the Hungarian algorithm with 0 for a missing edge
void bipartite_against_hungarian(){
    for (long long it = 0; it < stress::scaled(40); it++){
        int left = stress::rand_int(1, 60), right = stress::rand_int(1, 60), n = left + right;
        long long max_w = it % 2 ? 5 : 1000000;
        int density = stress::rand_int(5, 100);

        Matrix w(n, vector<long long>(n, 0));
        vector<vector<long long>> cost(left, vector<long long>(right, 0));
        for (int i = 0; i < left; i++){
            for (int j = 0; j < right; j++){
                if (stress::rand_int(0, 99) >= density) continue;
                cost[i][j] = w[i][left + j] = w[left + j][i] = stress::rand_int(1, max_w);
            }
        }

        WeightedBlossom wb = build(w);
        long long got = wb.maximum_matching();
        assert(got == hungarian(cost, false).first);
        check_mates(wb, w, got);
        check_duals(wb, w);
    }
}

/// Larger general graphs, certified by the duals: sparse odd cycles nest blossoms, dense ties expand them
void large_certified(){
    for (long long it = 0; it < stress::scaled(30); it++){
        int n = stress::rand_int(20, it % 5 == 0 ? 250 : 80);
        Matrix w;
        if (it % 3 == 0){
            w.assign(n, vector<long long>(n, 0));
            for (int start = 0; start < n;){
                int len = min(n - start, 2 * (int)stress::rand_int(1, 4) + 1);
                for (int i = 0; i + 1 < len; i++) w[start + i][start + i + 1] = w[start + i + 1][start + i] = 10;
                if (len >= 3) w[start][start + len - 1] = w[start + len - 1][start] = 10;
                if (start > 0){
                    int a = stress::rand_int(0, start - 1), b = start + stress::rand_int(0, len - 1);
                    w[a][b] = w[b][a] = stress::rand_int(1, 20);
                }
                start += len;
            }
        }
        else w = random_graph(n, it % 3 == 1 ? stress::rand_int(2, 10) : 100, it % 2 ? 4 : 1000000, false);

        WeightedBlossom wb = build(w);
        long long got = wb.maximum_matching();
        check_mates(wb, w, got);
        check_duals(wb, w);
    }
}

/// Around 40 vertices with weights up to 20, the regime where inner blossoms expand and relabel their odd path
void medium_certified(){
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(20, 45);
        Matrix w = random_graph(n, stress::rand_int(20, 70), it % 2 ? 5 : 20, false);
        WeightedBlossom wb = build(w);
        long long got = wb.maximum_matching();
        check_mates(wb, w, got);
        check_duals(wb, w);
    }
}

int main(){
    exhaustive_small();
    random_small();
    largest_size_recipe();
    bipartite_against_hungarian();
    medium_certified();
    large_certified();

    /// Edges arrive in two batches, the second call must rebuild from scratch and see the new edges
    for (long long it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(2, 10);
        Matrix w = random_graph(n, 60, 30, false), half = w;
        WeightedBlossom wb(n);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (w[u][v] > 0 && stress::rand_int(0, 1)) wb.add_edge(u, v, w[u][v]);
                else half[u][v] = half[v][u] = 0;
            }
        }
        assert(wb.maximum_matching() == brute(half, false).second);

        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) if (w[u][v] > 0 && half[u][v] == 0) wb.add_edge(v, u, w[u][v]);
        }
        long long got = wb.maximum_matching();
        assert(got == brute(w, false).second);
        check_mates(wb, w, got);
    }

    return 0;
}
