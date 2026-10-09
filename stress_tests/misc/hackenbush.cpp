#include "../common.h"

#define main library_main
#include "../../code_library/misc/hackenbush.cpp"
#undef main

typedef pair<int, int> Edge;

/// Edges of mask still connected to root through edges of mask
int grounded(int n, const vector<Edge>& edges, int mask, int root){
    vector<char> seen(n, 0);
    seen[root] = 1;
    for (bool changed = true; changed;){
        changed = false;
        for (int i = 0; i < (int)edges.size(); i++){
            if (!(mask >> i & 1)) continue;
            auto [u, v] = edges[i];
            if (seen[u] != seen[v]) seen[u] = seen[v] = 1, changed = true;
        }
    }

    int kept = 0;
    for (int i = 0; i < (int)edges.size(); i++){
        if ((mask >> i & 1) && seen[edges[i].first]) kept |= 1 << i;
    }
    return kept;
}

/// Green Hackenbush by the game rules: mex over every edge deletion, memoized on the set of remaining edges
int green_brute(int n, const vector<Edge>& edges, int root){
    int m = edges.size();
    vector<int> memo(1 << m, -1);
    function<int(int)> solve = [&](int mask){
        if (memo[mask] != -1) return memo[mask];
        vector<char> reached(m + 2, 0);
        for (int i = 0; i < m; i++){
            if (mask >> i & 1) reached[min(m + 1, solve(grounded(n, edges, mask ^ (1 << i), root)))] = 1;
        }
        int mex = 0;
        while (reached[mex]) mex++;
        return memo[mask] = mex;
    };
    return solve(grounded(n, edges, (1 << m) - 1, root));
}

/// Fusion and colon principles without lowlinks: bridges by deleting each edge, components by flood fill
int green_naive(int n, const vector<Edge>& edges, int root){
    int m = edges.size();
    auto connected = [&](int skip, int a, int b){
        vector<char> seen(n, 0);
        vector<int> todo = {a};
        seen[a] = 1;
        while (!todo.empty()){
            int u = todo.back();
            todo.pop_back();
            for (int i = 0; i < m; i++){
                if (i == skip) continue;
                auto [x, y] = edges[i];
                int w = x == u ? y : (y == u ? x : -1);
                if (w != -1 && !seen[w]) seen[w] = 1, todo.push_back(w);
            }
        }
        return (bool)seen[b];
    };

    vector<char> bridge(m, 0);
    for (int i = 0; i < m; i++) bridge[i] = edges[i].first != edges[i].second && !connected(i, edges[i].first, edges[i].second);

    vector<int> comp(n, -1);
    for (int s = 0; s < n; s++){
        if (comp[s] != -1) continue;
        comp[s] = s;
        for (bool changed = true; changed;){
            changed = false;
            for (int i = 0; i < m; i++){
                auto [u, v] = edges[i];
                if (!bridge[i] && (comp[u] == s) != (comp[v] == s)) comp[u] = comp[v] = s, changed = true;
            }
        }
    }

    function<int(int, int)> solve = [&](int c, int from){
        int value = 0;
        for (int i = 0; i < m; i++){
            auto [u, v] = edges[i];
            if (!bridge[i]){
                if (comp[u] == c) value ^= 1;
                continue;
            }
            if (i == from || (comp[u] != c && comp[v] != c)) continue;
            value ^= solve(comp[u] == c ? comp[v] : comp[u], i) + 1;
        }
        return value;
    };
    return solve(comp[root], -1);
}

vector<Edge> random_graph(int n, int m){
    vector<Edge> edges;
    for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
    return edges;
}

int green_lib(int n, const vector<Edge>& edges, int root){
    GreenHackenbush g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);
    return g.grundy(root);
}

const int FRAC = 20;
const long long ONE = 1LL << FRAC;

long long floor_div(long long a, long long b){
    return a / b - (a % b != 0 && (a < 0) != (b < 0));
}

/// Simplest dyadic strictly between lo and hi (FRAC fraction bits), the value of a game that is a number
long long simplest(long long lo, long long hi, bool has_lo, bool has_hi){
    if (has_lo && has_hi) assert(lo < hi);
    if ((!has_lo || lo < 0) && (!has_hi || hi > 0)) return 0;
    if (!has_hi || (has_lo && lo >= 0)){
        long long k = (floor_div(lo, ONE) + 1) * ONE;
        if (!has_hi || k < hi) return k;
    }
    else{
        long long k = (-floor_div(-hi, ONE) - 1) * ONE;
        if (!has_lo || k > lo) return k;
    }

    for (long long step = ONE / 2; step > 0; step /= 2){
        long long k = (floor_div(lo, step) + 1) * step;
        if (k < hi) return k;
    }
    assert(false);
    return 0;
}

/// Red-blue Hackenbush by the game rules: the simplest number between the best Left and best Right options
long long red_blue_brute(int n, const vector<Edge>& edges, const vector<int>& blue, int root){
    int m = edges.size();
    vector<long long> memo(1 << m);
    vector<char> done(1 << m, 0);
    function<long long(int)> solve = [&](int mask){
        if (done[mask]) return memo[mask];
        long long lo = 0, hi = 0;
        bool has_lo = false, has_hi = false;
        for (int i = 0; i < m; i++){
            if (!(mask >> i & 1)) continue;
            long long option = solve(grounded(n, edges, mask ^ (1 << i), root));
            if (blue[i]) lo = has_lo ? max(lo, option) : option, has_lo = true;
            else hi = has_hi ? min(hi, option) : option, has_hi = true;
        }
        done[mask] = 1;
        return memo[mask] = simplest(lo, hi, has_lo, has_hi);
    };
    return solve(grounded(n, edges, (1 << m) - 1, root));
}

Dyadic red_blue_lib(int n, const vector<Edge>& edges, const vector<int>& blue, int root){
    RedBlueHackenbush t(n);
    for (int i = 0; i < (int)edges.size(); i++) t.add_edge(edges[i].first, edges[i].second, blue[i]);
    return t.value(root);
}

long long to_fixed(const Dyadic& d){
    long long value = d.whole * ONE;
    for (int b : d.bits){
        assert(b <= FRAC);
        value += ONE >> b;
    }
    return value;
}

/// The edge rule in 64 fraction bits, exact while no subtree holds 64 or more edges
typedef __int128 Wide;
const Wide WIDE_ONE = (Wide)1 << 64;

Wide wide_floor(Wide x){
    Wide q = x / WIDE_ONE;
    return (x % WIDE_ONE != 0 && x < 0) ? q - 1 : q;
}

Wide wide_blue(Wide x){
    if (x >= 0) return x + WIDE_ONE;
    Wide whole = wide_floor(x), frac = x - whole * WIDE_ONE, top = WIDE_ONE + frac;
    int k = -(int)whole;
    assert((top & (((Wide)1 << k) - 1)) == 0);
    return top >> k;
}

Wide red_blue_wide(const vector<vector<pair<int, int>>>& children, int u){
    Wide value = 0;
    for (auto [v, b] : children[u]){
        Wide x = red_blue_wide(children, v);
        value += b ? wide_blue(x) : -wide_blue(-x);
    }
    return value;
}

Wide to_wide(const Dyadic& d){
    Wide value = (Wide)d.whole * WIDE_ONE;
    for (int b : d.bits){
        assert(1 <= b && b <= 64);
        value += WIDE_ONE >> b;
    }
    return value;
}

/// Fixed point with len fraction bits as a bit array, for exact values of long stalks
struct Fixed{
    long long whole = 0;
    vector<char> bits;

    Fixed(int len) : bits(len + 1, 0) {}

    Fixed(const Dyadic& d, int len) : whole(d.whole), bits(len + 1, 0){
        for (int b : d.bits) bits[b] = 1;
    }

    Fixed operator+(const Fixed& o) const{
        Fixed r(bits.size() - 1);
        int carry = 0;
        for (int i = (int)bits.size() - 1; i >= 1; i--){
            int s = bits[i] + o.bits[i] + carry;
            r.bits[i] = s & 1, carry = s >> 1;
        }
        r.whole = whole + o.whole + carry;
        return r;
    }

    Fixed operator-() const{
        Fixed r = *this;
        r.whole = -whole;
        if (count(bits.begin() + 1, bits.end(), 1) == 0) return r;

        r.whole--;
        for (int i = 1; i < (int)bits.size(); i++) r.bits[i] ^= 1;
        int i = bits.size() - 1;
        while (r.bits[i]) r.bits[i--] = 0;
        r.bits[i] = 1;
        return r;
    }

    Dyadic to_dyadic() const{
        Dyadic d;
        d.whole = whole;
        for (int i = 1; i < (int)bits.size(); i++){
            if (bits[i]) d.bits.push_back(i);
        }
        return d;
    }
};

/// Berlekamp's stalk rule: the first run is an integer, then the k-th edge after it is worth +-2^-k by its color
Fixed stalk_value(const vector<int>& blue, int len){
    Fixed positive(len), negative(len);
    if (blue.empty()) return positive;

    int run = 0;
    while (run < (int)blue.size() && blue[run] == blue[0]) run++;
    (blue[0] ? positive : negative).whole = run;
    for (int j = run; j < (int)blue.size(); j++) (blue[j] ? positive : negative).bits[j - run + 1] = 1;
    return positive + (-negative);
}

void check_green(int n, const vector<Edge>& edges, bool brute){
    for (int root = 0; root < n; root++){
        int expected = brute ? green_brute(n, edges, root) : green_naive(n, edges, root);
        assert(green_lib(n, edges, root) == expected);
    }
}

/// Recoloring every edge negates the value, and a tree beside its recolored copy is a second player win
void check_negation(int n, const vector<Edge>& edges, const vector<int>& blue){
    vector<int> flipped(blue.size());
    for (int i = 0; i < (int)blue.size(); i++) flipped[i] = !blue[i];

    Dyadic original = red_blue_lib(n, edges, blue, 0), negated = red_blue_lib(n, edges, flipped, 0);
    int len = 1;
    for (auto& d : {original, negated}){
        if (!d.bits.empty()) len = max(len, d.bits.back());
    }
    assert((-Fixed(original, len)).to_dyadic() == negated);

    vector<Edge> doubled = edges;
    vector<int> doubled_blue = blue;
    auto copy = [&](int v){ return v == 0 ? 0 : v + n - 1; };
    for (int i = 0; i < (int)edges.size(); i++){
        doubled.push_back({copy(edges[i].first), copy(edges[i].second)});
        doubled_blue.push_back(flipped[i]);
    }
    assert(red_blue_lib(2 * n - 1, doubled, doubled_blue, 0) == Dyadic());
}

Dyadic check_stalks(const vector<vector<int>>& stalks){
    int len = 1, n = 1;
    for (auto& s : stalks) len = max(len, (int)s.size()), n += s.size();

    RedBlueHackenbush t(n);
    Fixed expected(len);
    int next = 1;
    for (auto& s : stalks){
        for (int i = 0; i < (int)s.size(); i++, next++) t.add_edge(i == 0 ? 0 : next - 1, next, s[i]);
        expected = expected + stalk_value(s, len);
    }
    Dyadic got = t.value(0);
    assert(got == expected.to_dyadic());
    return got;
}

vector<int> random_stalk(int len){
    vector<int> s(len);
    int stay = stress::rand_int(0, 100);
    for (int i = 0; i < len; i++) s[i] = i == 0 || stress::rand_int(1, 100) > stay ? stress::rand_int(0, 1) : s[i - 1];
    return s;
}

/// Green: game search (exhaustive and random), then fusion without lowlinks on mid-size graphs, closed forms at 2e5
/// Red-blue: game search (also with extra components), the edge rule in 64-bit fixed point up to 63 edges,
/// Berlekamp's stalk rule at 2e5 edges, negation identities on 2e5-vertex random, binary, caterpillar and broom trees
int main(){
    for (int n = 1; n <= 4; n++){
        vector<Edge> kinds;
        for (int u = 0; u < n; u++){
            for (int v = u; v < n; v++) kinds.push_back({u, v});
        }

        int k = kinds.size(), max_m = n <= 3 ? 4 : 3;
        for (int m = 0; m <= max_m; m++){
            int total = 1;
            for (int i = 0; i < m; i++) total *= k;
            for (int code = 0; code < total; code++){
                vector<Edge> edges;
                for (int i = 0, c = code; i < m; i++, c /= k) edges.push_back(kinds[c % k]);
                check_green(n, edges, true);
            }
        }
    }

    for (long long it = 0; it < stress::scaled(600); it++){
        int n = stress::rand_int(1, 7), m = stress::rand_int(0, 10);
        check_green(n, random_graph(n, m), true);
    }

    for (long long it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(1, 32), m = stress::rand_int(0, 2 * n);
        auto edges = random_graph(n, m);
        if (it % 2){
            edges.clear();
            for (int v = 1; v < n; v++) edges.push_back({v, stress::rand_int(max(0, v - 3), v - 1)});
            for (int i = stress::rand_int(0, n / 3); i > 0; i--){
                int u = stress::rand_int(0, n - 1);
                edges.push_back({u, stress::rand_int(max(0, u - 4), min(n - 1, u + 4))});
            }
        }
        check_green(n, edges, false);
    }

    const int BIG = 200000;
    GreenHackenbush long_path(BIG + 1);
    for (int i = 0; i < BIG; i++) long_path.add_edge(i, i + 1);
    assert(long_path.grundy(0) == BIG && long_path.grundy(BIG / 2) == 0 && long_path.grundy(BIG / 2 + 1) == ((BIG / 2 + 1) ^ (BIG / 2 - 1)));

    for (int cycle : {BIG, BIG - 1}){
        GreenHackenbush ring(BIG + 1);
        for (int i = 0; i < BIG - cycle; i++) ring.add_edge(i, i + 1);
        for (int i = 0; i < cycle; i++) ring.add_edge(BIG - cycle + i, BIG - cycle + (i + 1) % cycle);
        assert(ring.grundy(0) == BIG - cycle + (cycle & 1));
    }

    GreenHackenbush lasso(BIG + 1);
    for (int i = 0; i < BIG / 2; i++) lasso.add_edge(i, i + 1);
    for (int i = BIG / 2; i < BIG; i++) lasso.add_edge(i, i + 1);
    lasso.add_edge(BIG, BIG / 2);
    assert(lasso.grundy(0) == BIG / 2 + 1);

    for (int n = 1; n <= 6; n++){
        int shapes = 1;
        for (int v = 1; v < n; v++) shapes *= v;
        for (int shape = 0; shape < shapes; shape++){
            vector<Edge> edges;
            for (int v = 1, c = shape; v < n; c /= v, v++) edges.push_back({c % v, v});
            for (int colors = 0; colors < 1 << (n - 1); colors++){
                vector<int> blue(n - 1);
                for (int i = 0; i < n - 1; i++) blue[i] = colors >> i & 1;
                assert(to_fixed(red_blue_lib(n, edges, blue, 0)) == red_blue_brute(n, edges, blue, 0));
            }
        }
    }

    for (long long it = 0; it < stress::scaled(600); it++){
        int n = stress::rand_int(1, 13);
        vector<int> label(n), blue(n - 1);
        iota(label.begin(), label.end(), 0);
        shuffle(label.begin(), label.end(), stress::rng());
        vector<Edge> edges;
        int span = stress::rand_int(1, n), p = stress::rand_int(0, 100);
        for (int v = 1; v < n; v++) edges.push_back({label[stress::rand_int(max(0, v - span), v - 1)], label[v]});
        for (auto& b : blue) b = stress::rand_int(1, 100) <= p;

        for (int root = 0; root < n; root++){
            Dyadic got = red_blue_lib(n, edges, blue, root);
            long long expected = red_blue_brute(n, edges, blue, root);
            assert(to_fixed(got) == expected && got.sign() == (expected > 0) - (expected < 0));
        }
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int tree = stress::rand_int(1, 9), extra = stress::rand_int(1, 4), n = tree + extra;
        vector<int> label(n);
        iota(label.begin(), label.end(), 0);
        shuffle(label.begin(), label.end(), stress::rng());
        vector<pair<Edge, int>> colored;
        for (int v = 1; v < tree; v++) colored.push_back({{label[stress::rand_int(0, v - 1)], label[v]}, (int)stress::rand_int(0, 1)});
        for (int i = stress::rand_int(0, 4); i > 0; i--){
            Edge e = {label[stress::rand_int(tree, n - 1)], label[stress::rand_int(tree, n - 1)]};
            colored.push_back({e, (int)stress::rand_int(0, 1)});
        }
        shuffle(colored.begin(), colored.end(), stress::rng());

        vector<Edge> edges;
        vector<int> blue;
        for (auto [e, b] : colored) edges.push_back(e), blue.push_back(b);
        for (int i = 0; i < tree; i++) assert(to_fixed(red_blue_lib(n, edges, blue, label[i])) == red_blue_brute(n, edges, blue, label[i]));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 64), span = stress::rand_int(1, n), p = stress::rand_int(0, 100);
        vector<vector<pair<int, int>>> children(n);
        vector<Edge> edges;
        vector<int> blue;
        for (int v = 1; v < n; v++){
            int u = stress::rand_int(max(0, v - span), v - 1), b = stress::rand_int(1, 100) <= p;
            children[u].push_back({v, b});
            edges.push_back({v, u});
            blue.push_back(b);
        }
        assert(to_wide(red_blue_lib(n, edges, blue, 0)) == red_blue_wide(children, 0));
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        vector<vector<int>> stalks(stress::rand_int(1, 12));
        int len = stress::rand_int(1, 400);
        for (auto& s : stalks) s = random_stalk(stress::rand_int(0, len));
        check_stalks(stalks);
    }

    vector<int> alternating(BIG);
    for (int i = 0; i < BIG; i++) alternating[i] = i % 2 == 0;
    assert(check_stalks({alternating}).bits.back() == BIG - 1);
    check_stalks({random_stalk(BIG / 2), random_stalk(BIG / 2)});
    check_stalks(vector<vector<int>>(50, random_stalk(BIG / 50)));

    vector<vector<int>> mixed;
    for (int i = 0; i < 40; i++) mixed.push_back(random_stalk(stress::rand_int(1, BIG / 40)));
    check_stalks(mixed);

    const int HALF = BIG / 2;
    for (int shape = 0; shape < 4; shape++){
        vector<Edge> edges;
        vector<int> blue;
        int p = stress::rand_int(20, 80);
        for (int v = 1; v < HALF; v++){
            int parent = 0;
            if (shape == 0) parent = stress::rand_int(0, v - 1);
            else if (shape == 1) parent = (v - 1) / 2;
            else if (shape == 2) parent = v < HALF / 2 ? v - 1 : stress::rand_int(0, HALF / 2 - 1);
            else parent = v < HALF / 2 ? v - 1 : HALF / 2 - 1;
            edges.push_back({parent, v});
            blue.push_back(stress::rand_int(1, 100) <= p);
        }
        check_negation(HALF, edges, blue);
    }

    return 0;
}
