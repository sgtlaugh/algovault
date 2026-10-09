#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/xor_basis.cpp"
#undef main

/// Distinct XORs of all subsets, built by doubling the set whenever a value falls outside it
vector<ull> brute_span(const vector<ull>& values){
    vector<ull> span = {0};
    set<ull> seen = {0};
    for (ull x : values){
        if (seen.count(x)) continue;
        int len = span.size();
        for (int i = 0; i < len; i++){
            span.push_back(span[i] ^ x);
            seen.insert(span[i] ^ x);
        }
    }

    sort(span.begin(), span.end());
    return span;
}

ull random_bits(int bits){
    return bits ? stress::rng()() >> (64 - bits) : 0;
}

ull random_probe(const vector<ull>& span){
    int kind = stress::rand_int(0, 4);
    ull x = span[stress::rand_int(0, span.size() - 1)];
    if (kind == 0) return x;
    if (kind == 1) return x + 1;
    if (kind == 2) return x - 1;
    if (kind == 3) return random_bits(stress::rand_int(0, 64));
    return x ^ (1ULL << stress::rand_int(0, 63));
}

/// Values drawn from the span of a few generators, so two sets share a nontrivial intersection
vector<ull> combos(const vector<ull>& gens, int n){
    vector<ull> values(n);
    for (auto& x : values){
        for (ull g : gens) if (stress::rand_int(0, 1)) x ^= g;
    }
    return values;
}

vector<ull> random_gens(int count, int bits){
    vector<ull> gens(count);
    for (auto& g : gens) g = random_bits(stress::rand_int(0, bits));
    if (count && stress::rand_int(0, 3) == 0) gens[0] = ULLONG_MAX;
    return gens;
}

void check_against_span(const XorBasis& b, const vector<ull>& span, int probes){
    assert(1ULL << b.rank == span.size());
    for (ull k = 0; k < span.size(); k++) assert(b.kth(k) == span[k]);
    for (ull x : span) assert(b.contains(x));

    for (int it = 0; it < probes; it++){
        ull x = random_probe(span);
        bool in = binary_search(span.begin(), span.end(), x);
        ull best = 0, worst = ULLONG_MAX;
        for (ull s : span) best = max(best, x ^ s), worst = min(worst, x ^ s);

        assert(b.contains(x) == in);
        assert(b.max_xor(x) == best && b.min_xor(x) == worst);
        assert(b.count_less(x) == (ull)(lower_bound(span.begin(), span.end(), x) - span.begin()));
    }
}

XorBasis build(const vector<ull>& values){
    XorBasis b;
    vector<ull> prefix;
    for (ull x : values){
        vector<ull> span = brute_span(prefix);
        assert(b.insert(x) == !binary_search(span.begin(), span.end(), x));
        prefix.push_back(x);
    }
    return b;
}

void check_intersection(const vector<ull>& p, const vector<ull>& q){
    XorBasis bp = build(p), bq = build(q);
    vector<ull> sp = brute_span(p), sq = brute_span(q), common;
    set_intersection(sp.begin(), sp.end(), sq.begin(), sq.end(), back_inserter(common));

    XorBasis both = bp.intersect(bq);
    vector<ull> gens;
    for (ull x : both.basis) if (x) gens.push_back(x);
    assert(brute_span(gens) == common);
    assert(1ULL << both.rank == common.size());
}

int main(){
    for (int mask = 0; mask < (1 << 8); mask++){
        vector<ull> values;
        for (int x = 0; x < 8; x++) if (mask >> x & 1) values.push_back(x);
        shuffle(values.begin(), values.end(), stress::rng());
        XorBasis b = build(values);
        check_against_span(b, brute_span(values), 20);

        int n = values.size();
        map<ull, int> hits;
        for (int s = 0; s < (1 << n); s++){
            ull x = 0;
            for (int i = 0; i < n; i++) if (s >> i & 1) x ^= values[i];
            hits[x]++;
        }
        assert(hits.size() == 1ULL << b.rank);
        for (auto [x, c] : hits) assert(c == 1 << (n - b.rank));
    }

    for (int p = 0; p < (1 << 7); p++){
        for (int q = 0; q < (1 << 7); q++){
            vector<ull> vp, vq;
            for (int x = 1; x < 8; x++){
                if (p >> (x - 1) & 1) vp.push_back(x);
                if (q >> (x - 1) & 1) vq.push_back(x);
            }
            check_intersection(vp, vq);
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 12), bits = stress::rand_int(0, it % 3 ? 6 : 64);
        vector<ull> values(n);
        for (auto& x : values) x = random_bits(stress::rand_int(0, bits));
        if (n && stress::rand_int(0, 3) == 0) values[stress::rand_int(0, n - 1)] = ULLONG_MAX;
        check_against_span(build(values), brute_span(values), 30);
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int bits = stress::rand_int(1, it % 2 ? 64 : 8);
        vector<ull> shared = random_gens(stress::rand_int(0, 4), bits);
        vector<ull> gp = shared, gq = shared;
        for (ull g : random_gens(stress::rand_int(0, 4), bits)) gp.push_back(g);
        for (ull g : random_gens(stress::rand_int(0, 4), bits)) gq.push_back(g);
        check_intersection(combos(gp, stress::rand_int(0, 10)), combos(gq, stress::rand_int(0, 10)));
    }

    for (long long it = 0; it < stress::scaled(30); it++){
        int n = stress::rand_int(500, 3000);
        vector<ull> values = combos(random_gens(stress::rand_int(0, 12), 64), n);
        XorBasis b;
        for (ull x : values) b.insert(x);
        check_against_span(b, brute_span(values), 200);
    }

    for (long long it = 0; it < stress::scaled(500); it++){
        int count = stress::rand_int(0, 70);
        vector<ull> gens = random_gens(count, stress::rand_int(1, 3) == 1 ? 40 : 64);

        /// Random bit lengths almost never span all 64 bits, so full-width generators are what reach rank 64
        if (count >= 64) for (auto& g : gens) g = stress::rng()();
        vector<ull> values = combos(gens, stress::rand_int(count >= 64 ? 64 : 0, 150));
        XorBasis b;
        ull plain[64] = {};
        int rank = 0;
        for (ull x : values){
            ull y = x;
            for (int i = 63; i >= 0 && y; i--){
                if (!(y >> i & 1)) continue;
                if (!plain[i]){
                    plain[i] = y, rank++;
                    break;
                }
                y ^= plain[i];
            }
            assert(b.insert(x) == (y != 0));
        }
        assert(b.rank == rank);

        ull last = rank == 64 ? ULLONG_MAX : (1ULL << rank) - 1;
        assert(b.kth(0) == 0 && b.count_less(b.kth(last)) == last);
        for (int q = 0; q < 50; q++){
            ull x = stress::rng()(), best = x, worst = x;
            for (int i = 63; i >= 0; i--) best = max(best, best ^ plain[i]), worst = min(worst, worst ^ plain[i]);
            assert(b.max_xor(x) == best && b.min_xor(x) == worst && b.contains(x) == (worst == 0));

            ull k = stress::rng()() & last, s = b.kth(k);
            assert(b.contains(s) && b.count_less(s) == k);
            assert(k == last || b.kth(k + 1) > s);
            ull c = b.count_less(x);
            assert(c == 0 || b.kth(c - 1) < x);
            assert((c > last && rank < 64) || (c <= last && b.kth(c) >= x));
        }
    }

    /// Rank 63 without bit 63: count_less returns early at the top bit with 2^63 values below it
    XorBasis low;
    for (int i = 0; i < 63; i++) low.insert(1ULL << i);
    assert(low.count_less(1ULL << 63) == 1ULL << 63 && low.count_less(ULLONG_MAX) == 1ULL << 63);

    return 0;
}
