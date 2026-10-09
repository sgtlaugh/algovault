#include "../common.h"

#define main library_main
#include "../../code_library/misc/simulated_annealing.cpp"
#undef main

using Matrix = vector<vector<long long>>;

long long tour_cost(const Matrix& w, const vector<int>& p){
    long long cost = 0;
    for (size_t i = 0; i < p.size(); i++) cost += w[p[i]][p[(i + 1) % p.size()]];
    return cost;
}

/// Every tour with city 0 first, the exact optimum
long long brute_tsp(const Matrix& w){
    int n = w.size();
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);

    long long best = LLONG_MAX;
    do best = min(best, tour_cost(w, p));
    while (next_permutation(p.begin() + 1, p.end()));
    return best;
}

vector<int> two_opt(vector<int> p, mt19937_64& rng){
    int i = rng() % p.size(), j = rng() % p.size();
    if (i > j) swap(i, j);
    reverse(p.begin() + i, p.begin() + j + 1);
    return p;
}

Matrix random_instance(int n, int max_w){
    Matrix w(n, vector<long long>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++) w[i][j] = w[j][i] = stress::rand_int(1, max_w);
    }
    return w;
}

/// TSP against all permutations, plus the bookkeeping the optimum alone cannot see:
/// the returned state is the lowest energy ever evaluated, and the Metropolis rule matches its stationary distribution
int main(){
    for (long long it = 0; it < stress::scaled(200); it++){
        int n = it < 40 ? it % 8 + 1 : stress::rand_int(5, 8);
        int max_w = it % 3 == 0 ? 1000000000 : 100;
        Matrix w = random_instance(n, max_w);
        vector<int> start(n);
        iota(start.begin(), start.end(), 0);
        shuffle(start.begin(), start.end(), stress::rng());

        long long lowest_seen = LLONG_MAX;
        auto energy = [&](const vector<int>& p){
            long long e = tour_cost(w, p);
            lowest_seen = min(lowest_seen, e);
            return e;
        };
        unsigned long long seed = stress::rng()();
        mt19937_64 rng(seed), replay(seed);
        AnnealingSchedule cooling{max_w * 0.5, max_w * 0.001, 20000, 0};
        auto best = simulated_annealing(start, energy, two_opt, cooling, rng);

        vector<int> sorted_best = best;
        sort(sorted_best.begin(), sorted_best.end());
        for (int i = 0; i < n; i++) assert(sorted_best[i] == i);
        assert(tour_cost(w, best) == lowest_seen);
        assert(tour_cost(w, best) == brute_tsp(w));
        assert(simulated_annealing(start, energy, two_opt, cooling, replay) == best);
    }

    /// Two states 0 and 1 with energies 0 and d at a fixed T: Metropolis visits state 1 a fraction e^(-d/T) / (1 + e^(-d/T))
    for (double ratio : {0.5, 1.0, 2.0}){
        long long d = 10, steps = stress::scaled(200000), visits_up = 0, calls = 0;
        double t = d / ratio;
        auto level = [&](int s){ return s * d; };
        auto flip = [&](int s, mt19937_64&){ calls++, visits_up += s; return 1 - s; };
        mt19937_64 rng(stress::rng()());
        simulated_annealing(0, level, flip, AnnealingSchedule{t, t, steps, 0}, rng);

        double expected = exp(-ratio) / (1 + exp(-ratio));
        assert(calls == steps);
        assert(abs((double)visits_up / steps - expected) < 0.01);
    }

    /// The same chain with floating point energies, which take the other branch of uphill_gap
    {
        double d = 0.37;
        long long steps = stress::scaled(200000), visits_up = 0;
        auto level = [&](int s){ return s * d; };
        auto flip = [&](int s, mt19937_64&){ visits_up += s; return 1 - s; };
        mt19937_64 rng(stress::rng()());
        simulated_annealing(0, level, flip, AnnealingSchedule{d, d, steps, 0}, rng);
        assert(abs((double)visits_up / steps - exp(-1) / (1 + exp(-1))) < 0.01);
    }

    /// Cooling from T = 100d to T = d / 1000: the walk sits in state 1 nearly half the time early on and never at the end
    for (int it = 0; it < 3; it++){
        long long d = 10, steps = 100000, early_up = 0, late_up = 0, calls = 0;
        auto level = [&](int s){ return s * d; };
        auto flip = [&](int s, mt19937_64&){
            if (calls < steps / 10) early_up += s;
            if (calls >= steps - steps / 10) late_up += s;
            calls++;
            return 1 - s;
        };
        mt19937_64 rng(stress::rng()());
        simulated_annealing(0, level, flip, AnnealingSchedule{100.0 * d, d / 1000.0, steps, 0}, rng);

        assert(early_up > steps / 10 * 0.4);
        assert(late_up == 0);
    }

    /// Near zero temperature only downhill or level steps pass, so the walk's energy never rises
    for (long long it = 0; it < stress::scaled(50); it++){
        Matrix w = random_instance(8, 100);
        long long last = LLONG_MAX;
        auto watched = [&](const vector<int>& p, mt19937_64& rng){
            long long e = tour_cost(w, p);
            assert(e <= last);
            last = e;
            return two_opt(p, rng);
        };
        vector<int> start(8);
        iota(start.begin(), start.end(), 0);
        mt19937_64 rng(stress::rng()());
        simulated_annealing(start, [&](const vector<int>& p){ return tour_cost(w, p); }, watched, AnnealingSchedule{1e-12, 1e-12, 2000, 0}, rng);
    }

    /// Energies spanning all of long long: every uphill gap is at least 2^63, far too steep to take at T <= 1e6
    for (long long it = 0; it < stress::scaled(20); it++){
        long long low = stress::rand_int(LLONG_MIN, LLONG_MIN / 2), high = stress::rand_int(LLONG_MAX / 2, LLONG_MAX);
        auto extremes = [&](int s){ return s ? low : high; };
        long long visits_up = 0;
        auto flip = [&](int s, mt19937_64&){ visits_up += s == 0; return 1 - s; };
        mt19937_64 rng(stress::rng()());
        assert(simulated_annealing(1, extremes, flip, AnnealingSchedule{1e6, 1, 2000, 0}, rng) == 1);
        assert(visits_up == 0);
    }

    /// The same two-state chain shifted far past 2^53, where a double cannot hold the energies apart:
    /// the stationary fraction depends only on the gap, so it must not move with the offset
    for (long long it = 0; it < stress::scaled(30); it++){
        long long d = stress::rand_int(1, 10), steps = stress::scaled(100000), visits_up = 0;
        long long off = it % 2 ? stress::rand_int(1LL << 53, LLONG_MAX - d) : stress::rand_int(LLONG_MIN, -(1LL << 53));
        auto level = [&](int s){ return off + s * d; };
        auto flip = [&](int s, mt19937_64&){ visits_up += s; return 1 - s; };
        mt19937_64 rng(stress::rng()());
        simulated_annealing(0, level, flip, AnnealingSchedule{(double)d, (double)d, steps, 0}, rng);
        assert(abs((double)visits_up / steps - exp(-1) / (1 + exp(-1))) < 0.01);
    }

    /// Time based: stops at the deadline, one cheap step past it at most, and still solves a small instance
    for (int it = 0; it < 3; it++){
        Matrix w = random_instance(6, 100);
        vector<int> start(6);
        iota(start.begin(), start.end(), 0);
        mt19937_64 rng(stress::rng()());
        auto t0 = chrono::steady_clock::now();
        auto best = simulated_annealing(start, [&](const vector<int>& p){ return tour_cost(w, p); }, two_opt, AnnealingSchedule{50, 0.1, 0, 0.2}, rng);
        double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        assert(secs >= 0.2 && secs < 0.3);
        assert(tour_cost(w, best) == brute_tsp(w));
    }

    return 0;
}
