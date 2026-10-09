/***
 *
 * Simulated Annealing
 * Generic minimizer: random walk over states that accepts worse neighbours with a probability that cools over time
 *
 * Complexity: O(steps * (neighbour + energy + state copy)), plus one steady_clock read per step in time mode
 *
 * simulated_annealing(initial, energy, neighbour, schedule, rng) returns the lowest energy state it visited
 *     energy(const State&): signed integer or floating point value to minimize, negate it to maximize, any range of the type is safe
 *     neighbour(const State&, mt19937_64&): a random nearby state, draw all randomness from the given rng
 *     a step from energy e to e' is always taken when e' <= e, otherwise with probability exp((e - e') / T)
 *
 * AnnealingSchedule{t_start, t_end, steps, seconds} cools geometrically, T = t_start * (t_end / t_start)^progress
 *     seconds > 0: time based, progress = elapsed / seconds, runs until the time is up
 *     seconds == 0: iteration based, progress = step / steps, exactly `steps` neighbours, same rng seed gives the same result
 *     requires t_start >= t_end > 0. A good t_start makes a typical uphill step pass about a third of the time,
 *     a good t_end makes the smallest uphill step almost never pass
 *
 * A heuristic: nothing guarantees the optimum, more steps and a neighbour that changes little help
 *
 * Example, TSP with 2-opt moves on a tour `p`:
 *     mt19937_64 rng(17);
 *     auto len = [&](const vector<int>& p){ ... };
 *     auto two_opt = [](vector<int> p, mt19937_64& rng){ int i = rng() % p.size(), j = rng() % p.size(); reverse(...); return p; };
 *     auto best = simulated_annealing(p, len, two_opt, AnnealingSchedule{100, 0.1, 200000, 0}, rng);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct AnnealingSchedule{
    double t_start, t_end;
    long long steps;
    double seconds;

    double temperature(double progress) const{
        return t_start * pow(t_end / t_start, progress);
    }
};

template<typename State, typename Energy, typename Neighbour>
State simulated_annealing(State cur, Energy energy, Neighbour neighbour, const AnnealingSchedule& schedule, mt19937_64& rng){
    auto cur_e = energy(cur), best_e = cur_e;
    State best = cur;
    uniform_real_distribution<double> unit(0.0, 1.0);
    auto start = chrono::steady_clock::now();

    for (long long step = 0; ; step++){
        double progress;
        if (schedule.seconds > 0){
            progress = chrono::duration<double>(chrono::steady_clock::now() - start).count() / schedule.seconds;
            if (progress >= 1) break;
        }
        else{
            if (step >= schedule.steps) break;
            progress = (double)step / schedule.steps;
        }

        State next = neighbour(cur, rng);
        auto next_e = energy(next);
        /// Subtracting in double: in the energy's own type a large signed integer gap overflows
        if (next_e <= cur_e || unit(rng) < exp(((double)cur_e - (double)next_e) / schedule.temperature(progress))){
            cur = move(next), cur_e = next_e;
            if (cur_e < best_e) best = cur, best_e = cur_e;
        }
    }

    return best;
}

int main(){
    AnnealingSchedule cooling{100, 1, 20000, 0};
    assert(abs(cooling.temperature(0) - 100) < 1e-9);
    assert(abs(cooling.temperature(0.5) - 10) < 1e-9);
    assert(abs(cooling.temperature(1) - 1) < 1e-9);

    mt19937_64 rng(17);
    auto parabola = [](long long x){ return (x - 37) * (x - 37); };
    auto step = [](long long x, mt19937_64& rng){ return x + (rng() & 1 ? 1 : -1); };
    assert(simulated_annealing(0LL, parabola, step, cooling, rng) == 37);
    assert(simulated_annealing(-5LL, parabola, step, AnnealingSchedule{100, 1, 0, 0}, rng) == -5);

    /// Energies at both ends of long long: the uphill gap of 2^64 - 1 has acceptance exp(-1.8e19 / 100) = 0
    auto extremes = [](int x){ return x ? LLONG_MIN : LLONG_MAX; };
    long long visits_up = 0;
    auto flip = [&](int x, mt19937_64&){ visits_up += x == 0; return 1 - x; };
    assert(simulated_annealing(1, extremes, flip, cooling, rng) == 1);
    assert(visits_up == 0);

    /// The 8 boundary points of a 3 x 3 grid: no two are closer than 1, so 8 edges cost at least 8, the perimeter
    vector<pair<int, int>> pts = {{0, 0}, {2, 2}, {1, 0}, {0, 2}, {2, 0}, {1, 2}, {0, 1}, {2, 1}};
    auto tour_len = [&](const vector<int>& p){
        double len = 0;
        for (size_t i = 0; i < p.size(); i++){
            auto [x1, y1] = pts[p[i]];
            auto [x2, y2] = pts[p[(i + 1) % p.size()]];
            len += hypot(x1 - x2, y1 - y2);
        }
        return len;
    };
    auto two_opt = [](vector<int> p, mt19937_64& rng){
        int i = rng() % p.size(), j = rng() % p.size();
        if (i > j) swap(i, j);
        reverse(p.begin() + i, p.begin() + j + 1);
        return p;
    };
    vector<int> tour = {0, 1, 2, 3, 4, 5, 6, 7};
    assert(abs(tour_len(simulated_annealing(tour, tour_len, two_opt, AnnealingSchedule{2, 0.01, 20000, 0}, rng)) - 8) < 1e-9);
    assert(abs(tour_len(simulated_annealing(tour, tour_len, two_opt, AnnealingSchedule{2, 0.01, 0, 0.05}, rng)) - 8) < 1e-9);

    return 0;
}
