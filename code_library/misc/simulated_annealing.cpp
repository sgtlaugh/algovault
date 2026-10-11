/***
 *
 * Simulated Annealing
 * Random walk over states that accepts worse moves with a probability that cools over time, as a loop driver
 * for in-place moves with an energy delta, and as a generic minimizer over whole states
 *
 * Complexity: O(steps * cost of one move), plus one steady_clock read per step in time mode
 *
 * Annealer::by_steps(t_start, t_end, steps, seed = 17): exactly `steps` moves, same seed gives the same result
 * Annealer::by_time(t_start, t_end, seconds, seed = 17): moves until the time is up, the clock starts at the first running()
 *     T = t_start * (t_end / t_start)^progress, progress = step / steps or elapsed / seconds
 *     requires t_start >= t_end > 0. A good t_start makes a typical uphill step pass about a third of the time,
 *     a good t_end makes the smallest uphill step almost never pass
 * sa.running() cools to the next step and returns false once the steps or the time are used up
 * sa.accept(delta) is true when delta <= 0, otherwise with probability exp(-delta / T); sa.rng is for the caller's moves
 *
 * simulated_annealing(initial, energy, neighbour, annealer) returns the lowest energy state it visited
 *     energy(const State&): signed integer or floating point value to minimize, negate it to maximize, any range of the type is safe
 *     neighbour(const State&, mt19937_64&): a random nearby state, draw all randomness from the given rng
 *     each step copies the state and evaluates it from scratch, use the loop driver when a move has a cheap delta
 *
 * A heuristic: nothing guarantees the optimum, more steps and a move that changes little help
 *
 * Example, TSP with in-place 2-opt moves on a tour `p` of length `cur`:
 *     auto sa = Annealer::by_time(2, 0.01, 1.9);
 *     while (sa.running()){
 *         int i = sa.rng() % n, j = sa.rng() % n;
 *         long long d = two_opt_delta(i, j);
 *         if (sa.accept(d)) apply_two_opt(i, j), cur += d, best = min(best, cur);
 *     }
 *     small states keep the one-liner: simulated_annealing(p, tour_len, two_opt, Annealer::by_steps(100, 0.1, 200000))
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Annealer{
    double t_start, t_end, seconds;
    long long steps, step = 0;
    double temperature = 0;
    mt19937_64 rng;
    chrono::steady_clock::time_point start;

    static Annealer by_steps(double t_start, double t_end, long long steps, unsigned long long seed = 17){
        return {t_start, t_end, 0, steps, 0, 0, mt19937_64(seed), {}};
    }

    static Annealer by_time(double t_start, double t_end, double seconds, unsigned long long seed = 17){
        return {t_start, t_end, seconds, 0, 0, 0, mt19937_64(seed), {}};
    }

    bool running(){
        if (step == 0) start = chrono::steady_clock::now();
        double progress = seconds > 0 ? chrono::duration<double>(chrono::steady_clock::now() - start).count() / seconds
                                      : step >= steps ? 1 : (double)step / steps;
        if (progress >= 1) return false;

        temperature = t_start * pow(t_end / t_start, progress), step++;
        return true;
    }

    bool accept(double delta){
        return delta <= 0 || uniform_real_distribution<double>(0, 1)(rng) < exp(-delta / temperature);
    }
};

/// Integers subtract in unsigned, exact and overflow free since hi > lo: casting each energy to double
/// first would round values past 2^53 together and turn small uphill gaps into free moves
template<typename E>
double uphill_gap(E lo, E hi){
    if constexpr (is_integral_v<E>) return (double)((make_unsigned_t<E>)hi - (make_unsigned_t<E>)lo);
    else return (double)(hi - lo);
}

template<typename State, typename Energy, typename Neighbour>
State simulated_annealing(State cur, Energy energy, Neighbour neighbour, Annealer sa){
    auto cur_e = energy(cur), best_e = cur_e;
    State best = cur;
    while (sa.running()){
        State next = neighbour(cur, sa.rng);
        auto next_e = energy(next);
        if (next_e <= cur_e || sa.accept(uphill_gap(cur_e, next_e))){
            cur = move(next), cur_e = next_e;
            if (cur_e < best_e) best = cur, best_e = cur_e;
        }
    }

    return best;
}

int main(){
    auto cooling = Annealer::by_steps(100, 1, 2);
    assert(cooling.running() && abs(cooling.temperature - 100) < 1e-9);
    assert(cooling.running() && abs(cooling.temperature - 10) < 1e-9);    /// geometric: 100 * (1 / 100)^0.5
    assert(!cooling.running());
    assert(!Annealer::by_steps(1, 1, 0).running() && !Annealer::by_time(1, 1, 0).running());

    auto parabola = [](long long x){ return (x - 37) * (x - 37); };
    auto step = [](long long x, mt19937_64& rng){ return x + (rng() & 1 ? 1 : -1); };
    assert(simulated_annealing(0LL, parabola, step, Annealer::by_steps(100, 1, 20000)) == 37);

    /// The 8 boundary points of a 3 x 3 grid: no two are closer than 1, so the best tour is the perimeter, 8
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
    auto best = simulated_annealing(tour, tour_len, two_opt, Annealer::by_time(2, 0.01, 0.05));  /// 50 ms
    assert(abs(tour_len(best) - 8) < 1e-9);

    /// Loop driver with in-place swaps minimizing sum |p[i] - i|, whose optimum is the identity, 0
    int n = 10;
    vector<int> p = {3, 9, 0, 7, 1, 8, 2, 6, 4, 5};
    long long cur = 0;
    for (int i = 0; i < n; i++) cur += abs(p[i] - i);

    long long lowest = cur;
    auto sa = Annealer::by_steps(5, 0.01, 20000);
    while (sa.running()){
        int i = sa.rng() % n, j = sa.rng() % n;
        long long d = abs(p[j] - i) + abs(p[i] - j) - abs(p[i] - i) - abs(p[j] - j);
        if (sa.accept(d)) swap(p[i], p[j]), cur += d, lowest = min(lowest, cur);
    }
    assert(lowest == 0);
    return 0;
}
