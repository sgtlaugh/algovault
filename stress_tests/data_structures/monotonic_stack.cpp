#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/monotonic_stack.cpp"
#undef main

template <typename Compare>
void check_nearest(const vector<long long>& a, Compare compare){
    int n = a.size();
    auto left = previous_smaller(a, compare), right = next_smaller(a, compare);
    for (int i = 0; i < n; i++){
        int l = i - 1, r = i + 1;
        while (l >= 0 && !compare(a[l], a[i])) l--;
        while (r < n && !compare(a[r], a[i])) r++;
        assert(left[i] == l && right[i] == r);
    }
}

/// Nearest elements under each of the four orderings by scanning outward, rectangle by trying every range
int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, 40);
        long long range = it % 2 ? 3 : 1000000000;
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(0, range);
        auto left = previous_smaller(a), right = next_smaller(a);

        for (int i = 0; i < n; i++){
            int l = i - 1, r = i + 1;
            while (l >= 0 && !(a[l] < a[i])) l--;
            while (r < n && !(a[r] < a[i])) r++;
            assert(left[i] == l && right[i] == r);
        }
        check_nearest(a, less_equal<long long>());
        check_nearest(a, greater<long long>());
        check_nearest(a, greater_equal<long long>());

        long long best = 0;
        for (int l = 0; l < n; l++){
            long long low = LLONG_MAX;
            for (int r = l; r < n; r++) low = min(low, a[r]), best = max(best, low * (r - l + 1));
        }
        assert(largest_rectangle(a) == best);
    }

    vector<int> up(300000);
    iota(up.begin(), up.end(), 1);
    assert(largest_rectangle(up) == 150000LL * 150001);

    return 0;
}
