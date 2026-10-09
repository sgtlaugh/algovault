#include "../common.h"

#define main library_main
#include "../../code_library/strings/de_bruijn_sequence.cpp"
#undef main

/// Length k^n, symbols in [0, k), every cyclic length n window distinct
bool is_de_bruijn(const vector<int>& s, int k, int n){
    long long len = 1;
    for (int i = 0; i < n; i++) len *= k;
    if ((long long)s.size() != len) return false;
    for (int x : s){
        if (x < 0 || x >= k) return false;
    }

    vector<bool> seen(len);
    for (long long i = 0; i < len; i++){
        long long code = 0;
        for (int j = 0; j < n; j++) code = code * k + s[(i + j) % len];
        if (seen[code]) return false;
        seen[code] = true;
    }

    return true;
}

/// The original recursive FKM (ShahjalalShohag code-library), with its arrays made local
vector<int> original(int k, int n){
    if (k == 1) return {0};
    vector<int> ans, aux(k * n + 1);
    function<void(int, int)> db = [&](int t, int p){
        if (t > n){
            if (n % p == 0){
                for (int i = 1; i <= p; i++) ans.push_back(aux[i]);
            }
            return;
        }
        aux[t] = aux[t - p];
        db(t + 1, p);
        for (int i = aux[t - p] + 1; i < k; i++){
            aux[t] = i;
            db(t + 1, t);
        }
    };

    db(1, 1);
    return ans;
}

/// First sequence in lexicographic order with all windows distinct, by backtracking
bool smallest(vector<int>& s, set<vector<int>>& seen, int k, int n, int len){
    if ((int)s.size() == len){
        set<vector<int>> all = seen;
        for (int i = max(0, len - n + 1); i < len; i++){
            vector<int> w;
            for (int j = 0; j < n; j++) w.push_back(s[(i + j) % len]);
            if (!all.insert(w).second) return false;
        }
        return true;
    }

    for (int c = 0; c < k; c++){
        s.push_back(c);
        vector<int> w;
        if ((int)s.size() >= n) w.assign(s.end() - n, s.end());
        if (w.empty() || seen.insert(w).second){
            if (smallest(s, seen, k, n, len)) return true;
            if (!w.empty()) seen.erase(w);
        }
        s.pop_back();
    }

    return false;
}

int main(){
    for (int k = 1; k <= 6; k++){
        long long len = 1;
        for (int n = 1; len * k <= 1000000; n++){
            len *= k;
            auto s = de_bruijn(k, n);
            assert(is_de_bruijn(s, k, n));
            assert(s == original(k, n));
            if (k == 1 && n > 8) break;
        }
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int k = stress::rand_int(1, 40), n = stress::rand_int(1, 3);
        if (k > 12 && n == 3) n = 2;
        auto s = de_bruijn(k, n);
        assert(is_de_bruijn(s, k, n));
        assert(s == original(k, n));
    }

    for (auto [k, n] : vector<pair<int, int>>{{1, 1}, {1, 3}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {3, 1}, {3, 2}, {4, 2}, {5, 1}, {2, 5}}){
        int len = 1;
        for (int i = 0; i < n; i++) len *= k;
        vector<int> s;
        set<vector<int>> seen;
        assert(smallest(s, seen, k, n, len));
        assert(de_bruijn(k, n) == s);
    }

    assert(is_de_bruijn(de_bruijn(2, 20), 2, 20));
    assert(is_de_bruijn(de_bruijn(3, 12), 3, 12));
    assert(is_de_bruijn(de_bruijn(1000, 2), 1000, 2));
    assert(is_de_bruijn(de_bruijn(1000000, 1), 1000000, 1));

    return 0;
}
