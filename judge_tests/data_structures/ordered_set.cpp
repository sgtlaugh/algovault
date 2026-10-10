// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/ordered_set
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/ordered_set.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

int read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    int x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

/// The problem's set never holds duplicates, so every OrderedMultiset answer must match the OrderedSet one
int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    OrderedSet<int> s;
    OrderedMultiset<int> m;
    for (int i = 0; i < n; i++){
        int x = read_int();
        s.insert(x);
        m.insert(x);
    }

    string out;
    while (q--){
        int t = read_int(), x = read_int(), res = 0;
        bool print = t >= 2;
        if (t == 0){
            if (s.find(x) == s.end()) s.insert(x), m.insert(x);
            assert(m.count(x) == 1);
        }
        else if (t == 1){
            bool erased = s.erase(x);
            assert(m.erase(x) == erased && m.count(x) == 0);
        }
        else if (t == 2){
            res = (int)s.size() < x ? -1 : *s.find_by_order(x - 1);
            assert(m.size() == (int)s.size() && (res == -1 || m.kth(x - 1) == res));
        }
        else if (t == 3){
            res = s.order_of_key(x + 1);
            assert(m.count_less(x + 1) == res);
        }
        else if (t == 4){
            auto it = s.upper_bound(x);
            res = it == s.begin() ? -1 : *prev(it);
        }
        else{
            auto it = s.lower_bound(x);
            res = it == s.end() ? -1 : *it;
        }

        if (print){
            out += to_string(res);
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
