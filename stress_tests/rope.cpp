#include "common.h"

#define main library_main
#include "../code_library/rope.cpp"
#undef main

/// Every version against a std::string copy kept alongside it
int main(){
    for (long long it = 0; it < stress::scaled(300); it++){
        VersionedText text;
        vector<string> ref = {""};
        int ops = stress::rand_int(1, 200);
        for (int op = 0; op < ops; op++){
            int v = stress::rand_int(0, ref.size() - 1), n = ref[v].size();
            if (n == 0 || stress::rand_int(0, 2)){
                int pos = stress::rand_int(0, n), len = stress::rand_int(0, 6);
                string s;
                for (int i = 0; i < len; i++) s += stress::rand_int(0, 9) ? char('a' + stress::rand_int(0, 25)) : '\0';
                assert(text.insert(v, pos, s) == (int)ref.size());
                ref.push_back(ref[v].substr(0, pos) + s + ref[v].substr(pos));
            }
            else{
                int pos = stress::rand_int(0, n - 1), len = stress::rand_int(0, n - pos);
                assert(text.erase(v, pos, len) == (int)ref.size());
                ref.push_back(ref[v].substr(0, pos) + ref[v].substr(pos + len));
            }

            int w = stress::rand_int(0, ref.size() - 1), m = ref[w].size();
            assert(text.size(w) == m);
            if (m){
                int pos = stress::rand_int(0, m - 1), len = stress::rand_int(0, m - pos);
                assert(text.substr(w, pos, len) == ref[w].substr(pos, len));
                assert(text.at(w, pos) == ref[w][pos]);
            }
        }
        for (int v = 0; v < (int)ref.size(); v++) assert(text.substr(v, 0, ref[v].size()) == ref[v]);
    }

    /// Many versions of a long text stay cheap because copies share structure
    VersionedText big;
    int v = big.insert(0, 0, string(200000, 'x'));
    for (int i = 0; i < 20000; i++) v = big.insert(v, stress::rand_int(0, big.size(v)), "ab");
    assert(big.size(v) == 240000 && big.size(1) == 200000);
    return 0;
}
