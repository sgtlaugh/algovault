/***
 *
 * Rope
 * GNU rope as a persistent text: every edit makes a new version and all old versions stay readable
 *
 * Complexity:
 *   - insert: O(log n + |s|)
 *   - erase, at: O(log n)
 *   - substr: O(log n + len)
 *   - copying a crope: O(1), versions share their unchanged pieces
 *
 * VersionedText text: version 0 is the empty text
 * text.insert(v, pos, s): id of a new version equal to version v with s inserted before position pos
 * text.erase(v, pos, len): id of a new version equal to version v without the len characters starting at pos
 * text.substr(v, pos, len): the len characters of version v starting at pos
 * text.size(v): length of version v
 * text.at(v, pos): the character at pos in version v
 *
 * Without the wrapper the same calls work on a crope directly, main() shows it
 *
***/

#include <bits/stdtr1c++.h>
#include <ext/rope>

using namespace std;
using namespace __gnu_cxx;

struct VersionedText{
    vector<crope> versions = {crope()};

    int insert(int v, int pos, const string& s){
        assert(0 <= pos && pos <= size(v));
        crope next = versions[v];
        next.insert(pos, s.data(), s.size());  /// with the length, so '\0' inside s is kept
        versions.push_back(next);
        return versions.size() - 1;
    }

    int erase(int v, int pos, int len){
        assert(0 <= pos && 0 <= len && pos + len <= size(v));
        crope next = versions[v];
        next.erase(pos, len);
        versions.push_back(next);
        return versions.size() - 1;
    }

    string substr(int v, int pos, int len) const{
        assert(0 <= pos && 0 <= len && pos + len <= size(v));
        crope part = versions[v].substr(pos, len);
        return string(part.begin(), part.end());
    }

    int size(int v) const{
        assert(0 <= v && v < (int)versions.size());
        return versions[v].size();
    }

    char at(int v, int pos) const{
        assert(0 <= pos && pos < size(v));
        return versions[v][pos];
    }
};

int main(){
    int n;

    scanf("%d", &n);

    VersionedText text;
    int now = 0, d = 0;
    while (n--){
        int op, p, c, v;
        scanf("%d", &op);
        if (op == 1){
            char s[128];
            scanf("%d %s", &p, s);
            now = text.insert(now, p - d, s);
        }
        else if (op == 2){
            scanf("%d %d", &p, &c);
            now = text.erase(now, p - d - 1, c - d);
        }
        else{
            scanf("%d %d %d", &v, &p, &c);
            string s = text.substr(v - d, p - d - 1, c - d);
            d += count(s.begin(), s.end(), 'c');
            puts(s.c_str());
        }
    }

    return 0;
}
