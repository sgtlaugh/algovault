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
 * text.erase(v, pos, len): id of a new version equal to version v without the len characters starting at pos
 * text.insert(v, pos, s): id of a new version equal to version v with s inserted before position pos
 * text.at(v, pos): the character at pos in version v
 * text.size(v): length of version v
 * text.substr(v, pos, len): the len characters of version v starting at pos
 *
 * Notes:
 *   - Without the wrapper the same calls work on a crope directly, main() shows it
 *
***/

#include <bits/stdtr1c++.h>
#include <ext/rope>

using namespace std;
using namespace __gnu_cxx;

struct VersionedText{
    vector<crope> versions = {crope()};

    int erase(int v, int pos, int len){
        assert(0 <= pos && 0 <= len && pos + len <= size(v));
        crope next = versions[v];
        next.erase(pos, len);
        versions.push_back(next);
        return versions.size() - 1;
    }

    int insert(int v, int pos, const string& s){
        assert(0 <= pos && pos <= size(v));
        crope next = versions[v];
        next.insert(pos, s.data(), s.size());  /// with the length, so '\0' inside s is kept
        versions.push_back(next);
        return versions.size() - 1;
    }

    char at(int v, int pos) const{
        assert(0 <= pos && pos < size(v));
        return versions[v][pos];
    }

    int size(int v) const{
        assert(0 <= v && v < (int)versions.size());
        return versions[v].size();
    }

    string substr(int v, int pos, int len) const{
        assert(0 <= pos && 0 <= len && pos + len <= size(v));
        crope part = versions[v].substr(pos, len);
        return string(part.begin(), part.end());
    }
};

int main(){
    /// Every edit is a new version, the old ones stay readable
    VersionedText text;
    int a = text.insert(0, 0, "hello");
    int b = text.insert(a, 5, " world");
    int c = text.erase(b, 0, 6);
    int d = text.insert(c, 0, "big ");

    assert(text.substr(a, 0, 5) == "hello");
    assert(text.substr(b, 0, 11) == "hello world");
    assert(text.substr(c, 0, 5) == "world");
    assert(text.substr(d, 0, 9) == "big world");
    assert(text.substr(b, 3, 4) == "lo w");
    assert(text.size(0) == 0);
    assert(text.size(b) == 11);
    assert(text.at(d, 4) == 'w');

    /// The same calls on a crope directly, a copy is an O(1) snapshot
    crope r("abc");
    r.insert(1, "xy");
    assert(string(r.c_str()) == "axybc");

    crope old = r;
    r.erase(0, 2);
    assert(string(r.c_str()) == "ybc");
    assert(string(r.substr(0, 2).c_str()) == "yb");
    assert(string(old.c_str()) == "axybc");

    return 0;
}
