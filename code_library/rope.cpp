/***
 *
 * Rope
 * GNU rope as a persistent text: every edit makes a new version and all old versions stay readable
 *
 * Complexity: O(log n) per insert / erase / character access, O(log n + k) for a substring of length k
 *             copying a rope is O(1), versions share their unchanged pieces
 *
 * VersionedText text;       version 0 is the empty text
 * text.insert(v, pos, s):   new version = version v with s inserted before position pos, returns its index
 * text.erase(v, pos, len):  new version = version v without the len characters starting at pos
 * text.substr(v, pos, len), text.size(v), text.at(v, pos)
 *
 * Without the wrapper the same calls work on a crope directly:
 *     crope r = "abc"; r.insert(1, "xy"); r.erase(0, 2); r.substr(0, 2); crope old = r;  (old is an O(1) snapshot)
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
        return versions[v].size();
    }

    char at(int v, int pos) const{
        assert(0 <= pos && pos < size(v));
        return versions[v][pos];
    }
};

int main(){
    VersionedText text;
    int a = text.insert(0, 0, "hello");
    int b = text.insert(a, 5, " world");
    int c = text.erase(b, 0, 6);
    int d = text.insert(c, 0, "big ");

    assert(text.substr(a, 0, 5) == "hello");
    assert(text.substr(b, 0, 11) == "hello world");
    assert(text.substr(c, 0, 5) == "world");
    assert(text.substr(d, 0, 9) == "big world");
    assert(text.size(0) == 0 && text.size(b) == 11 && text.at(d, 4) == 'w');
    assert(text.substr(b, 3, 4) == "lo w");
    return 0;
}
