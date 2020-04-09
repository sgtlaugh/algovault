#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/***
 *
 * compress values in input array in[] and return them in out[]
 * relative order of values are not necessarily preserved
 * 0 based indexing for arrays and compressed values
 *
***/

void compress(int n, int* in, int* out){ /// 0 based index
    unordered_map <int, int> mp;
    for (int i = 0; i < n; i++) out[i] = mp.emplace(in[i], mp.size()).first->second;
}


/***
 *
 * compress values in input array in[] and return them in out[]
 * relative order of values are preserved
 * 0 based indexing for arrays and compressed values
 *
***/

void compress_sorted(int n, int* in, int* out){
    int i, id = 0;
    vector <int> v;
    unordered_map <int, int> mp;

    for (i = 0; i < n; i++) v.push_back(in[i]);
    sort(v.begin(), v.end());
    if (n) v.push_back(v[n - 1] + 1);

    for (i = 0; i < n; i++){
        if (v[i] != v[i + 1]) mp[v[i]] = id++;
    }
    for (i = 0; i < n; i++) out[i] = mp[in[i]];
}

int main(){
    int ar[] = {2000000000, 1000000000, 2000000000, 1, 10, 5};
    int n = sizeof(ar) / sizeof(ar[0]);
    int i, out[100];

    compress(n, ar, out);
    for (i = 0; i < n; i++) printf("%d ", out[i]);  /// 0 1 0 2 3 4
    puts("");

    compress_sorted(n, ar, out);
    for (i = 0; i < n; i++) printf("%d ", out[i]);  /// 4 3 4 0 2 1
    puts("");

    return 0;
}
