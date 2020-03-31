#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Hashmap{
    int t, sz, hmod;
    vector <int> id;
    vector <long long> key, val;

    inline int randint(int a, int b){
		unsigned int x = rand(), y = rand();
		unsigned int num = (x << 15) ^ y;
		return num % (b - a + 1) + a;
	}

    inline unsigned long long smhash(unsigned long long h){
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccd;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53;
        h ^= h >> 33;
        return h;
    }

    inline int next_prime(int n){
        for (int i = n; ;i++){
            for (int j = 2; ;j = j + 1 + (j > 2)){
                if ((j * j) > i) return i;
                if ((i % j) == 0) break;
            }
        }
        return -1;
    }

    void clear(){t++;}

    inline int pos(unsigned long long x){
        int i = smhash(x) % hmod;
        while (id[i] == t && (unsigned long long)key[i] != x) i++;
        return i;
    }

    inline void erase(long long x){
        int i = pos(x);
        if (id[i] == t) key[i] = 0, val[i] = 0, id[i] = 0, sz--;
    }

    inline long long find(long long x){
        int i = pos(x);
        return (id[i] != t) ? -1 : val[i];
    }

    inline bool contains(long long x){
        int i = pos(x);
        return (id[i] == t);
    }

    inline void add(long long x, long long v){
        int i = pos(x);
        (id[i] == t) ? (val[i] += v) : (key[i] = x, val[i] = v, id[i] = t, sz++);
    }

    inline int size(){
        return sz;
    }

    Hashmap(int m=32){
        unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
		srand(seed);

        m = (m << 1) - randint(1, m);
        hmod = next_prime(max(100, m));

        sz = 0, t = 1;
        id.resize(hmod + 0x1FF, 0);
        key.resize(hmod + 0x1FF, 0), val.resize(hmod + 0x1FF, 0);
    }
};

int main(){
	Hashmap H = Hashmap();

	cout << H.find(5) << endl; 		/// -1
	H.add(5, 32);
	cout << H.find(5) << endl; 		/// 32

	H.add(5, 15);
	cout << H.find(5) << endl; 		/// 47
	H.add(5, -47);
	cout << H.find(5) << endl; 		/// 0

	cout << H.contains(5) << endl;  /// true
	H.erase(5);
	cout << H.contains(5) << endl;  /// false

	H.add(13, 3);
	H.add(666, 15);
	cout << H.find(666) << endl; 	/// 15
	cout << H.size() << endl;  		/// 2

	return 0;
}
