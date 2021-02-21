/***
 * Simple open addressing hash table with linear probing
 * Every operation is O(1) with very light constant factor
 *
 * Uses SM Hash for choosing the buckets randomly
 * Need to specify the maximum number of entries to be inserted on construction
 * This is because of better performance as keys are not redistributed
 *
 * Default values for keys is 0, like STL map
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename TKey, typename TValue>
class HashMap{
    int cur_id, _size, hash_mod;

    vector <int> id;
    vector <TKey> keys;
    vector <TValue> values;

    bool is_prime(int n){
        if (n <= 1) return false;
        for (int i = 2; i * i <= n; i++){
            if (n % i == 0) return false;
        }
        return true;
    }

    int next_prime(int n){
        while (!is_prime(++n)) {}
        return n;
    }

    inline unsigned long long sm_hash(unsigned long long h){
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccd;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53;
        h ^= h >> 33;

        return h;
    }

    inline int get_pos(TKey x){
        int i = sm_hash(x) % hash_mod;
        while (id[i] == cur_id && keys[i] != x) i++;
        return i;
    }

    void update(TKey x, TValue v, bool replace=false){
        int i = get_pos(x);

        if (id[i] == cur_id){
            if (replace) values[i] = 0;
            values[i] += v;
            return;
        }

        keys[i] = x, values[i] = v, id[i] = cur_id;
        _size++;
    }

    public:
        HashMap(int max_len){
            max_len = max(max_len, 16);
            mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

            int m = max_len / 2, s = sqrt(max_len);
            int salt = rng() % m;
            hash_mod = next_prime(2.5 * max_len - salt);

            _size = 0, cur_id = 1;
            id.resize(hash_mod + s, 0);
            keys.resize(hash_mod + s, 0);
            values.resize(hash_mod + s, 0);
        }

        void set(TKey x, TValue v){
            update(x, v, true);
        }

        void add(TKey x, TValue v){
            update(x, v, false);
        }

        TValue get(TKey x){
            int i = get_pos(x);
            return id[i] == cur_id ? values[i] : 0;
        }

        void erase(TKey x){
            int i = get_pos(x);

            if (id[i] == cur_id){
                keys[i] = 0, values[i] = 0, id[i] = 0;
                _size--;
            }
        }

        bool contains(TKey x){
            int i = get_pos(x);
            return id[i] == cur_id;
        }

        void clear(){
            _size = 0;
            cur_id++;
        }

        int size(){
            return _size;
        }

};

int main(){
    auto hashmap = HashMap<int, int>(100);

    assert(hashmap.get(5) == 0);

    hashmap.add(5, 10);
    assert(hashmap.get(5) == 10);

    hashmap.add(5, 5);
    assert(hashmap.get(5) == 15);

    hashmap.set(5, 32);
    assert(hashmap.get(5) == 32);

    assert(hashmap.size() == 1);
    hashmap.erase(5);
    assert(hashmap.size() == 0);

    hashmap.erase(13);
    assert(hashmap.size() == 0);

    hashmap.set(INT_MIN, INT_MIN);
    assert(hashmap.contains(INT_MIN));

    hashmap.add(INT_MIN, 1);
    assert(hashmap.get(INT_MIN) == (INT_MIN + 1));

    hashmap.add(INT_MAX, 1);
    assert(hashmap.contains(INT_MAX));

    assert(hashmap.size() == 2);
    hashmap.clear();
    assert(hashmap.size() == 0);

    return 0;
}
