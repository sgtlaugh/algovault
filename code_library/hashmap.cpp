/***
 * Simple open addressing hash table with linear probing
 * Every operation is O(1) with very light constant factor
 *
 * Uses a variation SM Hash for choosing the buckets randomly
 * Otherwise it'd be easy to generate counter cases leading to O(n) per operation
 * Need to specify the maximum number of entries to be inserted on construction
 * This is because of better performance as keys are not redistributed
 *
 * For maximal performance, declare once and clear for re-use
 *
 * Default values for keys is 0, like STL map
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename TKey, typename TValue>
class HashMap{
    private:
        const uint64_t seed = chrono::system_clock::now().time_since_epoch().count();

        uint64_t rnd = mt19937_64(seed)();
        const uint64_t salt = (rnd & 1) ? rnd : rnd--;

        int cur_id, _size, hash_mod;

        vector<int> id;
        vector<TKey> keys;
        vector<TValue> values;

        inline int get_pos(TKey x){
            int i = ((x ^ (x >> 33)) * salt) & hash_mod;
            while (id[i] == cur_id && keys[i] != x) i++;
            return i;
        }

    public:
        HashMap(int max_len){
            int buffer = 4 * sqrt(max_len);
            int max_bits = 32 - __builtin_clz(max_len);
            hash_mod = (1 << (max_bits + 1)) - 1;

            _size = 0, cur_id = 1;
            id.resize(hash_mod + buffer, 0);
            keys.resize(hash_mod + buffer, 0);
            values.resize(hash_mod + buffer, 0);
        }

        void set(TKey x, TValue v){
            int i = get_pos(x);
            if (id[i] == cur_id){
                values[i] = v;
                return;
            }

            keys[i] = x, values[i] = v, id[i] = cur_id;
            _size++;
        }

        void add(TKey x, TValue v){
            int i = get_pos(x);
            if (id[i] == cur_id){
                values[i] += v;
                return;
            }

            keys[i] = x, values[i] = v, id[i] = cur_id;
            _size++;
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
