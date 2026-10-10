/***
 * Open addressing hash table with linear probing for integral keys
 *
 * Complexity:
 *   - O(max_len) to construct, the table holds the smallest power of two >= 2 * max_len slots, so the load factor stays at most 1/2
 *   - Expected O(1) per set/add/get/contains/erase with a light constant factor
 *   - O(1) to clear
 *
 * Buckets come from splitmix64 of the key plus a salt chosen at runtime
 * Otherwise it'd be easy to generate counter cases leading to O(n) per operation
 * See https://codeforces.com/blog/entry/62393
 *
 * Need to specify the maximum number of entries to be inserted on construction, the table never grows
 * For maximal performance, declare once and clear for re-use
 *
 * Default value for missing keys is TValue(), like STL map
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename TKey, typename TValue>
class HashMap{
    static_assert(is_integral<TKey>::value, "HashMap needs an integral key type");

    private:
        /// home caches the bucket of key so erase never re-hashes, it fits in the padding next to id for 8 byte keys
        struct Slot{
            int id, home;
            TKey key;
            TValue value;
        };

        /// static so construction stays cheap for many small maps and copy assignment is not deleted
        static inline const uint64_t salt = mt19937_64(chrono::steady_clock::now().time_since_epoch().count())();

        int cur_id = 1, _size = 0, max_len, mask;
        vector<Slot> slots;

        static uint64_t splitmix64(uint64_t x){
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        inline int home(TKey x) const{
            return splitmix64((uint64_t)x + salt) & mask;
        }

        /// Returns the slot holding x, or the empty slot where x would be inserted
        inline int get_pos(TKey x, int i) const{
            while (slots[i].id == cur_id && slots[i].key != x) i = (i + 1) & mask;
            return i;
        }

        inline void insert_at(int i, int h, TKey x, TValue v){
            assert(_size < max_len);
            slots[i] = {cur_id, h, x, v};
            _size++;
        }

    public:
        HashMap(int max_len) : max_len(max(max_len, 1)){
            int bits = 1;
            while ((1 << bits) < 2 * this->max_len) bits++;  /// load factor stays at most 1/2
            mask = (1 << bits) - 1;
            slots.assign(mask + 1, Slot{0, 0, TKey(), TValue()});
        }

        void set(TKey x, TValue v){
            int h = home(x), i = get_pos(x, h);
            if (slots[i].id == cur_id) slots[i].value = v;
            else insert_at(i, h, x, v);
        }

        void add(TKey x, TValue v){
            int h = home(x), i = get_pos(x, h);
            if (slots[i].id == cur_id) slots[i].value += v;
            else insert_at(i, h, x, v);
        }

        TValue get(TKey x) const{
            const Slot& s = slots[get_pos(x, home(x))];
            return s.id == cur_id ? s.value : TValue();
        }

        bool contains(TKey x) const{
            return slots[get_pos(x, home(x))].id == cur_id;
        }

        void erase(TKey x){
            int i = get_pos(x, home(x));
            if (slots[i].id != cur_id) return;
            _size--;

            /// Backward shift instead of tombstones: an entry moves into the hole at i unless its home lies cyclically in (i, j]
            for (int j = (i + 1) & mask; slots[j].id == cur_id; j = (j + 1) & mask){
                if (((j - slots[j].home) & mask) >= ((j - i) & mask)){
                    slots[i] = slots[j];
                    i = j;
                }
            }
            slots[i].id = 0;
        }

        void clear(){
            _size = 0;
            cur_id++;
        }

        int size() const{
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
    assert(hashmap.size() == 0 && !hashmap.contains(INT_MIN) && !hashmap.contains(INT_MAX));

    /// Keys that differ only in their high bits, the classic way to break weak hashes
    auto spread = HashMap<long long, int>(100000);
    for (int i = 0; i < 100000; i++) spread.set((long long)i << 40, i);
    for (int i = 0; i < 100000; i++) assert(spread.get((long long)i << 40) == i);

    /// A tiny table with random inserts and erases keeps long, wrapping probe chains, compare against std::map
    mt19937 rng(20020523);
    auto small = HashMap<int, int>(40);
    map<int, int> expected;
    for (int op = 0; op < 200000; op++){
        int key = (int)(rng() % 120) - 60, value = rng() % 1000;
        if (rng() % 3 == 0) small.erase(key), expected.erase(key);
        else if (expected.count(key) || (int)expected.size() < 40) small.add(key, value), expected[key] += value;
        assert(small.size() == (int)expected.size());
        assert(small.get(key) == (expected.count(key) ? expected[key] : 0));
    }
    for (auto& kv : expected) assert(small.get(kv.first) == kv.second);

    return 0;
}
