/***
 *
 * Hash Map
 * Hash table for integral keys that crafted hack tests cannot slow down
 *
 * Complexity:
 *   - construction: O(max_len)
 *   - set, add, get, contains, erase: expected O(1)
 *   - clear: O(1)
 *
 * HashMap<TKey, TValue> mp(max_len): empty map for at most max_len distinct keys, the table never grows
 * mp.set(x, v): store v under x
 * mp.add(x, v): add v to the value under x, a missing x starts from TValue()
 * mp.get(x): the value under x, TValue() if x is missing, never inserts
 * mp.contains(x): true if x is in the map
 * mp.erase(x): remove x, nothing happens if x is missing
 * mp.clear(): remove every key
 * mp.size(): number of keys in the map
 *
 * Inserting a key beyond max_len fails an assert
 * Keys are hashed with splitmix64 plus a salt chosen at runtime, so fixed tests cannot force collisions, see https://codeforces.com/blog/entry/62393
 * Linear probing over the smallest power of two >= 2 * max_len slots, so the load factor stays at most 1/2
 * For speed, declare one map and clear() it between uses
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename TKey, typename TValue>
struct HashMap{
    static_assert(is_integral<TKey>::value, "HashMap needs an integral key type");

    /// home caches the bucket of key so erase never re-hashes, it fits in the padding next to id for 8 byte keys
    struct Slot{
        int id, home;
        TKey key;
        TValue value;
    };

    /// static so construction stays cheap for many small maps and copy assignment is not deleted
    static inline const uint64_t salt = mt19937_64(chrono::steady_clock::now().time_since_epoch().count())();

    /// A slot is live only while its id equals cur_id, so clear() empties every slot at once by bumping cur_id
    int cur_id = 1, _size = 0, max_len, mask;
    vector<Slot> slots;

    static uint64_t splitmix64(uint64_t x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    int home(TKey x) const{
        return splitmix64((uint64_t)x + salt) & mask;
    }

    /// Returns the slot holding x, or the empty slot where x would be inserted
    int get_pos(TKey x, int i) const{
        while (slots[i].id == cur_id && slots[i].key != x) i = (i + 1) & mask;
        return i;
    }

    void insert_at(int i, int h, TKey x, TValue v){
        assert(_size < max_len);
        slots[i] = {cur_id, h, x, v};
        _size++;
    }

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
    int q;
    scanf("%d", &q);

    HashMap<long long, long long> mp(q);
    while (q--){
        int t;
        long long k;
        scanf("%d %lld", &t, &k);
        if (t == 0){
            long long v;
            scanf("%lld", &v);
            mp.set(k, v);
        }
        else printf("%lld\n", mp.get(k));
    }
    return 0;
}
