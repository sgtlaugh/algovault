#include <stdio.h>
#include <bits/stdtr1c++.h>


#define MAX_LOG      20
#define MAX_LETTERS  26


using namespace std;


struct AhoCorasick{
    int id, edge[256];
    vector <long long> counter;
    vector <string> dictionary;
    vector <int> leaf, fail;

    vector <vector<int>> dp;
    vector <map<char, int>> trie;

    inline int node(){
        leaf.push_back(0);
        counter.push_back(0);
		dp.push_back(vector<int>(MAX_LETTERS, 0));
		trie.push_back(map<char, int>());
		return id++;
    }

    inline int size(){
        return dictionary.size();
    }

    void clear(){
        trie.clear();
        dictionary.clear();
        dp.clear(), fail.clear(), leaf.clear(), counter.clear();

        id = 0, node();
        for (int i = 'a'; i <= 'z'; i++) edge[i] = i - 'a'; /// change here if different character set
    }

    AhoCorasick(){
        clear();
    }

    inline void insert(const char* str){
        int j, x, cur = 0;
        for (j = 0; str[j] != 0; j++){
            x = edge[(int)str[j]];
            if (!trie[cur].count(x)){
                int next_node = node();
                trie[cur][x] = next_node;
            }
            cur = trie[cur][x];
        }

        leaf[cur]++;
        dictionary.push_back(str);
    }

	/// call build once after insertion is done
    inline void build(){
        vector <pair<int, pair<int, int> > > Q;
        fail.resize(id, 0);
        Q.push_back({0, {0, 0}});

        for (int i = 0; i < id; i++){
            for (int j = 0; j < MAX_LETTERS; j++){
                dp[i][j] = i;
            }
        }

        for(int i = 0; i < (int)Q.size(); i++){
            int u = Q[i].first;
            int p = Q[i].second.first;
            char c = Q[i].second.second;
            for(auto& it: trie[u]) Q.push_back({it.second, {u, it.first}});

            if (u){
                int f = fail[p];
                while (f && !trie[f].count(c)) f = fail[f];
                if(!trie[f].count(c) || trie[f][c] == u) fail[u] = 0;
                else fail[u] = trie[f][c];
                counter[u] = leaf[u] + counter[fail[u]];

                for (int j = 0; j < MAX_LETTERS; j++){
                    if (u && !trie[u].count(j)) dp[u][j] = dp[fail[u]][j];
                }
            }
        }
    }

    inline int next(int cur, char ch){
        int x = edge[(int)ch];
        cur = dp[cur][x];
        if (trie[cur].count(x)) cur = trie[cur][x];
        return cur;
    }

	/// total number of occurrences of all words from dictionary in str
    long long count(const char* str){
        int cur = 0;
        long long res = 0;

        for (int j = 0; str[j] && id > 1; j++){
            cur = next(cur, str[j]);
            res += counter[cur];
        }

        return res;
    }
};


/// dynamic aho corasick in N log N
/// no need to call build explicitly before count, unlike aho corasick

struct DynamicAhoCorasick{
    AhoCorasick ar[MAX_LOG];

    DynamicAhoCorasick(){
        for (int i = 0; i < MAX_LOG; i++) ar[i] = AhoCorasick();
    }

    inline void insert(const char* str){
        int i, k = 0;
        for (k = 0; k < MAX_LOG && ar[k].size(); k++){}

        ar[k].insert(str);
        for (i = 0; i < k; i++){
            for (auto s: ar[i].dictionary){
                ar[k].insert(s.c_str());
            }
            ar[i].clear();
        }
        ar[k].build();
    }

    long long count(const char* str){
        long long res = 0;
        for (int i = 0; i < MAX_LOG; i++) res += ar[i].count(str);
        return res;
    }
};


int main(){
	auto ac = DynamicAhoCorasick();

	ac.insert("hello");
	ac.insert("world");

	printf("%lld\n", ac.count("lol"));            									   /// 0
	ac.insert("lol");
	printf("%lld\n", ac.count("lol"));            									   /// 1
	ac.insert("lol");
	printf("%lld\n", ac.count("lol"));            									   /// 2

	ac.insert("abracadabra");
	ac.insert("abaababbaba");
	ac.insert("aaba");

	printf("%lld\n", ac.count("helloworldlol"));  									   /// 4
	printf("%lld\n", ac.count("abaababbaba"));   									   /// 2
	printf("%lld\n", ac.count("aba"));            									   /// 0
	printf("%lld\n", ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba"));  /// 7

	ac.insert("hello");
	ac.insert("world");
	ac.insert("lol");
	ac.insert("lol");

	ac.insert("a");
	ac.insert("baa");

	printf("%lld\n", ac.count("helloworldlol"));  									   /// 8
	printf("%lld\n", ac.count("abaababbaba"));   									   /// 9
	printf("%lld\n", ac.count("aba"));            									   /// 2
	printf("%lld\n", ac.count("baababaababbbabaabaabaababbabababbbbaaabababababba"));  /// 38

	return 0;
}
