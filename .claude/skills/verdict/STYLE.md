# Practice solution style

The code after the template copies reads like the user wrote it in a contest: short, plain, no comments. The excerpts at the end are the user's own accepted code; match them more than any rule here. Their `scanf` and `printf` predate the `fio` rule: read them as `read` and `write`.
For that code these rules replace CONTRIBUTING's `## API and code`, `## Style` and `main()` rules, which govern templates; CI's whitespace and brace checks still apply.

## Rules

- No comments at all after the template copies
- Input and output only through the `fio` copy, with `using namespace fio;` first in `main`: `read(n, m)`, `write(res)`, `write(fixed(x, 6))` for a float; never `scanf`, `printf`, `puts`, `cin` or `cout`, which would reorder the stream
- A line `write` cannot shape, like a field width or `Case 3: ` before the answer, is built with `snprintf(buf, sizeof(buf), "Case %d: %lld", ++T, res);` and written as `write(string(buf));`
- Declare the loop and scratch variables at the top of the function, `int t, i, j, res;`, and loop with `for (i = 0; ...)`; a variable used in one small spot may be declared there
- Never declare a variable that is not used
- Fixed-size input lives in globals sized by a macro, `#define MAX 100010` and `int n, ar[MAX];`; use a `vector` when a template takes one
- Small free functions above `main` for the steps, named for their job: `gen()`, `solve(k)`, `check()`, `is_equal()`
- `int main(){` reads input, calls the steps, prints, and ends with `return 0;`
- Test cases run as `read(t); while (t--){` or `while (read(n)){`, with a `Case %d: ` prefix only when the problem asks for it
- On UVa, POJ and HDU a test count is read until EOF, `while (read(t)){ while (t--){`, since their input files often repeat the block (UVa 1718 failed reading one)
- One-statement bodies go on the same line without braces: `for (i = 0; i < n; i++) read(ar[i]);`, `if (x > n) return false;`
- Related assignments share a line with commas: `low = 0, high = n - 1;`
- Parentheses make the order obvious: `if ((i + 1) < n)`, `res += (x * y)`; halving and parity use `>> 1` and `& 1`
- Binary search is `while ((low + 1) < high)`, then check `low` and `high`
- Blank lines separate the steps: input, setup, the main loop, output
- Braces as the repo: `){`, `else{`, an `else` goes on the line after the `}`
- `clr(ar)` and `dbg(x)` may be defined and used, other macros only when they shorten something used several times; macros go right after the template copies:
  ```cpp
  #define clr(ar) memset(ar, 0, sizeof(ar))
  #define dbg(x) cerr << #x << " = " << x << endl
  ```
- Nothing the excerpts never use: `auto [a, b]`, lambdas for code that runs once, `const auto&` everywhere, `std::` prefixes, `#define int long long`

## Excerpts

Codeforces 1547F, 2021, after its sparse table:

```cpp
int n, ar[MAX], dp[MAX];

bool is_equal(){
    for (int i = 1; i < n; i++){
        if (ar[i] != ar[0]) return false;
    }
    return true;
}

int main(){
    int t, i, res, low, high, mid;

    scanf("%d", &t);
    while (t--){
        scanf("%d", &n);
        for (i = 0; i < n; i++) scanf("%d", &ar[i]);

        if (is_equal()){
            puts("0");
            continue;
        }

        vector<uint32_t> v;
        for (i = 0; i < n; i++) v.push_back(ar[i]);
        for (i = 0; i < n; i++) v.push_back(ar[i]);

        auto T = SparseTable<uint32_t>(v);

        int g = T.query(0, n - 1);
        for (i = 0; i < n; i++){
            dp[i] = INT_MAX;

            low = i + 1, high = i + n - 1;
            while ((low + 1) < high){
                mid = (low + high) >> 1;
                if (T.query(i, mid) == g) high = mid;
                else low = mid;
            }

            res = -1;
            if (T.query(i, low) == g) res = low;
            else if (T.query(i, high) == g) res = high;

            if (res != -1) dp[i] = min(dp[i], (res - i));
        }

        res = 0;
        for (i = 0; i < n; i++) res = max(res, dp[i]);
        printf("%d\n", res);
    }
    return 0;
}
```

Codeforces 920F, 2021, the steps around its segment tree:

```cpp
inline void update(int l, int r){
    if (l > r) return;
    node v = query(1, 1, n, l, r);
    if (ar[v.pos_max] <= 2) return;

    update(1, 1, n, v.pos_max);
    update(l, v.pos_max - 1);
    update(v.pos_max + 1, r);
}

int main(){
    int q, i, j, l, r, flag;

    for (i = 1; i < MAXV; i++){
        for (j = i; j < MAXV; j += i){
            divcnt[j]++;
        }
    }

    scanf("%d %d", &n, &q);
    for (i = 1; i <= n; i++) scanf("%d", &ar[i]);
    build(1, 1, n);

    while (q--){
        scanf("%d %d %d", &flag, &l, &r);
        if (flag == 1) update(l, r);
        else{
            node v = query(1, 1, n, l, r);
            printf("%lld\n", v.sum);
        }
    }
    return 0;
}
```

UVa 12706 Zero-Knowledge Protocol:

```cpp
unordered_map<int, int> mp;
int n, m, r, C1[20010], C2[20010], ar[20010], ar2[20010];

bool check(){
    for (int i = 0; i <= r; i++){
        if (C1[i] != C2[i]) return false;
    }
    return true;
}

int main(){
    int t, i;

    scanf("%d", &t);
    while (t--){
        scanf("%d %d", &n, &m);
        for (i = 0; i < n; i++) scanf("%d", &ar[i]);
        for (i = 0; i < m; i++) scanf("%d", &ar2[i]);

        r = 0;
        mp.clear();
        bool flag = true;
        for (i = 0; i < n; i++){
            int x = ar[i];
            if (!mp[x]) mp[x] = ++r;
            ar[i] = mp[x];
        }

        for (i = 0; i < m; i++){
            int x = ar2[i];
            if (!mp[x]) flag = false;
            ar2[i] = mp[x];
        }

        long long res = 0;
        if (flag && m <= n){
            clr(C1), clr(C2);
            for (i = 0; i < m; i++) C1[ar[i]]++, C2[ar2[i]]++;

            for (i = 0; (i + m) <= n; i++){
                if (check()) res = res + ((i + 1) * (i + 1));
                C1[ar[i]]--;
                C1[ar[i + m]]++;
            }
        }

        printf("%lld\n", res);
    }
    return 0;
}
```
