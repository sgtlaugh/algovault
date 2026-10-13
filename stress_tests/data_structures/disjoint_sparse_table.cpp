#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/data_structures/disjoint_sparse_table.cpp"
#undef main

/// An out of range index also aborts under _GLIBCXX_ASSERTIONS but is undefined behaviour in a judge build, hence the message check
void check_rejects(int n, int l, int r){
    int fd[2];
    assert(pipe(fd) == 0);
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        dup2(fd[1], STDERR_FILENO);
        close(fd[0]), close(fd[1]);
        DisjointSparseTable<int> table(vector<int>(n, 1));
        table.query(l, r);
        _exit(0);
    }

    close(fd[1]);
    string err;
    char buf[4096];
    for (ssize_t got; (got = read(fd[0], buf, sizeof buf)) > 0;) err.append(buf, got);
    close(fd[0]);

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
    assert(err.find("code_library/data_structures/disjoint_sparse_table.cpp:") != string::npos && err.find("Assertion") != string::npos);
}

/// The default plus makes strings concatenate: associative but not commutative, operand order must be exact
/// A product mod p lambda checks a custom op against a direct scan
int main(){
    check_rejects(9, 3, 1);
    check_rejects(9, -1, 2);
    check_rejects(9, 0, 9);
    check_rejects(0, 0, 0);

    const long long MOD = 998244353;
    auto mul = [&](long long x, long long y){ return x * y % MOD; };

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 40 : 600);
        if (it % 7 == 0) n = 1 << stress::rand_int(0, 9);  /// block boundaries line up with n

        vector<string> words(n);
        vector<long long> nums(n), factors(n);
        for (int i = 0; i < n; i++){
            words[i] = string(1, 'a' + i % 26) + (i % 3 ? "" : to_string(i));
            nums[i] = stress::rand_int(-1000000000000LL, 1000000000000LL);
            factors[i] = stress::rand_int(0, it % 2 ? 3 : MOD - 1);
        }

        DisjointSparseTable<string> concat(words);
        DisjointSparseTable<long long> sum(nums);
        DisjointSparseTable prod(factors, mul);
        for (int q = 0; q < 50; q++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, min(n - 1, l + (q % 2 ? 5 : n)));
            string expected;
            for (int i = l; i <= r; i++) expected += words[i];
            assert(concat.query(l, r) == expected);
            assert(sum.query(l, r) == accumulate(nums.begin() + l, nums.begin() + r + 1, 0LL));

            long long product = 1;
            for (int i = l; i <= r; i++) product = product * factors[i] % MOD;
            assert(prod.query(l, r) == product);
        }
    }

    return 0;
}
