/***
 *
 * Shared helpers for the stress tests
 *
 * Each stress test includes this header, then the library file with its main renamed:
 *
 *     #include "common.h"
 *     #define main library_main
 *     #include "../code_library/fenwick_tree.cpp"
 *     #undef main
 *
 * STRESS_SEED fixes the random seed, the default is fixed so CI runs are reproducible
 * STRESS_SCALE multiplies iteration counts for longer runs, the default is 1
 * The seed is printed to stderr so any failing run can be replayed
 *
 * stress::expect_assert(library, call): call, run in a child, must abort on an assert in library, a path suffix like "data_structures/mo.cpp"
 *
***/

#include <bits/stdc++.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace std;

namespace stress{
    inline long long env_or(const char* name, long long fallback){
        const char* value = getenv(name);
        return (value && *value) ? atoll(value) : fallback;
    }

    inline unsigned long long seed(){
        static const unsigned long long s = env_or("STRESS_SEED", 20260105);
        return s;
    }

    inline long long scaled(long long iterations){
        static const long long k = max(1LL, env_or("STRESS_SCALE", 1));
        return iterations * k;
    }

    inline mt19937_64& rng(){
        static mt19937_64 generator(seed());
        return generator;
    }

    inline long long rand_int(long long lo, long long hi){
        return uniform_int_distribution<long long>(lo, hi)(rng());
    }

    /// _GLIBCXX_ASSERTIONS also aborts on a bad index, so only the message proves library's own assert fired; the alarm ends a child that loops
    template <typename Call>
    void expect_assert(const string& library, Call call){
        int fd[2];
        assert(pipe(fd) == 0);
        pid_t pid = fork();
        assert(pid >= 0);
        if (pid == 0){
            dup2(fd[1], STDERR_FILENO);
            close(fd[0]), close(fd[1]);
            alarm(10);
            call();
            _exit(0);
        }

        close(fd[1]);
        string err;
        char buf[4096];
        for (ssize_t got; (got = read(fd[0], buf, sizeof buf)) > 0;) err.append(buf, got);
        close(fd[0]);

        int status;
        assert(waitpid(pid, &status, 0) == pid);
        bool fired = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT && err.find("Assertion") != string::npos
            && err.find(library + ":") != string::npos;
        if (!fired) fprintf(stderr, "expected an assert in %s, the child printed:\n%s\n", library.c_str(), err.c_str());
        assert(fired);
    }

    static const bool announced = fprintf(stderr, "STRESS_SEED=%llu\n", seed()) > 0;
}
