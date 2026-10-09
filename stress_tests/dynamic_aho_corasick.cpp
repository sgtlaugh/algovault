#include "common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../code_library/dynamic_aho_corasick.cpp"
#undef main

string random_word(int len, int alphabet){
    string s(len, 'a');
    for (auto& c : s) c = 'a' + stress::rand_int(0, alphabet - 1);
    return s;
}

long long brute(const vector<string>& patterns, const string& text){
    long long res = 0;
    for (const auto& p : patterns){
        for (size_t i = 0; i + p.size() <= text.size(); i++) res += text.compare(i, p.size(), p) == 0;
    }
    return res;
}

/// With every slot full, the next insert must hit the capacity assert rather than touch ar[MAX_LOG]
/// Its message is checked because the overflow can also abort through a library assert on garbage memory
void check_full_insert_aborts(){
    int fd[2];
    assert(pipe(fd) == 0);
    pid_t pid = fork();
    if (pid == 0){
        dup2(fd[1], STDERR_FILENO);
        DynamicAhoCorasick ac;
        for (int i = 0; i < MAX_LOG; i++) ac.ar[i].insert("overbooked");
        ac.insert("onemore");
        _exit(0);
    }
    close(fd[1]);
    string err;
    char buf[4096];
    for (ssize_t len; (len = read(fd[0], buf, sizeof(buf))) > 0;) err.append(buf, len);
    close(fd[0]);

    int status;
    waitpid(pid, &status, 0);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
    assert(err.find("k < MAX_LOG") != string::npos);
}

/// Queried after every insert so that each merge of the binary decomposition is exercised
int main(){
    for (long long it = 0; it < stress::scaled(200); it++){
        int alphabet = stress::rand_int(1, it % 5 ? 3 : 26);
        vector<string> patterns;
        DynamicAhoCorasick ac;

        for (int ins = stress::rand_int(0, it % 10 ? 20 : 150); ins; ins--){
            patterns.push_back(random_word(stress::rand_int(1, 6), alphabet));
            ac.insert(patterns.back());

            string text = random_word(stress::rand_int(0, 60), alphabet);
            for (auto& c : text) if (stress::rand_int(0, 15) == 0) c = "A{ .\xe9\x80"[stress::rand_int(0, 5)];  /// outside the alphabet, including bytes above 0x7f
            assert(ac.count(text) == brute(patterns, text));
        }
    }

    check_full_insert_aborts();
    return 0;
}
