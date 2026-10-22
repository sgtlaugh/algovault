#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/misc/fast_io.cpp"
#undef main

/// Per process names, concurrent runs (seed sweeps, parallel CI) would otherwise overwrite each other's files
const string IN = "/tmp/algovault_fast_io_in_" + to_string(getpid()), OUT = "/tmp/algovault_fast_io_out_" + to_string(getpid());

/// Points stdin at fresh contents and resets the reader's buffer state
void set_input(const string& contents){
    ofstream(IN, ios::binary) << contents;
    assert(freopen(IN.c_str(), "rb", stdin));
    fio::inptr = fio::buf_len = 0;
}

string read_file(const string& path){
    ifstream f(path, ios::binary);
    return string(istreambuf_iterator<char>(f), istreambuf_iterator<char>());
}

string separator(){
    const string seps[] = {" ", "\t", "\n", "\r\n", "  \n ", "\n\n"};
    return seps[stress::rand_int(0, 5)];
}

template <typename T>
T random_value(){
    int mode = stress::rand_int(0, 4);
    if (mode == 0) return numeric_limits<T>::min();
    if (mode == 1) return numeric_limits<T>::max();
    if (mode == 2) return (T)stress::rand_int(-1000, 1000);
    return (T)stress::rng()();
}

int main(){
    using namespace fio;

    for (long long it = 0; it < stress::scaled(300); it++){
        /// Integers of three widths and string tokens, often longer than the 8 KB buffer
        int count = stress::rand_int(1, it % 5 ? 50 : 4000);
        vector<int> ints(count);
        vector<long long> lls(count);
        vector<unsigned long long> ulls(count);
        vector<string> strs(count);
        string input;
        for (int i = 0; i < count; i++){
            ints[i] = random_value<int>(), lls[i] = random_value<long long>(), ulls[i] = random_value<unsigned long long>();
            strs[i] = string(stress::rand_int(1, 12), 'a');
            for (auto& ch : strs[i]) ch = (char)stress::rand_int(33, 255);  /// printable ASCII and bytes 0x80 - 0xFF
            input += to_string(ints[i]) + separator() + to_string(lls[i]) + separator() + to_string(ulls[i]) + separator() + strs[i];
            if (i + 1 < count || stress::rand_int(0, 1)) input += separator();  /// often no whitespace after the last token
        }

        set_input(input);
        for (int i = 0; i < count; i++){
            int a; long long b; unsigned long long c; string d;
            assert(read(a, b, c, d) == 4);
            assert(a == ints[i] && b == lls[i] && c == ulls[i] && d == strs[i]);
        }
        int extra;
        assert(read(extra) == 0);

        /// Lines, including empty ones, CRLF endings and a last line without a newline
        vector<string> lines(stress::rand_int(0, it % 5 ? 20 : 2000));
        string text;
        bool crlf = stress::rand_int(0, 1);
        for (size_t i = 0; i < lines.size(); i++){
            lines[i] = string(stress::rand_int(0, 2) ? stress::rand_int(0, 30) : 0, 'x');
            for (auto& ch : lines[i]) ch = (char)stress::rand_int(32, 126);
            text += lines[i];
            /// An empty last line needs its newline, "a\n" means the single line "a"
            if (i + 1 < lines.size() || lines[i].empty() || stress::rand_int(0, 1)) text += crlf ? "\r\n" : "\n";
        }

        set_input(text);
        string line;
        for (auto& expected : lines) assert(read_line(line) && line == expected);
        assert(!read_line(line));

        /// Vectors read exactly their current size, nested ones fill a pre-sized grid, chars skip whitespace only
        int n = min(count, 50), rows = stress::rand_int(0, 4), cols = stress::rand_int(0, 4);
        vector<vector<long long>> grid(rows, vector<long long>(cols));
        string chars(stress::rand_int(1, 10), 'a');
        string vec_input = to_string(n);
        for (int i = 0; i < n; i++) vec_input += separator() + to_string(ints[i]);
        for (auto& row : grid){
            for (auto& x : row) x = random_value<long long>(), vec_input += separator() + to_string(x);
        }
        for (auto& ch : chars) ch = (char)stress::rand_int(33, 126);
        vec_input += separator() + chars + separator();

        set_input(vec_input);
        int len;
        assert(read(len) == 1 && len == n);
        vector<int> got_ints(len);
        vector<vector<long long>> got_grid(rows, vector<long long>(cols));
        assert(read(got_ints, got_grid) == 2);
        assert(got_ints == vector<int>(ints.begin(), ints.begin() + n) && got_grid == grid);
        for (char expected : chars){
            char ch;
            assert(read(ch) == 1 && ch == expected);
        }
        char past_end;
        assert(read(past_end) == 0);

        /// Output against to_string, including the minimum values
        assert(freopen(OUT.c_str(), "wb", stdout));
        string expected_out;
        for (int i = 0; i < min(count, 200); i++){
            write(ints[i], lls[i], ulls[i], strs[i]);
            expected_out += to_string(ints[i]) + " " + to_string(lls[i]) + " " + to_string(ulls[i]) + " " + strs[i] + "\n";
        }

        write(vector<int>(ints.begin(), ints.begin() + min(count, 5)));
        for (int i = 0; i < min(count, 5); i++) expected_out += (i ? " " : "") + to_string(ints[i]);
        expected_out += "\n";
        write(chars[0], ints[0], chars.back());
        expected_out += string(1, chars[0]) + " " + to_string(ints[0]) + " " + chars.back() + "\n";
        flush();
        fflush(stdout);
        assert(read_file(OUT) == expected_out);
    }

    /// Floating reads round-trip shortest exact forms, ints and floats mixed; fixed() output matches iostream's fixed
    for (long long it = 0; it < stress::scaled(300); it++){
        int count = stress::rand_int(1, it % 5 ? 50 : 3000);
        vector<double> ds(count);
        vector<long double> lds(count);
        vector<int> mixed(count);
        string input;
        char buf[64];
        for (int i = 0; i < count; i++){
            int mode = stress::rand_int(0, 4);
            double d = mode == 0 ? 0.0 : mode == 1 ? -0.0 : mode == 2 ? (double)stress::rand_int(-1000, 1000) / 8 : ldexp((double)(stress::rng()() >> 11), stress::rand_int(-80, 40)) * (stress::rand_int(0, 1) ? 1 : -1);
            ds[i] = d, lds[i] = (long double)d / 3, mixed[i] = random_value<int>();
            snprintf(buf, sizeof(buf), "%.17g", ds[i]);
            input += string(stress::rand_int(0, 1) && !signbit(ds[i]) ? "+" : "") + buf + separator();
            snprintf(buf, sizeof(buf), "%.21Lg", lds[i]);
            input += string(buf) + separator() + to_string(mixed[i]) + separator();
        }

        set_input(input);
        for (int i = 0; i < count; i++){
            double d;
            long double ld;
            int x;
            assert(read(d, ld, x) == 3);
            assert(d == ds[i] && signbit(d) == signbit(ds[i]));
            assert(ld == lds[i] && x == mixed[i]);
        }
        double past;
        assert(read(past) == 0);

        assert(freopen(OUT.c_str(), "wb", stdout));
        ostringstream expected;
        for (int i = 0; i < min(count, 300); i++){
            int digits = stress::rand_int(0, 18);
            long double big = ldexp((long double)1, stress::rand_int(300, 1000));
            write(fixed(ds[i], digits), fixed(lds[i], digits), fixed(big, digits % 4));
            expected << std::fixed << setprecision(digits) << (long double)ds[i] << ' ' << lds[i] << ' ' << setprecision(digits % 4) << big << '\n';
        }
        flush();
        fflush(stdout);
        assert(read_file(OUT) == expected.str());
    }

    /// Just above the midpoint of two doubles: strtold then a cast rounds to the even neighbor below, the correct double is above
    for (long long it = 0; it < stress::scaled(200); it++){
        double lo = ldexp((double)(((stress::rng()() >> 11) | (1ULL << 52)) & ~1ULL), -52), hi = nextafter(lo, 2.0);
        char buf[128];
        snprintf(buf, sizeof(buf), "%.60Le", ((long double)lo + hi) / 2);
        string token(buf);
        token.insert(token.find('e'), "1");

        set_input(token + "\n");
        double d;
        assert(read(d) == 1);
        assert(d == hi);
    }

    /// A forgotten flush() must not lose output, the child writes the pending buffer when it exits
    assert(freopen(OUT.c_str(), "wb", stdout));
    write("forgotten", 42);
    pid_t pid = fork();
    if (pid == 0) exit(0);
    outptr = 0;
    int status;
    assert(waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0);
    assert(read_file(OUT) == "forgotten 42\n");

    remove(IN.c_str()), remove(OUT.c_str());
    return 0;
}
