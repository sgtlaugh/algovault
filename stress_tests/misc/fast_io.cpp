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

string int128_to_string(__int128 v){
    unsigned __int128 u = v < 0 ? -(unsigned __int128)v : v;
    string s;
    do s += char('0' + u % 10), u /= 10; while (u);
    if (v < 0) s += '-';
    return string(s.rbegin(), s.rend());
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
            /// Often no whitespace after the last token
            if (i + 1 < count || stress::rand_int(0, 1)) input += separator();
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

        /// A vector argument ends its line even when more arguments follow
        write(vector<int>(ints.begin(), ints.begin() + min(count, 5)), lls[0]);
        for (int i = 0; i < min(count, 5); i++) expected_out += (i ? " " : "") + to_string(ints[i]);
        expected_out += "\n" + to_string(lls[0]) + "\n";
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

    /// Just above a midpoint of two doubles: strtold then a cast rounds to the even double below, the right one is up
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

    /// Double reads are bit-identical to strtod on the token's first 127 characters, around every fast path limit
    vector<string> fixed_tokens = {"9007199254740991", "9007199254740992", "9007199254740993", "1234567890123456789",
        "12345678901234567890", "18446744073709551621", "0000000000000000001.5", "00000000000000000001.5",
        "1e22", "1e23", "1e-22", "1e-23", "4.9e-324", "2.2250738585072014e-308", "1.7976931348623157e308", "1e309",
        "-1e-400", "inf", "-nan", "0x1p-3", "1e", "1e+", "2E-", ".", "-", "+.", "5.", ".5", "-.5e1", "e5", "1.5abc",
        "1e0000000000000000005", "-0", "+0.0e-0", string(200, '7'), "0." + string(150, '0') + "1e150"};
    for (long long it = 0; it < stress::scaled(300); it++){
        vector<string> tokens = fixed_tokens;
        for (int i = 0; i < 300; i++){
            const string signs[] = {"", "-", "+"};
            string token = signs[stress::rand_int(0, 2)];
            if (i < 60){
                /// Mantissas 2^53 - 1 to 2^53 + 2 times 10^-23 to 10^23, a double rounding past 2^53 misses by one ulp
                token += to_string((1LL << 53) + stress::rand_int(-1, 2));
                int e = stress::rand_int(-23, 23), point = stress::rand_int(0, min(16, max(0, -e)));
                token.insert(token.size() - point, ".");
                if (e + point) token += "e" + to_string(e + point);
            }
            else {
                for (int k = stress::rand_int(0, 25); k > 0; k--) token += char('0' + stress::rand_int(0, 9));
                if (stress::rand_int(0, 2)) token += ".";
                for (int k = stress::rand_int(0, 25); k > 0; k--) token += char('0' + stress::rand_int(0, 9));
                if (stress::rand_int(0, 1)){
                    token += (stress::rand_int(0, 1) ? "e" : "E") + signs[stress::rand_int(0, 2)];
                    token += to_string(stress::rand_int(0, 2) ? stress::rand_int(0, 30) : stress::rand_int(0, 400));
                }
                if (token.empty()) token = ".";
            }
            tokens.push_back(token);
        }
        shuffle(tokens.begin(), tokens.end(), stress::rng());

        string input;
        for (auto& token : tokens) input += token + separator();
        set_input(input);
        for (auto& token : tokens){
            double d, expected = strtod(token.substr(0, 127).c_str(), nullptr);
            assert(read(d) == 1);
            assert(memcmp(&d, &expected, sizeof(d)) == 0);
        }
    }

    /// __int128 extremes use the integer reader and writer, float parses at its own width, short as an integer
    for (long long it = 0; it < stress::scaled(300); it++){
        int count = stress::rand_int(1, 200);
        vector<__int128> bigs(count);
        vector<float> fs(count);
        vector<short> shorts(count);
        string input, expected_out;
        char buf[64];
        for (int i = 0; i < count; i++){
            int mode = stress::rand_int(0, 3);
            __int128 top = ~((unsigned __int128)1 << 127);
            bigs[i] = mode == 0 ? -top - 1 : mode == 1 ? top : (__int128)(((unsigned __int128)stress::rng()() << 64) | stress::rng()());
            fs[i] = ldexpf((float)stress::rand_int(-(1 << 23), 1 << 23), stress::rand_int(-60, 40));
            shorts[i] = random_value<short>();
            snprintf(buf, sizeof(buf), "%.9g", fs[i]);
            input += int128_to_string(bigs[i]) + separator() + buf + separator() + to_string(shorts[i]) + separator();
            expected_out += int128_to_string(bigs[i]) + "\n";
        }

        set_input(input);
        assert(freopen(OUT.c_str(), "wb", stdout));
        for (int i = 0; i < count; i++){
            __int128 big;
            float f;
            short sh;
            assert(read(big, f, sh) == 3);
            assert(big == bigs[i] && f == fs[i] && sh == shorts[i]);
            write(big);
        }
        flush();
        fflush(stdout);
        assert(read_file(OUT) == expected_out);
    }

    /// Every split of a number across a buffer refill, at the end of input too, with and without a trailing newline
    for (int pad = BUF_SIZE - 25; pad <= BUF_SIZE + 5; pad++){
        int small = random_value<int>();
        long long big = random_value<long long>();
        string tail = stress::rand_int(0, 1) ? "\n" : "";
        set_input(string(pad, ' ') + to_string(small) + " " + to_string(big) + tail);

        int a = 0;
        long long b = 0;
        assert(read(a, b) == 2);
        assert(a == small && b == big);
        assert(read(a) == 0);
    }

    /// Interactive, no flush() calls: a refill by fread or before our query leaves deadlocks, and the alarm fails it
    int to_judge[2], to_us[2];
    assert(!pipe(to_judge) && !pipe(to_us));
    long long rounds = stress::scaled(200);
    pid_t judge = fork();
    if (judge == 0){
        close(to_judge[1]), close(to_us[0]);
        FILE *in = fdopen(to_judge[0], "r"), *out = fdopen(to_us[1], "w");
        for (long long r = 0; r < rounds; r++){
            long long x;
            int k;
            if (fscanf(in, " ? %lld %d", &x, &k) != 2) _exit(1);

            string reply = "round " + to_string(r) + "\n";
            for (int j = 0; j < k; j++) reply += to_string(x + j) + (j + 1 < k ? " " : "\n");

            /// Two writes with a pause between, so a read can return a token cut in half
            size_t cut = stress::rand_int(0, reply.size());
            fwrite(reply.data(), 1, cut, out), fflush(out);
            if (r % 10 == 0) usleep(200);
            fwrite(reply.data() + cut, 1, reply.size() - cut, out), fflush(out);
        }
        _exit(0);
    }

    close(to_judge[0]), close(to_us[1]);
    int saved_in = dup(0), saved_out = dup(1);
    fflush(stdout);
    dup2(to_us[0], 0), dup2(to_judge[1], 1);
    close(to_us[0]), close(to_judge[1]);
    inptr = buf_len = 0;
    alarm(60);
    for (long long r = 0; r < rounds; r++){
        long long x = random_value<long long>() / 4;
        int k = stress::rand_int(1, r % 10 ? 5 : 3000);
        write('?', x, k);

        string line;
        assert(read_line(line) && line == "round " + to_string(r));
        vector<long long> got(k);
        assert(read(got) == 1);
        for (int j = 0; j < k; j++) assert(got[j] == x + j);
    }
    alarm(0);

    flush();
    dup2(saved_in, 0), dup2(saved_out, 1);
    close(saved_in), close(saved_out);
    inptr = buf_len = 0;
    int judge_status;
    assert(waitpid(judge, &judge_status, 0) == judge && WIFEXITED(judge_status) && WEXITSTATUS(judge_status) == 0);

    /// A forgotten flush() must not lose output, the child writes the pending buffer when it exits
    assert(freopen(OUT.c_str(), "wb", stdout));
    write("forgotten", 42);
    pid_t pid = fork();
    if (pid == 0) exit(0);
    outptr = 0;
    int status;
    assert(waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0);
    assert(read_file(OUT) == "forgotten 42\n");

    /// A failed read(), here on a closed stdin, is the end of input, never a replay of stale buffer bytes
    close(0);
    inptr = buf_len = 0, inbuf[0] = 'z';
    char stale;
    assert(read(stale) == 0);

    remove(IN.c_str()), remove(OUT.c_str());
    return 0;
}
