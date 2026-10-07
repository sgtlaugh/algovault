#include "common.h"

#define main library_main
#include "../code_library/fast_io.cpp"
#undef main

const string IN = "/tmp/algovault_fast_io_in.txt", OUT = "/tmp/algovault_fast_io_out.txt";

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
        flush();
        fflush(stdout);
        assert(read_file(OUT) == expected_out);
    }
    return 0;
}
