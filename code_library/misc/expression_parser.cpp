/***
 *
 * Expression Parser
 * Parses an infix expression once, then evaluates it over any number type or lists it in postfix order
 *
 * Complexity: O(|s|) to parse, O(|s|) operations on T per evaluate
 *
 * Expression e(s): parses s, asserting that it fits the grammar in the notes
 * e.evaluate(value_of): the value of the expression, value_of(token) giving each operand's value, all of one type T
 * e.postfix: the tokens in postfix order, unary minus as "~"
 *
 * Notes:
 *   - Precedence, lowest first: + and -, then * and /, then unary - and +, then ^, which is right-associative
 *   - So a-b-c is (a-b)-c, -2^2 is -4, 2^3^2 is 512, and 2^-1 is allowed
 *   - An operand is a maximal run of letters, digits, '_' and '.', like x, 42 or 3.5; value_of reads it
 *   - A number literal keeps the sign of its exponent, so 1e-5 and 2.5E+3 are single operands
 *   - Spaces and tabs are skipped anywhere, and there is no implicit product: 2x is the single operand "2x"
 *   - T needs +, -, *, / and unary -, with its own meaning: integer / truncates as in C++
 *   - ^ is pow for a floating T and exact squaring for an integer T or __int128, exponent >= 0; other T assert on ^
 *   - A struct a + b x with these operators evaluates a side of a linear equation in x
 *   - Assignments, comparisons, functions and side effects are not parsed: split them off first, as in main()
 *   - Each nesting level recurses: 2e4 levels of parentheses or unary minus fit an 8 MB stack at -O2
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Expression{
    vector<string> postfix;
    string s;
    int pos = 0;

    Expression(const string& text) : s(text){
        parse_sum();
        skip_spaces();
        assert(pos == (int)s.size());
    }

    template <typename F>
    auto evaluate(F value_of) const{
        using T = decltype(value_of(string()));
        vector<T> stack;
        for (const string& token : postfix){
            if (is_operand(token)){
                stack.push_back(value_of(token));
                continue;
            }
            if (token == "~"){
                stack.back() = -stack.back();
                continue;
            }

            T b = stack.back();
            stack.pop_back();
            T& a = stack.back();
            if (token == "+") a = a + b;
            else if (token == "-") a = a - b;
            else if (token == "*") a = a * b;
            else if (token == "/") a = a / b;
            else a = power(a, b);
        }

        return stack.back();
    }

    void parse_sum(){
        parse_product();
        while (next_is('+') || next_is('-')){
            char op = s[pos++];
            parse_product();
            postfix.push_back(string(1, op));
        }
    }

    void parse_product(){
        parse_unary();
        while (next_is('*') || next_is('/')){
            char op = s[pos++];
            parse_unary();
            postfix.push_back(string(1, op));
        }
    }

    void parse_unary(){
        if (next_is('-')){
            pos++;
            parse_unary();
            postfix.push_back("~");
        }
        else if (next_is('+')){
            pos++;
            parse_unary();
        }
        else parse_power();
    }

    /// The exponent is parsed as a unary term, which makes ^ right-associative and allows 2^-1
    void parse_power(){
        parse_primary();
        if (next_is('^')){
            pos++;
            parse_unary();
            postfix.push_back("^");
        }
    }

    void parse_primary(){
        if (next_is('(')){
            pos++;
            parse_sum();
            assert(next_is(')'));
            pos++;
            return;
        }

        int start = pos;
        while (pos < (int)s.size() && (is_operand_char(s[pos]) || is_exponent_sign(start))) pos++;
        assert(pos > start);
        postfix.push_back(s.substr(start, pos - start));
    }

    /// A sign is part of the operand only right after the e of a number literal and before a digit, as in 1e-5
    bool is_exponent_sign(int start) const{
        if ((s[pos] != '+' && s[pos] != '-') || pos + 1 == (int)s.size() || !isdigit((unsigned char)s[pos + 1])) return false;
        if (pos - start < 2 || (s[pos - 1] != 'e' && s[pos - 1] != 'E')) return false;
        for (int i = start; i < pos - 1; i++){
            if (!isdigit((unsigned char)s[i]) && s[i] != '.') return false;
        }

        return true;
    }

    bool next_is(char c){
        skip_spaces();
        return pos < (int)s.size() && s[pos] == c;
    }

    void skip_spaces(){
        while (pos < (int)s.size() && (s[pos] == ' ' || s[pos] == '\t')) pos++;
    }

    static bool is_operand_char(char c){
        return isalnum((unsigned char)c) || c == '_' || c == '.';
    }

    static bool is_operand(const string& token){
        return is_operand_char(token[0]);
    }

    /// Picks pow or exact squaring by type, so a double exponent is never truncated to an integer
    template <typename T>
    static T power(T base, T exponent){
        if constexpr (is_floating_point_v<T>) return pow(base, exponent);
        else if constexpr (is_integral_v<T> || is_same_v<T, __int128>){
            assert(exponent >= 0);
            T res = 1;
            while (exponent > 0){
                if (exponent & 1) res *= base;
                exponent >>= 1;
                if (exponent > 0) base *= base;  /// a square after the last bit would overflow for no reason
            }

            return res;
        }
        else{
            assert(!"^ needs an arithmetic T");
            return base;
        }
    }
};

int main(){
    /// 2 * (3 + 4) - 5 = 9, its operands read as integers
    Expression e("2 * (3 + 4) - 5");
    assert((e.postfix == vector<string>{"2", "3", "4", "+", "*", "5", "-"}));
    assert(e.evaluate([](const string& t){ return stoll(t); }) == 9);

    /// Unary minus binds looser than ^, and ^ groups to the right
    Expression unary("-2^2 + 2^3^2");
    assert(unary.evaluate([](const string& t){ return stoll(t); }) == 508);

    /// Variables through value_of, with doubles: x = 10, y = 4 gives 10 / 4 + 0.5 = 3
    map<string, double> memory = {{"x", 10}, {"y", 4}};
    Expression with_vars("x / y + 0.5");
    assert(with_vars.evaluate([&](const string& t){ return isdigit(t[0]) ? stod(t) : memory[t]; }) == 3);

    return 0;
}
