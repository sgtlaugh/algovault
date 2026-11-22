#include "../common.h"

#define main library_main
#include "../../code_library/misc/expression_parser.cpp"
#undef main

/// A random expression tree, printed with only the parentheses precedence needs plus some extra ones
struct Node{
    char op;  /// 'n' operand, '~' unary minus, else a binary operator
    string token;
    int left = -1, right = -1;
};

struct Generator{
    vector<Node> nodes;
    bool integer;
    vector<string> names = {"x", "y1", "_z", "var"};

    int leaf(bool small_exponent){
        if (small_exponent) return add({'n', to_string(stress::rand_int(0, 3))});
        int kind = stress::rand_int(0, 5);
        if (kind == 0 && !integer) return add({'n', names[stress::rand_int(0, 3)]});
        if (kind == 1 && !integer) return add({'n', to_string(stress::rand_int(0, 99)) + ".5"});
        if (kind == 2) return add({'n', to_string(stress::rand_int(10, 999))});
        if (kind == 3 && !integer){
            string mantissa = to_string(stress::rand_int(1, 9)) + (stress::rand_int(0, 1) ? ".5" : "");
            string sign = stress::rand_int(0, 1) ? "-" : "+";
            return add({'n', mantissa + (stress::rand_int(0, 1) ? "e" : "E") + sign + to_string(stress::rand_int(0, 3))});
        }

        return add({'n', to_string(stress::rand_int(0, 9))});
    }

    int build(int leaves){
        if (leaves == 1) return stress::rand_int(0, 4) ? leaf(false) : add({'~', "", build(1)});
        if (stress::rand_int(0, 5) == 0) return add({'~', "", build(leaves)});

        const char ops[] = "+-*/^";
        char op = ops[stress::rand_int(0, 4)];
        if (op == '^'){
            int base = build(leaves - 1);
            return add({'^', "", base, leaf(integer)});
        }

        int split = stress::rand_int(1, leaves - 1);
        int left = build(split);
        return add({op, "", left, build(leaves - split)});
    }

    int add(Node node){
        nodes.push_back(node);
        return nodes.size() - 1;
    }

    /// 1 for + -, 2 for * /, 3 for unary minus, 4 for ^, 5 for an operand
    int level(int v){
        char op = nodes[v].op;
        if (op == '+' || op == '-') return 1;
        if (op == '*' || op == '/') return 2;
        if (op == '~') return 3;
        if (op == '^') return 4;
        return 5;
    }

    string space(){
        int kind = stress::rand_int(0, 7);
        return kind == 0 ? " " : kind == 1 ? "\t" : "";
    }

    string print(int v, bool wrap){
        string body;
        const Node& node = nodes[v];
        if (node.op == 'n') body = node.token;
        else if (node.op == '~') body = "-" + space() + print(node.left, level(node.left) < 3);
        else if (node.op == '^') body = print(node.left, level(node.left) < 5) + space() + "^" + space() + print(node.right, level(node.right) < 3);
        else{
            int need = level(v);
            body = print(node.left, level(node.left) < need) + space() + node.op + space() + print(node.right, level(node.right) <= need);
        }

        if (!wrap && stress::rand_int(0, 7) == 0) wrap = true;
        if (!wrap) return body;
        return "(" + space() + (stress::rand_int(0, 5) == 0 ? "+" : "") + body + space() + ")";
    }

    void postfix(int v, vector<string>& out){
        const Node& node = nodes[v];
        if (node.left != -1) postfix(node.left, out);
        if (node.right != -1) postfix(node.right, out);
        out.push_back(node.op == 'n' ? node.token : string(1, node.op));
    }
};

/// The tree's value in __int128, false on division by zero or once a value leaves +-1e12, so long long cannot overflow
bool reference_integer(const Generator& g, int v, __int128& res){
    const Node& node = g.nodes[v];
    if (node.op == 'n'){
        res = stoll(node.token);
        return true;
    }

    __int128 a, b = 0;
    if (!reference_integer(g, node.left, a)) return false;
    if (node.right != -1 && !reference_integer(g, node.right, b)) return false;
    if (node.op == '~') res = -a;
    else if (node.op == '+') res = a + b;
    else if (node.op == '-') res = a - b;
    else if (node.op == '*') res = a * b;
    else if (node.op == '/'){
        if (b == 0) return false;
        res = a / b;
    }
    else{
        res = 1;
        for (int i = 0; i < b; i++) res *= a;
    }

    __int128 limit = 1000000000000LL;
    return -limit <= res && res <= limit;
}

double reference_double(const Generator& g, int v, const map<string, double>& memory){
    const Node& node = g.nodes[v];
    if (node.op == 'n') return isdigit(node.token[0]) ? stod(node.token) : memory.at(node.token);

    double a = reference_double(g, node.left, memory);
    if (node.op == '~') return -a;
    double b = reference_double(g, node.right, memory);
    if (node.op == '+') return a + b;
    if (node.op == '-') return a - b;
    if (node.op == '*') return a * b;
    if (node.op == '/') return a / b;
    return pow(a, b);
}

/// a + b x, closed under + and - always and under * and / by a constant, as a linear equation's sides are
struct Linear{
    long long a, b;

    Linear operator+(const Linear& o) const{ return {a + o.a, b + o.b}; }
    Linear operator-(const Linear& o) const{ return {a - o.a, b - o.b}; }
    Linear operator-() const{ return {-a, -b}; }
    Linear operator*(const Linear& o) const{
        assert(b == 0 || o.b == 0);
        return {a * o.a, a * o.b + b * o.a};
    }
    Linear operator/(const Linear& o) const{
        assert(o.b == 0 && o.a != 0);
        return {a / o.a, b / o.a};
    }
};

long long parse_integer(const string& t){
    return stoll(t);
}

void check_integer_trees(){
    int checked = 0;
    while (checked < stress::scaled(15000)){
        Generator g;
        g.integer = true;
        int root = g.build(stress::rand_int(1, 10));
        __int128 expected;
        if (!reference_integer(g, root, expected)) continue;

        string text = g.space() + g.print(root, false) + g.space();
        Expression e(text);
        vector<string> postfix;
        g.postfix(root, postfix);
        assert(e.postfix == postfix);
        assert(e.evaluate(parse_integer) == (long long)expected);
        checked++;
    }
}

void check_double_trees(){
    for (int it = 0; it < stress::scaled(15000); it++){
        Generator g;
        g.integer = false;
        int root = g.build(stress::rand_int(1, 10));
        map<string, double> memory;
        for (const string& name : g.names) memory[name] = stress::rand_int(-50, 50) / 4.0;

        Expression e(g.print(root, false));
        vector<string> postfix;
        g.postfix(root, postfix);
        assert(e.postfix == postfix);

        double expected = reference_double(g, root, memory);
        double got = e.evaluate([&](const string& t){ return isdigit(t[0]) ? stod(t) : memory.at(t); });
        assert((isnan(expected) && isnan(got)) || expected == got);
    }
}

/// x appears only where the expression stays linear, so evaluating each side gives a + b x exactly
void check_linear_type(){
    Expression left("3 * (x - 2) + -(4 - x) * 2");
    Expression right("(2 * x + 10) / 2 * 4 - x");
    auto value_of = [](const string& t){ return t == "x" ? Linear{0, 1} : Linear{stoll(t), 0}; };
    Linear l = left.evaluate(value_of), r = right.evaluate(value_of);
    assert(l.a == -14 && l.b == 5);
    assert(r.a == 20 && r.b == 3);
}

void check_fixed_cases(){
    auto value = [](const string& text){ return Expression(text).evaluate(parse_integer); };
    assert(value("10-4-3") == 3);
    assert(value("100/10/5") == 2);
    assert(value("2^3^2") == 512);
    assert(value("-2^2") == -4);
    assert(value("(-2)^2") == 4);
    assert(value("--2") == 2);
    assert(value("2*-3") == -6);
    assert(value("2- -3") == 5);
    assert(value("+7") == 7);
    assert(value("7 / -2") == -3);
    assert(value("0^0") == 1);
    assert(value("\t ( ( 42 ) ) ") == 42);

    double half = Expression("2^-1").evaluate([](const string& t){ return stod(t); });
    assert(half == 0.5);
    assert((Expression("x1_2.5").postfix == vector<string>{"x1_2.5"}));
    assert((Expression("2x").postfix == vector<string>{"2x"}));

    /// A sign after the e of a number literal, before a digit, stays in the operand; elsewhere it is an operator
    assert((Expression("1e-5*2").postfix == vector<string>{"1e-5", "2", "*"}));
    assert((Expression("2.5E+3").postfix == vector<string>{"2.5E+3"}));
    assert((Expression("2e-x").postfix == vector<string>{"2e", "x", "-"}));
    assert((Expression("x1e-5").postfix == vector<string>{"x1e", "5", "-"}));
    assert((Expression("1e5e-3").postfix == vector<string>{"1e5e", "3", "-"}));
    assert((Expression("e-5").postfix == vector<string>{"e", "5", "-"}));
    assert(Expression("1e-5 * 2").evaluate([](const string& t){ return stod(t); }) == 2e-5);

    __int128 nine = Expression("3^2").evaluate([](const string& t){ return (__int128)stoll(t); });
    assert(nine == 9);
}

/// Deep nesting through every parse level; 5000 is what fits 8 MB under the sanitizers, 2e4 fits at -O2
void check_deep_nesting(){
    int depth = 5000;
    string parens = string(depth, '(') + "5" + string(depth, ')');
    assert(Expression(parens).evaluate(parse_integer) == 5);

    string minuses = string(depth, '-') + "5";
    assert(Expression(minuses).evaluate(parse_integer) == 5);

    string chain = "1";
    for (int i = 0; i < 100000; i++) chain += "+1";
    assert(Expression(chain).evaluate(parse_integer) == 100001);
}

void check_asserts(){
    const char* bad[] = {"", "2+", "(2", "2)", "2 3", "*2", "()", "2^", "^2", "(2))", "2 + + "};
    for (const char* text : bad) stress::expect_assert("misc/expression_parser.cpp", [&]{ Expression e(text); });

    stress::expect_assert("misc/expression_parser.cpp", []{ Expression("2^-1").evaluate(parse_integer); });
    stress::expect_assert("misc/expression_parser.cpp", []{
        Expression("x^2").evaluate([](const string&){ return Linear{0, 1}; });
    });
}

int main(){
    check_fixed_cases();
    check_integer_trees();
    check_double_trees();
    check_linear_type();
    check_deep_nesting();
    check_asserts();
    return 0;
}
