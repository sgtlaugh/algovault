/***
 *
 * Next Palindrome
 * Smallest palindromic number strictly greater than a given number, as decimal strings of any length
 *
 * Complexity: O(n) for n digits
 *
 * next_palindrome(s): s is a non-negative integer without leading zeros ("0" allowed)
 *
 * Mirror the left half onto the right; if that is not larger than s, add one to the left half (middle digit
 * included) and mirror again; when the left half is all nines the answer is 1 0...0 1 with one more digit
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

string next_palindrome(const string& s){
    int n = s.size();
    assert(n > 0 && (s[0] != '0' || n == 1));

    string res = s;
    for (int i = 0; i < n / 2; i++) res[n - 1 - i] = res[i];
    if (res > s) return res;

    int i = (n - 1) / 2;
    while (i >= 0 && res[i] == '9') res[i--] = '0';
    if (i < 0) return "1" + string(n - 1, '0') + "1";
    res[i]++;
    for (int k = 0; k < n / 2; k++) res[n - 1 - k] = res[k];
    return res;
}

int main(){
    assert(next_palindrome("0") == "1");
    assert(next_palindrome("8") == "9");
    assert(next_palindrome("9") == "11");
    assert(next_palindrome("10") == "11");
    assert(next_palindrome("11") == "22");
    assert(next_palindrome("99") == "101");
    assert(next_palindrome("123") == "131");
    assert(next_palindrome("808") == "818");
    assert(next_palindrome("2133") == "2222");
    assert(next_palindrome("1991") == "2002");
    assert(next_palindrome("999") == "1001");
    assert(next_palindrome("12921") == "13031");
    return 0;
}
