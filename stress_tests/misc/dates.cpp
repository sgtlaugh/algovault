#include "../common.h"

#define main library_main
#include "../../code_library/misc/dates.cpp"
#undef main

bool leap_reference(long long y, Calendar cal){
    if (cal == JULIAN) return y % 4 == 0;
    if (y % 400 == 0) return true;
    if (y % 100 == 0) return false;
    return y % 4 == 0;
}

int month_length(long long y, int m, Calendar cal){
    if (m == 2) return leap_reference(y, cal) ? 29 : 28;
    if (m == 4 || m == 6 || m == 9 || m == 11) return 30;
    return 31;
}

/// Walks one day at a time from a known anchor, checking every call against the walked date and weekday
void walk(Calendar cal, Date start, long long start_number, long long days, int step){
    Date cur = start;
    long long n = start_number;
    int weekday = day_of_week(start_number);
    for (long long i = 0; i < days; i++){
        assert(day_number(cur.year, cur.month, cur.day, cal) == n);
        assert(day_number(cur, cal) == n);
        Date back = date_of(n, cal);
        assert(back.year == cur.year && back.month == cur.month && back.day == cur.day);
        assert(day_of_week(n) == weekday);
        assert(days_in_month(cur.year, cur.month, cal) == month_length(cur.year, cur.month, cal));
        assert(is_leap(cur.year, cal) == leap_reference(cur.year, cal));

        if (step > 0){
            if (++cur.day > month_length(cur.year, cur.month, cal)){
                cur.day = 1;
                if (++cur.month > 12) cur.month = 1, cur.year++;
            }
            n++, weekday = (weekday + 1) % 7;
        }
        else{
            if (--cur.day == 0){
                if (--cur.month == 0) cur.month = 12, cur.year--;
                cur.day = month_length(cur.year, cur.month, cal);
            }
            n--, weekday = (weekday + 6) % 7;
        }
    }
}

int main(){
    /// Day 0 is Monday, Julian -4712-01-01 and Gregorian -4713-11-24: forwards past year 5000, backwards past -10000
    assert(day_of_week(0) == 1);
    long long forward = stress::scaled(3500000), backward = stress::scaled(2200000);
    walk(JULIAN, {-4712, 1, 1}, 0, forward, 1);
    walk(JULIAN, {-4712, 1, 1}, 0, backward, -1);
    walk(GREGORIAN, {-4713, 11, 24}, 0, forward, 1);
    walk(GREGORIAN, {-4713, 11, 24}, 0, backward, -1);

    /// Huge years: round trips, and whole cycles of 400 Gregorian or 4 Julian years add a fixed day count
    for (long long it = 0; it < stress::scaled(200000); it++){
        Calendar cal = stress::rand_int(0, 1) ? GREGORIAN : JULIAN;
        long long y = stress::rand_int(-1000000000000000LL, 999999999000000LL);
        int m = stress::rand_int(1, 12), d = stress::rand_int(1, month_length(y, m, cal));
        long long n = day_number(y, m, d, cal);
        Date back = date_of(n, cal);
        assert(back.year == y && back.month == m && back.day == d);

        long long cycle_years = cal == GREGORIAN ? 400 : 4, cycle_days = cal == GREGORIAN ? 146097 : 1461;
        long long k = stress::rand_int(0, 1000);
        assert(day_number(y + k * cycle_years, m, d, cal) - n == k * cycle_days);
    }

    /// The year boundary and the leap day itself, around centuries that are and are not leap years
    for (long long y : {1600LL, 1700LL, 1900LL, 2000LL, 2100LL, -400LL, -100LL, 0LL}){
        assert(day_number(y, 3, 1) - day_number(y, 2, 28) == (leap_reference(y, GREGORIAN) ? 2 : 1));
        assert(day_number(y + 1, 1, 1) - day_number(y, 12, 31) == 1);
        assert(day_number(y, 3, 1, JULIAN) - day_number(y, 2, 28, JULIAN) == 2);
    }

    return 0;
}
