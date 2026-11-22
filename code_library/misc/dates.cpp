/***
 *
 * Dates
 * Calendar arithmetic: a date to its day number and back, weekdays, leap years, Gregorian or Julian calendar
 *
 * Complexity: O(1) per call
 *
 * Date: year, month 1..12, day 1..31
 * date_of(n, cal = GREGORIAN): the date of day number n
 * day_number(y, m, d, cal = GREGORIAN): the Julian Day Number of the date, consecutive dates differ by 1
 * day_number(t, cal = GREGORIAN): the same for a Date t
 * day_of_week(n): weekday of day number n, 0 = Sunday .. 6 = Saturday
 * days_in_month(y, m, cal = GREGORIAN): 28 to 31
 * is_leap(y, cal = GREGORIAN): Gregorian divisible by 4 and not by 100, or by 400; Julian divisible by 4
 *
 * Notes:
 *   - Days between two dates: day_number(b) - day_number(a); n days after a: date_of(day_number(a) + n)
 *   - Both calendars share one scale, day 0 is January 1, 4713 BC Julian, and 2000-01-01 Gregorian is 2451545
 *   - Calendars are proleptic, any year in [-1e15, 1e15], year 0 is 1 BC, -1 is 2 BC and so on
 *   - Rome switched after Julian 1582-10-04 to Gregorian 1582-10-15, Britain and the US after 1752-09-02 to 09-14
 *   - day_number checks nothing, a valid date has 1 <= m <= 12 and 1 <= d <= days_in_month(y, m)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

enum Calendar{GREGORIAN, JULIAN};

struct Date{
    long long year;
    int month, day;
};

/// Years start on March 1 here, so the leap day ends a year and month lengths follow the fixed 153-day pattern
Date date_of(long long n, Calendar cal = GREGORIAN){
    long long cycle = cal == GREGORIAN ? 146097 : 1461, z = n - (cal == GREGORIAN ? 1721120 : 1721118);
    long long era = (z >= 0 ? z : z - cycle + 1) / cycle, doe = z - era * cycle;

    long long yoe;
    if (cal == GREGORIAN) yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    else yoe = (doe - doe / 1460) / 365;
    long long doy = doe - 365 * yoe - yoe / 4 + (cal == GREGORIAN ? yoe / 100 : 0);

    int mp = (5 * doy + 2) / 153, d = doy - (153 * mp + 2) / 5 + 1, m = mp < 10 ? mp + 3 : mp - 9;
    long long y = yoe + era * (cal == GREGORIAN ? 400 : 4) + (m <= 2);
    return {y, m, d};
}

/// 146097 days per 400 Gregorian years, 1461 per 4 Julian years; 1721120 and 1721118 are March 1 of year 0 in each
long long day_number(long long y, int m, int d, Calendar cal = GREGORIAN){
    y -= m <= 2;
    long long years = cal == GREGORIAN ? 400 : 4;
    long long era = (y >= 0 ? y : y - years + 1) / years, yoe = y - era * years;
    long long doy = (153 * ((m + 9) % 12) + 2) / 5 + d - 1;

    if (cal == GREGORIAN) return era * 146097 + yoe * 365 + yoe / 4 - yoe / 100 + doy + 1721120;
    return era * 1461 + yoe * 365 + doy + 1721118;
}

long long day_number(const Date& t, Calendar cal = GREGORIAN){
    return day_number(t.year, t.month, t.day, cal);
}

int day_of_week(long long n){
    return ((n + 1) % 7 + 7) % 7;
}

int days_in_month(long long y, int m, Calendar cal = GREGORIAN){
    static const int days[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m == 2 && (cal == JULIAN ? y % 4 == 0 : (y % 4 == 0 && y % 100 != 0) || y % 400 == 0)) return 29;
    return days[m];
}

bool is_leap(long long y, Calendar cal = GREGORIAN){
    return days_in_month(y, 2, cal) == 29;
}

int main(){
    /// 2000-01-01 is Julian Day 2451545, a Saturday, and 60 days after 1999-12-31 is the leap day 2000-02-29
    long long n = day_number(2000, 1, 1);
    assert(n == 2451545);
    assert(day_of_week(n) == 6);

    Date leap_day = date_of(day_number(1999, 12, 31) + 60);
    assert(leap_day.year == 2000);
    assert(leap_day.month == 2);
    assert(leap_day.day == 29);

    assert(day_number(2004, 3, 1) - day_number(2004, 2, 1) == 29);

    /// Days between two Dates: 2024-02-28 to 2024-03-01 crosses the leap day
    Date start{2024, 2, 28}, end{2024, 3, 1};
    assert(day_number(end) - day_number(start) == 2);

    assert(!is_leap(1900));
    assert(is_leap(2000));
    assert(days_in_month(2023, 4) == 30);

    /// Julian calendar: 1900 is a leap year, and the day after Julian 1752-09-02 is Gregorian 1752-09-14
    assert(is_leap(1900, JULIAN));
    assert(day_number(1752, 9, 2, JULIAN) + 1 == day_number(1752, 9, 14));

    Date origin = date_of(0, JULIAN);
    assert(origin.year == -4712);
    assert(origin.month == 1);
    assert(origin.day == 1);

    return 0;
}
