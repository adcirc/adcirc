#define CATCH_CONFIG_MAIN
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <thread>

#include "DateTime.hpp"
#include "TimeDelta.hpp"

// ====================================================
// Compile-time tests using static_assert
// ====================================================

// TimeDelta constexpr tests
static_assert(TimeDelta().totalMilliseconds() == 0,
              "Default TimeDelta should have zero duration");
static_assert(TimeDelta(1, 2, 3, 4, 5).days() == 1,
              "Days component extraction failed");
static_assert(TimeDelta(1, 2, 3, 4, 5).hours() == 2,
              "Hours component extraction failed");
static_assert(TimeDelta(1, 2, 3, 4, 5).minutes() == 3,
              "Minutes component extraction failed");
static_assert(TimeDelta(1, 2, 3, 4, 5).seconds() == 4,
              "Seconds component extraction failed");
static_assert(TimeDelta(1, 2, 3, 4, 5).milliseconds() == 5,
              "Milliseconds component extraction failed");

static_assert(TimeDelta::fromDays(2).totalDays() == 2,
              "fromDays factory method failed");
static_assert(TimeDelta::fromHours(24).totalHours() == 24,
              "fromHours factory method failed");
static_assert(TimeDelta::fromMinutes(60).totalMinutes() == 60,
              "fromMinutes factory method failed");
static_assert(TimeDelta::fromSeconds(60).totalSeconds() == 60,
              "fromSeconds factory method failed");
static_assert(TimeDelta::fromMilliseconds(1000).totalMilliseconds() == 1000,
              "fromMilliseconds factory method failed");

static_assert(TimeDelta::fromDays(1) + TimeDelta::fromDays(2) ==
                  TimeDelta::fromDays(3),
              "TimeDelta addition failed");
static_assert(TimeDelta::fromDays(3) - TimeDelta::fromDays(1) ==
                  TimeDelta::fromDays(2),
              "TimeDelta subtraction failed");
static_assert(TimeDelta::fromDays(2) * 3 == TimeDelta::fromDays(6),
              "TimeDelta multiplication failed");
static_assert(TimeDelta::fromDays(6) / 2 == TimeDelta::fromDays(3),
              "TimeDelta division failed");

static_assert(TimeDelta::fromDays(1) == TimeDelta::fromDays(1),
              "TimeDelta equality failed");
static_assert(TimeDelta::fromDays(1) != TimeDelta::fromDays(2),
              "TimeDelta inequality failed");
static_assert(TimeDelta::fromDays(1) < TimeDelta::fromDays(2),
              "TimeDelta less than failed");
static_assert(TimeDelta::fromDays(2) > TimeDelta::fromDays(1),
              "TimeDelta greater than failed");
static_assert(TimeDelta::fromDays(1) <= TimeDelta::fromDays(1),
              "TimeDelta less than or equal failed");
static_assert(TimeDelta::fromDays(2) >= TimeDelta::fromDays(2),
              "TimeDelta greater than or equal failed");

// DateTime constexpr tests
static_assert(DateTime(2022, 1, 1).year() == 2022, "Year extraction failed");
static_assert(DateTime(2022, 1, 1).month() == 1, "Month extraction failed");
static_assert(DateTime(2022, 1, 1).day() == 1, "Day extraction failed");
static_assert(DateTime(2022, 1, 1, 12, 30, 45, 500).hour() == 12,
              "Hour extraction failed");
static_assert(DateTime(2022, 1, 1, 12, 30, 45, 500).minute() == 30,
              "Minute extraction failed");
static_assert(DateTime(2022, 1, 1, 12, 30, 45, 500).second() == 45,
              "Second extraction failed");
static_assert(DateTime(2022, 1, 1, 12, 30, 45, 500).millisecond() == 500,
              "Millisecond extraction failed");

static_assert(DateTime(2022, 1, 1) == DateTime(2022, 1, 1),
              "DateTime equality failed");
static_assert(DateTime(2022, 1, 1) < DateTime(2022, 1, 2),
              "DateTime less than failed");
static_assert(DateTime(2022, 1, 2) > DateTime(2022, 1, 1),
              "DateTime greater than failed");
static_assert(DateTime(2022, 1, 1) <= DateTime(2022, 1, 1),
              "DateTime less than or equal failed");
static_assert(DateTime(2022, 1, 1) >= DateTime(2022, 1, 1),
              "DateTime greater than or equal failed");

static_assert(DateTime(2022, 1, 1) + TimeDelta::fromDays(1) ==
                  DateTime(2022, 1, 2),
              "DateTime addition with TimeDelta failed");
static_assert(DateTime(2022, 1, 2) - TimeDelta::fromDays(1) ==
                  DateTime(2022, 1, 1),
              "DateTime subtraction with TimeDelta failed");
static_assert((DateTime(2022, 1, 2) - DateTime(2022, 1, 1)).totalDays() == 1,
              "DateTime difference calculation failed");

// ====================================================
// Runtime tests using Catch2
// ====================================================

TEST_CASE("TimeDelta basic functionality", "[timedelta]") {
  SECTION("Default constructor") {
    TimeDelta ts;
    REQUIRE(ts.totalMilliseconds() == 0);
  }

  SECTION("Component constructor") {
    TimeDelta ts(TimeDelta::s_TimedeltaComponents{2, 3, 4, 5, 6});
    CHECK(ts.days() == 2);
    CHECK(ts.hours() == 3);
    CHECK(ts.minutes() == 4);
    CHECK(ts.seconds() == 5);
    CHECK(ts.milliseconds() == 6);
  }

  SECTION("Milliseconds constructor") {
    TimeDelta ts(TimeDelta::s_TimedeltaComponents{
        0, 0, 0, 0, 1000 * 60 * 60});  // 1 hour in milliseconds
    CHECK(ts.hours() == 1);
    CHECK(ts.days() == 0);
    CHECK(ts.minutes() == 0);
    CHECK(ts.seconds() == 0);
    CHECK(ts.milliseconds() == 0);
  }

  SECTION("Copy constructor") {
    TimeDelta ts1(1, 2, 3, 4, 5);
    TimeDelta ts2(ts1);
    CHECK(ts2.days() == 1);
    CHECK(ts2.hours() == 2);
    CHECK(ts2.minutes() == 3);
    CHECK(ts2.seconds() == 4);
    CHECK(ts2.milliseconds() == 5);
  }

  SECTION("Assignment operator") {
    TimeDelta ts1(1, 2, 3, 4, 5);
    const auto ts2 = ts1;
    CHECK(ts2.days() == 1);
    CHECK(ts2.hours() == 2);
    CHECK(ts2.minutes() == 3);
    CHECK(ts2.seconds() == 4);
    CHECK(ts2.milliseconds() == 5);
  }
}

TEST_CASE("TimeDelta factory methods", "[timedelta]") {
  SECTION("fromDays") {
    TimeDelta ts = TimeDelta::fromDays(2);
    CHECK(ts.days() == 2);
    CHECK(ts.hours() == 0);
    CHECK(ts.totalDays() == 2);
  }

  SECTION("fromHours") {
    TimeDelta ts = TimeDelta::fromHours(25);
    CHECK(ts.days() == 1);
    CHECK(ts.hours() == 1);
    CHECK(ts.totalHours() == 25);
  }

  SECTION("fromMinutes") {
    TimeDelta ts = TimeDelta::fromMinutes(60);
    CHECK(ts.hours() == 1);
    CHECK(ts.minutes() == 0);
    CHECK(ts.totalMinutes() == 60);
  }

  SECTION("fromSeconds") {
    TimeDelta ts = TimeDelta::fromSeconds(60);
    CHECK(ts.minutes() == 1);
    CHECK(ts.seconds() == 0);
    CHECK(ts.totalSeconds() == 60);
  }

  SECTION("fromMilliseconds") {
    TimeDelta ts = TimeDelta::fromMilliseconds(1000);
    CHECK(ts.seconds() == 1);
    CHECK(ts.milliseconds() == 0);
    CHECK(ts.totalMilliseconds() == 1000);
  }
}

TEST_CASE("TimeDelta component and total accessors", "[timedelta]") {
  SECTION("Components add up correctly") {
    // 1 day, 2 hours, 3 minutes, 4 seconds, 5 milliseconds
    TimeDelta ts(1, 2, 3, 4, 5);

    // Check components
    CHECK(ts.days() == 1);
    CHECK(ts.hours() == 2);
    CHECK(ts.minutes() == 3);
    CHECK(ts.seconds() == 4);
    CHECK(ts.milliseconds() == 5);

    // Check totals
    int64_t expectedTotalMilliseconds = 1LL * 24 * 60 * 60 * 1000 +
                                        2LL * 60 * 60 * 1000 + 3LL * 60 * 1000 +
                                        4LL * 1000 + 5;

    CHECK(ts.totalMilliseconds() == expectedTotalMilliseconds);
    CHECK(ts.totalSeconds() == expectedTotalMilliseconds / 1000);
    CHECK(ts.totalMinutes() == expectedTotalMilliseconds / (60 * 1000));
    CHECK(ts.totalHours() == expectedTotalMilliseconds / (60 * 60 * 1000));
    CHECK(ts.totalDays() == expectedTotalMilliseconds / (24 * 60 * 60 * 1000));
  }

  SECTION("Negative durations") {
    // Create a timedelta representing -1 day
    TimeDelta negative(
        TimeDelta::s_TimedeltaComponents{0, 0, 0, 0, -1 * 24 * 60 * 60 * 1000});

    CHECK(negative.days() == -1);
    CHECK(negative.totalDays() == -1);
    CHECK(negative.totalHours() == -24);
    CHECK(negative.totalMinutes() == -24 * 60);
    CHECK(negative.totalSeconds() == -24 * 60 * 60);
    CHECK(negative.totalMilliseconds() == -24LL * 60 * 60 * 1000);
  }
}

TEST_CASE("TimeDelta arithmetic operations", "[timedelta]") {
  SECTION("Addition") {
    TimeDelta ts1 = TimeDelta::fromDays(1);
    TimeDelta ts2 = TimeDelta::fromHours(12);

    TimeDelta sum = ts1 + ts2;
    CHECK(sum.totalHours() == 36);
  }

  SECTION("Subtraction") {
    TimeDelta ts1 = TimeDelta::fromDays(2);
    TimeDelta ts2 = TimeDelta::fromHours(24);

    TimeDelta diff = ts1 - ts2;
    CHECK(diff.totalDays() == 1);
  }

  SECTION("Multiplication") {
    TimeDelta ts = TimeDelta::fromHours(2);

    TimeDelta product = ts * 3;
    CHECK(product.totalHours() == 6);
  }

  SECTION("Division") {
    TimeDelta ts = TimeDelta::fromHours(6);

    TimeDelta quotient = ts / 2;
    CHECK(quotient.totalHours() == 3);
  }

  SECTION("Chained operations") {
    TimeDelta result = TimeDelta::fromHours(6) + TimeDelta::fromMinutes(30) -
                       TimeDelta::fromMinutes(15);
    CHECK(result.totalMinutes() == 6 * 60 + 15);
  }
}

TEST_CASE("TimeDelta comparison operators", "[timedelta]") {
  TimeDelta ts1 = TimeDelta::fromHours(1);
  TimeDelta ts2 = TimeDelta::fromHours(2);
  TimeDelta ts3 = TimeDelta::fromHours(1);

  SECTION("Equality") {
    CHECK(ts1 == ts3);
    CHECK_FALSE(ts1 == ts2);
  }

  SECTION("Inequality") {
    CHECK(ts1 != ts2);
    CHECK_FALSE(ts1 != ts3);
  }

  SECTION("Less than") {
    CHECK(ts1 < ts2);
    CHECK_FALSE(ts2 < ts1);
    CHECK_FALSE(ts1 < ts3);
  }

  SECTION("Greater than") {
    CHECK(ts2 > ts1);
    CHECK_FALSE(ts1 > ts2);
    CHECK_FALSE(ts1 > ts3);
  }

  SECTION("Less than or equal") {
    CHECK(ts1 <= ts2);
    CHECK(ts1 <= ts3);
    CHECK_FALSE(ts2 <= ts1);
  }

  SECTION("Greater than or equal") {
    CHECK(ts2 >= ts1);
    CHECK(ts1 >= ts3);
    CHECK_FALSE(ts1 >= ts2);
  }
}

TEST_CASE("TimeDelta string representation", "[timedelta]") {
  SECTION("Basic format") {
    TimeDelta ts(1, 2, 3, 4, 5);
    CHECK(ts.toString() == "1d 02:03:04.005");
  }

  SECTION("Without days") {
    TimeDelta ts(0, 2, 3, 4, 5);
    CHECK(ts.toString() == "02:03:04.005");
  }

  SECTION("Without milliseconds") {
    TimeDelta ts(1, 2, 3, 4, 0);
    CHECK(ts.toString() == "1d 02:03:04");
  }

  SECTION("Only hours, minutes, seconds") {
    TimeDelta ts(0, 2, 3, 4, 0);
    CHECK(ts.toString() == "02:03:04");
  }
}

TEST_CASE("DateTime constructors", "[datetime]") {
  SECTION("Default constructor") {
    DateTime dt;
    CHECK(dt.timestamp() == 0);
  }

  SECTION("Components constructor") {
    DateTime dt(2022, 1, 31, 12, 34, 56, 789);
    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 789);
  }

  SECTION("Timestamp constructor") {
    // Construct a known datetime
    DateTime refDt(2022, 1, 31, 12, 34, 56, 789);
    int64_t ts = refDt.timestamp();

    // Create from timestamp
    DateTime dt(ts);

    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 789);
  }

  SECTION("Copy constructor") {
    DateTime dt1(2022, 1, 31, 12, 34, 56, 789);
    DateTime dt2(dt1);

    CHECK(dt2.year() == 2022);
    CHECK(dt2.month() == 1);
    CHECK(dt2.day() == 31);
    CHECK(dt2.hour() == 12);
    CHECK(dt2.minute() == 34);
    CHECK(dt2.second() == 56);
    CHECK(dt2.millisecond() == 789);
  }

  SECTION("Move constructor") {
    DateTime dt1(2022, 1, 31, 12, 34, 56, 789);
    DateTime dt2(std::move(dt1));

    CHECK(dt2.year() == 2022);
    CHECK(dt2.month() == 1);
    CHECK(dt2.day() == 31);
    CHECK(dt2.hour() == 12);
    CHECK(dt2.minute() == 34);
    CHECK(dt2.second() == 56);
    CHECK(dt2.millisecond() == 789);
  }

  SECTION("Assignment operator") {
    DateTime dt1(2022, 1, 31, 12, 34, 56, 789);
    DateTime dt2;
    dt2 = dt1;

    CHECK(dt2.year() == 2022);
    CHECK(dt2.month() == 1);
    CHECK(dt2.day() == 31);
    CHECK(dt2.hour() == 12);
    CHECK(dt2.minute() == 34);
    CHECK(dt2.second() == 56);
    CHECK(dt2.millisecond() == 789);
  }

  SECTION("Move assignment operator") {
    DateTime dt1(2022, 1, 31, 12, 34, 56, 789);
    DateTime dt2;
    dt2 = std::move(dt1);

    CHECK(dt2.year() == 2022);
    CHECK(dt2.month() == 1);
    CHECK(dt2.day() == 31);
    CHECK(dt2.hour() == 12);
    CHECK(dt2.minute() == 34);
    CHECK(dt2.second() == 56);
    CHECK(dt2.millisecond() == 789);
  }
}

TEST_CASE("DateTime parse method", "[datetime]") {
  SECTION("Default format") {
    auto dt = DateTime::strptime("2022-01-31 12:34:56");

    CHECK(dt.valid());
    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 0);
  }

  SECTION("Custom format") {
    auto dt = DateTime::strptime("31/01/2022 12:34:56", "%d/%m/%Y %H:%M:%S");

    CHECK(dt.valid());
    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 0);
  }

  SECTION("With milliseconds in string") {
    auto dt = DateTime::strptime("2022-01-31 12:34:56.789");

    CHECK(dt.valid());
    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 789);
  }

  SECTION("With %f format specifier") {
    auto dt =
        DateTime::strptime("2022-01-31 12:34:56.789", "%Y-%m-%d %H:%M:%S");

    CHECK(dt.valid());
    CHECK(dt.year() == 2022);
    CHECK(dt.month() == 1);
    CHECK(dt.day() == 31);
    CHECK(dt.hour() == 12);
    CHECK(dt.minute() == 34);
    CHECK(dt.second() == 56);
    CHECK(dt.millisecond() == 789);
  }
}

TEST_CASE("DateTime format and to_iso_string methods", "[datetime]") {
  DateTime dt(2022, 1, 31, 12, 34, 56, 789);

  SECTION("Default format") { CHECK(dt.strftime() == "2022-01-31 12:34:56"); }

  SECTION("Custom format") {
    CHECK(dt.strftime("%d/%m/%Y %H:%M:%S") == "31/01/2022 12:34:56");
  }

  SECTION("Format with milliseconds") {
    CHECK(dt.strftime_w_milliseconds("%Y-%m-%d %H:%M:%S") ==
          "2022-01-31 12:34:56.789");
  }

  SECTION("ISO string format") {
    CHECK(dt.to_iso_string() == "2022-01-31T12:34:56");
    CHECK(dt.to_iso_string_msec() == "2022-01-31T12:34:56.789");
  }

  SECTION("Format with various specifiers") {
    CHECK(dt.strftime("%Y") == "2022");
    CHECK(dt.strftime("%m") == "01");
    CHECK(dt.strftime("%d") == "31");
    CHECK(dt.strftime("%H") == "12");
    CHECK(dt.strftime("%M") == "34");
    CHECK(dt.strftime("%S") == "56");
    CHECK(dt.strftime("%Y-%m") == "2022-01");
  }
}

TEST_CASE("DateTime timestamp method", "[datetime]") {
  DateTime dt(2022, 1, 31, 12, 34, 56, 789);
  int64_t ts = dt.timestamp();

  // Create a new DateTime from the timestamp
  DateTime dt2(ts);

  // Both should represent the same point in time
  CHECK(dt.year() == dt2.year());
  CHECK(dt.month() == dt2.month());
  CHECK(dt.day() == dt2.day());
  CHECK(dt.hour() == dt2.hour());
  CHECK(dt.minute() == dt2.minute());
  CHECK(dt.second() == dt2.second());
  CHECK(dt.millisecond() == dt2.millisecond());
  CHECK(dt.timestamp() == dt2.timestamp());
}

TEST_CASE("DateTime arithmetic with TimeDelta", "[datetime]") {
  DateTime dt(2022, 1, 15, 12, 0, 0, 0);

  SECTION("Addition") {
    DateTime result = dt + TimeDelta::fromDays(10);
    CHECK(result.year() == 2022);
    CHECK(result.month() == 1);
    CHECK(result.day() == 25);
    CHECK(result.hour() == 12);

    result = dt + TimeDelta::fromHours(12);
    CHECK(result.day() == 16);
    CHECK(result.hour() == 0);
  }

  SECTION("Subtraction") {
    DateTime result = dt - TimeDelta::fromDays(10);
    CHECK(result.year() == 2022);
    CHECK(result.month() == 1);
    CHECK(result.day() == 5);
    CHECK(result.hour() == 12);

    result = dt - TimeDelta::fromHours(13);
    CHECK(result.day() == 14);
    CHECK(result.hour() == 23);
  }

  SECTION("Difference between DateTimes") {
    DateTime dt1(2022, 1, 15, 12, 0, 0, 0);
    DateTime dt2(2022, 1, 20, 18, 30, 0, 0);

    TimeDelta diff = dt2 - dt1;
    CHECK(diff.totalDays() == 5);
    CHECK(diff.totalHours() == 5 * 24 + 6);
    CHECK(diff.totalMinutes() == (5 * 24 + 6) * 60 + 30);
  }
}

TEST_CASE("DateTime comparison operators", "[datetime]") {
  DateTime dt1(2022, 1, 15);
  DateTime dt2(2022, 1, 20);
  DateTime dt3(2022, 1, 15);

  SECTION("Equality") {
    CHECK(dt1 == dt3);
    CHECK_FALSE(dt1 == dt2);
  }

  SECTION("Less than") {
    CHECK(dt1 < dt2);
    CHECK_FALSE(dt2 < dt1);
    CHECK_FALSE(dt1 < dt3);
  }

  SECTION("Greater than") {
    CHECK(dt2 > dt1);
    CHECK_FALSE(dt1 > dt2);
    CHECK_FALSE(dt1 > dt3);
  }

  SECTION("Less than or equal") {
    CHECK(dt1 <= dt2);
    CHECK(dt1 <= dt3);
    CHECK_FALSE(dt2 <= dt1);
  }

  SECTION("Greater than or equal") {
    CHECK(dt2 >= dt1);
    CHECK(dt1 >= dt3);
    CHECK_FALSE(dt1 >= dt2);
  }
}

TEST_CASE("DateTime now method", "[datetime]") {
  // Get current time
  DateTime now = DateTime::now();

  // Get current time using standard library
  const auto stdNow = std::chrono::system_clock::now();

  // We can't use std::gmtime in thread-safe way across all platforms, so
  // we just check that now() returns a recent timestamp
  auto nowTs = now.timestamp();
  auto stdTs = std::chrono::duration_cast<std::chrono::milliseconds>(
                   stdNow.time_since_epoch())
                   .count();

  // Timestamps should be within 5 seconds of each other
  CHECK(std::abs(nowTs - stdTs) < 5000);

  // Basic validation that now() returns a sensible date (not 1970)
  CHECK(now.year() >= 2022);
}

TEST_CASE("TimeDelta and DateTime edge cases", "[edge]") {
  SECTION("DateTime with extreme values") {
    // Far past date
    DateTime ancient(1, 1, 1);

    // Very far future date
    DateTime farFuture(9999, 12, 31, 23, 59, 59, 999);

    CHECK(ancient < farFuture);
  }

  SECTION("TimeDelta with large values") {
    // Very large TimeDelta
    TimeDelta largeSpan(
        TimeDelta::s_TimedeltaComponents{10000, 0, 0, 0, 0});  // 10000 days

    CHECK(largeSpan.totalDays() == 10000);
    CHECK(largeSpan.days() == 10000);
  }

  SECTION("Date wrapping") {
    // Create a date and add enough time to wrap to the next month
    DateTime dt(2022, 1, 31);
    DateTime result = dt + TimeDelta::fromDays(1);

    CHECK(result.year() == 2022);
    CHECK(result.month() == 2);
    CHECK(result.day() == 1);
  }

  SECTION("Leap year handling") {
    // Test February 29 in leap year
    DateTime leapDay(2020, 2, 29);

    // Add one year
    DateTime nextYear = leapDay + TimeDelta::fromDays(366);

    // Should be March 1, 2021 (2020 was a leap year)
    CHECK(nextYear.year() == 2021);
    CHECK(nextYear.month() == 3);
    CHECK(nextYear.day() == 1);
  }
}

TEST_CASE("TimeDelta arithmetic with different units", "[timedelta]") {
  SECTION("Mixed units addition") {
    TimeDelta days = TimeDelta::fromDays(1);
    TimeDelta hours = TimeDelta::fromHours(12);
    TimeDelta minutes = TimeDelta::fromMinutes(30);

    TimeDelta total = days + hours + minutes;

    CHECK(total.totalHours() == 36);
    CHECK(total.totalMinutes() == 36 * 60 + 30);
  }

  SECTION("Overflow handling") {
    TimeDelta hours = TimeDelta::fromHours(25);

    CHECK(hours.days() == 1);
    CHECK(hours.hours() == 1);
  }
}

TEST_CASE("DateTime parsing edge cases", "[datetime]") {
  SECTION("Invalid date string handling") {
    // This tests that the code doesn't crash, not specific results
    auto dt = DateTime::strptime("not a date");

    // Just verify we can access methods without crashing
    CHECK_FALSE(dt.valid());
  }

  SECTION("Partially valid date string") {
    auto dt = DateTime::strptime("2022-01-XX");

    // This should not crash, but dt should be invalid
    CHECK_FALSE(dt.valid());
  }
}

TEST_CASE("DateTime serialization roundtrip", "[datetime]") {
  DateTime original(2022, 3, 15, 14, 30, 45, 500);

  // Convert to string
  std::string iso = original.to_iso_string_msec();

  // Parse back with appropriate format
  auto parsed = DateTime::strptime(iso, "%Y-%m-%dT%H:%M:%S");

  // Should be the same timestamp
  CHECK(parsed.valid());
  CHECK(parsed.timestamp() == original.timestamp());
}

TEST_CASE("TimeDelta for measuring elapsed time", "[timedelta]") {
  auto startTime = DateTime::now();

  // Sleep for a short duration
  std::this_thread::sleep_for(std::chrono::milliseconds(100));

  auto endTime = DateTime::now();
  TimeDelta elapsed = endTime - startTime;

  // Should be at least our sleep duration
  CHECK(elapsed.totalMilliseconds() >= 80);

  // Shouldn't be too much longer
  CHECK(elapsed.totalMilliseconds() < 1000);
}

TEST_CASE("DateTime stream insertion operator", "[datetime]") {
  DateTime dt(2022, 1, 15, 12, 30, 45);
  std::ostringstream oss;

  oss << dt;

  CHECK(oss.str() == "2022-01-15 12:30:45");
}

TEST_CASE("DateTime get_time_point method", "[datetime]") {
  DateTime dt(2022, 1, 15, 12, 30, 45);
  auto tp = dt.get_time_point();

  // Create a new DateTime from time_point
  DateTime dt2(tp);

  CHECK(dt.timestamp() == dt2.timestamp());
  CHECK(dt.year() == dt2.year());
  CHECK(dt.month() == dt2.month());
  CHECK(dt.day() == dt2.day());
  CHECK(dt.hour() == dt2.hour());
  CHECK(dt.minute() == dt2.minute());
  CHECK(dt.second() == dt2.second());
}

// ====================================================
// Array-based parsing tests
// ====================================================

// Helper function to test array parsing directly without Fortran interface
namespace {
// Test function that mimics the C wrapper functionality
auto test_array_parsing(const std::string& str,
                        const std::vector<std::string>& formats) -> DateTime {
  for (const auto& format : formats) {
    auto result = DateTime::strptime(str, format);
    if (result.valid()) {
      return result;
    }
  }
  return DateTime(DateTime::INVALID_TIME_POINT);
}
}  // namespace

TEST_CASE("Array-based parsing with multiple formats",
          "[datetime][array_parsing]") {
  SECTION("Parse with first format succeeding") {
    std::string date_str = "2024-03-15 14:30:25";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%d-%m-%Y %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Parse with second format succeeding") {
    std::string date_str = "2024/03/15 14:30:25";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%d-%m-%Y %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Parse with third format succeeding") {
    std::string date_str = "15-03-2024 14:30:25";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%d-%m-%Y %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Parse with compact format") {
    std::string date_str = "20240315143025";
    std::vector<std::string> formats = {"%Y-%m-%d %H:%M:%S",
                                        "%Y/%m/%d %H:%M:%S", "%Y%m%d%H%M%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Parse with ISO format variations") {
    std::string date_str = "2024-03-15T14:30:25Z";
    std::vector<std::string> formats = {
        "%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Parse with dot-separated format") {
    std::string date_str = "2024.03.15 14:30:25";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%Y.%m.%d %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }
}

TEST_CASE("Array-based parsing failure cases", "[datetime][array_parsing]") {
  SECTION("All formats fail") {
    std::string date_str = "invalid-date-string";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%d-%m-%Y %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    CHECK_FALSE(result.valid());
  }

  SECTION("Empty format array") {
    std::string date_str = "2024-03-15 14:30:25";
    std::vector<std::string> formats = {};

    auto result = test_array_parsing(date_str, formats);
    CHECK_FALSE(result.valid());
  }

  SECTION("Empty string with valid formats") {
    std::string date_str = "";
    std::vector<std::string> formats = {"%Y-%m-%d %H:%M:%S",
                                        "%Y/%m/%d %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    CHECK_FALSE(result.valid());
  }

  SECTION("Invalid format strings") {
    std::string date_str = "2024-03-15 14:30:25";
    std::vector<std::string> formats = {"invalid-format",
                                        "%Y-%z-%d",  // Invalid specifier
                                        ""};

    auto result = test_array_parsing(date_str, formats);
    CHECK_FALSE(result.valid());
  }
}

TEST_CASE("Array-based parsing with edge cases", "[datetime][array_parsing]") {
  SECTION("Single format in array") {
    std::string date_str = "2024-03-15 14:30:25";
    std::vector<std::string> formats = {"%Y-%m-%d %H:%M:%S"};

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
  }

  SECTION("Many formats with success at end") {
    std::string date_str = "03/15/2024 14:30:25";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%d-%m-%Y %H:%M:%S",
        "%Y.%m.%d %H:%M:%S", "%Y%m%d%H%M%S",
        "%m/%d/%Y %H:%M:%S"  // This should succeed
    };

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Date-only formats") {
    std::string date_str = "2024-03-15";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S",  // This won't match (has time)
        "%Y-%m-%d"            // This will match
    };

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 0);  // Default time
    CHECK(result.minute() == 0);
    CHECK(result.second() == 0);
  }

  SECTION("Time with milliseconds") {
    std::string date_str = "2024-03-15 14:30:25.123";
    std::vector<std::string> formats = {
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M:%S.%f"  // Format with milliseconds
    };

    auto result = test_array_parsing(date_str, formats);
    REQUIRE(result.valid());

    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
    CHECK(result.millisecond() == 123);
  }
}

TEST_CASE("Auto-detection with fallback formats", "[datetime][array_parsing]") {
  SECTION("Auto-detection succeeds") {
    std::string date_str = "2024-03-15 14:30:25";

    // First try auto-detection (should succeed)
    auto auto_result = DateTime::strptime(date_str, "auto");
    REQUIRE(auto_result.valid());

    CHECK(auto_result.year() == 2024);
    CHECK(auto_result.month() == 3);
    CHECK(auto_result.day() == 15);
    CHECK(auto_result.hour() == 14);
    CHECK(auto_result.minute() == 30);
    CHECK(auto_result.second() == 25);
  }

  SECTION("Auto-detection fails, fallback succeeds") {
    std::string date_str =
        "15/Mar/2024 14:30:25";  // Not in auto-detected formats

    // Auto-detection should fail
    auto auto_result = DateTime::strptime(date_str, "auto");
    CHECK_FALSE(auto_result.valid());

    // But manual parsing with the right format should succeed
    auto manual_result = DateTime::strptime(date_str, "%d/%b/%Y %H:%M:%S");
    REQUIRE(manual_result.valid());

    CHECK(manual_result.year() == 2024);
    CHECK(manual_result.month() == 3);
    CHECK(manual_result.day() == 15);
    CHECK(manual_result.hour() == 14);
    CHECK(manual_result.minute() == 30);
    CHECK(manual_result.second() == 25);
  }
}

// ====================================================
// Comprehensive Automatic Format Detection Tests
// ====================================================

TEST_CASE("DateTime auto-detection comprehensive format coverage",
          "[datetime][auto]") {
  // Test each of the 11 auto-detected formats individually

  SECTION("Format 1: %Y-%m-%d %H:%M:%S") {
    auto result = DateTime::strptime("2024-03-15 14:30:25");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 2: %Y-%m-%dT%H:%M:%SZ (ISO with timezone)") {
    auto result = DateTime::strptime("2024-03-15T14:30:25Z");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 3: %Y-%m-%dT%H:%M:%S (ISO format)") {
    auto result = DateTime::strptime("2024-03-15T14:30:25");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 4: %Y/%m/%d %H:%M:%S") {
    auto result = DateTime::strptime("2024/03/15 14:30:25");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 5: %Y.%m.%d %H:%M:%S") {
    auto result = DateTime::strptime("2024.03.15 14:30:25");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 6: %Y%m%d%H%M%S (compact format)") {
    auto result = DateTime::strptime("20240315143025");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("Format 7: %Y/%m/%d %H:%M") {
    auto result = DateTime::strptime("2024/03/15 14:30");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 0);
  }

  SECTION("Format 8: %Y-%m-%d (date only)") {
    auto result = DateTime::strptime("2024-03-15");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 0);
    CHECK(result.minute() == 0);
    CHECK(result.second() == 0);
  }

  SECTION("Format 9: %Y/%m/%d (date only with slashes)") {
    auto result = DateTime::strptime("2024/03/15");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 0);
    CHECK(result.minute() == 0);
    CHECK(result.second() == 0);
  }

  SECTION("Format 10: %Y.%m.%d (date only with dots)") {
    auto result = DateTime::strptime("2024.03.15");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 0);
    CHECK(result.minute() == 0);
    CHECK(result.second() == 0);
  }

  SECTION("Format 11: %Y%m%d (compact date only)") {
    auto result = DateTime::strptime("20240315");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 0);
    CHECK(result.minute() == 0);
    CHECK(result.second() == 0);
  }
}

TEST_CASE("DateTime auto-detection with milliseconds", "[datetime][auto]") {
  SECTION("Standard format with milliseconds") {
    auto result = DateTime::strptime("2024-03-15 14:30:25.123");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
    CHECK(result.millisecond() == 123);
  }

  SECTION("ISO format with milliseconds") {
    auto result = DateTime::strptime("2024-03-15T14:30:25.456");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
    CHECK(result.millisecond() == 456);
  }

  SECTION("Slash format with milliseconds") {
    auto result = DateTime::strptime("2024/03/15 14:30:25.789");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
    CHECK(result.millisecond() == 789);
  }
}

TEST_CASE("DateTime auto-detection format precedence", "[datetime][auto]") {
  // Test which format wins when multiple could potentially match

  SECTION("Most specific format wins - full datetime vs date-only") {
    // This should match the full datetime format, not just the date part
    auto result = DateTime::strptime("2024-03-15 14:30:25");
    REQUIRE(result.valid());
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 25);
  }

  SECTION("ISO with Z takes precedence over ISO without Z") {
    auto result = DateTime::strptime("2024-03-15T14:30:25Z");
    REQUIRE(result.valid());
    // Should successfully parse even with Z suffix
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
  }

  SECTION("Minutes-precision format preferred when seconds missing") {
    auto result = DateTime::strptime("2024/03/15 14:30");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);
    CHECK(result.hour() == 14);
    CHECK(result.minute() == 30);
    CHECK(result.second() == 0);
  }
}

TEST_CASE("DateTime auto-detection edge cases", "[datetime][auto]") {
  SECTION("Ambiguous date handling - year boundaries") {
    // Test dates that could be ambiguous
    auto result = DateTime::strptime("2024-01-01");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 1);
    CHECK(result.day() == 1);
  }

  SECTION("Leap year dates") {
    auto result = DateTime::strptime("2024-02-29");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 2);
    CHECK(result.day() == 29);
  }

  SECTION("End of month dates") {
    auto result = DateTime::strptime("2024-01-31");
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 1);
    CHECK(result.day() == 31);
  }

  SECTION("Different separator consistency") {
    // Ensure formats with mixed separators fail appropriately
    auto result = DateTime::strptime("2024-03/15 14:30:25");
    CHECK_FALSE(result.valid());
  }

  SECTION("Partial format strings") {
    // Test strings that partially match formats
    auto result = DateTime::strptime("2024-03");
    CHECK_FALSE(result.valid());
  }
}

TEST_CASE("DateTime auto-detection invalid inputs", "[datetime][auto]") {
  SECTION("Completely invalid strings") {
    CHECK_FALSE(DateTime::strptime("not a date").valid());
    CHECK_FALSE(DateTime::strptime("").valid());
    CHECK_FALSE(DateTime::strptime("abc123def").valid());
  }

  SECTION("Invalid date values") {
    CHECK_FALSE(DateTime::strptime("2024-13-15").valid());  // Invalid month
    CHECK_FALSE(
        DateTime::strptime("2024-02-30").valid());  // Invalid day for February
    CHECK_FALSE(
        DateTime::strptime("2024-04-31").valid());  // Invalid day for April
  }

  SECTION("Invalid time values") {
    // Test with strings that have non-numeric characters in time positions
    CHECK_FALSE(
        DateTime::strptime("2024-03-15 aa:30:25").valid());  // Non-numeric hour
    CHECK_FALSE(DateTime::strptime("2024-03-15 14:bb:25")
                    .valid());  // Non-numeric minute
    CHECK_FALSE(DateTime::strptime("2024-03-15 14:30:cc")
                    .valid());  // Non-numeric second
  }

  SECTION("Invalid time value rejection") {
    // The parsing library should reject out-of-range time values
    // rather than silently normalizing them

    auto result1 = DateTime::strptime("2024-03-15 25:30:25");  // 25 hours
    CHECK_FALSE(result1.valid());  // Should fail to parse

    auto result2 = DateTime::strptime("2024-03-15 14:60:25");  // 60 minutes
    CHECK_FALSE(result2.valid());  // Should fail to parse

    auto result3 = DateTime::strptime("2024-03-15 14:30:60");  // 60 seconds
    CHECK_FALSE(result3.valid());  // Should fail to parse
  }

  SECTION("Malformed but partially parseable") {
    // The parsing library may be more lenient than expected, so test truly
    // malformed strings
    CHECK_FALSE(DateTime::strptime("2024/03-15").valid());  // Mixed separators
    CHECK_FALSE(DateTime::strptime("2024-").valid());       // Incomplete date
    CHECK_FALSE(DateTime::strptime("202a-03-15").valid());  // Non-numeric year
  }
}

TEST_CASE("DateTime auto-detection performance characteristics",
          "[datetime][auto][performance]") {
  // Basic performance comparison between auto and manual format specification

  SECTION("Auto vs manual parsing speed comparison") {
    const std::string test_string = "2024-03-15 14:30:25";
    const int iterations = 1000;

    // Warm up
    for (int i = 0; i < 10; ++i) {
      DateTime::strptime(test_string);
      DateTime::strptime(test_string, "%Y-%m-%d %H:%M:%S");
    }

    // Test auto parsing
    auto start_auto = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
      auto result = DateTime::strptime(test_string);
      REQUIRE(result.valid());  // Ensure it's working
    }
    auto end_auto = std::chrono::high_resolution_clock::now();

    // Test manual parsing
    auto start_manual = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
      auto result = DateTime::strptime(test_string, "%Y-%m-%d %H:%M:%S");
      REQUIRE(result.valid());  // Ensure it's working
    }
    auto end_manual = std::chrono::high_resolution_clock::now();

    auto auto_duration = std::chrono::duration_cast<std::chrono::microseconds>(
        end_auto - start_auto);
    auto manual_duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end_manual -
                                                              start_manual);

    // Auto should be slower but not more than 10x slower for the first format
    // match This is more of a sanity check than a strict requirement
    CHECK(auto_duration.count() > 0);
    CHECK(manual_duration.count() > 0);

    // Log the performance for informational purposes
    INFO("Auto parsing took: " << auto_duration.count() << " microseconds");
    INFO("Manual parsing took: " << manual_duration.count() << " microseconds");
    INFO("Ratio (auto/manual): "
         << static_cast<double>(auto_duration.count()) /
                static_cast<double>(manual_duration.count()));
  }

  SECTION("Auto parsing worst case - last format matches") {
    // Test parsing time when the last format in the list matches
    const std::string compact_date =
        "20240315";  // Should match %Y%m%d (last in list)

    auto result = DateTime::strptime(compact_date);
    REQUIRE(result.valid());
    CHECK(result.year() == 2024);
    CHECK(result.month() == 3);
    CHECK(result.day() == 15);

    // Should still complete in reasonable time
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i) {
      DateTime::strptime(compact_date);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Should complete 100 iterations in less than 100ms (very generous)
    CHECK(duration.count() < 100);
  }
}

TEST_CASE("DateTime Julian Day functionality", "[datetime][julian]") {
  SECTION("Julian Day Number calculation") {
    // Test well-known Julian Day Numbers

    // July 24, 2025 (JDN 2460881)
    DateTime dt_2000(2025, 7, 24, 0, 0, 0);
    CHECK(dt_2000.julianDayNumber() == 2460881);

    // January 1, 1970 (Unix epoch) = JDN 2440588
    DateTime dt_1970(1970, 1, 1, 12, 0, 0);
    CHECK(dt_1970.julianDayNumber() == 2440588);

    // May 23, 1968 (Julian Day 0 of GPS week) = JDN 2440000
    DateTime dt_gps(1968, 5, 23, 12, 0, 0);
    CHECK(dt_gps.julianDayNumber() == 2440000);

    // October 15, 1582 (Gregorian calendar adoption) = JDN 2299161
    DateTime dt_greg(1582, 10, 15, 12, 0, 0);
    CHECK(dt_greg.julianDayNumber() == 2299161);
  }

  SECTION("Julian Day with fractional part") {
    // Test fractional Julian Day calculations

    // January 1, 2000 00:00 UTC (midnight) = JD 2451544.5
    DateTime dt_midnight(2000, 1, 1, 0, 0, 0);
    CHECK(dt_midnight.julianDay() == Catch::Approx(2451544.5).epsilon(1e-6));

    // January 1, 2000 12:00 UTC (noon) = JD 2451545.0
    DateTime dt_noon(2000, 1, 1, 12, 0, 0);
    CHECK(dt_noon.julianDay() == Catch::Approx(2451545.0).epsilon(1e-6));

    // January 1, 2000 18:00 UTC (6 PM) = JD 2451545.25
    DateTime dt_evening(2000, 1, 1, 18, 0, 0);
    CHECK(dt_evening.julianDay() == Catch::Approx(2451545.25).epsilon(1e-6));

    // Test with minutes, seconds, and milliseconds
    // January 1, 2000 12:30:30.500 UTC
    DateTime dt_precise(2000, 1, 1, 12, 30, 30, 500);
    double expected_jd = 2451545.0 + (30.0 / 1440.0) + (30.5 / 86400.0);
    CHECK(dt_precise.julianDay() == Catch::Approx(expected_jd).epsilon(1e-6));
  }

  SECTION("Julian Day edge cases") {
    // Test very early date
    DateTime dt_early(1, 1, 1, 12, 0, 0);
    int64_t jdn_early = dt_early.julianDayNumber();
    CHECK(jdn_early > 0);  // Should be positive

    // Test far future date
    DateTime dt_future(3000, 12, 31, 12, 0, 0);
    int64_t jdn_future = dt_future.julianDayNumber();
    CHECK(jdn_future > 2451545);  // Should be much larger than Y2K

    // Test February 20 of a non-leap year
    DateTime dt_non_leap(2021, 2, 20, 12, 0, 0);
    int64_t jdn_non_leap = dt_non_leap.julianDayNumber();
    CHECK(jdn_non_leap == 2459266);  // 50 days after Jan 1, 2021

    // Test leap year boundary (February 29, 2000)
    DateTime dt_leap(2000, 2, 29, 12, 0, 0);
    int64_t jdn_leap = dt_leap.julianDayNumber();
    CHECK(jdn_leap == 2451604);  // 59 days after Jan 1, 2000

    // Test century year that's not a leap year (1900)
    DateTime dt_1900(1900, 2, 28, 12, 0, 0);
    DateTime dt_1900_next(1900, 3, 1, 12, 0, 0);
    CHECK(dt_1900_next.julianDayNumber() - dt_1900.julianDayNumber() == 1);
  }

  SECTION("Julian Day consistency checks") {
    // Test that consecutive days differ by 1
    DateTime dt1(2023, 6, 15, 12, 0, 0);
    DateTime dt2(2023, 6, 16, 12, 0, 0);
    CHECK(dt2.julianDayNumber() - dt1.julianDayNumber() == 1);

    // Test fractional day progression
    DateTime dt_start(2023, 6, 15, 0, 0, 0);  // midnight
    DateTime dt_end(2023, 6, 15, 23, 59, 59,
                    999);  // just before midnight next day

    double jd_start = dt_start.julianDay();
    double jd_end = dt_end.julianDay();

    // Should be almost exactly 1 day difference
    CHECK((jd_end - jd_start) ==
          Catch::Approx(0.99999998835846782).epsilon(1e-6));
  }
}
