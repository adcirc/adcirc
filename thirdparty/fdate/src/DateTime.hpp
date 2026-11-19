/*
 * FDate - A Fortran Date and Time Library based on C++
 * Copyright (C) 2025 Zach Cobell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <limits>
#include <sstream>
#include <string>

#include "TimeDelta.hpp"
#include "date_hh.h"

/**
 * @brief A high-precision DateTime class with millisecond accuracy
 *
 * The DateTime class provides a comprehensive interface for date and time
 * operations with millisecond precision. It supports parsing from multiple
 * string formats, formatting to various output formats, and arithmetic
 * operations with TimeDelta objects.
 *
 * The class uses std::chrono internally for time calculations and the Howard
 * Hinnant date library for parsing and formatting operations.
 *
 * @note All timestamps are stored as milliseconds since the Unix epoch
 * (1970-01-01 00:00:00 UTC)
 * @see TimeDelta for time duration operations
 *
 * @note Thread Safety: All DateTime operations are thread-safe for read-only
 * operations. Parsing operations are thread-safe as they don't modify the
 * global state. Error handling writes to stderr which may interleave in
 * multithreaded environments.
 */
class DateTime {
  /** @brief Internal time point type with millisecond precision */
  using t_time_point = std::chrono::time_point<std::chrono::system_clock,
                                               std::chrono::milliseconds>;

  /** @brief Internal time point storage */
  t_time_point m_tp;

  /**
   * @brief Parses a DateTime from a string using a specific format
   *
   * Uses the Howard Hinnant date library for parsing. Automatically detects
   * and handles millisecond precision based on the presence of a decimal point
   * in the appropriate position within the input string.
   *
   * @param str The string to parse
   * @param format The strftime-style format string (default: "%Y-%m-%d
   * %H:%M:%S")
   * @return DateTime containing the parsed DateTime, or an invalid DateTime
   * if parsing failed (check with .valid())
   *
   * @note If the string has a '.' at position (length-4), millisecond parsing
   * is attempted
   * @see parse() for automatic format detection
   */
  static auto parse_string(const std::string& str,
                           const std::string& format = "%Y-%m-%d %H:%M:%S")
      -> DateTime {
    std::istringstream input_stream(str);
    date::sys_time<std::chrono::milliseconds> this_time_point;

    // If the user has provided a string that has milliseconds (i.e., position 4
    // from the end has a period)
    if (str.size() >= 4 && str[str.size() - 4] == '.') {
      // Parse using milliseconds
      date::from_stream(input_stream, format.c_str(), this_time_point);
    } else {
      // Parse without milliseconds
      auto tp_temp = date::sys_time<std::chrono::seconds>{};
      date::from_stream(input_stream, format.c_str(), tp_temp);
      this_time_point =
          std::chrono::time_point_cast<std::chrono::milliseconds>(tp_temp);
    }

    // Check if parsing was successful
    if (input_stream.fail() || input_stream.bad()) {
      return DateTime(INVALID_TIME_POINT);  // Parsing failed
    }

    // Additional validation: ensure the entire string was consumed
    // This prevents partial matches where invalid time components are ignored
    if (char remaining_char; input_stream >> remaining_char) {
      // There are unconsumed non-whitespace characters, which suggests
      // the format didn't match the full input string
      return DateTime(INVALID_TIME_POINT);
    }

    return DateTime(this_time_point);
  }

 public:
  /** @brief Constant representing an invalid timestamp value */
  static constexpr auto INVALID_TIMESTAMP =
      -std::numeric_limits<int64_t>::max();

  /** @brief Constant representing an invalid DateTime object */
  static constexpr auto INVALID_TIME_POINT =
      t_time_point(std::chrono::milliseconds(INVALID_TIMESTAMP));

  static constexpr int DATETIME_MIN_YEAR = 0;
  static constexpr int DATETIME_MIN_MONTH = 1;
  static constexpr int DATETIME_MAX_MONTH = 12;
  static constexpr int DATETIME_MAX_DAYS = 31;
  static constexpr int DATETIME_MIN_DAYS = 1;
  static constexpr int DATETIME_MAX_HOURS = 23;
  static constexpr int DATETIME_MIN_HOURS = 0;
  static constexpr int DATETIME_MAX_MINUTES = 59;
  static constexpr int DATETIME_MIN_MINUTES = 0;
  static constexpr int DATETIME_MAX_SECONDS = 59;
  static constexpr int DATETIME_MIN_SECONDS = 0;
  static constexpr int DATETIME_MAX_MILLISECONDS = 999;
  static constexpr int DATETIME_MIN_MILLISECONDS = 0;

  // Constructors

  /**
   * @brief Default constructor creating a DateTime at Unix epoch (1970-01-01
   * 00:00:00.000)
   */
  constexpr DateTime() noexcept : m_tp(std::chrono::milliseconds(0)) {};

  /**
   * @brief Default destructor
   */
  ~DateTime() = default;

  /**
   * @brief Constructs a DateTime from individual date and time components
   *
   * @param year The year (e.g., 2023)
   * @param month The month (1-12)
   * @param day The day of month (1-31)
   * @param hour The hour (0-23, default: 0)
   * @param minute The minute (0-59, default: 0)
   * @param second The second (0-59, default: 0)
   * @param millisecond The millisecond (0-999, default: 0)
   *
   * @note No validation is performed on input values. Invalid dates may produce
   * unexpected results.
   */
  constexpr DateTime(int year, unsigned month, unsigned day, unsigned hour = 0,
                     unsigned minute = 0, unsigned second = 0,
                     unsigned millisecond = 0) noexcept
      : m_tp(std::chrono::time_point_cast<std::chrono::milliseconds>(
            date::sys_days{date::year{year} / date::month{month} /
                           date::day{day}} +
            std::chrono::hours{hour} + std::chrono::minutes{minute} +
            std::chrono::seconds{second} +
            std::chrono::milliseconds{millisecond})) {
#ifndef NDEBUG
    assert(month >= DATETIME_MIN_MONTH && month <= DATETIME_MAX_MONTH);
    assert(day >= DATETIME_MIN_DAYS && day <= DATETIME_MAX_DAYS);
    assert(hour <= DATETIME_MAX_HOURS);
    assert(minute <= DATETIME_MAX_MINUTES);
    assert(second <= DATETIME_MAX_SECONDS);
    assert(millisecond <= DATETIME_MAX_MILLISECONDS);
#endif
  }

  /**
   * @brief Constructs a DateTime from a timestamp
   *
   * @param timestamp Milliseconds since Unix epoch (1970-01-01 00:00:00.000
   * UTC)
   */
  explicit constexpr DateTime(int64_t timestamp) noexcept
      : m_tp(t_time_point(std::chrono::milliseconds(timestamp))) {}

  /**
   * @brief Constructs a DateTime from a time_point
   *
   * @param milliseconds A std::chrono::time_point with millisecond precision
   */
  constexpr explicit DateTime(const t_time_point& milliseconds) noexcept
      : m_tp(milliseconds) {}

  /**
   * @brief Parses a DateTime from a string using automatic or specified format
   * detection
   *
   * When format is "auto", the function attempts to parse the string using a
   * predefined list of common date/time formats. All auto-detected formats use
   * YYYY-MM-DD date component ordering to avoid ambiguity.
   *
   * Supported auto-detection formats:
   * - "%Y-%m-%d %H:%M:%S" (2023-12-01 14:30:25)
   * - "%Y-%m-%dT%H:%M:%SZ" (2023-12-01T14:30:25Z - ISO with timezone)
   * - "%Y-%m-%dT%H:%M:%S" (2023-12-01T14:30:25 - ISO format)
   * - "%Y/%m/%d %H:%M:%S" (2023/12/01 14:30:25)
   * - "%Y.%m.%d %H:%M:%S" (2023.12.01 14:30:25)
   * - "%Y%m%d%H%M%S" (20231201143025 - compact format)
   * - "%Y/%m/%d %H:%M" (2023/12/01 14:30)
   * - "%Y-%m-%d" (2023-12-01)
   * - "%Y/%m/%d" (2023/12/01)
   * - "%Y.%m.%d" (2023.12.01)
   * - "%Y%m%d" (20231201 - compact date only)
   *
   * @param str The string to parse
   * @param format The format string to use, or "auto" for automatic detection
   * (default: "auto")
   * @return DateTime containing the parsed DateTime, or an invalid DateTime
   * if parsing failed (check with .valid())
   *
   * @note Formats are tried in order from most specific to least specific
   * @see parse_string() for single format parsing
   */
  static auto strptime(const std::string& str,
                       const std::string& format = "auto") -> DateTime {
    if (format == "auto") {
      // All formats use YYYY-MM-DD ordering to avoid ambiguity
      // Ordered from the most specific to the least specific
      static const std::array<std::string, 11> format_options = {
          "%Y-%m-%d %H:%M:%S",   // 2023-12-01 14:30:25
          "%Y-%m-%dT%H:%M:%SZ",  // 2023-12-01T14:30:25Z (ISO with timezone)
          "%Y-%m-%dT%H:%M:%S",   // 2023-12-01T14:30:25 (ISO)
          "%Y/%m/%d %H:%M:%S",   // 2023/12/01 14:30:25
          "%Y.%m.%d %H:%M:%S",   // 2023.12.01 14:30:25
          "%Y%m%d%H%M%S",        // 20231201143025
          "%Y/%m/%d %H:%M",      // 2023/12/01 14:30
          "%Y-%m-%d",            // 2023-12-01
          "%Y/%m/%d",            // 2023/12/01
          "%Y.%m.%d",            // 2023.12.01
          "%Y%m%d"               // 20231201
      };

      for (const auto& fmt : format_options) {
        const auto result = parse_string(str, fmt);
        if (result.valid()) {
          return result;
        }
      }
      return DateTime(INVALID_TIME_POINT);
    } else {
      return parse_string(str, format);
    }
  }

  /**
   * @brief Gets the timestamp as milliseconds since Unix epoch
   *
   * @return int64_t The number of milliseconds since 1970-01-01 00:00:00.000
   * UTC
   */
  [[nodiscard]] constexpr auto timestamp() const noexcept -> int64_t {
    return m_tp.time_since_epoch().count();
  }

  /**
   * @brief Formats the DateTime as a string without milliseconds
   *
   * @param fmt The strftime-style format string (default: "%Y-%m-%d %H:%M:%S")
   * @return std::string The formatted date/time string
   *
   * @note This method truncates to seconds precision
   * @see format_w_milliseconds() for millisecond precision formatting
   */
  [[nodiscard]] auto strftime(
      const std::string& fmt = "%Y-%m-%d %H:%M:%S") const -> std::string {
    auto tp_sec = std::chrono::time_point_cast<std::chrono::seconds>(m_tp);
    std::ostringstream oss;
    date::to_stream(oss, fmt.c_str(), tp_sec);
    return oss.str();
  }

  /**
   * @brief Formats the DateTime as a string with millisecond precision
   *
   * @param fmt The strftime-style format string (default: "%Y-%m-%d %H:%M:%S")
   * @return std::string The formatted date/time string including milliseconds
   *
   * @see format() for seconds precision formatting
   */
  [[nodiscard]] auto strftime_w_milliseconds(
      const std::string& fmt = "%Y-%m-%d %H:%M:%S") const -> std::string {
    std::ostringstream oss;
    date::to_stream(oss, fmt.c_str(), m_tp);
    return oss.str();
  }

  /**
   * @brief Converts to ISO 8601 format string without milliseconds
   *
   * @return std::string The DateTime in format "YYYY-MM-DDTHH:MM:SS"
   *
   * @see toISOStringMsec() for ISO format with milliseconds
   */
  [[nodiscard]] auto to_iso_string() const -> std::string {
    return strftime("%Y-%m-%dT%H:%M:%S");
  }

  /**
   * @brief Converts to ISO 8601 format string with milliseconds
   *
   * @return std::string The DateTime in format "YYYY-MM-DDTHH:MM:SS.mmm"
   *
   * @see toISOString() for ISO format without milliseconds
   */
  [[nodiscard]] auto to_iso_string_msec() const -> std::string {
    return strftime_w_milliseconds("%Y-%m-%dT%H:%M:%S");
  }

  // Getters

  /**
   * @brief Gets the year component
   *
   * @return int64_t The year (e.g., 2023)
   */
  [[nodiscard]] constexpr auto year() const noexcept -> int64_t {
    return static_cast<int>(
        date::year_month_day{date::floor<date::days>(m_tp)}.year());
  }

  /**
   * @brief Gets the month component
   *
   * @return unsigned The month (1-12, where 1 = January)
   */
  [[nodiscard]] constexpr auto month() const noexcept -> unsigned {
    return static_cast<unsigned>(
        date::year_month_day{date::floor<date::days>(m_tp)}.month());
  }

  /**
   * @brief Gets the day of month component
   *
   * @return unsigned The day (1-31)
   */
  [[nodiscard]] constexpr auto day() const noexcept -> unsigned {
    return static_cast<unsigned>(
        date::year_month_day{date::floor<date::days>(m_tp)}.day());
  }

  /**
   * @brief Gets the hour component
   *
   * @return unsigned The hour (0-23)
   */
  [[nodiscard]] constexpr auto hour() const noexcept -> unsigned {
    return static_cast<unsigned>((m_tp - date::floor<date::days>(m_tp)) /
                                 std::chrono::hours(1));
  }

  /**
   * @brief Gets the minute component
   *
   * @return unsigned The minute (0-59)
   */
  [[nodiscard]] constexpr auto minute() const noexcept -> unsigned {
    const auto time_since_midnight = m_tp - date::floor<date::days>(m_tp);
    return static_cast<unsigned>(
        (time_since_midnight - std::chrono::hours(hour())) /
        std::chrono::minutes(1));
  }

  /**
   * @brief Gets the second component
   *
   * @return unsigned The second (0-59)
   */
  [[nodiscard]] constexpr auto second() const noexcept -> unsigned {
    const auto time_since_midnight = m_tp - date::floor<date::days>(m_tp);
    return static_cast<unsigned>((time_since_midnight -
                                  std::chrono::hours(hour()) -
                                  std::chrono::minutes(minute())) /
                                 std::chrono::seconds(1));
  }

  /**
   * @brief Gets the millisecond component
   *
   * @return unsigned The millisecond (0-999)
   */
  [[nodiscard]] constexpr auto millisecond() const noexcept -> unsigned {
    const auto time_since_midnight = m_tp - date::floor<date::days>(m_tp);
    return static_cast<unsigned>(
        (time_since_midnight - std::chrono::hours(hour()) -
         std::chrono::minutes(minute()) - std::chrono::seconds(second())) /
        std::chrono::milliseconds(1));
  }

  /**
   * @brief Gets the Julian Day Number (JDN) for this date.
   *
   * The Julian Day Number is the integer number of days that have elapsed
   * since the beginning of the Julian Period on January 1, 4713 BCE
   * (proleptic Julian calendar) which corresponds to November 24, 4714 BCE
   * in the proleptic Gregorian calendar.
   *
   * @return int64_t The Julian Day Number (integer days since JD epoch)
   */
  // cppcheck-suppress functionStatic
  [[nodiscard]] constexpr auto julianDayNumber() const noexcept -> int64_t {
    const auto y = this->year();
    const auto m = static_cast<int64_t>(this->month());
    const auto d = static_cast<int64_t>(this->day());

    // Julian Day Number formula for Gregorian calendar
    const auto a = (14 - m) / 12;
    const auto y_adj = y + 4800 - a;
    const auto m_adj = m + 12 * a - 3;

    return d + (153 * m_adj + 2) / 5 + 365 * y_adj + y_adj / 4 - y_adj / 100 +
           y_adj / 400 - 32045;
  }

  /**
   * @brief Gets the Julian Day (JD) with fractional day for this date and time
   *
   * The Julian Day is the continuous count of days since the beginning of
   * the Julian Period, including the fractional part representing the time
   * of day. Julian Day 0.0 begins at noon (12:00 UTC) on January 1, 4713 BCE.
   *
   * @return double The Julian Day with fractional part (days.fraction since JD
   * epoch)
   */
  [[nodiscard]] constexpr auto julianDay() const noexcept -> double {
    const auto jdn = static_cast<double>(julianDayNumber());
    const auto h = static_cast<double>(hour());
    const auto min = static_cast<double>(minute());
    const auto s = static_cast<double>(second());
    const auto ms = static_cast<double>(millisecond());

    // Convert time to fractional day (Julian Day starts at noon, so subtract 12
    // hours)
    const auto fractional_day =
        (h - 12.0) / 24.0 + min / 1440.0 + s / 86400.0 + ms / 86400000.0;

    return jdn + fractional_day;
  }

  /**
   *  @brief Computes the Julian century (JC) for this DateTime
   *
   *  The Julian century is a time unit used in astronomy, defined as
   *  100 Julian years, where one Julian year is exactly 365.25 days.
   *  It is used to measure time intervals in astronomical calculations.
   */
  [[nodiscard]] constexpr auto julianCentury() const noexcept -> double {
    // Julian Day 2451545.0 corresponds to J2000.0 epoch
    const auto jd = julianDay();
    return (jd - 2451545.0) / 36525.0;  // Convert to Julian centuries
  }

  /**
   * @brief Gets the internal time_point representation
   *
   * @return time_point The internal std::chrono::time_point with millisecond
   * precision
   *
   * @note This method provides direct access to the internal representation
   */
  [[nodiscard]] constexpr auto get_time_point() const noexcept -> t_time_point {
    return m_tp;
  }

  // TimeDelta operations

  /**
   * @brief Adds a TimeDelta to this DateTime
   *
   * @param span The TimeDelta to add
   * @return DateTime A new DateTime representing this time plus the span
   *
   * @see operator-(const TimeDelta&) for subtraction
   */
  constexpr auto operator+(const TimeDelta& span) const noexcept -> DateTime {
    return DateTime(m_tp + span.duration());
  }

  /**
   * @brief Subtracts a TimeDelta from this DateTime
   *
   * @param span The TimeDelta to subtract
   * @return DateTime A new DateTime representing this time minus the span
   *
   * @see operator+(const TimeDelta&) for addition
   */
  constexpr auto operator-(const TimeDelta& span) const noexcept -> DateTime {
    return DateTime(m_tp - span.duration());
  }

  /**
   * @brief Calculates the time difference between two DateTimes
   *
   * @param other The DateTime to subtract from this one
   * @return TimeDelta The time difference (this - other)
   *
   * @note If other is later than this DateTime, the result will be negative
   */
  constexpr auto operator-(const DateTime& other) const noexcept -> TimeDelta {
    const auto dt_diff = m_tp - other.m_tp;
    return TimeDelta::fromMilliseconds(dt_diff.count());
  }

  // Comparison operators

  /**
   * @brief Equality comparison operator
   *
   * @param other The DateTime to compare with
   * @return bool True if both DateTimes represent the exact same point in time
   */
  constexpr auto operator==(const DateTime& other) const noexcept -> bool {
    return m_tp == other.m_tp;
  }

  /**
   * @brief Less than comparison operator
   *
   * @param other The DateTime to compare with
   * @return bool True if this DateTime is earlier than other
   */
  constexpr auto operator<(const DateTime& other) const noexcept -> bool {
    return m_tp < other.m_tp;
  }

  /**
   * @brief Greater than comparison operator
   *
   * @param other The DateTime to compare with
   * @return bool True if this DateTime is later than other
   */
  constexpr auto operator>(const DateTime& other) const noexcept -> bool {
    return m_tp > other.m_tp;
  }

  /**
   * @brief Less than or equal comparison operator
   *
   * @param other The DateTime to compare with
   * @return bool True if this DateTime is earlier than or equal to other
   */
  constexpr auto operator<=(const DateTime& other) const noexcept -> bool {
    return m_tp <= other.m_tp;
  }

  /**
   * @brief Greater than or equal comparison operator
   *
   * @param other The DateTime to compare with
   * @return bool True if this DateTime is later than or equal to other
   */
  constexpr auto operator>=(const DateTime& other) const noexcept -> bool {
    return m_tp >= other.m_tp;
  }

  /**
   * @brief Stream insertion operator for easy printing
   *
   * @param output_stream The output stream to write to
   * @param date_time The DateTime to output
   * @return std::ostream& Reference to the output stream for chaining
   */
  friend auto operator<<(std::ostream& output_stream, const DateTime& date_time)
      -> std::ostream& {
    return output_stream << date_time.strftime();
  }

  /**
   * @brief Creates a DateTime representing the current system time
   *
   * @return DateTime A new DateTime object set to the current system time with
   * millisecond precision
   *
   * @note Uses std::chrono::system_clock::now() internally
   */
  [[nodiscard]] static auto now() -> DateTime {
    const auto now_tp = std::chrono::time_point_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now());
    return DateTime(now_tp);
  }

  /**
   * @brief Checks if the DateTime is valid (not equal to INVALID_TIME_POINT)
   *
   * @return bool True if the DateTime is valid, false otherwise
   */
  [[nodiscard]] constexpr auto valid() const noexcept -> bool {
    return m_tp != INVALID_TIME_POINT;
  }
};
