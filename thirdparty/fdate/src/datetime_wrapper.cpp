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

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

#include "DateTime.hpp"

extern "C" {

//=============================================================================
// TimeDelta functions
//=============================================================================

/**
 * @brief Create a TimeDelta from days, hours, minutes, seconds, and
 * milliseconds
 *
 * @param days Number of days
 * @param hours Number of hours
 * @param minutes Number of minutes
 * @param seconds Number of seconds
 * @param milliseconds Number of milliseconds
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_create(int days, int hours, int minutes, int seconds,
                        int milliseconds) -> int64_t {
  const TimeDelta time_delta(days, hours, minutes, seconds, milliseconds);
  return time_delta.totalMilliseconds();
}

/**
 * @brief Create a TimeDelta from a number of days
 *
 * @param days Number of days
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_from_days(const int64_t days) -> int64_t {
  return TimeDelta::fromDays(days).totalMilliseconds();
}

/**
 * @brief Create a TimeDelta from a number of hours
 *
 * @param hours Number of hours
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_from_hours(const int64_t hours) -> int64_t {
  return TimeDelta::fromHours(hours).totalMilliseconds();
}

/**
 * @brief Create a TimeDelta from a number of minutes
 *
 * @param minutes Number of minutes
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_from_minutes(const int64_t minutes) -> int64_t {
  return TimeDelta::fromMinutes(minutes).totalMilliseconds();
}

/**
 * @brief Create a TimeDelta from a number of seconds
 *
 * @param seconds Number of seconds
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_from_seconds(const int64_t seconds) -> int64_t {
  return TimeDelta::fromSeconds(seconds).totalMilliseconds();
}

/**
 * @brief Create a TimeDelta from a number of milliseconds
 *
 * @param milliseconds Number of milliseconds
 * @return int64_t Total milliseconds representation of the TimeDelta
 */
auto f_timedelta_from_milliseconds(const int64_t milliseconds) -> int64_t {
  return TimeDelta::fromMilliseconds(milliseconds).totalMilliseconds();
}

/**
 * @brief Get the days component from a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int Days component
 */
auto f_timedelta_get_days(const int64_t ts_ms) -> int64_t {
  return TimeDelta::to_components(ts_ms).days;
}

/**
 * @brief Get the hours component from a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int Hours component
 */
auto f_timedelta_get_hours(const int64_t ts_ms) -> int64_t {
  return TimeDelta::to_components(ts_ms).hours;
  ;
}

/**
 * @brief Get the minute component from a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int Minutes component
 */
auto f_timedelta_get_minutes(const int64_t ts_ms) -> int64_t {
  return TimeDelta::to_components(ts_ms).minutes;
}

/**
 * @brief Get the second component from a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int Seconds component
 */
auto f_timedelta_get_seconds(const int64_t ts_ms) -> int64_t {
  return TimeDelta::to_components(ts_ms).seconds;
}

/**
 * @brief Get the millisecond component from a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int Milliseconds component
 */
auto f_timedelta_get_milliseconds(const int64_t ts_ms) -> int64_t {
  return TimeDelta::to_components(ts_ms).milliseconds;
}

/**
 * @brief Get the total days representation of a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Total days
 */
auto f_timedelta_get_total_days(const int64_t ts_ms) -> int64_t {
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return time_delta.totalDays();
}

/**
 * @brief Get the total hours representation of a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Total hours
 */
auto f_timedelta_get_total_hours(const int64_t ts_ms) -> int64_t {
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return time_delta.totalHours();
}

/**
 * @brief Get the total minutes representation of a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Total minutes
 */
auto f_timedelta_get_total_minutes(const int64_t ts_ms) -> int64_t {
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return time_delta.totalMinutes();
}

/**
 * @brief Get the total seconds representation of a TimeDelta
 *
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Total seconds
 */
auto f_timedelta_get_total_seconds(const int64_t ts_ms) -> int64_t {
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return time_delta.totalSeconds();
}

/**
 * @brief Add two TimeDeltas together
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return int64_t Resulting TimeDelta as milliseconds
 */
auto f_timedelta_add(const int64_t ts1_ms, const int64_t ts2_ms) -> int64_t {
  return ts1_ms + ts2_ms;
}

/**
 * @brief Subtract one TimeDelta from another
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return int64_t Resulting TimeDelta as milliseconds
 */
auto f_timedelta_subtract(const int64_t ts1_ms, const int64_t ts2_ms)
    -> int64_t {
  return ts1_ms - ts2_ms;
}

/**
 * @brief Multiply a TimeDelta by a factor
 *
 * @param ts_ms TimeDelta as milliseconds
 * @param factor Multiplication factor
 * @return int64_t Resulting TimeDelta as milliseconds
 */
auto f_timedelta_multiply(const int64_t ts_ms, const int factor) -> int64_t {
  return ts_ms * static_cast<int64_t>(factor);
}

/**
 * @brief Divide a TimeDelta by a divisor
 *
 * @param ts_ms TimeDelta as milliseconds
 * @param divisor Division factor
 * @return int64_t Resulting TimeDelta as milliseconds
 */
auto f_timedelta_divide(const int64_t ts_ms, const int divisor) -> int64_t {
  return ts_ms / static_cast<int64_t>(divisor);
}

/**
 * @brief Convert a TimeDelta to a string representation
 *
 * @param ts_ms TimeDelta as milliseconds
 * @param buffer Output buffer for the string
 * @param buffer_size Size of the output buffer
 */
void f_timedelta_to_string(const int64_t ts_ms, char* buffer,
                           const int buffer_size) {
  if (buffer_size <= 0 || buffer == nullptr) {
    std::cerr << "[ERROR]: Invalid buffer size or null buffer\n";
    return;
  }

  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  const std::string str = time_delta.toString();

  const size_t copy_length =
      std::min(static_cast<size_t>(buffer_size - 1), str.length());

  strncpy(buffer, str.c_str(), copy_length);
  buffer[copy_length] = '\0';
}

/**
 * @brief Compare two TimeDeltas for equality
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return true if equal, false otherwise
 */
auto f_timedelta_equals(const int64_t ts1_ms, const int64_t ts2_ms) -> bool {
  return ts1_ms == ts2_ms;
}

/**
 * @brief Check if one TimeDelta is less than another
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return true if ts1 < ts2, false otherwise
 */
auto f_timedelta_less_than(const int64_t ts1_ms, const int64_t ts2_ms) -> bool {
  return ts1_ms < ts2_ms;
}

/**
 * @brief Check if one TimeDelta is greater than another
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return true if ts1 > ts2, false otherwise
 */
auto f_timedelta_greater_than(const int64_t ts1_ms, const int64_t ts2_ms)
    -> bool {
  return ts1_ms > ts2_ms;
}

/**
 * @brief Check if one TimeDelta is less than or equal to another
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return true if ts1 <= ts2, false otherwise
 */
auto f_timedelta_less_equal(const int64_t ts1_ms, const int64_t ts2_ms)
    -> bool {
  return ts1_ms <= ts2_ms;
}

/**
 * @brief Check if one TimeDelta is greater than or equal to another
 *
 * @param ts1_ms First TimeDelta as milliseconds
 * @param ts2_ms Second TimeDelta as milliseconds
 * @return true if ts1 >= ts2, false otherwise
 */
auto f_timedelta_greater_equal(const int64_t ts1_ms, const int64_t ts2_ms)
    -> bool {
  return ts1_ms >= ts2_ms;
}

//=============================================================================
// DateTime functions
//=============================================================================

/**
 * @brief Create a DateTime from year, month, day, etc.
 *
 * @param year Year
 * @param month Month (1-12)
 * @param day Day (1-31)
 * @param hour Hour (0-23)
 * @param minute Minute (0-59)
 * @param second Second (0-59)
 * @param millisecond Millisecond (0-999)
 * @return int64_t DateTime as milliseconds since epoch
 */
auto f_datetime_create(const int year, const int month, const int day,
                       const int hour, const int minute, const int second,
                       const int millisecond) -> int64_t {
  if (year < DateTime::DATETIME_MIN_YEAR ||
      month < DateTime::DATETIME_MIN_MONTH ||
      month > DateTime::DATETIME_MAX_MONTH ||
      day < DateTime::DATETIME_MIN_DAYS || day > DateTime::DATETIME_MAX_DAYS ||
      hour < DateTime::DATETIME_MIN_HOURS ||
      hour > DateTime::DATETIME_MAX_HOURS ||
      minute < DateTime::DATETIME_MIN_MINUTES ||
      minute > DateTime::DATETIME_MAX_MINUTES ||
      second < DateTime::DATETIME_MIN_SECONDS ||
      second > DateTime::DATETIME_MAX_SECONDS ||
      millisecond < DateTime::DATETIME_MIN_MILLISECONDS ||
      millisecond > DateTime::DATETIME_MAX_MILLISECONDS) {
    return DateTime::INVALID_TIMESTAMP;
  } else {
    const auto u_month = static_cast<unsigned>(month);
    const auto u_day = static_cast<unsigned>(day);
    const auto u_hour = static_cast<unsigned>(hour);
    const auto u_minute = static_cast<unsigned>(minute);
    const auto u_second = static_cast<unsigned>(second);
    const auto u_millisecond = static_cast<unsigned>(millisecond);
    const DateTime date(year, u_month, u_day, u_hour, u_minute, u_second,
                        u_millisecond);
    return date.timestamp();
  }
}

/**
 * @brief Get the current DateTime
 *
 * @return int64_t Current DateTime as milliseconds since epoch
 */
auto f_datetime_now() -> int64_t { return DateTime::now().timestamp(); }

/**
 * @brief Parse a DateTime from a string
 *
 * @param str String representation of a DateTime
 * @param format Format string (similar to strftime)
 * @param str_len Length of the string
 * @param format_len Length of the format string
 * @return int64_t DateTime as milliseconds since epoch
 */
auto f_datetime_strptime(const char* str, const char* format, const int str_len,
                         const int format_len) -> int64_t {
  // Cast the lengths coming from fortran over to size_t
  if (str_len <= 0 || format == nullptr) {
    return DateTime::INVALID_TIMESTAMP;
  }

  if (format_len <= 0 || str == nullptr) {
    return DateTime::INVALID_TIMESTAMP;
  }

  try {
    const auto str_len_t = static_cast<size_t>(str_len);
    const auto format_len_t = static_cast<size_t>(format_len);

    // Create strings from Fortran character arrays
    const std::string str_cpp(str, str_len_t);
    const std::string format_cpp(format, format_len_t);

    const auto date_time = DateTime::strptime(str_cpp, format_cpp);
    if (date_time.valid()) {
      return date_time.timestamp();
    } else {
      return DateTime::INVALID_TIMESTAMP;
    }
  } catch (...) {
    return DateTime::INVALID_TIMESTAMP;
  }
}

/**
 * @brief Get the year component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Year component
 */
auto f_datetime_get_year(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.year();
}

/**
 * @brief Get the month component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Month component (1-12)
 */
auto f_datetime_get_month(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.month();
}

/**
 * @brief Get the day component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Day component (1-31)
 */
auto f_datetime_get_day(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.day();
}

/**
 * @brief Get the hour component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Hour component (0-23)
 */
auto f_datetime_get_hour(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.hour();
}

/**
 * @brief Get the minute component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Minute component (0-59)
 */
auto f_datetime_get_minute(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.minute();
}

/**
 * @brief Get the second component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Second component (0-59)
 */
auto f_datetime_get_second(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.second();
}

/**
 * @brief Get the millisecond component from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int Millisecond component (0-999)
 */
auto f_datetime_get_millisecond(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.millisecond();
}

/**
 * @brief Get the Julian Day Number (JDN) for a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return int64_t Julian Day Number (integer days since JD epoch)
 */
auto f_datetime_get_julian_day_number(const int64_t dt_ms) -> int64_t {
  const DateTime date(dt_ms);
  return date.julianDayNumber();
}

/**
 * @brief Get the Julian Day (JD) with fractional day for a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return double Julian Day with fractional part (days.fraction since JD epoch)
 */
auto f_datetime_get_julian_day(const int64_t dt_ms) -> double {
  const DateTime date(dt_ms);
  return date.julianDay();
}

/**
 * @brief Get the Julian Century (JC) for a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return double Julian Century (time unit used in astronomy)
 */
auto f_datetime_get_julian_century(const int64_t dt_ms) -> double {
  const DateTime date(dt_ms);
  return date.julianCentury();
}

/**
 * @brief Add a TimeDelta to a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Resulting DateTime as milliseconds since epoch
 */
auto f_datetime_add_timedelta(const int64_t dt_ms, const int64_t ts_ms)
    -> int64_t {
  const DateTime date(dt_ms);
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return (date + time_delta).timestamp();
}

/**
 * @brief Subtract a TimeDelta from a DateTime
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @param ts_ms TimeDelta as milliseconds
 * @return int64_t Resulting DateTime as milliseconds since epoch
 */
auto f_datetime_subtract_timedelta(const int64_t dt_ms, const int64_t ts_ms)
    -> int64_t {
  const DateTime date(dt_ms);
  const TimeDelta time_delta(TimeDelta::to_components(ts_ms));
  return (date - time_delta).timestamp();
}

/**
 * @brief Calculate the difference between two DateTimes
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return int64_t Resulting TimeDelta as milliseconds
 */
auto f_datetime_difference(const int64_t dt1_ms, const int64_t dt2_ms)
    -> int64_t {
  const DateTime dt1(dt1_ms);
  const DateTime dt2(dt2_ms);
  return (dt1 - dt2).totalMilliseconds();
}

/**
 * @brief Format a DateTime to a string
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @param format Format string (similar to strftime)
 * @param buffer Output buffer for the string
 * @param format_len Length of the format string
 * @param buffer_size Size of the output buffer
 */
void f_datetime_strftime(const int64_t dt_ms, const char* format, char* buffer,
                         const int format_len, const int buffer_size) {
  if (format_len <= 0 || buffer_size <= 0 || format == nullptr ||
      buffer == nullptr) {
    std::cerr << "[ERROR]: Invalid format/buffer size or null buffer/format\n";
    return;
  }

  const auto format_len_t = static_cast<size_t>(format_len);
  const auto buffer_size_t = static_cast<size_t>(buffer_size);

  const DateTime date(dt_ms);
  const std::string format_cpp(format, format_len_t);
  const std::string str = date.strftime(format_cpp);
  strncpy(buffer, str.c_str(), buffer_size_t - 1);
  buffer[buffer_size_t - 1] = '\0';
}

/**
 * @brief Format a DateTime to a string with milliseconds
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @param format Format string (similar to strftime)
 * @param buffer Output buffer for the string
 * @param format_len Length of the format string
 * @param buffer_size Size of the output buffer
 */
void f_datetime_strftime_milliseconds(const int64_t dt_ms, const char* format,
                                      char* buffer, const int format_len,
                                      const int buffer_size) {
  if (format_len <= 0 || buffer_size <= 0 || format == nullptr ||
      buffer == nullptr) {
    std::cerr << "[ERROR]: Invalid format/buffer size or null format/buffer\n";
    return;
  }

  const auto format_len_t = static_cast<size_t>(format_len);
  const auto buffer_size_t = static_cast<size_t>(buffer_size);

  const DateTime date(dt_ms);
  const std::string format_cpp(format, format_len_t);
  const std::string str = date.strftime_w_milliseconds(format_cpp);
  strncpy(buffer, str.c_str(), buffer_size_t - 1);
  buffer[buffer_size_t - 1] = '\0';
}

/**
 * @brief Convert a DateTime to ISO string format
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @param buffer Output buffer for the string
 * @param buffer_size Size of the output buffer
 */
void f_datetime_to_iso_string(const int64_t dt_ms, char* buffer,
                              const int buffer_size) {
  if (buffer_size <= 0 || buffer == nullptr) {
    std::cerr << "[ERROR]: Invalid buffer size or null buffer\n";
    return;
  }

  const auto buffer_size_t = static_cast<size_t>(buffer_size);

  const DateTime date(dt_ms);
  const std::string str = date.to_iso_string();
  strncpy(buffer, str.c_str(), buffer_size_t - 1);
  buffer[buffer_size_t - 1] = '\0';
}

/**
 * @brief Compare two DateTimes for equality
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return true if equal, false otherwise
 */
auto f_datetime_equals(const int64_t dt1_ms, const int64_t dt2_ms) -> bool {
  return dt1_ms == dt2_ms;
}

/**
 * @brief Check if one DateTime is less than another
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return true if dt1 < dt2, false otherwise
 */
auto f_datetime_less_than(const int64_t dt1_ms, const int64_t dt2_ms) -> bool {
  return dt1_ms < dt2_ms;
}

/**
 * @brief Check if one DateTime is greater than another
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return true if dt1 > dt2, false otherwise
 */
auto f_datetime_greater_than(const int64_t dt1_ms, const int64_t dt2_ms)
    -> bool {
  return dt1_ms > dt2_ms;
}

/**
 * @brief Check if one DateTime is less than or equal to another
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return true if dt1 <= dt2, false otherwise
 */
auto f_datetime_less_equal(const int64_t dt1_ms, const int64_t dt2_ms) -> bool {
  return dt1_ms <= dt2_ms;
}

/**
 * @brief Check if one DateTime is greater than or equal to another
 *
 * @param dt1_ms First DateTime as milliseconds since epoch
 * @param dt2_ms Second DateTime as milliseconds since epoch
 * @return true if dt1 >= dt2, false otherwise
 */
auto f_datetime_greater_equal(const int64_t dt1_ms, const int64_t dt2_ms)
    -> bool {
  return dt1_ms >= dt2_ms;
}

/**
 * @brief Check if a DateTime is valid
 *
 * @param dt_ms DateTime as milliseconds since epoch
 * @return true if valid, false otherwise
 */
auto f_datetime_is_valid(const int64_t dt_ms) -> bool {
  return dt_ms != DateTime::INVALID_TIMESTAMP;
}

/**
 * Returns the invalid timestamp constant.
 * @return int64_t The invalid timestamp constant.
 */
auto f_datetime_invalid_timestamp() -> int64_t {
  return DateTime::INVALID_TIMESTAMP;
}

/**
 * @brief Parse a DateTime from a string using an array of format options
 *
 * @param str String representation of a DateTime
 * @param formats Array of format strings (similar to strftime)
 * @param str_len Length of the string
 * @param formats_len Array of lengths for each format string
 * @param num_formats Number of format strings in the array
 * @return int64_t DateTime as milliseconds since epoch, or INVALID_TIMESTAMP if
 * all formats fail
 */
auto f_datetime_strptime_with_formats(const char* str, const char** formats,
                                      const int str_len, const int* formats_len,
                                      const int num_formats) -> int64_t {
  // Enhanced input validation
  if (str_len <= 0 || str == nullptr || formats == nullptr ||
      formats_len == nullptr || num_formats <= 0) {
    return DateTime::INVALID_TIMESTAMP;
  }

  // Prevent excessive format counts that could cause performance issues
  if (num_formats > 1000) {
    return DateTime::INVALID_TIMESTAMP;
  }

  // Validate string length to prevent buffer issues
  if (str_len > 10000) {
    return DateTime::INVALID_TIMESTAMP;
  }

  try {
    const auto str_len_t = static_cast<size_t>(str_len);
    const std::string str_cpp(str, str_len_t);

    // Try each format until one succeeds
    for (int i = 0; i < num_formats; ++i) {
      if (formats[i] == nullptr || formats_len[i] <= 0) {
        continue;
      }

      // Validate individual format length to prevent excessive memory usage
      if (formats_len[i] > 1000) {
        continue;
      }

      const auto format_len_t = static_cast<size_t>(formats_len[i]);
      const std::string format_cpp(formats[i], format_len_t);

      const auto date_time = DateTime::strptime(str_cpp, format_cpp);
      if (date_time.valid()) {
        return date_time.timestamp();
      }
    }

    return DateTime::INVALID_TIMESTAMP;
  } catch (...) {
    return DateTime::INVALID_TIMESTAMP;
  }
}

/**
 * @brief Parse a DateTime with automatic format detection, falling back to
 * custom formats
 *
 * @param str String representation of a DateTime
 * @param fallback_formats Array of fallback format strings to try after
 * auto-detection fails
 * @param str_len Length of the string
 * @param formats_len Array of lengths for each fallback format string
 * @param num_formats Number of fallback format strings in the array
 * @return int64_t DateTime as milliseconds since epoch, or INVALID_TIMESTAMP if
 * all attempts fail
 */
auto f_datetime_strptime_auto_with_fallback(const char* str,
                                            const char** fallback_formats,
                                            const int str_len,
                                            const int* formats_len,
                                            const int num_formats) -> int64_t {
  if (str_len <= 0 || str == nullptr) {
    return DateTime::INVALID_TIMESTAMP;
  }

  try {
    const auto str_len_t = static_cast<size_t>(str_len);
    const std::string str_cpp(str, str_len_t);

    // First try automatic format detection
    const auto date_time = DateTime::strptime(str_cpp, "auto");
    if (date_time.valid()) {
      return date_time.timestamp();
    }

    // If auto-detection failed and we have fallback formats, try them
    if (fallback_formats != nullptr && formats_len != nullptr &&
        num_formats > 0) {
      for (int i = 0; i < num_formats; ++i) {
        if (fallback_formats[i] == nullptr || formats_len[i] <= 0) {
          continue;
        }

        const auto format_len_t = static_cast<size_t>(formats_len[i]);
        const std::string format_cpp(fallback_formats[i], format_len_t);

        const auto fallback_date_time = DateTime::strptime(str_cpp, format_cpp);
        if (fallback_date_time.valid()) {
          return fallback_date_time.timestamp();
        }
      }
    }

    return DateTime::INVALID_TIMESTAMP;
  } catch (...) {
    return DateTime::INVALID_TIMESTAMP;
  }
}

}  // extern "C"
