/*
 * FDate - A Fortran Date and Time Library based on C++
 * Copyright (C) 2024 Zach Cobell (zcobell@gmail.com)
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

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

#include "date_hh.h"

/**
 * @brief A high-precision time duration class with millisecond accuracy
 *
 * The TimeDelta class represents time intervals or durations with millisecond
 * precision. It can be used for time arithmetic operations and provides
 * convenient methods for creating time spans from various time units (days,
 * hours, minutes, seconds, milliseconds).
 *
 * TimeDelta objects can be added to or subtracted from DateTime objects, and
 * can be combined with other TimeDelta objects using standard arithmetic
 * operations.
 *
 * The class provides both individual component access (e.g., days(), hours())
 * and total duration access (e.g., totalDays(), totalHours()) for flexible time
 * calculations.
 *
 * @note All internal calculations use std::chrono::milliseconds for precision
 * @note Negative time spans are supported and handled correctly
 * @see DateTime for date and time point operations
 */
class TimeDelta {
  /** @brief Internal duration storage with millisecond precision */
  std::chrono::milliseconds m_duration;

  /**
   * @brief Private constructor from milliseconds to initialize the TimeDelta
   *
   * @param milliseconds The duration in milliseconds
   *
   * @note This constructor is private to enforce the use of factory methods and
   * public constructors
   */
  explicit constexpr TimeDelta(std::chrono::milliseconds milliseconds) noexcept
      : m_duration(milliseconds) {}

 public:
  /**
   * @brief Structure to hold the decomposed components of a TimeDelta
   *
   * This structure represents a TimeDelta broken down into its constituent
   * parts. All components can be negative if the overall TimeDelta is negative.
   */
  struct s_TimedeltaComponents {
    int64_t days;          ///< Number of days component
    int64_t hours;         ///< Number of hours component (0-23)
    int64_t minutes;       ///< Number of minutes component (0-59)
    int64_t seconds;       ///< Number of seconds component (0-59)
    int64_t milliseconds;  ///< Number of milliseconds component (0-999)
  };

  // Constructors

  /**
   * @brief Default constructor initializes the TimeDelta to zero duration
   */
  constexpr TimeDelta() noexcept : m_duration(std::chrono::milliseconds{0}) {}

  /**
   * @brief Constructs a TimeDelta from a components structure
   *
   * @param components A structure containing the individual time components
   *
   * @note All components are summed to create the total duration
   * @note Negative components are handled correctly
   */
  constexpr explicit TimeDelta(const s_TimedeltaComponents& components) noexcept
      : m_duration(std::chrono::milliseconds{
            date::days(components.days) + std::chrono::hours(components.hours) +
            std::chrono::minutes(components.minutes) +
            std::chrono::seconds(components.seconds) +
            std::chrono::milliseconds(components.milliseconds)}) {}

  /**
   * @brief Constructs a TimeDelta from individual time components
   *
   * @param days Number of days
   * @param hours Number of hours
   * @param minutes Number of minutes
   * @param seconds Number of seconds
   * @param milliseconds Number of milliseconds
   *
   * @note All components are summed to create the total duration
   * @note Negative values are allowed and handled correctly
   * @note No validation is performed; values outside typical ranges (e.g.,
   * hours > 23) are accepted
   */
  explicit constexpr TimeDelta(const int days, const int hours,
                               const int minutes, const int seconds,
                               const int milliseconds) noexcept
      : TimeDelta(s_TimedeltaComponents{days, hours, minutes, seconds,
                                        milliseconds}) {}

  // Static factory methods

  /**
   * @brief Creates a TimeDelta from a number of days
   *
   * @param days Number of days (can be negative)
   * @return TimeDelta object representing the specified duration
   *
   * @see fromHours(), fromMinutes(), fromSeconds(), fromMilliseconds()
   */
  [[nodiscard]] static constexpr auto fromDays(const int64_t days) noexcept
      -> TimeDelta {
    return TimeDelta(s_TimedeltaComponents{days, 0, 0, 0, 0});
  }

  /**
   * @brief Creates a TimeDelta from a number of hours
   *
   * @param hours Number of hours (can be negative)
   * @return TimeDelta object representing the specified duration
   *
   * @see fromDays(), fromMinutes(), fromSeconds(), fromMilliseconds()
   */
  [[nodiscard]] static constexpr auto fromHours(const int64_t hours) noexcept
      -> TimeDelta {
    return TimeDelta(s_TimedeltaComponents{0, hours, 0, 0, 0});
  }

  /**
   * @brief Creates a TimeDelta from a number of minutes
   *
   * @param minutes Number of minutes (can be negative)
   * @return TimeDelta object representing the specified duration
   *
   * @see fromDays(), fromHours(), fromSeconds(), fromMilliseconds()
   */
  [[nodiscard]] static constexpr auto fromMinutes(
      const int64_t minutes) noexcept -> TimeDelta {
    return TimeDelta(s_TimedeltaComponents{0, 0, minutes, 0, 0});
  }

  /**
   * @brief Creates a TimeDelta from a number of seconds
   *
   * @param seconds Number of seconds (can be negative)
   * @return TimeDelta object representing the specified duration
   *
   * @see fromDays(), fromHours(), fromMinutes(), fromMilliseconds()
   */
  [[nodiscard]] static constexpr auto fromSeconds(
      const int64_t seconds) noexcept -> TimeDelta {
    return TimeDelta(s_TimedeltaComponents{0, 0, 0, seconds, 0});
  }

  /**
   * @brief Creates a TimeDelta from a number of milliseconds
   *
   * @param milliseconds Number of milliseconds (can be negative)
   * @return TimeDelta object representing the specified duration
   *
   * @see fromDays(), fromHours(), fromMinutes(), fromSeconds()
   */
  [[nodiscard]] static constexpr auto fromMilliseconds(
      const int64_t milliseconds) noexcept -> TimeDelta {
    return TimeDelta(s_TimedeltaComponents{0, 0, 0, 0, milliseconds});
  }

  /**
   * @brief Gets the TimeDelta decomposed into its constituent components
   *
   * @return s_TimedeltaComponents structure containing days, hours, minutes,
   * seconds, and milliseconds
   *
   * @note For negative TimeDeltas, all components will be negative
   * @note The sum of all components equals the total duration
   * @see to_components() for static version
   */
  [[nodiscard]] constexpr auto components() const noexcept {
    return to_components(m_duration.count());
  }

  /**
   * @brief Static method to convert a duration in milliseconds to components
   *
   * Decomposes a duration given in milliseconds into days, hours, minutes,
   * seconds, and milliseconds components. This method handles negative
   * durations correctly by applying the negative sign to all components.
   *
   * @param duration_ms Duration in milliseconds (can be negative)
   * @return s_TimedeltaComponents structure containing the decomposed time
   * components
   *
   * @note For negative durations, all returned components will be negative
   * @note The method uses std::chrono::duration_cast for accurate conversions
   * @see components() for instance method version
   */
  static constexpr auto to_components(const int64_t duration_ms) noexcept
      -> s_TimedeltaComponents {
    const bool negative = duration_ms < 0;
    int64_t this_duration = negative ? -duration_ms : duration_ms;

    const auto n_days = std::chrono::duration_cast<date::days>(
        std::chrono::milliseconds(this_duration));
    this_duration -=
        std::chrono::duration_cast<std::chrono::milliseconds>(n_days).count();
    const auto n_hours = std::chrono::duration_cast<std::chrono::hours>(
        std::chrono::milliseconds(this_duration));
    this_duration -=
        std::chrono::duration_cast<std::chrono::milliseconds>(n_hours).count();
    const auto n_minutes = std::chrono::duration_cast<std::chrono::minutes>(
        std::chrono::milliseconds(this_duration));
    this_duration -=
        std::chrono::duration_cast<std::chrono::milliseconds>(n_minutes)
            .count();
    const auto n_seconds = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::milliseconds(this_duration));
    this_duration -=
        std::chrono::duration_cast<std::chrono::milliseconds>(n_seconds)
            .count();
    const auto n_milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::milliseconds(this_duration));
    this_duration -=
        std::chrono::duration_cast<std::chrono::milliseconds>(n_milliseconds)
            .count();

    assert(this_duration == 0);

    s_TimedeltaComponents this_components{n_days.count(), n_hours.count(),
                                          n_minutes.count(), n_seconds.count(),
                                          n_milliseconds.count()};

    // Apply negative sign if needed
    if (negative) {
      this_components.days = -this_components.days;
      this_components.hours = -this_components.hours;
      this_components.minutes = -this_components.minutes;
      this_components.seconds = -this_components.seconds;
      this_components.milliseconds = -this_components.milliseconds;
    }

    return this_components;
  }

  // Component Accessors

  /**
   * @brief Gets the days component of the TimeDelta
   *
   * @return int64_t The number of days (can be negative)
   *
   * @note This returns only the days component, not the total number of days
   * @see totalDays() for the total duration expressed in days
   */
  [[nodiscard]] constexpr auto days() const noexcept -> int64_t {
    return components().days;
  }

  /**
   * @brief Gets the hours component of the TimeDelta
   *
   * @return int64_t The number of hours component (0-23, or negative for
   * negative TimeDeltas)
   *
   * @note This returns only the hours component, not the total number of hours
   * @see totalHours() for the total duration expressed in hours
   */
  [[nodiscard]] constexpr auto hours() const noexcept -> int64_t {
    return components().hours;
  }

  /**
   * @brief Gets the minutes component of the TimeDelta
   *
   * @return int64_t The number of minutes component (0-59, or negative for
   * negative TimeDeltas)
   *
   * @note This returns only the minutes component, not the total number of
   * minutes
   * @see totalMinutes() for the total duration expressed in minutes
   */
  [[nodiscard]] constexpr auto minutes() const noexcept -> int64_t {
    return components().minutes;
  }

  /**
   * @brief Gets the seconds component of the TimeDelta
   *
   * @return int64_t The number of seconds component (0-59, or negative for
   * negative TimeDeltas)
   *
   * @note This returns only the seconds component, not the total number of
   * seconds
   * @see totalSeconds() for the total duration expressed in seconds
   */
  [[nodiscard]] constexpr auto seconds() const noexcept -> int64_t {
    return components().seconds;
  }

  /**
   * @brief Gets the milliseconds component of the TimeDelta
   *
   * @return int64_t The number of milliseconds component (0-999, or negative
   * for negative TimeDeltas)
   *
   * @note This returns only the milliseconds component, not the total number of
   * milliseconds
   * @see totalMilliseconds() for the total duration expressed in milliseconds
   */
  [[nodiscard]] constexpr auto milliseconds() const noexcept -> int64_t {
    return components().milliseconds;
  }

  // Total Duration Accessors

  /**
   * @brief Gets the total duration expressed as days
   *
   * @return int64_t The total number of days in the TimeDelta (can be negative)
   *
   * @note This converts the entire duration to days, truncating fractional
   * parts
   * @see days() for only the days component
   */
  [[nodiscard]] constexpr auto totalDays() const noexcept -> int64_t {
    return std::chrono::duration_cast<date::days>(m_duration).count();
  }

  /**
   * @brief Gets the total duration expressed as hours
   *
   * @return int64_t The total number of hours in the TimeDelta (can be
   * negative)
   *
   * @note This converts the entire duration to hours, truncating fractional
   * parts
   * @see hours() for only the hours component
   */
  [[nodiscard]] constexpr auto totalHours() const noexcept -> int64_t {
    return std::chrono::duration_cast<std::chrono::hours>(m_duration).count();
  }

  /**
   * @brief Gets the total duration expressed as minutes
   *
   * @return int64_t The total number of minutes in the TimeDelta (can be
   * negative)
   *
   * @note This converts the entire duration to minutes, truncating fractional
   * parts
   * @see minutes() for only the minutes component
   */
  [[nodiscard]] constexpr auto totalMinutes() const noexcept -> int64_t {
    return std::chrono::duration_cast<std::chrono::minutes>(m_duration).count();
  }

  /**
   * @brief Gets the total duration expressed as seconds
   *
   * @return int64_t The total number of seconds in the TimeDelta (can be
   * negative)
   *
   * @note This converts the entire duration to seconds, truncating fractional
   * parts
   * @see seconds() for only the seconds component
   */
  [[nodiscard]] constexpr auto totalSeconds() const noexcept -> int64_t {
    return std::chrono::duration_cast<std::chrono::seconds>(m_duration).count();
  }

  /**
   * @brief Gets the total duration expressed as milliseconds
   *
   * @return int64_t The total number of milliseconds in the TimeDelta (can be
   * negative)
   *
   * @note This is the most precise representation of the TimeDelta's duration
   * @see milliseconds() for only the milliseconds component
   */
  [[nodiscard]] constexpr auto totalMilliseconds() const noexcept -> int64_t {
    return m_duration.count();
  }

  // Internal access

  /**
   * @brief Gets the internal std::chrono::milliseconds representation
   *
   * @return std::chrono::milliseconds The internal duration object
   *
   * @note This method provides direct access to the internal representation
   * @note Primarily used for integration with other chrono-based operations
   */
  [[nodiscard]] constexpr auto duration() const noexcept
      -> std::chrono::milliseconds {
    return m_duration;
  }

  // Arithmetic Operators

  /**
   * @brief Adds two TimeDelta objects
   *
   * @param other The TimeDelta to add to this one
   * @return TimeDelta A new TimeDelta representing the sum of both durations
   *
   * @see operator-() for subtraction
   */
  [[nodiscard]] constexpr auto operator+(const TimeDelta& other) const noexcept
      -> TimeDelta {
    return TimeDelta(
        to_components(m_duration.count() + other.m_duration.count()));
  }

  /**
   * @brief Subtracts two TimeDelta objects
   *
   * @param other The TimeDelta to subtract from this one
   * @return TimeDelta A new TimeDelta representing the difference (this -
   * other)
   *
   * @see operator+() for addition
   */
  [[nodiscard]] constexpr auto operator-(const TimeDelta& other) const noexcept
      -> TimeDelta {
    return TimeDelta(
        to_components(m_duration.count() - other.m_duration.count()));
  }

  /**
   * @brief Multiplies the TimeDelta by a scalar factor
   *
   * @param factor The factor to multiply by (can be negative)
   * @return TimeDelta A new TimeDelta representing the product
   *
   * @note Multiplying by a negative factor reverses the TimeDelta's direction
   * @see operator/() for division
   */
  [[nodiscard]] constexpr auto operator*(const int64_t factor) const noexcept
      -> TimeDelta {
    return TimeDelta(to_components(m_duration.count() * factor));
  }

  /**
   * @brief Divides the TimeDelta by a scalar divisor
   *
   * @param divisor The divisor to divide by
   * @return TimeDelta A new TimeDelta representing the quotient
   *
   * @note Division by zero results in a zero TimeDelta being returned
   * @note Integer division truncates fractional parts
   * @see operator*() for multiplication
   */
  [[nodiscard]] constexpr auto operator/(const int64_t divisor) const noexcept
      -> TimeDelta {
    return divisor == 0
               ? fromMilliseconds(0)
               : TimeDelta(to_components(m_duration.count() / divisor));
  }

  // Comparison Operators

  /**
   * @brief Equality comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if both TimeDeltas represent the same duration
   */
  [[nodiscard]] constexpr auto operator==(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration == other.m_duration;
  }

  /**
   * @brief Inequality comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if the TimeDeltas represent different durations
   */
  [[nodiscard]] constexpr auto operator!=(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration != other.m_duration;
  }

  /**
   * @brief Less than comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if this TimeDelta is shorter than the other
   */
  [[nodiscard]] constexpr auto operator<(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration < other.m_duration;
  }

  /**
   * @brief Greater than comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if this TimeDelta is longer than the other
   */
  [[nodiscard]] constexpr auto operator>(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration > other.m_duration;
  }

  /**
   * @brief Less than or equal comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if this TimeDelta is shorter than or equal to the other
   */
  [[nodiscard]] constexpr auto operator<=(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration <= other.m_duration;
  }

  /**
   * @brief Greater than or equal comparison operator
   *
   * @param other The TimeDelta to compare with
   * @return bool True if this TimeDelta is longer than or equal to the other
   */
  [[nodiscard]] constexpr auto operator>=(const TimeDelta& other) const noexcept
      -> bool {
    return m_duration >= other.m_duration;
  }

  // String representation

  /**
   * @brief Converts the TimeDelta to a human-readable string representation
   *
   * The format depends on the duration:
   * - If days are present: "Nd HH:MM:SS" or "Nd HH:MM:SS.mmm"
   * - If no days: "HH:MM:SS" or "HH:MM:SS.mmm"
   * - Milliseconds are only included if non-zero
   * - Negative TimeDeltas are indicated by negative component values
   *
   * Examples:
   * - "02:30:45" (2 hours, 30 minutes, 45 seconds)
   * - "1d 14:30:25.123" (1 day, 14 hours, 30 minutes, 25 seconds, 123
   * milliseconds)
   * - "00:00:00.500" (500 milliseconds)
   *
   * @return std::string A formatted string representation of the TimeDelta
   *
   * @note All components are zero-padded for consistent formatting
   * @note Days are shown without padding when present
   */
  [[nodiscard]] auto toString() const -> std::string {
    std::ostringstream oss;
    if (days() != 0) {
      oss << days() << "d ";
    }
    oss << std::setfill('0') << std::setw(2) << hours() << ":"
        << std::setfill('0') << std::setw(2) << minutes() << ":"
        << std::setfill('0') << std::setw(2) << seconds();

    if (milliseconds() != 0) {
      oss << "." << std::setfill('0') << std::setw(3) << milliseconds();
    }

    return oss.str();
  }
};
