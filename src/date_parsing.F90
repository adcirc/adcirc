!-----------------------------------------------------------------------
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------
!> @file date_parsing.F90
!> @brief Lenient parsing of the date strings found in ADCIRC input files
!>
!> Dates in fort.15 metadata, namelists, and forcing files come in many
!> slightly different spellings, and older ADCIRC versions read most of them
!> with fixed-column reads that ignored the separators. A single strptime
!> format cannot express "seconds are optional" or "any of -, / or .", so
!> this module accepts the family of spellings that ADCIRC has historically
!> accepted while still rejecting anything ambiguous or that would be
!> misinterpreted (e.g. day-first dates or non-UTC offsets).
!>
!> This module intentionally depends only on mod_datetime so that it can be
!> compiled before every other ADCIRC module.
!-----------------------------------------------------------------------
module mod_date_parsing
   use mod_datetime, only: t_datetime, null_datetime
   implicit none

   private

   public :: parse_datetime_lenient

contains

   !> @brief Parse a date/time string, accepting the spellings used in ADCIRC inputs
   !>
   !> Accepted forms (leading/trailing whitespace ignored):
   !>  - Delimited: YYYY{-/.}MM{-/.}DD[{ |T}hh[:mm[:ss[{.|,}fff]]]]
   !>    Month, day and time fields may have one or two digits. Missing time
   !>    fields are zero.
   !>  - Compact:   YYYYMMDD[hhmmss] or YYYYMMDD{.| |T}hhmmss
   !>    Time digits are read two at a time; a final single digit is taken as
   !>    its own value. This matches the legacy (I4,I2,I2,1X,I2,I2,I2) read
   !>    used for NWS13ColdStartString, so "20080826.0000" is 00:00:00.
   !>  - Optionally followed by a UTC designator: Z, UTC, GMT, +00, +0000,
   !>    +00:00 (or with a minus sign). Any other trailing text, including a
   !>    non-zero offset, makes the string invalid because ADCIRC treats all
   !>    times as UTC.
   !>
   !> The year must have exactly four digits so that day-first strings such as
   !> 08/26/2008 are rejected instead of being misread.
   !>
   !> @param[in] str string to parse
   !> @return parsed datetime, or an invalid datetime if the string is not recognized
   function parse_datetime_lenient(str) result(dt)
      implicit none
      character(len=*), intent(in) :: str
      type(t_datetime) :: dt

      character(len=len(str)) :: s
      integer :: n, pos, ndig
      integer :: year, month, day, hour, minute, second, millisecond
      logical :: ok

      dt = null_datetime()
      s = adjustl(str)
      n = len_trim(s)
      if (n == 0) return

      hour = 0
      minute = 0
      second = 0
      millisecond = 0

      pos = 1
      ndig = count_digits(s, pos, n)

      if (ndig >= 8) then
         ! Compact form: YYYYMMDD with any time digits run on or after one separator
         if (ndig > 14) return
         read (s(1:4), '(I4)') year
         read (s(5:6), '(I2)') month
         read (s(7:8), '(I2)') day
         if (ndig > 8) then
            call read_compact_time(s(9:ndig), hour, minute, second, ok)
            if (.not. ok) return
            pos = ndig + 1
         else
            pos = 9
            if (pos < n .and. index(". T", s(pos:pos)) > 0 .and. is_digit(s(pos + 1:pos + 1))) then
               ndig = count_digits(s, pos + 1, n)
               if (ndig > 6) return
               call read_compact_time(s(pos + 1:pos + ndig), hour, minute, second, ok)
               if (.not. ok) return
               pos = pos + 1 + ndig
            end if
         end if
      elseif (ndig == 4) then
         ! Delimited form
         read (s(1:4), '(I4)') year
         pos = 5
         if (.not. read_separated_field(s, n, pos, "-/.", 2, month)) return
         if (.not. read_separated_field(s, n, pos, "-/.", 2, day)) return
         if (pos <= n) then
            if (index(" T", s(pos:pos)) > 0) then
               ! Skip the date/time separator (and any extra spaces) only when a time follows
               call skip_date_time_separator(s, n, pos)
               if (pos <= n .and. is_digit(s(pos:pos))) then
                  if (.not. read_field(s, n, pos, 2, hour)) return
                  if (read_separated_field(s, n, pos, ":", 2, minute)) then
                     if (read_separated_field(s, n, pos, ":", 2, second)) then
                        call read_fraction(s, n, pos, millisecond)
                     end if
                  end if
               end if
            end if
         end if
      else
         return
      end if

      if (.not. is_utc_designator(s(min(pos, n + 1):n))) return

      dt = t_datetime(year, month, day, hour, minute, second, millisecond)

   end function parse_datetime_lenient

   !> Number of consecutive digits in s starting at pos (not past n)
   pure integer function count_digits(s, pos, n) result(cnt)
      implicit none
      character(len=*), intent(in) :: s
      integer, intent(in) :: pos, n
      cnt = 0
      do while (pos + cnt <= n)
         if (.not. is_digit(s(pos + cnt:pos + cnt))) exit
         cnt = cnt + 1
      end do
   end function count_digits

   pure logical function is_digit(c)
      implicit none
      character(len=1), intent(in) :: c
      is_digit = (c >= '0' .and. c <= '9')
   end function is_digit

   !> Read an unsigned integer of 1..maxdig digits at pos and advance pos
   logical function read_field(s, n, pos, maxdig, val) result(ok)
      implicit none
      character(len=*), intent(in) :: s
      integer, intent(in) :: n, maxdig
      integer, intent(inout) :: pos
      integer, intent(out) :: val
      integer :: ndig
      ok = .false.
      val = 0
      if (pos > n) return
      ndig = count_digits(s, pos, n)
      if (ndig < 1 .or. ndig > maxdig) return
      read (s(pos:pos + ndig - 1), *) val
      pos = pos + ndig
      ok = .true.
   end function read_field

   !> Read a one-character separator from seps followed by a field. On failure pos is unchanged.
   logical function read_separated_field(s, n, pos, seps, maxdig, val) result(ok)
      implicit none
      character(len=*), intent(in) :: s, seps
      integer, intent(in) :: n, maxdig
      integer, intent(inout) :: pos
      integer, intent(inout) :: val
      integer :: trial, v
      ok = .false.
      if (pos >= n) return
      if (index(seps, s(pos:pos)) == 0) return
      trial = pos + 1
      if (.not. read_field(s, n, trial, maxdig, v)) return
      val = v
      pos = trial
      ok = .true.
   end function read_separated_field

   !> Skip a 'T' or a run of spaces between the date and time
   subroutine skip_date_time_separator(s, n, pos)
      implicit none
      character(len=*), intent(in) :: s
      integer, intent(in) :: n
      integer, intent(inout) :: pos
      integer :: trial
      trial = pos
      if (s(trial:trial) == 'T') then
         trial = trial + 1
      else
         do while (trial <= n)
            if (s(trial:trial) /= ' ') exit
            trial = trial + 1
         end do
      end if
      ! Only consume the separator when a time actually follows; otherwise the
      ! remaining text is checked as a timezone designator
      if (trial <= n) then
         if (is_digit(s(trial:trial))) pos = trial
      end if
   end subroutine skip_date_time_separator

   !> Read optional fractional seconds; only the first three digits are kept
   subroutine read_fraction(s, n, pos, millisecond)
      implicit none
      character(len=*), intent(in) :: s
      integer, intent(in) :: n
      integer, intent(inout) :: pos
      integer, intent(out) :: millisecond
      character(len=3) :: frac
      integer :: ndig
      millisecond = 0
      if (pos >= n) return
      if (index(".,", s(pos:pos)) == 0) return
      ndig = count_digits(s, pos + 1, n)
      if (ndig == 0) return
      frac = '000'
      frac(1:min(ndig, 3)) = s(pos + 1:pos + min(ndig, 3))
      read (frac, '(I3)') millisecond
      pos = pos + 1 + ndig
   end subroutine read_fraction

   !> Interpret compact time digits as hh, mm, ss pairs (legacy fixed-column semantics)
   subroutine read_compact_time(digits, hour, minute, second, ok)
      implicit none
      character(len=*), intent(in) :: digits
      integer, intent(out) :: hour, minute, second
      logical, intent(out) :: ok
      integer :: vals(3), i, first
      vals = 0
      ok = len(digits) <= 6
      if (.not. ok) return
      do i = 1, 3
         first = 2*i - 1
         if (first > len(digits)) exit
         read (digits(first:min(first + 1, len(digits))), *) vals(i)
      end do
      hour = vals(1)
      minute = vals(2)
      second = vals(3)
   end subroutine read_compact_time

   !> True if the text is empty or designates UTC
   logical function is_utc_designator(text) result(is_utc)
      implicit none
      character(len=*), intent(in) :: text
      character(len=len(text)) :: t
      integer :: i
      t = adjustl(text)
      do i = 1, len_trim(t)
         if (t(i:i) >= 'a' .and. t(i:i) <= 'z') t(i:i) = achar(iachar(t(i:i)) - 32)
      end do
      select case (trim(t))
      case ("", "Z", "UTC", "GMT", "+00", "+0000", "+00:00", "-00", "-0000", "-00:00")
         is_utc = .true.
      case default
         is_utc = .false.
      end select
   end function is_utc_designator

end module mod_date_parsing
