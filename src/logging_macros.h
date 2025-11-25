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
!  LOGGING MACROS
!-----------------------------------------------------------------------
!> @file logging_macros.h
!> @brief Preprocessor macros for RAII-based logging
!>
!> This header file provides the LOG_SCOPE macro for automatic message
!> source management. Include this file in any Fortran source that uses
!> the logging system with: #include "logging_macros.h"
!-----------------------------------------------------------------------

!> @brief LOG_SCOPE macro for automatic message source management
!>
!> This macro creates a log scope guard that automatically manages the
!> message source stack. The cleanup happens automatically when the
!> scope exits, even with early returns.
!>
!> Usage:
!>   subroutine my_function()
!>      use mod_logging, only: allMessage, INFO, t_log_scope, init_log_scope
!>      #include "logging_macros.h"
!>
!>      LOG_SCOPE("my_function")
!>
!>      ! ... function body ...
!>      if (error) return  ! cleanup happens automatically
!>   end subroutine
!>
!> Note: The macro expands to two statements - a declaration and a subroutine call.
!> This avoids temporary object creation issues with Intel compilers that can
!> cause premature finalization when using function-based constructors.
!-----------------------------------------------------------------------
#ifndef LOG_SCOPE
#define LOG_SCOPE(fname) type(t_log_scope) :: ls__ ;\
call init_log_scope(ls__, fname)
#endif

!-----------------------------------------------------------------------
!> @brief LOG_SCOPE_TRACED macro for automatic message source management with tracing
!>
!> This macro creates a log scope guard that automatically manages the
!> message source stack AND emits Enter/Exit debug messages when a trace
!> flag is enabled.
!>
!> Usage:
!>   subroutine my_function()
!>      use mod_logging, only: allMessage, INFO, t_log_scope, init_log_scope, DEBUG
!>      #include "logging_macros.h"
!>
!>      LOG_SCOPE_TRACED("my_function", COUPLE2SWAN_TRACING)
!>
!>      ! When COUPLE2SWAN_TRACE or ALL_TRACE is defined:
!>      ! - "Enter." is logged at DEBUG level on entry
!>      ! - "Exit." is logged at DEBUG level on exit (even with early returns)
!>   end subroutine
!-----------------------------------------------------------------------
#ifndef LOG_SCOPE_TRACED
#define LOG_SCOPE_TRACED(fname, trace_flag) type(t_log_scope) :: ls__ ;\
call init_log_scope(ls__, fname, trace_flag .or. ALL_TRACING)
#endif

!-----------------------------------------------------------------------
! Global trace flag - enables tracing for all modules
!-----------------------------------------------------------------------
#ifdef ALL_TRACE
#define ALL_TRACING .true.
#else
#define ALL_TRACING .false.
#endif

!-----------------------------------------------------------------------
! Module-specific trace flags
!-----------------------------------------------------------------------
#ifdef ADCIRC_TRACE
#define ADCIRC_TRACING .true.
#else
#define ADCIRC_TRACING .false.
#endif

#ifdef COLDSTART_TRACE
#define COLDSTART_TRACING .true.
#else
#define COLDSTART_TRACING .false.
#endif

#ifdef COUPLE2SWAN_TRACE
#define COUPLE2SWAN_TRACING .true.
#else
#define COUPLE2SWAN_TRACING .false.
#endif

#ifdef GLOBALIO_TRACE
#define GLOBALIO_TRACING .true.
#else
#define GLOBALIO_TRACING .false.
#endif

#ifdef GLOBAL_TRACE
#define GLOBAL_TRACING .true.
#else
#define GLOBAL_TRACING .false.
#endif

#ifdef GLOBAL_3DVS_TRACE
#define GLOBAL_3DVS_TRACING .true.
#else
#define GLOBAL_3DVS_TRACING .false.
#endif

#ifdef HARM_TRACE
#define HARM_TRACING .true.
#else
#define HARM_TRACING .false.
#endif

#ifdef HOTSTART_TRACE
#define HOTSTART_TRACING .true.
#else
#define HOTSTART_TRACING .false.
#endif

#ifdef MESH_TRACE
#define MESH_TRACING .true.
#else
#define MESH_TRACING .false.
#endif

#ifdef MESSENGER_TRACE
#define MESSENGER_TRACING .true.
#else
#define MESSENGER_TRACING .false.
#endif

#ifdef NETCDF_TRACE
#define NETCDF_TRACING .true.
#else
#define NETCDF_TRACING .false.
#endif

#ifdef NODALATTR_TRACE
#define NODALATTR_TRACING .true.
#else
#define NODALATTR_TRACING .false.
#endif

#ifdef OWIWIND_TRACE
#define OWIWIND_TRACING .true.
#else
#define OWIWIND_TRACING .false.
#endif

#ifdef READ_INPUT_TRACE
#define READ_INPUT_TRACING .true.
#else
#define READ_INPUT_TRACING .false.
#endif

#ifdef PREPINPUT_TRACE
#define PREPINPUT_TRACING .true.
#else
#define PREPINPUT_TRACING .false.
#endif

#ifdef READ_GLOBAL_TRACE
#define READ_GLOBAL_TRACING .true.
#else
#define READ_GLOBAL_TRACING .false.
#endif

#ifdef SUBDOMAIN_TRACE
#define SUBDOMAIN_TRACING .true.
#else
#define SUBDOMAIN_TRACING .false.
#endif

#ifdef SUBPREP_TRACE
#define SUBPREP_TRACING .true.
#else
#define SUBPREP_TRACING .false.
#endif

#ifdef SWAN_TRACE
#define SWAN_TRACING .true.
#else
#define SWAN_TRACING .false.
#endif

#ifdef TIMESTEP_TRACE
#define TIMESTEP_TRACING .true.
#else
#define TIMESTEP_TRACING .false.
#endif

#ifdef TVW_TRACE
#define TVW_TRACING .true.
#else
#define TVW_TRACING .false.
#endif

#ifdef VSMY_TRACE
#define VSMY_TRACING .true.
#else
#define VSMY_TRACING .false.
#endif

#ifdef WEIR_TRACE
#define WEIR_TRACING .true.
#else
#define WEIR_TRACING .false.
#endif

#ifdef WIND_TRACE
#define WIND_TRACING .true.
#else
#define WIND_TRACING .false.
#endif

#ifdef WRITER_TRACE
#define WRITER_TRACING .true.
#else
#define WRITER_TRACING .false.
#endif

#ifdef WRITE_OUTPUT_TRACE
#define WRITE_OUTPUT_TRACING .true.
#else
#define WRITE_OUTPUT_TRACING .false.
#endif

#ifdef GWCE_TRACE
#define GWCE_TRACING .true.
#else
#define GWCE_TRACING .false.
#endif

#ifdef CSTART_TRACE
#define CSTART_TRACING .true.
#else
#define CSTART_TRACING .false.
#endif

#ifdef HSTART_TRACE
#define HSTART_TRACING .true.
#else
#define HSTART_TRACING .false.
#endif

#ifdef NETCDFIO_TRACE
#define NETCDFIO_TRACING .true.
#else
#define NETCDFIO_TRACING .false.
#endif

#ifdef TIDALPOTENTIAL_TRACE
#define TIDALPOTENTIAL_TRACING .true.
#else
#define TIDALPOTENTIAL_TRACING .false.
#endif

#ifdef TRANSPORT_TRACE
#define TRANSPORT_TRACING .true.
#else
#define TRANSPORT_TRACING .false.
#endif

#ifdef WEIR_BOUNDARY_TRACE
#define WEIR_BOUNDARY_TRACING .true.
#else
#define WEIR_BOUNDARY_TRACING .false.
#endif

#ifdef IO_TRACE
#define IO_TRACING .true.
#else
#define IO_TRACING .false.
#endif

#ifdef COUPLE2BAROCLINIC3D_TRACE
#define COUPLE2BAROCLINIC3D_TRACING .true.
#else
#define COUPLE2BAROCLINIC3D_TRACING .false.
#endif
