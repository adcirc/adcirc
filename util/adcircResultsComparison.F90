!-----------------------------------------------------------------------
! adcircResultCompare.F90
!   Written by Zachary Cobell
!              The Water Institute
!              zcobell@thewaterinstitute.org
!
! This code will comapre two sets of ADCIRC files to see if the
! solution is different between them. The solution can be compared
! across formats as well. A return code of 0 is returned if the files
! are considered to match. A return code of 1 is returned otherwise.
!-----------------------------------------------------------------------
! @jasonfleming 20181129: Sample command line for compiling this program:
!
! gfortran -o adcircResultsComparison -I/usr/include -L/usr/lib adcircResultsComparison.F90 -lnetcdf -lnetcdff
!
! RENCI hatteras:
! ifort -o adcircResultsComparison -I/usr/share/Modules/software/CentOS-7/netcdf-Fortran/4.4.0_intel-18.0.0/include -L/usr/share/Modules/software/CentOS-7/netcdf-Fortran/4.4.0_intel-18.0.0/lib adcircResultsComparison.F90 -lnetcdf -lnetcdff
!-----------------------------------------------------------------------

module adcircCompare_module

#ifdef ADCNETCDF
   use NETCDF, only: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_INQUIRE, &
                     NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, NF90_GET_VAR, &
                     NF90_INQ_VARID, NF90_STRERROR, NF90_INQUIRE_VARIABLE, NF90_NOWRITE
#endif
   implicit none

#ifdef ADCNETCDF
   integer, parameter   :: nNetCDFVariables = 40
   character(200), save :: netcdf_types(nNetCDFVariables)
   character(200), save :: nc_longname(nNetCDFVariables)
   character(200), save :: nc_stdname(nNetCDFVariables)
#endif

   real(8), parameter   :: eps = epsilon(1d0)

   public

contains

   logical function findFile(filename, returnfalse)
      character(*), intent(IN) :: Filename
      logical, intent(IN), optional :: returnfalse
      logical                 :: exists

      inquire (FILE=trim(filename), EXIST=exists)

      FindFile = .true.
      if (.not. exists) then
         if (present(RETURNFALSE)) then
            if (returnfalse) then
               FindFile = .false.
               return
            else
               write (*, '(3A)') "Specified file ", trim(filename), &
                  " does not exist."
               call exit(1)
            end if
         else
            write (*, '(3A)') "Specified file ", trim(filename), &
               " does not exist."
            call exit(1)
         end if
      end if
   end function findFile

   integer function getFreeUnit()
      integer :: I
      logical :: ISOPEN
      I = 0
      do
         I = I + 1
         inquire (UNIT=I, OPENED=ISOPEN)
         if (ISOPEN) cycle
         GETFREEUNIT = I
         exit
      end do
      return
   end function getFreeUnit

#ifdef ADCNETCDF
   subroutine CHECK(stat, error, fatal)
      implicit none
      integer, intent(IN)             :: stat
      logical, intent(OUT), optional   :: error
      logical, intent(IN), optional    :: fatal
      logical                        :: fatal_local
#ifdef EBUG
      integer, allocatable            :: Dmy(:)
#endif

      if (.not. present(FATAL)) then
         FATAL_LOCAL = .true.
      else
         FATAL_LOCAL = FATAL
      end if

      if (fatal_local) then
         if (stat /= NF90_NOERR) then
            write (*, '(A,A)') "FATAL ERROR from ", trim(NF90_STRERROR(stat))
#ifdef EBUG
            !...Intentional segfault to trigger stack trace
            allocate (Dmy(1))
            Dmy(2) = 1
#endif
            call exit(1)
         else
            if (present(error)) error = .false.
         end if
      else
         if (stat /= NF90_NOERR) then
            if (present(error)) error = .true.
         else
            if (present(error)) error = .false.
         end if
      end if

   end subroutine CHECK
#endif

#ifdef ADCNETCDF
   subroutine initializeNetcdf()
      implicit none

      NETCDF_TYPES(:) = ""
      NETCDF_TYPES(1) = "sigmat" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(2) = "salinity" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(3) = "temperature" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(4) = "u-vel3D" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(5) = "v-vel3D" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(6) = "w-vel3D" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(7) = "q20" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(8) = "l" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(9) = "ev" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(10) = "qsurfkp1" !...Not implemented, placed here just to follow ADCIRC
      NETCDF_TYPES(11) = "zeta"
      NETCDF_TYPES(12) = "zeta_max"
      NETCDF_TYPES(13) = "u-vel"
      NETCDF_TYPES(14) = "v-vel"
      NETCDF_TYPES(15) = "vel_max"
      NETCDF_TYPES(16) = "pressure"
      NETCDF_TYPES(17) = "pressure_min"
      NETCDF_TYPES(18) = "windx"
      NETCDF_TYPES(19) = "windy"
      NETCDF_TYPES(20) = "wind_max"
      NETCDF_TYPES(21) = "radstress_x"
      NETCDF_TYPES(22) = "radstress_y"
      NETCDF_TYPES(23) = "radstress_max"
      NETCDF_TYPES(24) = "swan_HS"
      NETCDF_TYPES(25) = "swan_HS_max"
      NETCDF_TYPES(26) = "swan_DIR"
      NETCDF_TYPES(27) = "swan_DIR_max"
      NETCDF_TYPES(28) = "swan_TM01"
      NETCDF_TYPES(29) = "swan_TM01_max"
      NETCDF_TYPES(30) = "swan_TPS"
      NETCDF_TYPES(31) = "swan_TPS_max"
      NETCDF_TYPES(32) = "swan_windx"
      NETCDF_TYPES(33) = "swan_windy"
      NETCDF_TYPES(34) = "swan_wind_max"
      NETCDF_TYPES(35) = "swan_TM02"
      NETCDF_TYPES(36) = "swan_TM02_max"
      NETCDF_TYPES(37) = "swan_TMM10"
      NETCDF_TYPES(38) = "swan_TMM10_max"
      NETCDF_TYPES(39) = "winddrag"
      NETCDF_TYPES(40) = "dynamicWaterlevelCorrection"

      NC_LONGNAME(:) = ""
      NC_LONGNAME(1) = "water column vertically varying density"
      NC_LONGNAME(2) = "water column vertically varying salinity"
      NC_LONGNAME(3) = "water column vertically varying temperature"
      NC_LONGNAME(4) = "water column vertically varying east/west velocity"
      NC_LONGNAME(5) = "water column vertically varying north/south velocity"
      NC_LONGNAME(6) = "water column vertically varying up/down velocity"
      NC_LONGNAME(7) = "water column vertically varying turbulent kinetic energy"
      NC_LONGNAME(8) = "water column vertically varying mixing length"
      NC_LONGNAME(9) = "water column vertically varying eddy viscosity"
      NC_LONGNAME(10) = "sea surface temperature at the k+1 time level"
      NC_LONGNAME(11) = "water surface elevation above geoid"
      NC_LONGNAME(12) = "maximum water surface elevation above geoid"
      NC_LONGNAME(13) = "water column vertically averaged east/west velocity"
      NC_LONGNAME(14) = "water column vertically averaged north/south velocity"
      NC_LONGNAME(15) = "maximum water column vertically averaged velocity"
      NC_LONGNAME(16) = "air pressure at sea level"
      NC_LONGNAME(17) = "minimum air pressure at sea level"
      NC_LONGNAME(18) = "wind velocity in x-direction"
      NC_LONGNAME(19) = "wind velocity in y-direction"
      NC_LONGNAME(20) = "maximum wind velocity"
      NC_LONGNAME(21) = "radiation stress gradient x component"
      NC_LONGNAME(22) = "radiation stress gradient y component"
      NC_LONGNAME(23) = "maximum radiation stress gradient"
      NC_LONGNAME(24) = "significant wave height"
      NC_LONGNAME(25) = "maximum significant wave height"
      NC_LONGNAME(26) = "mean wave direction"
      NC_LONGNAME(27) = "maximum mean wave direction"
      NC_LONGNAME(28) = "mean absolute wave period"
      NC_LONGNAME(29) = "maximum TM01 mean wave period"
      NC_LONGNAME(30) = "smoothed peak period"
      NC_LONGNAME(31) = "maximum smoothed peak period"
      NC_LONGNAME(32) = "wind velocity in x-direction"
      NC_LONGNAME(33) = "wind velocity in y-direction"
      NC_LONGNAME(34) = "maximum wind stress"
      NC_LONGNAME(35) = "mean absoloute zero crossing period"
      NC_LONGNAME(36) = "maximum TM02 mean wave period"
      NC_LONGNAME(37) = "mean absolute wave period"
      NC_LONGNAME(38) = "maximum TMM10 mean wave period"
      NC_LONGNAME(39) = "wind drag"
      NC_LONGNAME(40) = "dynamic water surface correction above water level"

      NC_STDNAME(:) = ""
      NC_STDNAME(1) = "water_density_vertically_varying"
      NC_STDNAME(2) = "water_salinity_vertically_varying"
      NC_STDNAME(3) = "water_temperature_vertically_varying"
      NC_STDNAME(4) = "eastward_water_velocity_vertically_varying"
      NC_STDNAME(5) = "northward_water_velocity_vertically_varying"
      NC_STDNAME(6) = "upward_water_velocity_vertically_varying"
      NC_STDNAME(7) = "turbulent_kinetic_energy_vertically_varying"
      NC_STDNAME(8) = "water_mixing_length_vertically_varying"
      NC_STDNAME(9) = "water_eddy_viscosity_vertically_varying"
      NC_STDNAME(10) = "future sea surface temperature"
      NC_STDNAME(11) = "sea_surface_height_above_geoid"
      NC_STDNAME(12) = "maximum_sea_surface_height_above_geoid"
      NC_STDNAME(13) = "x_water_velocity_depth_averaged"
      NC_STDNAME(14) = "y_water_velocity_depth_averaged"
      NC_STDNAME(15) = "maximum_water_velocity_depth_averaged"
      NC_STDNAME(16) = "air_pressure_at_sea_level"
      NC_STDNAME(17) = "minimum_air_pressure_at_sea_level"
      NC_STDNAME(18) = "x_wind"
      NC_STDNAME(19) = "y_wind"
      NC_STDNAME(20) = "maximum_wind"
      NC_STDNAME(21) = "radiation_stress_gradient_x"
      NC_STDNAME(22) = "radiation_stress_gradient_y"
      NC_STDNAME(23) = "maximum_radiation_stress"
      NC_STDNAME(24) = "sea_surface_wave_significant_height"
      NC_STDNAME(25) = "maximum_sea_surface_wave_significant_height"
      NC_STDNAME(26) = "sea_surface_wave_to_direction"
      NC_STDNAME(27) = "maximum_sea_surface_wave_to_direction"
      NC_STDNAME(28) = "sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"
      NC_STDNAME(29) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"
      NC_STDNAME(30) = "sea_surface_wave_period_at_variance_spectral_density_maximum"
      NC_STDNAME(31) = "maximum_sea_surface_wave_period_at_variance_spectral_density_maximum"
      NC_STDNAME(32) = "x_wind"
      NC_STDNAME(33) = "y_wind"
      NC_STDNAME(34) = "maximum_wind"
      NC_STDNAME(35) = "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment"
      NC_STDNAME(36) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment"
      NC_STDNAME(37) = "sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"
      NC_STDNAME(38) = "maximum_sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"
      NC_STDNAME(39) = "wind_drag"
      NC_STDNAME(40) = "dynamic_sea_surface_correction_above_water_level"

      return

   end subroutine initializeNetcdf
#endif

#ifdef ADCNETCDF
   subroutine findMyNetCDFVariable(NCID)
      implicit none
      integer, intent(IN)  :: NCID
      integer :: I
      integer :: J
      integer :: NVAR
      character(200) :: NC_NAME

      call CHECK(NF90_INQUIRE(NCID, NVARIABLES=NVAR))

      call initializeNetcdf()
      do I = 1, NVAR
         call CHECK(NF90_INQUIRE_VARIABLE(NCID, I, &
                                          NAME=NC_NAME))
         do J = 1, size(NETCDF_TYPES)
            if (trim(adjustl(NC_NAME)) == &
                trim(adjustl(NETCDF_TYPES(J)))) then
               return
            end if
         end do
         if (I == NVAR) then
            write (*, '(A)') &
               "ADCIRC NetCDF Variable not found in file."
            call exit(1)
         end if
      end do

   end subroutine findMyNetCDFVariable
#endif

#ifdef ADCNETCDF
   subroutine getNetCDFVarId(NCID, VARID1, VARID2, NCOLS, VarName1, VarName2)
      implicit none
      integer, intent(IN)  :: NCID
      integer, intent(OUT) :: VARID1
      integer, intent(OUT) :: VARID2
      integer, intent(OUT) :: NCOLS
      character(256), intent(OUT), optional :: VarName1, VarName2
      integer             :: I
      integer             :: J
      integer             :: VARPOS2
      integer             :: NVAR
      character(200)      :: NC_NAME

      call CHECK(NF90_INQUIRE(NCID, NVARIABLES=NVAR))

      do I = 1, NVAR
         call CHECK(NF90_INQUIRE_VARIABLE(NCID, I, &
                                          NAME=NC_NAME))
         do J = 1, size(NETCDF_TYPES)
            if (trim(NC_NAME) == trim(NETCDF_TYPES(J))) then
               VARID1 = I
               if (present(VarName1)) then
                  VarName1 = NC_NAME
               end if
               select case (J)
               case (13, 18, 21, 33)
                  NCOLS = 2
                  VARPOS2 = J + 1
                  call CHECK(NF90_INQ_VARID(NCID, &
                                            trim(NETCDF_TYPES(VARPOS2)), VARID2))
                  if (present(VarName2)) then
                     VarName2 = NETCDF_TYPES(VARPOS2)
                  end if
               case DEFAULT
                  NCOLS = 1
                  VARID2 = -1
               end select
               return
            end if
         end do
      end do
      call exit(1)

   end subroutine getNetCDFVarId
#endif

   subroutine showHelp()
      implicit none

      write (*, '(A)') "Usage: adcircCompare [OPTIONS]"
      write (*, '(A)') ""
      write (*, '(A)') "Options:"
      write (*, '(A)') "  -h,--help          Displays this help"
      write (*, '(A)') "  -f1,--file1        File 1 for comparison"
      write (*, '(A)') "  -f2,--file2        File 2 for comparison"
      write (*, '(A)') "  -t,--tolerance     Maximum allowable difference in results"
      write (*, '(A)') "  -w,--wetdry        Ignore wet/dry differences"
      write (*, '(A)') "  -v,--verbose       Verbose output"
      write (*, '(A)') "  -c,--continue      Continue even when output is considered to have differences"
   end subroutine showHelp

   subroutine processCommandLineArgs(file1, file2, tolerance, wetdry, minmax, verbose, cont)
      implicit none

      character(256), intent(OUT)   :: file1, file2
      real(8), intent(OUT)        :: tolerance
      logical, intent(OUT)        :: wetdry
      logical, intent(OUT)        :: minmax !.true. if this is a min or max file (like maxvel.63)
      logical, intent(OUT)        :: verbose
      logical, intent(OUT)        :: cont
      integer                    :: iargc
      integer                    :: I
      character(2000)            :: CMD
      logical                    :: foundFile1, foundFile2, exists

      tolerance = 0.00001d0
      wetDry = .false.
      minmax = .false.
      verbose = .false.
      cont = .false.
      foundFile1 = .false.
      foundFile2 = .false.

      if (command_argument_count() == 0) then
         write (*, '(A)') "ERROR: No command line arguments specified."
         call showHelp()
         call exit(1)
      end if

      I = 0
      do while (I < command_argument_count())

         I = I + 1
         call get_command_argument(I, CMD)

         select case (trim(CMD))
         case ("-f1", "--file1")
            I = I + 1
            call get_command_argument(I, file1)
            foundFile1 = .true.
         case ("-f2", "--file2")
            I = I + 1
            call get_command_argument(I, file2)
            foundFile2 = .true.
         case ("-t", "--tolerance")
            I = I + 1
            call get_command_argument(I, CMD)
            read (CMD, *) tolerance
         case ("-w", "--wetdry")
            wetdry = .true.
         case ("-m", "--minmax")
            minmax = .true.
         case ("-h", "--help")
            call showHelp()
            call exit(0)
         case ("-v", "--verbose")
            verbose = .true.
         case ("-c", "--continue")
            cont = .true.
         case DEFAULT
            write (*, '(2A)') "ERROR: Unrecognized argument ", trim(CMD)
            call showHelp()
            call exit(1)
         end select
      end do

      if (.not. foundFile1) then
         write (*, '(A)') "ERROR: File 1 not specified."
         call showHelp()
         call exit(1)
      end if

      if (.not. foundFile2) then
         write (*, '(A)') "ERROR: File 2 not specified."
         call showHelp()
         call exit(1)
      end if

      exists = FindFile(file1)

      exists = FindFile(file2)

      return

   end subroutine processCommandLineArgs

   integer function determineFileType(filename)
      implicit none
      character(*), intent(IN) :: filename

      !...Check the various output formats to determine
      !   the type of file that is specified
#ifdef ADCNETCDF
      !...Check 1: netcdf
      if (checkFormat_netcdf(filename)) then
         determineFileType = 3
         return
      end if

      !...Check 2: xdmf
      if (checkFormat_xdmf(filename)) then
         determineFileType = 7
         return
      end if
#endif
      !...Check 3: full format ascii
      if (checkFormat_fullFormatAscii(filename)) then
         determineFileType = 1
         return
      end if

      !...Check 4: sparse format ascii
      if (checkFormat_sparseFormatAscii(filename)) then
         determineFileType = 4
         return
      end if

      determineFileType = -9999
      return

   end function determineFiletype

#ifdef ADCNETCDF
   logical function checkFormat_netcdf(filename)
      implicit none
      character(*), intent(IN) :: filename
      integer                 :: ncid
      integer                 :: ierr

      ierr = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
      if (ierr == NF90_NOERR) then
         checkFormat_netCDF = .true.
         return
      else
         checkFormat_netCDF = .false.
      end if
      return

   end function checkFormat_netcdf

   logical function checkFormat_xdmf(filename)
      implicit none
      character(*), intent(IN) :: filename

      checkFormat_xdmf = .false.
      write (*, '(A)') "WARNING: Check for XDMF not implemented: "//trim(filename)

   end function checkFormat_xdmf
#endif

   logical function checkFormat_fullFormatAscii(filename)
      implicit none
      character(*), intent(IN) :: filename
      integer                 :: io
      character(400)          :: header
      integer                 :: JunkI
      real(8)                 :: JunkR

      io = getFreeUnit()
      open (FILE=trim(filename), UNIT=io, ACTION="READ")
      read (io, '(A)') header
      read (io, '(A)') header
      read (io, '(A)') header
      close (io)

      read (header, *, ERR=100, end=100) JunkR, JunkI, JunkI, JunkR
      checkFormat_fullFormatAscii = .false.
      return
100   continue
      checkFormat_fullFormatAscii = .true.
      return
   end function checkFormat_fullFormatAscii

   logical function checkFormat_sparseFormatAscii(filename)
      implicit none
      character(*), intent(IN) :: filename
      integer                 :: io
      character(400)          :: header
      integer                 :: JunkI
      real(8)                 :: JunkR

      io = getFreeUnit()
      open (FILE=trim(filename), UNIT=io, ACTION="READ")
      read (io, '(A)') header
      read (io, '(A)') header
      read (io, '(A)') header
      close (io)

      read (header, *, ERR=100, end=100) JunkR, JunkI, JunkI, JunkR
      checkFormat_sparseFormatAscii = .true.
      return
100   continue
      checkFormat_sparseFormatAscii = .false.
      return
   end function checkFormat_sparseFormatAscii

   subroutine getMetadata(filename, noutfile, numsnaps, numnodes, ncol)
      implicit none
      integer, intent(IN)      :: noutfile
      character(*), intent(IN) :: filename
      integer, intent(OUT)     :: numsnaps, numnodes, ncol

      ncol = 0

      if (noutfile == 1 .or. noutfile == 4) then
         call getMetadataAscii(filename, numsnaps, numnodes, ncol)
      elseif (noutfile == 3) then
#ifdef ADCNETCDF
         call getMetadataNetCDF(filename, numsnaps, numnodes)
#else
         write (*, '(A)') "ERROR: must compile with netCDF support"
         call exit(1)
#endif
      else
         numsnaps = -9999
         numnodes = -9999
      end if

      return

   end subroutine getMetadata

   subroutine getMetadataAscii(filename, numsnaps, numnodes, numvalues)
      implicit none
      character(*), intent(IN) :: filename
      integer, intent(OUT)     :: numsnaps, numnodes, numvalues
      integer                 :: io
      integer                 :: junki
      real(8)                 :: junkr
      character(200)          :: header
      io = getFreeUnit()
      open (FILE=trim(filename), UNIT=io, ACTION="READ")
      read (io, '(A)') header
      read (io, '(A)') header
      close (io)
      ! FIXME: the numsnaps in the ascii output file header
      ! will be wrong if the file has been appended after
      ! hotstart (numsnaps is computed and written at cold
      ! start before any output data are added to the file;
      ! it is not updated/corrected if the run ends unexpectedly
      ! or if the run is hotstarted and more data are appended
      read (header, *) numsnaps, numnodes, junkR, junkI, numvalues
      return
   end subroutine getMetadataAscii

#ifdef ADCNETCDF
   subroutine getMetadataNetCDF(filename, numsnaps, numnodes)
      implicit none
      character(*), intent(IN) :: filename
      integer, intent(OUT)     :: numsnaps, numnodes
      integer                 :: ncid
      integer                 :: dimid_time
      integer                 :: dimid_node
      integer                 :: ierr

      call CHECK(NF90_OPEN(trim(filename), NF90_NOWRITE, ncid))
      call CHECK(NF90_INQ_DIMID(ncid, "time", dimid_time))

      ierr = NF90_INQ_DIMID(ncid, "node", dimid_node)
      if (ierr /= NF90_NOERR) then
         ierr = NF90_INQ_DIMID(ncid, "station", dimid_node)
      end if
      call CHECK(NF90_INQUIRE_DIMENSION(ncid, dimid_time, LEN=numsnaps))
      call CHECK(NF90_INQUIRE_DIMENSION(ncid, dimid_node, LEN=numnodes))
      call FindMyNetCDFVariable(ncid)
      call CHECK(NF90_CLOSE(ncid))

      return
   end subroutine getMetadataNetCDF
#endif

   subroutine openFile(filename, filetype, fileunit)
      implicit none
      character(*), intent(IN) :: filename
      integer, intent(IN)      :: filetype
      integer, intent(OUT)     :: fileunit
      character(200)          :: header

      if (filetype == 1 .or. filetype == 4) then
         fileunit = getFreeUnit()
         open (FILE=trim(filename), UNIT=fileunit, ACTION="READ")
         read (fileunit, *) header
         read (fileunit, *) header
      elseif (filetype == 3) then
#ifdef ADCNETCDF
         call CHECK(NF90_OPEN(trim(filename), NF90_NOWRITE, fileunit))
#endif
      else
         fileunit = -9999
      end if

      return
   end subroutine openFile

   subroutine closeFile(fileunit, filetype)
      implicit none
      integer, intent(IN) :: fileunit
      integer, intent(IN) :: filetype

      if (filetype == 1 .or. filetype == 4) then
         close (fileunit)
      elseif (filetype == 3) then
         !CALL CHECK(NF90_CLOSE(fileunit))
      else
         call exit(1)
      end if

      return
   end subroutine closeFile

   subroutine readNextSnapFullAscii(io_unit, nnodes, nvalues, nodaldata, time)
      implicit none
      integer, intent(IN)    :: io_unit
      integer, intent(IN)    :: nnodes
      integer, intent(IN)    :: nvalues
      real(8), intent(INOUT) :: nodaldata(:, :)
      real(8), intent(OUT)   :: time
      integer               :: timestep
      integer               :: I, J
      integer               :: junki

      read (io_unit, *) time, timestep
      do I = 1, nnodes
         read (io_unit, *) junkI, (nodalData(I, J), J=1, nvalues)
      end do

      return
   end subroutine readNextSnapFullAscii

   subroutine readNextSnapSparseAscii(io_unit, nvalues, nodaldata, time)
      implicit none
      integer, intent(IN)    :: io_unit
      integer, intent(IN)    :: nvalues
      real(8), intent(INOUT) :: nodaldata(:, :)
      real(8), intent(OUT)   :: time
      integer               :: node
      integer               :: timestep
      integer               :: nnondefault
      integer               :: I, J
      real(8)               :: defaultvalue

      read (io_unit, *) time, timestep, nnondefault, defaultvalue
      nodaldata(:, :) = defaultvalue
      do I = 1, nnondefault
         read (io_unit, *) node, (nodalData(node, J), J=1, nvalues)
      end do

      return
   end subroutine readNextSnapSparseAscii

#ifdef ADCNETCDF
   subroutine readNextSnapNetCDF(io_unit, snap, nnodes, nvalues, varid1, varid2, nodaldata)
      implicit none
      integer, intent(IN)    :: io_unit
      integer, intent(IN)    :: nnodes
      integer, intent(IN)    :: nvalues
      integer, intent(IN)    :: varid1, varid2
      integer, intent(IN)    :: snap
      real(8), intent(INOUT) :: nodaldata(:, :)

      call CHECK(NF90_GET_VAR(io_unit, varid1, nodaldata(:, 1), START=[1, snap], COUNT=[nnodes, 1]))
      if (nvalues == 2) then
         call CHECK(NF90_GET_VAR(io_unit, varid2, nodaldata(:, 2), START=[1, snap], COUNT=[nnodes, 1]))
      end if

      return
   end subroutine readNextSnapNetCDF
#endif

   subroutine compareData(nnodes, data1, data2, tolerance, wetdry, ndiff, avgdiff, maxdiff, nerror)
      implicit none
      integer, intent(IN)  :: nnodes
      real(8), intent(IN)  :: data1(:)
      real(8), intent(IN)  :: data2(:)
      real(8), intent(IN)  :: tolerance
      logical, intent(IN)  :: wetdry
      integer, intent(OUT) :: ndiff
      real(8), intent(OUT) :: avgdiff
      real(8), intent(OUT) :: maxdiff
      integer, intent(OUT) :: nerror
      integer             :: I
      real(8)             :: d

      ndiff = 0
      nerror = 0
      avgdiff = 0d0
      maxdiff = 0d0

      do I = 1, nnodes
         if (wetdry) then
            if (data1(I) <= -9999d0 .or. data2(I) <= -9999d0) then
               cycle
            end if
         end if
         d = data1(I) - data2(I)
         if (abs(d) > eps) then
            ndiff = ndiff + 1
         end if
         if (abs(d) > tolerance) then
            nerror = nerror + 1
         end if
         if (abs(d) > abs(maxdiff) .and. abs(d) > eps) then
            maxdiff = d
         end if
         avgdiff = avgdiff + d
      end do
      avgdiff = avgdiff/dble(nnodes)

      return
   end subroutine compareData

   subroutine displaySummary(nsnaps, nvalues, tolerance, wetdry, ndiff, avgdiff, maxdiff, nerror)
      implicit none
      integer, intent(IN) :: nsnaps
      integer, intent(IN) :: nvalues
      real(8), intent(IN) :: tolerance
      logical, intent(IN) :: wetdry
      real(8), intent(IN) :: avgdiff(:, :)
      real(8), intent(IN) :: maxdiff(:, :)
      integer, intent(IN) :: nerror(:, :)
      integer, intent(IN) :: ndiff(:, :)
      integer            ::I

      write (*, '(A)') "Summary of results file comparison"
      write (*, '(A)') ""
      write (*, '(A)') "|-----------------------------------------------------------------------|"
      write (*, '(A)') "|  TIMESNAP  |  NUMDIFF  |  AVERAGE DIFF  |  MAXIMUM DIFF  |  NUMERROR  |"
      write (*, '(A)') "|-----------------------------------------------------------------------|"
      do I = 1, nsnaps
         if (nvalues == 1) then
            if (ndiff(I, 1) > 0) then
               write (*, '(A,I9,A,I9,A,1PE12.4E3,A,1PE12.4E3,A,I9,A)') &
                  "|", I, "   | ", ndiff(I, 1), " |  ", avgdiff(I, 1), "  |  ", maxdiff(I, 1), "  | ", &
                  nerror(I, 1), "  |"
            end if
         elseif (nvalues == 2) then
            if (ndiff(I, 1) > 0 .or. ndiff(I, 2) > 0) then
               write (*, '(A,I9,A,I9,A,1PE12.4E3,A,1PE12.4E3,A,I9,A)') &
                  "|", I, "   | ", ndiff(I, 1) + ndiff(I, 2), " |  ", (avgdiff(I, 1) + avgdiff(I, 2))/2d0, &
                  "  |  ", max(abs(maxdiff(I, 1)), abs(maxdiff(I, 2))), "  | ", &
                  nerror(I, 1), "  |"
            end if
         end if
      end do
      write (*, '(A)') "|-----------------------------------------------------------------------|"
      write (*, '(A)') ""
      write (*, '(A,1PE12.4E3)') "  -->       Error Tolerance: ", tolerance
      if (wetdry) then
         write (*, '(A)') "  --> Wet/Dry Diffs Skipped: Yes"
      else
         write (*, '(A)') "  --> Wet/Dry Diffs Skipped: No"
      end if

      return
   end subroutine displaySummary

end module adcircCompare_module

program adcircResultCompare
   use adcircCompare_module, only: &
      processCommandLineArgs, determineFileType, &
      getMetadata, openFile, closeFile, &
      readNextSnapFullAscii, readNextSnapSparseAscii, &
      compareData, displaySummary

#ifdef ADCNETCDF
   use adcircCompare_module, only: readNextSnapNetCDF, getNetCDFVarId
#endif

   implicit none

   character(2000)     :: file1, file2
   integer             :: noutfile1, noutfile2
   integer             :: nsnapfile1, nsnapfile2
   integer             :: nnodefile1, nnodefile2
   integer             :: ncol1, ncol2
   integer             :: io_unit1, io_unit2
   integer             :: I
   integer, allocatable :: nerror(:, :), ndiff(:, :)
   real(8), allocatable :: maxdiff(:, :), avgdiff(:, :)
   real(8)             :: tolerance
   real(8)             :: time1, time2
   real(8), allocatable :: nodaldata1(:, :), nodaldata2(:, :)
   logical             :: wetdry, verbose, cont
   logical :: minmax ! .true. if this is a min or max file (e.g., maxvel.63)

#ifdef ADCNETCDF
   integer             :: varid11, varid12, varid21, varid22
#endif

   !...Process the command line arguments from the user
   call processCommandLineArgs(file1, file2, tolerance, wetdry, minmax, verbose, cont)

   !...Determine the file type of each file
   noutfile1 = determineFileType(file1)
   noutfile2 = determineFileType(file2)

   !...Pull the number of time snaps and nodes from each file
   call getMetadata(file1, noutfile1, nsnapfile1, nnodefile1, ncol1)
   call getMetadata(file2, noutfile2, nsnapfile2, nnodefile2, ncol2)

   !...Check for the number of snaps in each file
   if (nsnapfile1 /= nsnapfile2) then
      write (*, '(A)') "ERROR: There are different numbers of output snaps in the files."
      call exit(1)
   end if

   !...Check for the number of nodes in each file
   if (nnodefile1 /= nnodefile2) then
      write (*, '(A)') "ERROR: There are different numbers of nodes in the files."
      call exit(1)
   end if

   !...If this is an ascii minmax file, the 2nd dataset is the
   ! time of occurrence, which can have much larger differences
   ! than are allowed by the error tolerance, and since we don't
   ! have a good way to set the error tolerance fortime-of-peak
   ! in a way that will avoid false positives (i.e., tests
   ! failing unnecessarily) let's avoid the comparison of
   ! time-of-peak values for now.
   if (minmax .eqv. .true.) then
      nsnapfile1 = 1
      nsnapfile2 = 1
      write (*, '(A)') "INFO: The time-of-peak-occurrence values will not be compared, only the peak values themselves."
   end if

   io_unit1 = 0
   io_unit2 = 0
   call openFile(file1, noutfile1, io_unit1)
   call openFile(file2, noutfile2, io_unit2)

   !...For netCDF files, locate the netCDF variable ID's for later use
   if (noutfile1 == 3) then
#ifdef ADCNETCDF
      call getNetCDFVarId(io_unit1, varid11, varid12, ncol1)
#else
      write (*, '(A)') "ERROR: must compile with netCDF support"
      call exit(1)
#endif
   end if
   if (noutfile2 == 3) then
#ifdef ADCNETCDF
      call getNetCDFVarId(io_unit2, varid21, varid22, ncol2)
#else
      write (*, '(A)') "ERROR: must compile with netCDF support"
      call exit(1)
#endif
   end if

   !...Check if comparing like files
   if (ncol1 /= ncol2) then
      write (*, '(A)') "ERROR: The number of values in the files is different."
      call exit(1)
   end if

   allocate (nodaldata1(1:nnodefile1, 1:ncol1))
   allocate (nodaldata2(1:nnodefile2, 1:ncol2))
   allocate (ndiff(1:nsnapfile1, 1:ncol1))
   allocate (nerror(1:nsnapfile1, 1:ncol1))
   allocate (avgdiff(1:nsnapfile1, 1:ncol1))
   allocate (maxdiff(1:nsnapfile1, 1:ncol1))

   !...Time comparison loop
   do I = 1, nsnapfile1

      if (verbose) then
         write (*, '(A,I6,A,I6)') "Processing snap ", I, " of ", nsnapfile1
      end if

      !...Read data from file 1
      if (noutfile1 == 1) then
         call readNextSnapFullAscii(io_unit1, nnodefile1, ncol1, nodaldata1, time1)
      elseif (noutfile1 == 4) then
         call readNextSnapSparseAscii(io_unit1, ncol1, nodaldata1, time1)
      elseif (noutfile1 == 3) then
#ifdef ADCNETCDF
         call readNextSnapNetCDF(io_unit1, I, nnodefile1, ncol1, varid11, varid12, nodaldata1)
#else
         write (*, '(A)') "ERROR: must compile with netCDF support"
         call exit(1)
#endif
      else
         call exit(1)
      end if

      !...Read data from file 2
      if (noutfile2 == 1) then
         call readNextSnapFullAscii(io_unit2, nnodefile2, ncol2, nodaldata2, time2)
      elseif (noutfile2 == 4) then
         call readNextSnapSparseAscii(io_unit2, ncol2, nodaldata2, time2)
      elseif (noutfile1 == 3) then
#ifdef ADCNETCDF
         call readNextSnapNetCDF(io_unit2, I, nnodefile2, ncol2, varid21, varid22, nodaldata2)
#else
         write (*, '(A)') "ERROR: must compile with netCDF support"
         call exit(1)
#endif
      else
         call exit(1)
      end if

      !...Perform comparison
      call compareData(nnodefile1, nodalData1(:, 1), nodalData2(:, 1), tolerance, wetdry, &
                       ndiff(I, 1), avgdiff(I, 1), maxdiff(I, 1), nerror(I, 1))
      if (ncol1 == 2) then
         call compareData(nnodefile1, nodalData1(:, 2), nodalData2(:, 2), tolerance, wetdry, &
                          ndiff(I, 2), avgdiff(I, 2), maxdiff(I, 2), nerror(I, 2))
      end if

      !...Check to stop code on any error
      if (.not. cont) then
         if (ncol1 == 2) then
            if (nerror(I, 1) > 0 .or. nerror(I, 2) > 0) then
               write (*, '(A,I0,A)') "Code stopped at snap ", I, "."
               call displaySummary(I, ncol1, tolerance, wetdry, ndiff, avgdiff, maxdiff, nerror)
               call closeFile(io_unit1, noutfile1)
               call closeFile(io_unit2, noutfile2)
               write (*, '(A)') "Results do not match within specified tolerance."
               call exit(1)
            end if
         else
            if (nerror(I, 1) > 0) then
               write (*, '(A,I0,A)') "Code stopped at snap ", I, "."
               call displaySummary(I, ncol1, tolerance, wetdry, ndiff, avgdiff, maxdiff, nerror)
               call closeFile(io_unit1, noutfile1)
               call closeFile(io_unit2, noutfile2)
               write (*, '(A)') "Results do not match within specified tolerance."
               call exit(1)
            end if
         end if
      end if

   end do

   call closeFile(io_unit1, noutfile1)
   call closeFile(io_unit2, noutfile2)

   call displaySummary(nsnapfile1, ncol1, tolerance, wetdry, ndiff, avgdiff, maxdiff, nerror)
   do I = 1, nsnapfile1
      if (nerror(I, 1) > 0) then
         write (*, '(A)') "Results do not match within specified tolerance."
         call exit(1)
      end if
      if (ncol1 == 2) then
         if (nerror(I, 2) > 0) then
            write (*, '(A)') "Results do not match within specified tolerance."
            call exit(1)
         end if
      end if
   end do
   write (*, '(A)') "Results match within specified tolerance."

   deallocate (nodaldata1)
   deallocate (nodaldata2)
   deallocate (ndiff)
   deallocate (nerror)
   deallocate (avgdiff)
   deallocate (maxdiff)

end program adcircResultCompare
