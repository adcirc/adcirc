!-------------------------------------------------------------------------------!
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
!-------------------------------------------------------------------------------!
<<<<<<< HEAD
#include "logging_macros.h"
C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C                                                                             *
C This module handles the model initialization for a cold start.  The primary *
C 2D coldstart initialization is handled in SUBROUTINE COLDSTART.  The primary*
C 3D initialization is handled in SUBROUTINE COLDSTART_3D.                    *
C                                                                             *
C******************************************************************************

C
      SUBROUTINE COLDSTART()
C
      USE SIZES
      USE GLOBAL, only: DP1, DP2, eta1, qx2, qy2, tk, etas, eta2, ch1,
=======
#include "logging_macros.h"
!-------------------------------------------------------------------------
!> PADCIRC VERSION 45.12 03/17/2006
!> last changes in this file VERSION 45.12
!>
!> @brief This module handles adcirc initialization for a cold start.
!>
!> @todo This file has no module structure.
!>       Also need to update global file init to move from open_gbl_file
!>       and headerXX to the newer writeDomainHeader subroutine.
!-------------------------------------------------------------------------



!-------------------------------------------------------------------------
!> @brief: Handles all subroutine calls for a 2d coldstart initialization.
!-------------------------------------------------------------------------
      subroutine coldstart()
!-------------------------------------------------------------------------

      use sizes
      use global, only: DP1, DP2, eta1, qx2, qy2, tk, etas, eta2, ch1,
>>>>>>> Subject: Add additional source terms and update style
     &    en2, nddt, btime1, statim, btime2, btiminc, btime_end, qn1,
     &    qn0, nffr, ftiminc, qnin1, enin1, qnin2, enin2, qtime1, nws,
     &    qtime2, noute, noutv, noutc, noutge, noutgv, noutgc, nrs, qn2,
     &    noute, noutm, noutgw, metonly, rstime1, rstime2, rstiminc,
     &    rsnx1, rsny1, rsnx2, rsny2, ncice, cice1, cice2, cice_time1,
     &    cice_time2, cice_timinc, ctip, tip2, nscoue, nstae_g,
     &    nstae, rdes4, rid4, rdes8, rid8, aid8, esbin1, vv2, nolifa,
<<<<<<< HEAD
     &    uu2, inundationoutput, esbin2, nspoolv, nscouc, icstp, nstac,
     &    nstac_g, ntrspc, nspoolc, nscoum, bchgtiminc,
     &    h0, ntrspe, dt, nspoole, nscouv, np_g,
     &    outputsponge, ne_g, usingDynamicWaterLevelCorrection, ndsetse,
=======
     &    uu2, inundationOutput, esbin2, nspoolv, nscouc, nstac,
     &    nstac_g, ntrspc, nspoolc, nscoum, bchgtiminc, h0, ntrspe, dt, 
     &    nspoole, nscouv, np_g, outputSponge, ne_g,  ndsetse, 
     &    usingDynamicWaterLevelCorrection, nspoolge, nscougv, ndsetsv,
>>>>>>> Subject: Add additional source terms and update style
     &    nspoolgv, use_tvw, nscou_tvw, nout_tvw, ndsets_tvw, nspool_tvw,
     &    nscougw, ndsetsw, nspoolgw, nscougc, c2ddi, iden, 
<<<<<<< HEAD
     &    rundes, runid, dtdp, ntrspm, nspoolm, nscoub, ibstp, nscouge,
     &    ndsetsc, nspoolgc, nscougb, igbp, c3d, nstav_g, nstav,
     &    ntrspv, iicestp, nstam_g, nstam, cbaroclinic
      use mod_io, only: openFileForRead
      use boundaries, only: csii, siii, nbv
      use mesh, only: sfcx, sfcy, ycsfac, sfmx
      use mod_normal_flow_boundary, only: rotate_normal_flux
      use mod_logging, only: t_log_scope, init_log_scope,
     &    allMessage, DEBUG, screenUnit, nscreen
      USE GLOBAL_3DVS, ONLY: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT
      USE MESH, ONLY : NE, NP, X, Y, SLAM, SFEA, ICS, DP, NM, MJU, AID4,
=======
     &    rundes, runid, dtdp, ntrspm, nspoolm, nscoub, nscouge,
     &    ndsetsc, nspoolgc, nscougb, igbp, c3d, nstav_g, nstav,
     &    ntrspv, nstam_g, nstam, openfileforread, cbaroclinic
      use mod_io, only: openFileForRead
      use boundaries, only: csii, siii, nbv, neta, nfluxf, nope, nvel,
     &    lbcodei, npebc
      use mod_normal_flow_boundary, only: rotate_normal_flux
      use mod_logging, only: t_log_scope, init_log_scope,
     &    allMessage, DEBUG, screenUnit, nscreen
      use global_3Dvs, only: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT
      use mesh, only : ne, np, x, y, slam, sfea, ics, dp, nm, mju, aid4,
<<<<<<< HEAD
>>>>>>> Subject: Add additional source terms and update style
     &    agrid
=======
     &    agrid, sfcx, sfcy, ycsfac, sfmx, totalareacalc
>>>>>>> Subject: changes necessary for compilation after merging source terms
      use harm
      use wetdry, only : initializeWettingAndDrying
      use wind, only : coldStartMeteorologicalForcing, rsget
      use global_io, only: header63, header64, header77, header73,
     &    header74, header93, header91, header61, header62, header71,
     &    header72, header81, writedomainheader, open_gbl_file, nddt2get
<<<<<<< HEAD
      USE SUBDOMAIN, ONLY : subdomainOn, readFort015
      USE mod_terminate, ONLY : terminate, ADCIRC_EXIT_FAILURE
      USE OWI_ICE, ONLY : NCICE1_INIT, NCICE1_GET  !tcm v49.64.01

      USE RS2,ONLY : RS2INIT,RS2GET            ! sb added 09/xx/2006
      USE NODALATTRIBUTES, ONLY : outputTau0
      USE WRITE_OUTPUT, only : outputNodeCode, outputNOFF,
=======
      use subdomain, only : subdomainOn, readFort015
      USE mod_terminate, ONLY : terminate, ADCIRC_EXIT_FAILURE
      use owi_ice, only : ncice1_init, ncice1_get  !tcm v49.64.01
      use rain, only : useRain, coldstartRainForcing  !cms 11/20/25
      use hydrology, only : useInteriorHydrology,
     &    coldstartInteriorHydrology    ! cms 11/20/25
      use rs2,only : rs2init,rs2get            ! sb added 09/xx/2006
      use nodalattributes, only : outputTau0
      use write_output, only : outputNodeCode, outputNoff,
>>>>>>> Subject: Add additional source terms and update style
     &    initInundationOutput
      use time_varying_weir_boundary,only : weir_setup
      use gwce, only : etime1, etime2, etiminc
#ifdef adcnetcdf
      use couple2bc3d, only: initial_bc3d_netcdf
#endif
<<<<<<< HEAD
      IMPLICIT NONE

      INTEGER I, J,IOS, IJK
      INTEGER NM1, NM2, NM3
      INTEGER NC1, NC2, NC3, NCELE

      REAL(8) :: HTOT
      REAL(8) :: QNIN_TMP
!kmd - Evan's updates for rivers above MSL
!      REAL(8),ALLOCATABLE :: et_WSE(:)
      LOGICAL FOUND

C     jgf48.4628 Add capability to turn off solution if only met output
C     was requested.
      LOGICAL metOutputOn     ! .true. if met output was requested
      LOGICAL nonMetOutputOff ! .true. no other output was requested

c.... tcm v50.66.02 additions for time varying bathymetry
      CHARACTER(80) :: CDUM80
C

      LOG_SCOPE_TRACED("coldstart", CSTART_TRACING)
      metOutputOn = .FALSE.
      nonMetOutputOff = .FALSE.
      !
=======
      implicit none

      integer i, j,ios, ijk
      integer nm1, nm2, nm3
      integer nc1, nc2, nc3, ncele

      real(8) :: htot
      real(8) :: qnin_tmp
      logical found

      ! jgf48.4628 add capability to turn off solution if only met output
      ! was requested.
      logical metOutputOn     ! .true. if met output was requested
      logical nonMetOutputOff ! .true. no other output was requested

      ! tcm v50.66.02 additions for time varying bathymetry
      character(80) :: cdum80

      LOG_SCOPE_TRACED("coldstart", CSTART_TRACING)
      metOutputOn = .false.
      nonMetOutputOff = .false.

>>>>>>> Subject: Add additional source terms and update style
      call initializeWettingAndDrying()

      ! set up time varying weirs if necessary
      call weir_setup()

      if (subdomainOn) then
         call readFort015() ! ncsu subdomain modeling
      endif

<<<<<<< HEAD
C------------------------------------------------------------------------------
C... tcm v50.66.01 Time Varying Bathymetry
C------------------------------------------------------------------------------
C
      IF (NDDT.NE.0) then
C...    Set the first two bathymetry values, these arrays are
C...    allocated in alloc_main13 called during read_input
         DP1(:) = DP(:)  !set dp1 to initial depth from fort.14 grid
         DP2(:) = DP(:)  !set dp2 to the same values

         BTIME1 = STATIM*86400.D0
         BTIME2=BTIME2+BTIMINC  !new ending time for this record
         BTIME_END = BTIME1 + BCHGTIMINC

C...    read in the second depth dp2 for btime2
         IOS = 0
         CALL openFileForRead(141,TRIM(INPUTDIR)//'/'//'fort.141',IOS)
         IF (IOS.GT.0) THEN
            CALL terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &         message='Unable to open fort.141 for time '
     &                //'varying bathymetry')
         ENDIF

         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1976)
         WRITE(16,1976)
 1976    FORMAT(/,1X,'TIME VARYING BATHYMETRY INFORMATION ',
     &        'READ FROM UNIT 141',/)

         IF (ABS(NDDT).EQ.1) THEN
            DO I=1,NP
               READ(141,*) IJK,DP2(IJK)
            ENDDO
         ENDIF
         IF (ABS(NDDT).EQ.2) THEN
!!!        go get new record for only some nodes, all
!!!        other nodes keep their current value
            CALL NDDT2GET(141,DP2(:),-99999.d0)
         ENDIF

C...     IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...     DEPTHS ARE > OR = TO H0.
         IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
            DO I=1,NP
               IF (DP2(I).LT.H0) DP2(I) = H0
            ENDDO
         ENDIF
      ENDIF  !end test for time varying bathymetry
C
C------------------------------------------------------------------------------
C
C...  SET AT REST INITIAL CONDITION OVER WHOLE DOMAIN
C...  IF BOTTOM IS ABOVE THE GEIOD -> DRY NODE
C...  IF BOTTOM IS INITIALLY BELOW THE GEOID AND STARTDRY=1 -> DRY NODE
C...
      UU2(:) = 0.D0
      VV2(:) = 0.D0
      QX2(:) = 0.D0
      QY2(:) = 0.D0
      TK(:) = 0.D0
      ETA1(:) = 0.D0
      ETA2(:) = 0.D0
      ETAS(:) = 0.D0
      CH1(:) = 0.d0
      EN2(:) = 0.d0
      !
      ! Set the initial wet/dry state and water surface elevation
=======
      ! tcm v50.66.01 additions for time varying bathymetry
      if (nddt.ne.0) then
         ! set the first two bathymetry values, these arrays are
         ! allocated in alloc_main13 called during read_input
         dp1(:) = dp(:)  !set dp1 to initial depth from fort.14 grid
         dp2(:) = dp(:)  !set dp2 to the same values

         btime1 = statim*86400.d0
         btime2=btime2+btiminc  !new ending time for this record
         btime_end = btime1 + bchgtiminc

         ! read in the second depth dp2 for btime2
         ios = 0
         call openfileforread(141,trim(inputdir)//'/'//'fort.141',ios)
         if (ios.gt.0) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &         message='Unable to open fort.141 for time '
     &                //'varying bathymetry')
         endif

         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1976)
         write(16,1976)
 1976    format(/,1x,'Time varying bathymetry INFOrmation ',
     &        'read from unit 141',/)

         if (abs(nddt).eq.1) then
            do i=1,np
               read(141,*) ijk,dp2(ijk)
            enddo
         endif
         if (abs(nddt).eq.2) then
            ! go get new record for only some nodes, all
            ! other nodes keep their current value
            call nddt2get(141,dp2(:),-99999.d0)
         endif

         !if wetting and drying will not be used, make sure all
         ! depths are > or = to h0.
         if ((nolifa.eq.0).or.(nolifa.eq.1)) then
            do i=1,np
               if (dp2(i).lt.h0) dp2(i) = h0
            enddo
         endif
      endif  !end test for time varying bathymetry


      ! set at rest initial condition over whole domain
      ! if bottom is above the geoid -> dry node
      ! if bottom is initially below the geoid and startdry=1 -> dry node
      uu2(:) = 0.d0
      vv2(:) = 0.d0
      qx2(:) = 0.d0
      qy2(:) = 0.d0
      tk(:) = 0.d0
      eta1(:) = 0.d0
      eta2(:) = 0.d0
      etas(:) = 0.d0
      ch1(:) = 0.d0
      en2(:) = 0.d0
 
      ! Compute initial wet/dry state and water surface elevation:
      ! Checks for nodes that should be dry and updates nnodecode
      !cms: mvd nodecode = nnodecode here from totalAreaCalc
>>>>>>> Subject: Add additional source terms and update style
      call setColdWetDryStateAndWaterSurfaceElevation()

      ! cms 11/20/25: Compute active elements and total area for each
      ! node, mvd call from adcirc_init since needed within c/hstart for
      ! hydrology. nodecode reset mvd from totalAreaCalc to setColdWet...
      call totalareacalc()

      ! jgf52.08.01: initialize inundation arrays if necessary. this
      ! can only be called after nnodecode has been set.
      if (inundationOutput.eqv..true.) then
         call initInundationOutput()
      endif

      ! Initialize the elevation specified boundary condition
      ! if it requires the use of the unit 19 file.
      if (npebc) then
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1977)
         write(16,1977)
 1977    format(/,1x,'Elevation specified INFOrmation read from unit ',
     &        '19',/)
         ! NCSU Subdomain Modeling
         if (subdomainOn.eqv..false.) then
            open(19,file=trim(inputdir)//'/'//'fort.19')
            read(19,*) etiminc
            do j=1,neta
               read(19,*) esbin1(j)
            end do
            do j=1,neta
               read(19,*) esbin2(j)
            end do
            etime1 = statim*86400.d0
            etime2 = etime1 + etiminc
         endif
      endif

      ! Initialize the Normal Flow Boundary Condition
      qn2=0.d0
      qn1=0.d0
      qn0=0.d0

      if(nfluxf.eq.1.and.nffr.le.0) then
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1979)
         write(16,1979)
 1979    format(/,1x,'Normal flow INFOrmation read from unit 20',/)
         open(20,file=trim(inputdir)//'/'//'fort.20')
         read(20,*) ftiminc

         qnin1 = 0d0
         qnin2 = 0d0
         do j=1,nvel
            if ((lbcodei(j).eq.2).or.(lbcodei(j).eq.12).or.(lbcodei(j).eq.22)) then
               read(20,*) qnin_tmp
               qnin1(j) = rotate_normal_flux(ics, j, qnin_tmp)
            elseif (lbcodei(j).eq.32) then
               read(20,*) qnin_tmp,enin1(j)
               qnin1(j) = rotate_normal_flux(ics, j, qnin_tmp)
            endif
         end do
         do j=1,nvel
            if ((lbcodei(j).eq.2).or.(lbcodei(j).eq.12).or.(lbcodei(j).eq.22)) then
               read(20,*) qnin_tmp
               qnin2(j) = rotate_normal_flux(ics, j, qnin_tmp)
            elseif (lbcodei(j).eq.32) then
               read(20,*) qnin_tmp,enin2(j)
               qnin2(j) = rotate_normal_flux(ics, j, qnin_tmp)
            endif
         end do
         qtime1 = statim*86400.d0
         qtime2 = qtime1 + ftiminc
      endif


      ! Input meteorological Information from unit 22 or unit 200 series.
      ! If fleet numeric wind data is used, find beginning time in file,
      ! NOTE: can't deal with wind that starts after wreftim!!!!!!!!!!!!
      ! Read in and interpolate in space onto the adcirc grid both the
      ! time level 1 and level 2 wind fields.

      if(nws.ne.0) then
         ! jgf48.4627 If user only requests meteorological output, then
         ! set a flag that turns off the calculation of the gwce and
         ! momentum equations.
         if ((noute.eq.0) .and. (noutv.eq.0) .and. (noutc.eq.0) .and.
     &       (noutge.eq.0) .and. (noutgv.eq.0) .and. (noutgc.eq.0 ).and.
     &       (i3dsd.eq.0) .and. (i3dsv.eq.0) .and. (i3dst.eq.0) .and.
     &       (i3dgd.eq.0) .and. (i3dgv.eq.0) .and. (i3dgt.eq.0) .and.
     &       (nhase.eq.0) .and. (nhasv.eq.0) .and.
     &       (nhage.eq.0) .and. (nhagv.eq.0) )
     &       then
            nonMetOutputOff = .true.
         endif
         if ( (noutm.ne.0) .or. (noutgw.ne.0) ) then
            metOutputOn = .true.
         endif
         if ( metOutputOn .and. nonMetOutputOff ) then
            if ((nscreen.ne.0).and.myproc.eq.0) then
               write(screenunit,*)
     &      'INFO: only meteorological output was requested. ADCIRC'
               write(screenunit,*)
     &      'will not solve the gwce or momentum equations.'
             endif
             write(16,*)
     &      'INFO: only meteorological output was requested. ADCIRC'
             write(16,*)
     &      'will not solve the gwce or momentum equations.'
            metOnly = .true. ! set flag
         endif

         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1980)
         write(16,1980)
 1980    format(/,1x,'wind (and pressure) INFOrmation read.',/)

         call coldStartMeteorologicalForcing(nws)

      endif


      ! Input Radiation Stress Information from unit 23
      !   read in the time level 1 and level 2 fields
      if(nrs.ge.1) then ! sb46.28sb03
         rstime1 = statim*86400.d0
         rstime2 = rstime1+rstiminc
         if(nrs.eq.1) then
           open(23,file=trim(inputdir)//'/'//'fort.23')
           call rsget(rsnx1,rsny1)
           call rsget(rsnx2,rsny2)
         !endif  !cms see comment below
         ! nrs=2 was added. 09/xx/2006 sb
         elseif(nrs.eq.2) then !cms 11/20/25 changed to elseif
           call rs2init(rsnx1,rsny1,np)
           call rs2get(rsnx1,rsny1,np)
           call rs2get(rsnx2,rsny2,np)
         endif
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1981)
         write(16,1981)
 1981    format(/,1x,'Radiation stress INFOrmation read.',/)
      endif

      ! tcm: v49.64.01
      ! Input Ice Concentration Information from unit 25
      !   read in the time level 1 and level 2 fields
      if (ncice.ge.1) then
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,1112)
         write(16,1112)
         if(nscreen.ne.0.and.myproc.eq.0) write(screenunit,2980)
         write(16,2980)
 2980    format(/,1x,'ice concentration INFOrmation read.',/)

         if(ncice.eq.12) then
            call ncice1_init(cice1,np)
            call ncice1_get(cice1,np)
            call ncice1_get(cice2,np)
            cice_time1 = statim*86400.d0
            cice_time2 = cice_time1 + cice_timinc
         endif
      endif

      ! cms: v56 11/20/2025
      ! Additional Source Terms (Rain, Interior Hydrology)
      if (useRain) then
         if (nscreen.ne.0.and.myproc.eq.0) write(Screenunit,1112)
         write(16,1112)
         if (nscreen.ne.0.and.myproc.eq.0) write(Screenunit,2990)
         write(16,2990)
 2990    format(/,1X,'Reading rainfall intensities.',/)
         call coldStartRainForcing()
      end if

      if (useInteriorHydrology) then
         if (nscreen.ne.0.and.myproc.eq.0) write(Screenunit,1112)
         write(16,1112)
         if (nscreen.ne.0.and.myproc.eq.0) write(Screenunit,2991)
         write(16,2991)
 2991    format(/,1X,'Reading Interior Hydrology Discharges.',/)
         call coldStartInteriorHydrology()
      END IF



#ifdef adcnetcdf
      ! Coupling to 3d Baroclinic Mode
      if (c2ddi.and.cbaroclinic.and.abs(iden).ge.5) then
         call initial_bc3d_netcdf(statim*86400.d0)
      endif
#endif


      ! Initialize Tidal Potential Forcing
      if (ctip) then
         do i=1,np
            tip2(i)=0.0
         end do
      endif


      !------------------------------------------------------------------
      !  BEGIN OUTPUT INITIALIZATION
      !
      !  cms: 11/20/25 for both station and global output
      !       * Grouped together by control
      !       * Simplified logic for station and global
      !         since all logical values are either T or F
      !           if(variable.eqv..false) -> if(.not. variable)
      !           if(variable.eqv..true) -> if(variable)
      !------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Initialize Elevation Station Output Spool
      !-----------------------------------------------------------------
      ! Multiple files are controlled by the elevation output, for each:
      !   Open station output file
      !   Write out header INFOrmation including:
      !      ntrspe (no. of data pts. at each elevation station)
      !      nstae, dt*nspoole, nspoole, irtype
      ! Variables include:
      !   elevation (fort.61)
      !   bathymetry (fort.75) -- tcm v50.66.01
      !   rain.61 (lun=94) -- cms v56.02 
      !   dynamicWaterLevelCorrection.61  (lun=109) -- jgf
      !-----------------------------------------------------------------

      nscoue=0
      nscoub=0          !tcm: add for bathymetry

      if(abs(noute).eq.1) then
         if (.not. write_local_files) then
            call open_gbl_file(61, trim(globaldir)//'/'//'fort.61',
     &       nstae_g, nstae, header61)
            if (nddt.ne.0) then
               call open_gbl_file(75, trim(globaldir)//'/'//'fort.75',
     &            nstae_g, nstae, header61)
            end if
            if (useRain) then
               call open_gbl_file(94, trim(globaldir)//'/'//'rain.61',
     &         nstae_g, nstae, header61)
            end if
            if (usingDynamicWaterLevelCorrection) then
               call open_gbl_file(109, trim(globaldir)//'/'//
     &         'dynamicWaterLevelCorrection.61', nstae_g, nstae, header61)
            end if
         else  ! Write Local Files
            call open_gbl_file(61, trim(localdir)//'/'//'fort.61',
     &       nstae_g, nstae, header61)
            if (nddt.ne.0) then
               call open_gbl_file(75, trim(localdir)//'/'//'fort.75',
     &            nstae_g, nstae, header61)
            end if
            if (useRain) then
               call open_gbl_file(94, trim(localdir)//'/'//'fort.94',
     &         nstae_g, nstae, header61)
            end if
            if (usingDynamicWaterLevelCorrection) then
               call open_gbl_file(109, trim(localdir)//'/'//
     &         'dynamicWaterLevelCorrection.61', nstae_g, nstae, header61)
            end if
         endif
      endif


      ! Initialize Velocity Station Spool Counter
      !   open velocity station output file
      !   write out header INFOrmation including ntrspv (no. of data pts.
      !      at each velocity station), nstav, dt*nspoolv,
      !      nspoolv, irtype
      !-----------------------------------------------------------------

      nscouv=0

      if(abs(noutv).eq.1) then
         if (.not. write_local_files) then
            call open_gbl_file(62, trim(globaldir)//'/'//'fort.62',
     &          nstav_g, nstav, header62)
         else
            call open_gbl_file(62, trim(localdir)//'/'//'fort.62',
     &          nstav_g, nstav, header62)
         endif
      endif


      ! Initialize Concentration Station Spool Counter
      !   open concentration station output file
      !   write out header INFOrmation including ntrspc (no. of data pts.
      !      at each concentration station), nstac, dt*nspoolc,
      !      nspoolc, irtype
      !-----------------------------------------------------------------

      nscouc=0

      if(abs(noutc).eq.1) then
         if (.not. write_local_files) then
            call open_gbl_file(81, trim(globaldir)//'/'//'fort.81',
     &          nstac_g, nstac, header81)
         else
            call open_gbl_file(81, trim(localdir)//'/'//'fort.81',
     &          nstac_g, nstac, header81)
         endif
      endif


      !-----------------------------------------------------------------
      ! Initialize Meteorological Station Spool Output
      !-----------------------------------------------------------------
      ! Multiple files controlled by same output controls, for each
      !   Open meteorological station output file
      !   Write out header INFOrmation including:
      !      ntrspm (no. of data pts. at each meteorological station),
      !      nstam, dt*nspoolm, nspoolm, irtype
      ! Variables include:
      !   wind speed (fort.71)
      !   pressure (fort.72)
      !   ice (fort.91)
      !-----------------------------------------------------------------

      nscoum=0

      if(abs(noutm).eq.1) then
         if (.not. write_local_files) then
            call open_gbl_file(71, trim(globaldir)//'/'//'fort.71',
     &          nstam_g, nstam, header71)
            call open_gbl_file(72, trim(globaldir)//'/'//'fort.72',
     &         nstam_g, nstam, header72)
            if(ncice.ne.0) then
               call open_gbl_file(91, trim(globaldir)//'/'//'fort.91',
     &             nstam_g, nstam, header91)
            endif
         else
            call open_gbl_file(71, trim(localdir)//'/'//'fort.71',
     &          nstam_g, nstam, header71)
            call open_gbl_file(72, trim(localdir)//'/'//'fort.72',
     &         nstam_g, nstam, header72)
            if(ncice.ne.0) then
               call open_gbl_file(91, trim(localdir)//'/'//'fort.91',
     &             nstam_g, nstam, header91)
            endif
         endif
      endif


      ! --------------------------------------------------------------
      ! Initialize Global Elevation Output
      ! --------------------------------------------------------------
      !   Many files are produced on same schedule as global elevation
      !   For each:
      !      Open global output file
      !      Write out header INFOrmation including:
      !           ndsetse (no. of global data sets to write)
      !           np, dt*nspoolg*, nspoolg*, irtype
      !   Initialize a counter for spooling output
      !   Variables controlled by the global elevation include:
      !      elevation:          fort.63
      !      variable bathy:     fort.76
      !      tau0:               fort.90
      !      sponge layer:       fort.92
      !      nodecode (wetdry):  nodecode.63 (lun=101)
      !      noff (wetdry):      noff.100 (lun=102)
      !      dynamicWaterLevelCorrection.63  (lun=108)
      ! --------------------------------------------------------------

      nscouge=0   ! initialize elevation spool counter
      nscougb=0   ! tcm: variable bathy

      if((abs(noutge).eq.1).or.(abs(noutge).eq.4)) then

         ! Write Global Files
         if (.not. write_local_files) then
            call open_gbl_file(63, trim(globaldir)//'/'//'fort.63',
     &       np_g, np, header63)

            ! tcm v50.66.01 additions for time varying bathymetry
            if (nddt.ne.0) then
               call open_gbl_file(76, trim(globaldir)//'/'//'fort.76',
     &         np_g, np, header63)
            endif

            ! jgf47.06 tau0
            !casey 090302: corrected the unit number from 10 to 90.
            if (outputTau0) then
               call writedomainheader(90,
     &              trim(globaldir)//'/'//'fort.90',
     &              np_g, np, 'tau0      ')
            endif

             !kmd48.33 sponge layer
             if (outputSponge) then
                 call writedomainheader(92,
     &             trim(globaldir)//'/'//'fort.92',
     &             np_g, np, 'sponge     ')
             end if

             ! jgf: dynamicWaterLevelCorrection
             if (usingDynamicWaterLevelCorrection) then
                 call writedomainheader(108,trim(globaldir)//'/'//
     &            'dynamicwaterlevelcorrection.63',np_g, np,
     &            'dynamicwaterlevelcorrection')
             end if

         else ! Write Local Files
            call open_gbl_file(63, trim(localdir)//'/'//'fort.63',
     &         np_g, np, header63)

            !tcm variable bathy
            if (nddt .ne. 0) then
               call open_gbl_file(76, trim(localdir)//'/'//'fort.76',
     &         np_g, np, header63)
            endif

            if (outputTau0) then
               call writedomainheader(90,
     &           trim(globaldir)//'/'//'fort.90',
     &           np_g, np, 'tau0      ')
            endif

            if (outputSponge) then
               call writedomainheader(92,
     &           trim(localdir)//'/'//'fort.92',
     &           np_g, np, 'sponge     ')
            end if

             ! jgf: dynamicWaterLevelCorrection
             if (usingDynamicWaterLevelCorrection) then
                 call writedomainheader(108,trim(localdir)//'/'//
     &            'dynamicWaterLevelCorrection.63',np_g, np,
     &            'dynamicWaterLevelCorrection')
             end if

         endif
      endif


      ! jgf52.04.03: create files for nodecode and noff output
      ! (wet/dry state of nodes and wet/dry state of elements)
      if (abs(noutge).ne.0) then
         if (outputNodeCode) then
            if (write_local_files.eqv..true.) then
               call writeDomainHeader(101,trim(localdir)//'/'//
     &            'nodecode.63',np_g, np, 'nodecode    ')
            else
               call writedomainheader(101,trim(globaldir)//'/'//
     &            'nodecode.63',np_g, np, 'nodecode    ')
            endif
         endif
         if (outputNoff) then
            if (write_local_files.eqv..true.) then
               call writedomainheader(102,trim(localdir)//'/'//
     &            'noff.100',ne_g, ne, 'noff       ')
            else
               call writedomainheader(102,trim(globaldir)//'/'//
     &            'noff.100',ne_g, ne, 'noff       ')
            endif
         endif
      endif


      ! Initialize Global Velocity
      !   fort.64
      !---------------------------
      nscougv=0

      !casey 090828: added the possibility of noutgv = +/- 4.
      if((abs(noutgv).eq.1).or.(abs(noutgv).eq.4)) then
         if (.not. write_local_files) then
            call open_gbl_file(64, trim(globaldir)//'/'//'fort.64',
     &        np_g, np, header64)
         else
            call open_gbl_file(64, trim(localdir)//'/'//'fort.64',
     &        np_g, np, header64)
         endif
      endif

      ! Initialize Global Weir Elevation Change -- cobell
      !    fort.77
      !--------------------------------------------------
      if(use_tvw)then
         nscou_tvw = 0
         if((abs(nout_tvw).eq.1).or.(abs(nout_tvw).eq.4))then
            if (.not. write_local_files) then
               call open_gbl_file(77,trim(globaldir)//'/'//'fort.77',
     &           np_g,np,header77)
            else
               call open_gbl_file(77,trim(localdir)//'/'//'fort.77',
     &           np_g,np,header77)
            end if
         endif
      endif


      ! Initialize Global Meteorology Spool Counter
      !------------------------------------------------------------------
      ! Same control INFOrmation and spool increment for all
      !
      !   global radiation stress file (rads.64)
      !   global wind and pressure output files (fort.73/74)
      !   global ice output file (fort.93)
      !
      ! Write out header INFOrmation including ndsetsw (no. of global
      !      meteorology sets to be spooled), np, dt*nspoolgw,
      !      nspoolgw, irtype
      !------------------------------------------------------------------

      nscougw=0

      ! casey 090302: added the possibility of noutgw = +/- 4.
      if((abs(noutgw).eq.1).or.(abs(noutgw).eq.4)) then
         if (.not. write_local_files) then
            call open_gbl_file(73, trim(globaldir)//'/'//'fort.73',
     &         np_g, np, header73)
            call open_gbl_file(74, trim(globaldir)//'/'//'fort.74',
     &        np_g, np, header74)
      ! casey 090302: added radiation stress gradients.
      ! tcm v56 2024: added nrs=1 and nrs=5
            if ( (abs(nrs).eq.1).or.(abs(nrs).eq.3).or.(abs(nrs).eq.4).or.
     &       (abs(nrs).eq.5) ) then
               call open_gbl_file(164, trim(globaldir)//'/'//'rads.64',
     &             np_g, np, header74)
            end if
      ! tcm v49.64.01 added for ice
            if(ncice.ne.0) then
               call open_gbl_file(93, trim(globaldir)//'/'//'fort.93',
     &           np_g, np, header93)
            endif
         else
            call open_gbl_file(73, trim(localdir)//'/'//'fort.73',
     &            np_g, np, header73)
            call open_gbl_file(74, trim(localdir)//'/'//'fort.74',
     &        np_g, np, header74)
      ! casey 090302: added radiation stress gradients.
      ! tcm v56 2024: added nrs=1 and nrs=5
            if ( (abs(nrs).eq.1).or.(abs(nrs).eq.3).or.(abs(nrs).eq.4).or.
     &       (abs(nrs).eq.5) ) then
               call open_gbl_file(164, trim(localdir)//'/'//'rads.64',
     &             np_g, np, header74)
            end if
      ! tcm v49.64.01 added for ice
            if(ncice.ne.0) then
               call open_gbl_file(93, trim(localdir)//'/'//'fort.93',
     &           np_g, np, header93)
            endif
         endif
      endif


      ! Initialize Global Concentration Spool Counter
      !   Open global concentration output file (fort.83)
      !   Write out header INFOrmation including ndsetsc (no. of global
      !      concentration data sets to be spooled), np, dt*nspoolgc,
      !      nspoolgc, irtype

      nscougc=0

      if(abs(noutgc).eq.1) then
         open(83,file=trim(inputdir)//'/'//'fort.83')
         write(83,3220) rundes,runid,agrid
         write(83,3645) ndsetsc,np,dtdp*nspoolgc,nspoolgc,1
         close(83)
      endif

 3220 format(1x,a32,2x,a24,2x,a24)
 3645 format(1x,i10,1x,i10,1x,e15.7,1x,i5,1x,i5)



      !Initialize 3D Solution

c...  lines to run the code in 3d vs mode.

      if (c3d) then
         call coldstart_3d()
      endif

c...lines to run the code in 3d dss mode

c     if (c3ddss) then
c       call dssstup(dt,nt)
c     endif
c...
      ! Formatting for screen unit Messages during file reads
 1112 format(/,1x,73('_'))

<<<<<<< HEAD
      RETURN
      END SUBROUTINE COLDSTART


C******************************************************************************
C   Subroutine to initialize the 3D routines for a cold start including       *
C   reading in an initial density field.                                      *
C                                                                             *
C                                                                             *
c******************************************************************************

      SUBROUTINE COLDSTART_3D()
      USE SIZES, ONLY : WRITE_LOCAL_FILES, MNPROC, GLOBALDIR
      USE WRITE_OUTPUT, ONLY : initOutput3D
      USE GLOBAL_3DVS
      USE GLOBAL, ONLY : BCFLAG_TEMP, RES_BC_FLAG, BCFLAG_LNM,
     &                   RBCTIME1, BCSTATIM, RBCTIME2,
     &                   RBCTIMEINC, SBCTIME1,
     &                   SBCTIME2, SBCSTATIM, TTBCTIME1,
     &                   TTBCTIME2, TBCSTATIM, SBCTIMEINC,
     &                   TBCTIMEINC, TTBCTIMEINC, TTBCSTATIM,
     &                   TBCTIME1, TBCTIME2, q_heat1, q_heat2,
     &                   LNM_BC, LNM_BC1, LNM_BC2, RIVBCTIMINC,
     &                   RIVBCTIME1, RIVBCTIME2,
     &                   RIVBCSTATIM, NP_G
      USE mod_logging, ONLY : t_log_scope, init_log_scope, allMessage, DEBUG
      USE MESH, ONLY : NP, AGRID
      USE BOUNDARIES, ONLY : NETA, NOPE, NBD,
     &  totalbcrivernodes, bndbcriver
      USE mod_terminate, ONLY : terminate, ADCIRC_EXIT_FAILURE
      IMPLICIT NONE

      INTEGER :: IRType
      INTEGER :: N                       !loop counter
      INTEGER :: NH, IHN, IVN            !horizontal & vertical loop counters
      INTEGER :: NHNN, NVNN              !horizontal and vertical node numbers
      INTEGER :: NVN
      INTEGER NVP                        ! number of horizontal nodes

ckmd  additional variable for baroclinic boundary conditions
      CHARACTER(80) :: CDUM80
      INTEGER :: NumofBCNodes
      INTEGER :: NOD, J, NumofNodes, K
      INTEGER :: IINDX
      REAL(8), ALLOCATABLE :: TMP(:,:)
C
      LOG_SCOPE_TRACED("coldstart_3D", CSTART_TRACING)
C...
C...  Define format statements used to initialize 3D output files
C...
 499  FORMAT(1X,A32,2X,A24,2X,A24)

C...  IF A BAROCLINIC RUN, READ IN INITIAL DENSITY FIELD

      IF (CBaroclinic.and..not.C3DVS) THEN
         WRITE(16,*) "CBaroclinic is ",CBaroclinic
         WRITE(16,*) "CB3DVS is ",C3DVS
         CALL terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &      message='2DDI baroclinic runs not currently '
     &             //'supported.')
      ENDIF
C
      IF(CBaroclinic .and. C3DVS) THEN
         WRITE(16,424)
 424     FORMAT(/,5X,'INITIAL DENSITY FIELD READ IN FROM UNIT 11',/)
         OPEN(11,FILE=TRIM(INPUTDIR)//'/'//'fort.11')
         READ(11,*)             !skip over header line
         READ(11,*)             !skip over header line
         READ(11,*) NVN, NVP
         IF(NVN.NE.NFEN) THEN
            IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
               WRITE(ScreenUnit,351) NVN,NFEN
            ENDIF
            WRITE(16,351) NVN,NFEN
 351        FORMAT(/,2X,'***** INVALID INPUT IN THE DENSITY INITIAL ',
     &           'CONDITION FILE (UNIT 11) *****',
     &           /,2X,'***** NVN = ',I4,' MUST MATCH NFEN = ',I4,
     &           ' *****')
            CALL terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &         message='Invalid input in density initial '
     &                //'condition file (unit 11). NVN must match NFEN.')
         ENDIF
         DO IHN=1,NP
            DO IVN=1,NFEN
               IF(ABS(IDen).EQ.1)
     &              READ(11,*) NHNN,NVNN,SigT(NHNN,NVNN)
               IF(ABS(IDen).EQ.2)
     &              READ(11,*) NHNN,NVNN,Sal(NHNN,NVNN)
               IF(ABS(IDen).EQ.3)
     &              READ(11,*) NHNN,NVNN,Temp(NHNN,NVNN)
               IF(ABS(IDen).EQ.4)
     &              READ(11,*) NHNN,NVNN,Temp(NHNN,NVNN),Sal(NHNN,NVNN)
            ENDDO
         ENDDO
         CLOSE(11)
      ENDIF

      IF (CBaroclinic) THEN
!     Kendra45.12: Read in the temperature boundary condition for the
!     temperature field. Note, we will add the salinity condition
!     later. Currently, we are adopting the sign convention that the
!     flux into the domain is consider positive and the flux out of the
!     domain is negative. However, remember the constitutive law changes
!     the direction of the flux.
         IF ((IDEN.EQ.3).OR.(IDEN.EQ.-3)) THEN
            IF (NTF.EQ.0) THEN
               DO IHN=1,NP
                  qsurfkp1(IHN)=0.d0
                  qsurf(IHN)=0.d0
               END DO
            ELSE IF (NTF.EQ.1) THEN
               OPEN(111,FILE=TRIM(INPUTDIR)//'/'//'fort.111')
               DO IHN=1,NP
                  READ(111,*) NHNN,qsurfkp1(IHN),qsurf(IHN)
               END DO
               CLOSE(111)
            END IF
         END IF
      END IF

ckmd  Read in boundary condition information from files if
C...  we are using the transport portion of the code (IDEN.GE.1)
C...  and baroclinic.

      IF ((CBAROCLINIC).AND.(IDEN.GE.1)) THEN
         IF ((ABS(RES_BC_FLAG).GE.1).AND.(NOPE.GT.0)) THEN
           IF (BCFLAG_LNM.EQ.1) THEN
             OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
             RBCTIME1=BCSTATIM*86400.d0
             RBCTIME2=RBCTIME1+RBCTIMEINC
             READ(35,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC1(NumofBCNodes)
             END DO
             READ(35,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC2(NumofBCNodes)
             END DO
           ELSE IF (BCFLAG_LNM.EQ.2) THEN
             ! to be added later but currently information
             ! is taken in two ways from the HYCOM results
           ELSE IF (BCFLAG_LNM.EQ.3) THEN
             DO NumofBCNodes=1,NETA
                IINDX=NBD(NumofBCNodes)
                LNM_BC1(NumofBCNodes)=ETA2(IINDX)
                LNM_BC2(NumofBCNodes)=LNM_BC1(NumofBCNodes)
             END DO
           END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.2).OR.(ABS(RES_BC_FLAG).EQ.4)) THEN
           IF (NOPE.GT.0) THEN
             OPEN(36,FILE=TRIM(INPUTDIR)//'/'//'fort.36')
             SBCTIME1=SBCSTATIM*86400.D0
             SBCTIME2=SBCTIME1+SBCTIMEINC
             READ(36,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(36,*) NOD,(RESSAL1(NumofBCNodes,J),J=1,NFEN)
             END DO
             READ(36,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
             END DO
           END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &            (ABS(RES_BC_FLAG).EQ.4)) THEN
           IF (NOPE.GT.0) THEN
             OPEN(37,FILE=TRIM(INPUTDIR)//'/'//'fort.37')
             TBCTIME1=TBCSTATIM*86400.D0
             TBCTIME2=TBCTIME1+TBCTIMEINC
             READ(37,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(37,*) NOD,(RESTEMP1(NumofBCNodes,J),J=1,NFEN)
             END DO
             READ(37,'(A)') CDUM80
             DO NumofBCNodes=1,NETA
                READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
             END DO
           END IF
           IF (BCFLAG_TEMP.EQ.1) THEN ! read in file with multiple values but no calcs
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              DO NumofNodes=1,NP
                 READ(38,*) NOD, q_heat1(Numofnodes)
              END DO
              DO NumofNodes=1,NP
                 READ(38,*) NOD, q_heat2(NumofNodes)
              END DO
           ELSE IF (BCFLAG_TEMP.EQ.2) THEN ! read in from file with 6 components
              ALLOCATE (TMP(NP,6))
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=-TMP(NumofNodes,1)-
     &              TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &              TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &              TMP(NumofNodes,6)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &              TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &              TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &              TMP(NumofNodes,6)
              END DO
              DEALLOCATE(TMP)
           ELSE IF (BCFLAG_TEMP.EQ.3) THEN ! read in from file with 4 components
              ALLOCATE (TMP(NP,4))
              OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
              TTBCTIME1=TTBCSTATIM*86400.D0
              TTBCTIME2=TTBCTIME1+TTBCTIMEINC
              READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat1(NumofNodes)=TMP(NumofNodes,4)+
     &              TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &              TMP(NumofNodes,2)
              END DO
              READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
              DO NumofNodes=1,NP
                 q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &              TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &              TMP(NumofNodes,2)
              END DO
              DEALLOCATE(TMP)
           END IF ! end if for heat flux file
         END IF ! end if for temperature boundary and heat flux
      ELSE IF ((CBAROCLINIC).AND.(IDEN.LE.-1)) THEN
         IF ((RES_BC_FLAG.LT.0).AND.(NOPE.GT.0)) THEN
          IF (BCFLAG_LNM.EQ.1) THEN
            OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
            READ(35,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC(NumofBCNodes)
            END DO
          ELSE IF (BCFLAG_LNM.EQ.2) THEN
=======

      return
!-------------------------------------------------------------------------
      end subroutine coldstart
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!> @brief: Handles all subroutine calls for a 3D coldstart initialization.
!>
!> @details: Intializes all 3D routines for a coldstart, including
!>           reading the initial density fields and opening output.
!-------------------------------------------------------------------------
      subroutine coldstart_3D()
!-------------------------------------------------------------------------
      use sizes, only : write_local_files, mnproc, globaldir
      use write_output, only : initOutput3D
      use global_3dvs
      use global, only : bcflag_temp, res_bc_flag, bcflag_lnm,
     &                   rbctime1, bcstatim, rbctime2,
     &                   rbctimeinc, sbctime1,
     &                   sbctime2, sbcstatim, ttbctime1,
     &                   ttbctime2, tbcstatim, sbctimeinc,
     &                   tbctimeinc, ttbctimeinc, ttbcstatim,
     &                   tbctime1, tbctime2, q_heat1, q_heat2,
     &                   lnm_bc, lnm_bc1, lnm_bc2, rivbctiminc,
     &                   rivbctime1, rivbctime2,
     &                   rivbcstatim, np_g
      USE mod_logging, ONLY : t_log_scope, init_log_scope, allMessage, DEBUG
      USE mod_terminate, ONLY : terminate, ADCIRC_EXIT_FAILURE
      use mesh, only : np, agrid
      use boundaries, only : neta, nope, nbd,
     &  totalbcrivernodes, bndbcriver

      implicit none

      integer :: irtype
      integer :: n                !loop counter
      integer :: nh, ihn, ivn     !horizontal & vertical loop counters
      integer :: nhnn, nvnn       !horizontal and vertical node numbers
      integer :: nvn
      integer nvp                 ! number of horizontal nodes

      !kmd: additional variable for baroclinic boundary conditions
      character(80) :: cdum80
      integer :: numOfBCNodes
      integer :: nod, j, numofnodes, k
      integer :: iindx
      real(8), allocatable :: tmp(:,:)

      LOG_SCOPE_TRACED("coldstart_3D", CSTART_TRACING)

      ! Define format statements used to initialize 3D output files
 499  format(1x,a32,2x,a24,2x,a24)

      ! If a baroclinic run, read in initial density field
      if (cbaroclinic.and..not.c3dvs) then
         write(16,*) "CBaroclinic is ",CBaroclinic
         write(16,*) "CB3DVS is ",c3dvs
         call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &      message='2DDI baroclinic runs not currently supported.')
      endif

      if(cbaroclinic .and. c3dvs) then
         write(16,424)
 424     format(/,5x,'Initial density field read in from unit 11',/)
         open(11,file=trim(inputdir)//'/'//'fort.11')
         read(11,*)             !skip over header line
         read(11,*)             !skip over header line
         read(11,*) nvn, nvp
         if(nvn.ne.nfen) then
            if((nscreen.ne.0).and.(myproc.eq.0)) then
               write(screenunit,351) nvn,nfen
            endif
            write(16,351) nvn,nfen
 351        format(/,2x,'***** Invalid input in the density initial ',
     &           'condition file (unit 11) *****',
     &           /,2x,'***** nvn = ',i4,' must match nfen = ',i4,
     &           ' *****',
     &           /,10x,'****** RUN TERMINATED ******')
            call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &         message='Invalid input in density initial '
     &                //'condition file (unit 11). NVN must match NFEN.')
         endif
         do ihn=1,np
            do ivn=1,nfen
               if(abs(iden).eq.1)
     &              read(11,*) nhnn,nvnn,sigt(nhnn,nvnn)
               if(abs(iden).eq.2)
     &              read(11,*) nhnn,nvnn,sal(nhnn,nvnn)
               if(abs(iden).eq.3)
     &              read(11,*) nhnn,nvnn,temp(nhnn,nvnn)
               if(abs(iden).eq.4)
     &              read(11,*) nhnn,nvnn,temp(nhnn,nvnn),sal(nhnn,nvnn)
            enddo
         enddo
         close(11)
      endif

      if (CBaroclinic) then
      ! kmd 45.12: Read in the temperature boundary condition for the
      ! temperature field. Note, we will add the salinity condition
      ! later. Currently, we are adopting the sign convention that the
      ! flux into the domain is consider positive and the flux out of the
      ! domain is negative. However, remember the constitutive law changes
      ! the direction of the flux.
         if ((iden.eq.3).or.(iden.eq.-3)) then
            if (ntf.eq.0) then
               do ihn=1,np
                  qsurfkp1(ihn)=0.d0
                  qsurf(ihn)=0.d0
               end do
            else if (ntf.eq.1) then
               open(111,file=trim(inputdir)//'/'//'fort.111')
               do ihn=1,np
                  read(111,*) nhnn,qsurfkp1(ihn),qsurf(ihn)
               end do
               close(111)
            end if
         end if
      end if

      ! kmd: Read in boundary condition INFOrmation from files if
      ! we are using the transport portion of the code (iden.ge.1)
      ! and baroclinic.
      if ((CBaroclinic).and.(iden.ge.1)) then
         if ((abs(res_bc_flag).ge.1).and.(nope.gt.0)) then
           if (bcflag_lnm.eq.1) then
             open(35,file=trim(inputdir)//'/'//'fort.35')
             rbctime1=bcstatim*86400.d0
             rbctime2=rbctime1+rbctimeinc
             read(35,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(35,*) nod,lnm_bc1(numOfBCNodes)
             end do
             read(35,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(35,*) nod,lnm_bc2(numOfBCNodes)
             end do
           else if (bcflag_lnm.eq.2) then
             ! to be added later but currently INFOrmation
             ! is taken in two ways from the hycom results
           else if (bcflag_lnm.eq.3) then
             do numOfBCNodes=1,neta
                iindx=nbd(numOfBCNodes)
                lnm_bc1(numOfBCNodes)=eta2(iindx)
                lnm_bc2(numOfBCNodes)=lnm_bc1(numOfBCNodes)
             end do
           end if
         end if
         if ((abs(res_bc_flag).eq.2).or.(abs(res_bc_flag).eq.4)) then
           if (nope.gt.0) then
             open(36,file=trim(inputdir)//'/'//'fort.36')
             sbctime1=sbcstatim*86400.d0
             sbctime2=sbctime1+sbctimeinc
             read(36,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(36,*) nod,(ressal1(numOfBCNodes,j),j=1,nfen)
             end do
             read(36,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(36,*) nod,(ressal2(numOfBCNodes,j),j=1,nfen)
             end do
           end if
         end if
         if ((abs(res_bc_flag).eq.3).or.
     &            (abs(res_bc_flag).eq.4)) then
           if (nope.gt.0) then
             open(37,file=trim(inputdir)//'/'//'fort.37')
             tbctime1=tbcstatim*86400.d0
             tbctime2=tbctime1+tbctimeinc
             read(37,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(37,*) nod,(restemp1(numOfBCNodes,j),j=1,nfen)
             end do
             read(37,'(a)') cdum80
             do numOfBCNodes=1,neta
                read(37,*) nod,(restemp2(numOfBCNodes,j),j=1,nfen)
             end do
           end if
           if (bcflag_temp.eq.1) then ! read in file with multiple values but no calcs
              open(38,file=trim(inputdir)//'/'//'fort.38')
              ttbctime1=ttbcstatim*86400.d0
              ttbctime2=ttbctime1+ttbctimeinc
              do numofnodes=1,np
                 read(38,*) nod, q_heat1(numofnodes)
              end do
              do numofnodes=1,np
                 read(38,*) nod, q_heat2(numofnodes)
              end do
           else if (bcflag_temp.eq.2) then ! read in from file with 6 components
              allocate (tmp(np,6))
              open(38,file=trim(inputdir)//'/'//'fort.38')
              ttbctime1=ttbcstatim*86400.d0
              ttbctime2=ttbctime1+ttbctimeinc
              read(38,*) (nod,(tmp(k,j),j=1,6),k=1,np)
              do numofnodes=1,np
                 q_heat1(numofnodes)=-tmp(numofnodes,1)-
     &              tmp(numofnodes,2)+tmp(numofnodes,3)-
     &              tmp(numofnodes,5)+tmp(numofnodes,4)-
     &              tmp(numofnodes,6)
              end do
              read(38,*) (nod,(tmp(k,j),j=1,6),k=1,np)
              do numofnodes=1,np
                 q_heat2(numofnodes)=-tmp(numofnodes,1)-
     &              tmp(numofnodes,2)+tmp(numofnodes,3)-
     &              tmp(numofnodes,5)+tmp(numofnodes,4)-
     &              tmp(numofnodes,6)
              end do
              deallocate(tmp)
           else if (bcflag_temp.eq.3) then ! read in from file with 4 components
              allocate (tmp(np,4))
              open(38,file=trim(inputdir)//'/'//'fort.38')
              ttbctime1=ttbcstatim*86400.d0
              ttbctime2=ttbctime1+ttbctimeinc
              read(38,*) (nod,(tmp(k,j),j=1,4),k=1,np)
              do numofnodes=1,np
                 q_heat1(numofnodes)=tmp(numofnodes,4)+
     &              tmp(numofnodes,3)-tmp(numofnodes,1)+
     &              tmp(numofnodes,2)
              end do
              read(38,*) (nod,(tmp(k,j),j=1,4),k=1,np)
              do numofnodes=1,np
                 q_heat2(numofnodes)=tmp(numofnodes,4)+
     &              tmp(numofnodes,3)-tmp(numofnodes,1)+
     &              tmp(numofnodes,2)
              end do
              deallocate(tmp)
           end if ! end if for heat flux file
         end if ! end if for temperature boundary and heat flux
      else if ((cbaroclinic).and.(iden.le.-1)) then
         if ((res_bc_flag.lt.0).and.(nope.gt.0)) then
          if (bcflag_lnm.eq.1) then
            open(35,file=trim(inputdir)//'/'//'fort.35')
            read(35,'(a)') cdum80
            do numOfBCNodes=1,neta
                read(35,*) nod,lnm_bc(numOfBCNodes)
            end do
          else if (bcflag_lnm.eq.2) then
>>>>>>> Subject: Add additional source terms and update style
            ! add later
          else if (bcflag_lnm.eq.3) then
            do numOfBCNodes=1,neta
               iindx=nbd(numOfBCNodes)
               lnm_bc(numOfBCNodes)=eta2(iindx)
            end do
          end if
         end if
      end if

       ! kmd - river INFOrmation for boundary conditions in baroclinic simulation
      if (bndbcriver) then
        open(39,file=trim(inputdir)//'/'//'fort.39')
        read(39,*) rivbctiminc, rivbcstatim
        rivbctime1=rivbcstatim*86400.d0
        rivbctime2=rivbctime1+rivbctiminc
        do j=1,totalbcrivernodes
         if (iden.eq.2) then
           read(39,*) (bcrivsaln1(j,k),k=1,nfen)
         else if (iden.eq.3) then
           read(39,*) (bcrivtempn1(j,k),k=1,nfen)
         else if (iden.eq.4) then
           read(39,*) (bcrivsaln1(j,k), bcrivtempn1(j,k),k=1,nfen)
         end if
        end do
        do j=1,totalbcrivernodes
         if (iden.eq.2) then
           read(39,*) (bcrivsaln2(j,k),k=1,nfen)
         else if (iden.eq.3) then
           read(39,*) (bcrivtempn2(j,k),k=1,nfen)
         else if (iden.eq.4) then
           read(39,*) (bcrivsaln2(j,k), bcrivtempn2(j,k),k=1,nfen)
         end if
        end do
      end if


      ! Zero out stuff passed from 3D solution to external mode
      do nh=1,np
         duu(nh)=0.d0
         duv(nh)=0.d0
         dvv(nh)=0.d0
         uu(nh)=0.d0
         vv(nh)=0.d0
         dafluxx(nh)=0.d0
         dafluxy(nh)=0.d0
         bsx(nh)=0.d0
         bsy(nh)=0.d0
      enddo


      ! Initialize 3D Velocity and Turbulence Solution
      ! kmd45.12 need to initialize the new qkp1 and wzkp1 arrays
      do nh=1,np
         do n=1,nfen
            q(nh,n)=(0.d0,0.d0)
            qkp1(nh,n)=(0.d0,0.d0)
            ! Same problem as with l
            q20(nh,n)=-9999d0

            ! rjw bug fix: Cannot initialize l to zero, must initialize
            ! to lmin or non-zero value, so when elements wet, the first
            ! turb length scale is nonzero
            l(nh,n)=-9999d0
            wz(nh,n)=0.d0
            wzkp1(nh,n)=0.d0
         enddo
      enddo
      bpg(:,:)=(0.d0,0.d0) ! jgf50.44: was not initialized anywhere else

<<<<<<< HEAD
C...  INITIALIZE 3D OUTPUT FILES
C
C.... Initialize the 3D density station output file (Unit 41)
      IF((IDen.EQ.1).OR.(IDen.EQ.0)) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
      CALL initOutput3D(41, I3DSD, NDSET3DSD, NSTA3DD, NSTA3DD_G,
     &                  NSPO3DSD, IRType, I3DSDRec)
C
C.... Initialize the 3D velocity station output file (Unit 42)
      IRType=3
      CALL initOutput3D(42, I3DSV, NDSET3DSV, NSTA3DV, NSTA3DV_G,
     &                  NSPO3DSV, IRType, I3DSVRec)
C
C.... Initialize the 3D turbulence station output file (Unit 43)
      IRType=3
      CALL initOutput3D(43, I3DST, NDSET3DST, NSTA3DT, NSTA3DT_G,
     &                  NSPO3DST, IRType, I3DSTRec)
C.... Initialize the global 3D density output file (Unit 44)

      IF((IDen.EQ.1).OR.(IDen.EQ.0)) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
C     jgf46.27 Replaced IRType with IDen in header.
!kmd48.33bc - add in infomration for file containing the
!            top temperature boundary condition information
!      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_G,
!     &                  NSPO3DGD, IRType, I3DGDRec)
! arash
      CALL initOutput3D(44, I3DGD, NDSET3DGD, NP, NP_G,
     &                  NSPO3DGD, IDen, I3DGDRec)

      IF (ABS(I3DGD).EQ.1) THEN
         IF ((MNPROC.EQ.1.).OR.(WRITE_LOCAL_FILES.eqv..true.)) THEN
            IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
               IF (BCFLAG_TEMP.NE.0) THEN
                  OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47')
                  WRITE(47,499) RUNDES,RUNID,AGRID
                  WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
                  CLOSE(47)
               END IF
            END IF
         ENDIF
         IF ((MNPROC.GT.1.).AND.(WRITE_LOCAL_FILES.eqv..false.).and.
     &       (myProc.eq.0)) THEN
            IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
               IF (BCFLAG_TEMP.NE.0) THEN
                  OPEN(47,FILE=TRIM(GLOBALDIR)//'/'//'fort.47')
                  WRITE(47,499) RUNDES,RUNID,AGRID
                  WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
                  CLOSE(47)
               END IF
            END IF
         ENDIF
      ENDIF
C..RJW bug fix in 498 (kendra found this)
 498  FORMAT(1X,I10,1X,I10,1X,E15.7,I10,1X,I10,1X,I3)

C
c.... Initialize the global 3D velocity output file (Unit 45)
      CALL initOutput3D(45, I3DGV, NDSET3DGV, NP, NP_G,
     &                  NSPO3DGV, IRType, I3DGVRec)
C
c.... Initialize the global 3D turbulence output file (Unit 46)
      CALL initOutput3D(46, I3DGT, NDSET3DGT, NP, NP_G,
     &                  NSPO3DGT, IRType, I3DGTRec)
C
C.... Set up a few final odds and ends for a 3D run
      CALL VSSTUP ()

      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE COLDSTART_3D
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C                      S U B R O U T I N E
C           S E T   C O L D   W E T   D R Y   S T A T E
C       A N D   W A T E R   S U R F A C E   E L E V A T I O N
C----------------------------------------------------------------------
C     jgf52.08.01: Sets value of nnodecode at cold start. The values
C     of nnodecode are copied to the nodecode array in totalAreaCalc().
C     This is called from subroutine hotstart() for use in producing
C     detailed inundation output only. When called from hotstart(),
C     the hotstart wet/dry state and water surface elevation ends up
C     overwriting the values computed here.
C----------------------------------------------------------------------
=======
      !---------------------------
      ! Initialize 3D Output Files
      !---------------------------

      ! Initialize the 3D Density Station Output (unit 41)
      if((iden.eq.1).or.(iden.eq.0)) irtype=1
      if((iden.eq.2).or.(iden.eq.3)) irtype=2
      if(iden.eq.4) irtype=3
      call initOutput3D(41, i3dsd, ndset3dsd, nsta3dd, nsta3dd_g,
     &                  nspo3dsd, irtype, i3dsdrec)

      ! Initialize the 3D Velocity Station Output (unit 42)
      irtype=3
      call initOutput3D(42, i3dsv, ndset3dsv, nsta3dv, nsta3dv_g,
     &                  nspo3dsv, irtype, i3dsvrec)

      ! Initialize the 3D Turbulence Station Output (unit 43)
      irtype=3
      call initOutput3D(43, i3dst, ndset3dst, nsta3dt, nsta3dt_g,
     &                  nspo3dst, irtype, i3dstrec)

      ! Initialize the Global 3D Density Output (unit 44)
      if((iden.eq.1).or.(iden.eq.0)) irtype=1
      if((iden.eq.2).or.(iden.eq.3)) irtype=2
      if(iden.eq.4) irtype=3
      ! jgf46.27 replaced irtype with iden in header.
      ! kmd48.33bc - add information for file containing the
      !    top temperature boundary condition information
      ! arash
      call initOutput3D(44, i3dgd, ndset3dgd, np, np_g,
     &                  nspo3dgd, iden, i3dgdrec)

      if (abs(i3dgd).eq.1) then
         if ((mnproc.eq.1.).or.(write_local_files.eqv..true.)) then
            if((iden.eq.3).or.(iden.eq.4)) then
               if (bcflag_temp.ne.0) then
                  open(47,file=trim(localdir)//'/'//'fort.47')
                  write(47,499) rundes,runid,agrid
                  write(47,498) ndset3dgd,np,dtdp*nspo3dgd,nspo3dgd,1
                  close(47)
               end if
            end if
         endif
         if ((mnproc.gt.1.).and.(write_local_files.eqv..false.).and.
     &       (myproc.eq.0)) then
            if((iden.eq.3).or.(iden.eq.4)) then
               if (bcflag_temp.ne.0) then
                  open(47,file=trim(globaldir)//'/'//'fort.47')
                  write(47,499) rundes,runid,agrid
                  write(47,498) ndset3dgd,np,dtdp*nspo3dgd,nspo3dgd,1
                  close(47)
               end if
            end if
         endif
      endif
      !..rjw bug fix in 498 (kendra found this)
 498  format(1x,i10,1x,i10,1x,e15.7,i10,1x,i10,1x,i3)


      ! Initialize the Global 3D Velocity Output (unit 45)
      call initOutput3D(45, i3dgv, ndset3dgv, np, np_g,
     &                  nspo3dgv, irtype, i3dgvrec)

      ! Initialize the Global 3D Turbulence Output (unit 46)
      call initOutput3D(46, i3dgt, ndset3dgt, np, np_g,
     &                  nspo3dgt, irtype, i3dgtrec)

      ! Set up a few final odds and ends for a 3D run
      call vsstup ()

      return
!----------------------------------------------------------------------
      end subroutine coldstart_3d
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> @brief Sets initial wet/dry state and WSE values for a coldstart.
!>
!> @author jgf52.08.01
!>
!> @details Sets value of nnodecode at cold start by checking for
!>          landlocked nodes: must be connected to at least one
!>          functioning  element. The values of nnodecode are then
!>          copied to the nodecode array at the end of the subroutine.
!>          Note: this occurred in the call to totalAreaCalc from
!>          adcirc_init previously, but was moved here (cms 11/20/25)
!>          to allow for reuse of the totalAreaCalc subroutine.
!>
!>          Also computes the initial water surface elevation after
!>          checking for geoid offsets and River_et_WSE nodal
!>          attributes.
!>
!>          This subroutine is also called from hotstart() for use in
!>          producing initial conditions for detailed inundation output.
!>          When called from hotstart(), the hotstart wet/dry state and
!>          water surface elevation overwrites the values computed here.
!-------------------------------------------------------------------------
>>>>>>> Subject: Add additional source terms and update style
      subroutine setColdWetDryStateAndWaterSurfaceElevation()
!-------------------------------------------------------------------------
      use sizes, only : inputdir
<<<<<<< HEAD
      use mod_terminate, only : terminate, ADCIRC_EXIT_FAILURE
      use mesh, only : dp, np, mju, ne, nm
      use boundaries, only : nvel, lbcodei, nbv
      use global, only : eta2, nnodecode, noff, h0,
     &   eta1, nolifa, EtaDisc, River_above_MSL
      use mod_logging, only : logMessage, allMessage, t_log_scope, init_log_scope,
     &   DEBUG, WARNING, INFO, ERROR
      use nodalattributes, only : loadStartDry, STARTDRY, GeoidOffset,
     &   LoadGeoidOffset, OutputTau0, LoadRiver_et_WSE, River_et_WSE,
     &   BFCdlLimit, fric, fgamma, ftheta, manningsn, hbreak,
     &   loadManningsN, iflinbf, ifnlbf, ifhybf
!JLW: adding subgrid
      use subgrid, only: level0, readsubgridlookup
     &                   , getvertlookup
=======
      use mod_terminate, only : terminate, ADCIRC_EXIT_FAILURE
      use mesh, only : dp, np, mju, ne, nm
      use boundaries, only : nvel, lbcodei, nbv
      use global, only : nodecode, nnodecode, noff, h0, nolifa,
     &   eta1, eta2, EtaDisc, River_above_MSL
      use mod_logging, only : logMessage, allMessage, t_log_scope,
     &   init_log_scope, DEBUG, WARNING, INFO, ERROR
      use nodalattributes, only : loadstartdry, startdry, GeoidOffset,
     &   LoadGeoidOffset, OutputTau0, LoadRiver_et_WSE, River_et_WSE,
     &   BFCdlLimit, fric, fgamma, ftheta, manningsn, hbreak,
     &   loadManningsN, iflinbf, ifnlbf, ifhybf
      ! jlw: adding subgrid
      use subgrid, only: level0, readSubgridLookup, getVertLookup
>>>>>>> Subject: Add additional source terms and update style
     &                   , initVar
     &                   , gridDepthVertETA2
     &                   , subgridVertList
     &                   , wetFracVertTab
     &                   , numPhi
     &                   , wetDepthVertETA2
     &                   , gridDepthVertTab
     &                   , phi_time_derivative, dphidt
      implicit none
      integer :: nm1, nm2, nm3  ! node numbers around the element
      integer :: nc1, nc2, nc3  ! nodecodes around the element
      integer :: ncele          ! wet dry state of the element
      integer :: i
<<<<<<< HEAD
      character(1024) :: scratchMessage
      real(8) :: H1
      real(8) :: HABSMIN
      integer :: numGreater
=======
      character(1024) :: scratchMessage
      real(8) :: h1
      real(8) :: habsmin
      integer :: numGreater
>>>>>>> Subject: Add additional source terms and update style
      integer :: j
      ! jgf52.08.24: investigate reasons why a node was dried by the
      ! setcoldwetdrystate() subroutine.
      character(len=1000), allocatable :: dryReason(:)
      logical, allocatable :: dryReasonInitialized(:)
      logical :: found ! .true. if find obsolete init river elevation file

 
<<<<<<< HEAD
      LOG_SCOPE_TRACED("setColdWDSAndWSE", CSTART_TRACING)
=======
      LOG_SCOPE_TRACED("setColdWDSAndWSE", CSTART_TRACING)
>>>>>>> Subject: Add additional source terms and update style

      habsmin = h0*0.8d0

      ! tcm v51.27 20140502 changed from using a fort.88 file
      ! to a nodal attribute
      if (LoadRiver_et_WSE.eqv..true.) then
         river_above_msl=.true.
         call allMessage(INFO,
     &      'initial river elevations were specified as a nodal '
     &       //'attribute and will be used.')
      endif

      ! tcm v51.27 check to see if a fort.88 file is supplied
      ! if so issue a WARNING.
      inquire(file=trim(inputdir)//'/'//'fort.88',exist=found)
      if (found.eqv..true.) then
         call allMessage(WARNING,
     &    'a fort.88 file was found but is no longer supported. '
     &   //'initial river elevations must be set as a nodal attribute.')
      endif

      ! jgf46.01 added the ability to include steric effects.
      if (loadgeoidoffset.eqv..true.) then
         eta2(:)=geoidoffset(:)
      endif

      ! tcm v51.27 20140502 Made the geoidoffset and initial river
      ! elevation an additive quantity instead of the max of the two.
      !
      ! jgf52.08.01: Simplified to just add the initial river elevation
      ! to eta2 wherever the initial river elevation is greater than zero.
      if (LoadRiver_et_WSE.eqv..true.) then
         where (river_et_wse.gt.0.d0)
            eta2 = eta2 + river_et_wse
         end where
         ! tcmv51.27 20140502 added to remove this nodal attribute
         ! as it is only used at cold start
         deallocate (river_et_wse)
      endif

      ! Set the Radiation BC Information
      do i=1,nvel
         if(lbcodei(i).eq.52) then
            etadisc(nbv(i))=eta2(nbv(i))
         endif
      end do

      ! jlw: Adding subgrid lookup table read in
      if(level0.eqv..true.)then
#ifdef adcnetcdf
         ! cms 11/20/25: These are already initialized in call
         ! to initializeWettingAndDrying() at beginning of subroutine.
         ! Set all nodes and elements to wet by default.
         !nnodecode(:)=1
         !noff(:)=1

         ! jlw: set dphidt term
         if(phi_time_derivative)then
            dphidt = 1
         else
            dphidt = 0
         endif
         ! jlw: call subgrid lookup table read in
         call readSubgridLookup()
         if(nolifa.eq.2)then
            ! jlw: initialize subgrid variables
            call initVar()
            ! jlw: call the vertex lookup subroutine
            call getVertLookup()
            ! jgf52.08.23: collect, store, and report the reasons why each
            ! node was dried.
            allocate(dryReason(np))
            allocate(dryReasonInitialized(np))
            dryReasonInitialized = .false.
            do i=1,np
               h1 = gridDepthVertETA2(i)
               if (h1.le.h0) then
                  nnodecode(i) = 0
                  write(dryReason(i),'("node=",i0)') i
                  write(scratchMessage,
     &               '(" dp+eta2<=h0: ",f10.3," + ",f10.3," <= ",f6.3)')
     &               dp(i), eta2(i), h0
                  dryReason(i) = trim(dryReason(i)) // trim(scratchMessage)
                  dryReasonInitialized(i) = .true.
               endif
               if ((loadstartdry.eqv..true.).and.(startdry(i).eq.1)) then
                  nnodecode(i) = 0
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),'("node=",i0)') i
                     dryReasonInitialized(i) = .true.
                  else
                     dryReason(i) = trim(dryReason(i)) // ", "
                  endif
                  dryReason(i) = trim(dryReason(i)) // "startdry=1"
               endif
            end do

            ! Dry any landlocked nodes by checking that they are connected
            ! to at least 1 functioning element.
            mju(:)=0
            do i=1,ne
               nm1=nm(i,1)
               nm2=nm(i,2)
               nm3=nm(i,3)
               nc1=nnodecode(nm1)
               nc2=nnodecode(nm2)
               nc3=nnodecode(nm3)
               ncele=nc1*nc2*nc3*noff(i)
               mju(nm1)=mju(nm1)+ncele
               mju(nm2)=mju(nm2)+ncele
               mju(nm3)=mju(nm3)+ncele
            enddo
            do i=1,np
               if((nnodecode(i).eq.1).and.(mju(i).eq.0)) then
                  nnodecode(i)=0
                  write(scratchMessage,9883) i
                  call allMessage(INFO,scratchMessage)
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),'("node=",i0)') i
                     dryReasonInitialized(i) = .true.
                  else
                     dryReason(i) = trim(dryReason(i)) // ", "
                  endif
                  dryReason(i) = trim(dryReason(i)) // " landlocked"
               endif
            enddo
! 9883       format('node ',i0,' dried (landlocking).')

            !jgf52.08.23: report reasons why a node was initially dried.
            call logMessage(DEBUG,
     &         "nodes were set to an initially dry state "
     &         // "for the following reasons:")
            do i=1,np
               if (nnodecode(i).eq.0) then
                  if (dryReasonInitialized(i).eqv..false.) then
                     write(dryReason(i),
     &            '("node=",i0," this node is dry for no good reason.")') 
                  dryReasonInitialized(i) = .true.
                  endif
                  call logMessage(DEBUG,trim(dryReason(i)))
               endif
            end do
            ! release memory
            deallocate(dryReason)
            deallocate(dryReasonInitialized)
            !
            ! finish initializing water surface elevation (eta2 and eta1)
            do i=1,np
               if(nnodecode(i).eq.0)then
                  if(gridDepthVertETA2(i).lt.habsmin)then
                     gridDepthVertETA2(i) = habsmin
                     wetDepthVertETA2(i) = habsmin
                     if(subgridVertList(i).eq.1)then
                        numgreater = 0
                        do j=1,numPhi
                           if(gridDepthVertETA2(i).lt.
     &                        gridDepthVertTab(i,j))then
                              numgreater = numgreater + 1
                           endif
                        enddo
                        if(numgreater.eq.numPhi)then
                           eta2(i) = wetFracVertTab(i,1)
                        elseif(numgreater.eq.0)then
                           eta2(i) = wetFracVertTab(i,numPhi)
     &                        + (habsmin-gridDepthVertTab(i,numPhi))
                        else
                        !jlw: now back calculate eta2 based on this change
                           eta2(i) = ((gridDepthVertETA2(i)
     &                     - gridDepthVertTab(i,numPhi-numgreater))
     &                     / (gridDepthVertTab(i,numPhi-numgreater+1)
     &                     - gridDepthVertTab(i,numPhi-numgreater))
     &                     * (wetFracVertTab(i,numPhi-numgreater+1)
     &                     - wetFracVertTab(i,numPhi-numgreater))
     &                     + wetFracVertTab(i,numPhi-numgreater))
                        endif
                     else
                           eta2(i)=habsmin-dp(i)
                     endif
                  endif
               endif
            enddo
         endif
#else
<<<<<<< HEAD
         call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &      message='Subgrid calculation requires building '
     &             //'with netcdf enabled')
=======
         call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &      message='Subgrid calculation requires building '
     &             //'with netcdf enabled')
>>>>>>> Subject: Add additional source terms and update style
#endif
      else

      ! cms 11/20/25: These are already initialized in call
      ! to initializeWettingAndDrying() at beginning of subroutine.
      ! Set all nodes and elements to wet by default.
      !nnodecode(:)=1
      !noff(:)=1

      ! if wetting and drying is active, set the nodal wet dry state
      if (nolifa.eq.2) then
         !
         ! jgf52.08.23: collect, store, and report the reasons why each
         ! node was dried.
         allocate(dryReason(np))
         allocate(dryReasonInitialized(np))
         dryReasonInitialized = .false.
         do i=1,np
            ! set dry areas based on eta2
            if (dp(i)+eta2(i).le.h0) then
               nnodecode(i) = 0
               write(dryReason(i),'("node=",i0)') i
               write(scratchMessage,
     &            '(" dp+eta2<=h0: ",f10.3," + ",f10.3," <= ",f6.3)')
     &            dp(i), eta2(i), h0
               dryReason(i) = trim(dryReason(i)) // trim(scratchMessage)
               dryReasonInitialized(i) = .true.
            endif
            if ((loadstartdry.eqv..true.).and.(startdry(i).eq.1)) then
               nnodecode(i) = 0
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),'("node=",i0)') i
                  dryReasonInitialized(i) = .true.
               else
                  dryReason(i) = trim(dryReason(i)) // ", "
               endif
               dryReason(i) = trim(dryReason(i)) // "startdry=1"
            endif
         end do

         ! Dry any landlocked nodes by checking that they are connected
         ! to at least 1 functioning element.
         mju(:)=0
         do i=1,ne
            nm1=nm(i,1)
            nm2=nm(i,2)
            nm3=nm(i,3)
            nc1=nnodecode(nm1)
            nc2=nnodecode(nm2)
            nc3=nnodecode(nm3)
            ncele=nc1*nc2*nc3*noff(i)
            mju(nm1)=mju(nm1)+ncele
            mju(nm2)=mju(nm2)+ncele
            mju(nm3)=mju(nm3)+ncele
         enddo
         do i=1,np
            if((nnodecode(i).eq.1).and.(mju(i).eq.0)) then
               nnodecode(i)=0
               write(scratchMessage,9883) i
               call allMessage(INFO,scratchMessage)
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),'("node=",i0)') i
                  dryReasonInitialized(i) = .true.
               else
                  dryReason(i) = trim(dryReason(i)) // ", "
               endif
               dryReason(i) = trim(dryReason(i)) // " landlocked"
            endif
         enddo
 9883    format('node ',i0,' dried (landlocking).')
         !
         !jgf52.08.23: report the reasons why a node was initially dried.
         call logMessage(DEBUG,
     &      "nodes were set to an initially dry state "
     &      // "for the following reasons:")
         do i=1,np
            if (nnodecode(i).eq.0) then
               if (dryReasonInitialized(i).eqv..false.) then
                  write(dryReason(i),
     &            '("node=",i0," this node is dry for no good reason.")') i
                  dryReasonInitialized(i) = .true.
               endif
               call logMessage(DEBUG,trim(dryReason(i)))
            endif
         end do
         ! release memory
         deallocate(dryReason)
         deallocate(dryReasonInitialized)
         !
         ! finish initializing water surface elevation (eta2 and eta1)
         where (nnodecode.eq.0)
            eta2=h0-dp
         end where
      endif
<<<<<<< HEAD
!JLW: end subgrid if statement
      ENDIF
      eta1(:) = eta2(:)

C----------------------------------------------------------------------
=======
!jlw: end subgrid if statement
      endif

      !cms v56 moved update from totalAreaCalc here
      nodecode(1:np) = nnodecode(1:np)
      eta1(:) = eta2(:)

!----------------------------------------------------------------------
>>>>>>> Subject: Add additional source terms and update style
      end subroutine setColdWetDryStateAndWaterSurfaceElevation
!----------------------------------------------------------------------

