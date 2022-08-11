C----------------------------------------------------------------------------
C
C                           MODULE READ_GLOBAL
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 2.1 ( 03/17/2006 )
C
C                     current for ADCIRC v43.03   5/20/2003
C                 jgf Updated for ADCIRC v45.07  11/08/2005
C                 jgf Updated for ADCIRC v45.09  01/09/2006
C                 jgf Updated for ADCIRC v45.11  01/09/2006
C                 jgf Updated for ADCIRC v45.12  03/17/2006
C----------------------------------------------------------------------------

C
C  version 1.2 vjp 12/7/99
c  version 1.3 ral 10/10/01 revisions
C  version 1.6 meb 3/03 & ral 5/21/03
C  version 1.7 meb 4/04 added by jgf
C
C
C----------------------------------------------------------------------
C         S U B R O U T I N E     P R E P   R E A D   1 4
C----------------------------------------------------------------------
C     Reads the global ADCIRC grid data file for data decomposition
C     program ADCPP.  This version compatible with ADCIRC_v34.03
C     vjp  2/28/98
C---------------------------------------------------------------------
      subroutine prepRead14()
      use pre_global, only : x, y, dp, nneg, ibtypee, nvdllmsg, nvdll,
     &   nbdv, nvellmsg, nvell, lbcode, ibtype, nbvv, ibconnr,
     &   bar1, bar2, bar3, bar4, bar5, bar6, weir, weird, flux14_ary, 
     &   bcts14_ary, agrid, aperiodic_bc_ts, bndbcriver, exist_bc_ts, 
     &   exist_flux, iden, mnbou, mnvel, nbou, nboumsg, nelg, neta, 
     &   netamsg, nfover, nnodg, nope, nopemsg, nvel, nvelmsg, nweir, 
     &   river_above_msl, numLBCodeValues, nfluxf, node_dict, labels,
     &   NPERBC, IPERCONN, IPERPTR 
      use hashtable
      use memory_usage
      implicit none

      integer :: nbytes = 0
      INTEGER I,J,JW,K,ITEMP,ITYPE
      INTEGER DISC,BBN,IBP,I1
      INTEGER MAXNEIGH
      INTEGER NUMLINES
      INTEGER :: n1,n2,n3
c
c Dummy variables
c
      INTEGER DNBOU,DNVEL,DITYPE
      INTEGER,ALLOCATABLE :: DNVELL(:),DNBVV(:,:)
      INTEGER,ALLOCATABLE :: DIBTYPE(:)
      CHARACTER(80) DNBOUMSG,DNVELMSG
      CHARACTER(80),ALLOCATABLE ::  DNVELLMSG(:)
      CHARACTER(80) DUM
! kmd - variable for river in baroclinic simulation
      LOGICAL :: BC_Rivers = .False.

      ALLOCATE ( DNVELL(MNBOU),DNBVV(MNBOU,0:MNVEL),DIBTYPE(MNBOU) )
      ALLOCATE ( DNVELLMSG(MNBOU+1) )
      nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
      call memory_alloc(nbytes)
C
C--The Grid file was opened in SIZEUP performed a rewind, and is ready here
C
C--Read Grid Title
C
      EXIST_FLUX = 0
      READ(14,80) AGRID
C
C--Read Total Number of Elements and Nodes
C
      READ(14,*) NELG,NNODG
c     READ(14,80) NSIZES
c     READ(NSIZES,*) NELG,NNODG
c     CALL GETMSG(NSIZES,SIZEMSG)

C--Initialize hashtable
      ALLOCATE ( labels(1:NNODG) )
      CALL dict(node_dict,NNODG)
C
C--Read Nodal Coordinates and Bathymetry
C  If ICS=2 Will Convert later in read15
C
      NUMLINES=2
      DO I = 1,NNODG
         READ(14,*) labels(I),X(I),Y(I),DP(I)
         IF (labels(I).NE.I) THEN
            print *, I,J
            WRITE(*,'(A)') 'Error: Node Numbering not in Sequential Order'
            CALL EXIT(1) 
         ENDIF
         call add_ipair(node_dict,labels(I),I)
      ENDDO
C
C--Read Element Connectivity Table
C
       DO I = 1,NELG
          READ(14,*) J,ITEMP,n1,n2,n3
          IF (J.NE.I) THEN
            print *, I,J
            WRITE(*,'(A)') 'Error: Element Numbering not in Sequential Order'
            CALL EXIT(1) 
          ENDIF
          NNEG(1,I) = find(node_dict,n1)
          NNEG(2,I) = find(node_dict,n2)
          NNEG(3,I) = find(node_dict,n3)
       ENDDO
      NUMLINES=NUMLINES+NNODG+NELG
C...
C...Read Total Number of Open Boundary Segments
C...
      READ(14,80) NOPEMSG
      READ(NOPEMSG,*) NOPE
      NUMLINES=NUMLINES+1
C...
C...Read Total Number of Open Boundary Forcing Nodes
C...
      READ(14,80) NETAMSG
      READ(NETAMSG,*) NETA
      NUMLINES=NUMLINES+1
C...
C...Read Number of Nodes on Open Boundary Segment
C...and Segment Nodes Numbers
C...
      J=0
      IBTYPEE(:) = 0 !jgf49.32 initialize all elements of the array to zero
      DO K=1,NOPE
         READ(14,80) NVDLLMSG(K) !jgf49.32: TODO: IBTYPEE should be read here
         NUMLINES=NUMLINES+1
         READ(NVDLLMSG(K),*) NVDLL(K)
         DO I=1,NVDLL(K)
            READ(14,*) n1
            NBDV(K,I) = find(node_dict,n1)
            NUMLINES=NUMLINES+1
         ENDDO
         J=J+NVDLL(K)
      ENDDO
C
      IF (NETA.NE.J) THEN
         print *, "Total Number of Boundary Nodes = ",J
         print *, "This exceeds NETA = ",NETA
         IF (NFOVER.EQ.1) THEN  !jgf51.21.27: Has NFOVER been read at this point?
            NETA = J
            print *, "ADCPP corrected this error"
         ELSE
            CALL EXIT(1)
         ENDIF
      ENDIF
C
C--Read Total Number of Land Boundary Segments
C
      READ(14,80) NBOUMSG
      READ(NBOUMSG,*) NBOU
C
C--Read Total of Land Boundary Nodes
C
      READ(14,80) NVELMSG
      READ(NVELMSG,*) NVEL
C
C--Read Number of Nodes in the Land Boundary Segment and Boundary Type
C  and construct LBCODE array for read15 routine
C
      J=0
      NWEIR = 0
C.... DW
      NPERBC = 0 ;      
C.... DW
      ! jgf51.52.23: The LBCODE array is dimensioned by MNVEL which is 
      ! larger than the number of values in the LBCODE array. So we
      ! need to initialize to a value that will hopefully cause obvious
      ! isses if an initialized value is accidentally used.
      LBCODE(:) = -99
      DO K = 1,NBOU
C
         READ(14,80) NVELLMSG(K)
         READ(NVELLMSG(K),*) NVELL(K),IBTYPE(K)
         ITYPE = IBTYPE(K)
! kmd - added for rivers in a baroclinic simulation
         IF (ABS(ITYPE/100).EQ.1) THEN
            IF (IDEN.LE.0) THEN
               ITYPE=(ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
            ELSE IF (IDEN.GT.0) THEN
               APERIODIC_BC_TS=.TRUE.
               BndBCRiver=.TRUE.
               ITYPE=(ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
            END IF
         END IF
C
         DO I=1, NVELL(K)
            J = J+1
            LBCODE(J) = ITYPE
         ENDDO
         !
         ! jgf51.21: Cleaned up reading of boundary arrays by replacing
         ! a complex set of if/thens with a simple select case.
         select case(ITYPE)
         case(0,10,20,30,40,1,11,21,41)
            DO I=1,NVELL(K)
               READ(14,*) n1
               NBVV(K,I) = find(node_dict,n1)
              IBCONNR(K,I) = 0
            ENDDO
         case(3,13,23)
            DO I=1,NVELL(K)
               READ(14,*) n1,BAR1(K,I),BAR2(K,I)
               NBVV(K,I) = find(node_dict,n1)
              IBCONNR(K,I) = 0
            ENDDO
         case(4,24,64)
            DO I=1,NVELL(K)
               READ(14,*) n1,IBCONNR(K,I),
     &                    BAR1(K,I),BAR2(K,I),BAR3(K,I)
               NBVV(K,I) = find(node_dict,n1)
               !--Construct List of WEIR nodes and their duals
               NWEIR = NWEIR+1
               WEIR(NWEIR) = NBVV(K,I)
               WEIRD(NWEIR) = IBCONNR(K,I)
            ENDDO
         case(5,25)
            DO I =1,NVELL(K)
               READ(14,*) n1,n2,
     &                   BAR1(K,I),BAR2(K,I),BAR3(K,I),
     &                   BAR4(K,I),BAR5(K,I),BAR6(K,I)
               NBVV(K,I) = find(node_dict,n1)
               IBCONNR(K,I) = find(node_dict,n2)
               !--Construct List of WEIR nodes and their duals
               NWEIR = NWEIR+1
               WEIR(NWEIR) = NBVV(K,I)
               WEIRD(NWEIR) = IBCONNR(K,I)
            ENDDO
         case(2,12,22,32,52)
           NFLUXF=1 !jgf51.52.23
           DO I=1,NVELL(K)
              EXIST_FLUX = EXIST_FLUX + 1
              READ(14,*) n1
              NBVV(K,I) = find(node_dict,n1)
              IBCONNR(K,I) = 0
!kmd - added river boundary conditions for baroclinic
              IF (BndBCRiver) THEN
                 EXIST_BC_TS=EXIST_BC_TS+1
              END IF
           ENDDO
           BndBCRiver=.FALSE.
         case(94)
C..... DW, periodic bc
           DO I=1,NVELL(K)
              IF ( .NOT. ALLOCATED(IPERCONN) ) THEN
                 ALLOCATE( IPERCONN(MNVEL,2) ) ;
                 IPERCONN = 0 ;

                 ALLOCATE( IPERPTR(NBOU+1) ) ;
                 IPERPTR = 1 ; 
              END IF
              NPERBC = NPERBC+1 ;
C
C--   Construct List of periodic nodes: primary and their slaves
C
              READ(14,*) IPERCONN(NPERBC,1), IPERCONN(NPERBC,2) ;
              NBVV(K,I) = IPERCONN(NPERBC,1) ;
C
              IPERPTR(K+1:NBOU+1) = (IPERPTR(K)) + NVELL(K) ;      
           ENDDO
C..... DW

         case default
            print *, "IBTYPE not set correctly for segment ",K
            CALL EXIT(1)
         end select
      ENDDO
      !jgf51.52.23: Record the actual number of values stored in the
      !LBCODE array.
      numLBCodeValues = J


C     jgf45.06 added the following section
C MEB --------------------------------------------------
C MEB This entire section written 04/01/04 to accomodate
C MEB division of fort.20 files between PE directories.
C MEB This must be modified if fort.14 format changes.
C MEB --------------------------------------------------
      IF (EXIST_FLUX.GT.0) THEN
         REWIND(14)
         ALLOCATE ( FLUX14_ARY(EXIST_FLUX) )
!kmd - added river boundary conditions for baroclinic
         ALLOCATE ( BCTS14_ARY(EXIST_BC_TS) )
         EXIST_FLUX=0
         EXIST_BC_TS=0
c
C  SKIP OVER UNNEEDED LINES
c
         DO I=1,NUMLINES
            READ(14,*) DUM
         ENDDO
C
         READ(14,80) DNBOUMSG
         READ(NBOUMSG,*) DNBOU
C
C--Read Total of Land Boundary Nodes
C
         READ(14,80) DNVELMSG
         READ(NVELMSG,*) DNVEL
c
         DO K=1,DNBOU
            READ(14,80) DNVELLMSG(K)
            READ(DNVELLMSG(K),*) DNVELL(K), DIBTYPE(K)
            DITYPE=DIBTYPE(K)
! kmd - added in boundary conditions for rivers in a baroclinic simulation
            IF (ABS(DITYPE/100).EQ.1) THEN
               IF (IDEN.LE.0) THEN
                  DITYPE=(ABS(DITYPE)-100)*(DITYPE/ABS(DITYPE))
               ELSE IF (IDEN.GT.0) THEN
                  BC_Rivers=.TRUE.
                  DITYPE=(ABS(DITYPE)-100)*(DITYPE/ABS(DITYPE))
               END IF
            END IF
            ! jgf50.21 Use select case and add IBTYPE 32.
            select case(DITYPE)
            case(2,12,22,32,52)
               DO I=1, DNVELL(K)
                  EXIST_FLUX = EXIST_FLUX + 1
                  READ(14,*) n1
                  DNBVV(K,I) = find(node_dict,n1)
                  FLUX14_ARY(EXIST_FLUX) = DNBVV(K,I)
! kmd - added in boundary conditions for rivers in a baroclinic simulation
                  IF (BC_Rivers) THEN
                     EXIST_BC_TS=EXIST_BC_TS+1
                     BCTS14_ARY(EXIST_BC_TS)=DNBVV(K,I)
                  END IF
                  IF (DP(DNBVV(K,I)).LT.0.d0) THEN
                     River_above_MSL=.TRUE.  ! kmd - Evan's changes for rivers above MSL
                  END IF
               ENDDO
            case default
               DO I=1, DNVELL(K)
                  READ(14,*) n1
                  DNBVV(K,I) = find(node_dict,n1)
               ENDDO
            end select
            BC_Rivers = .FALSE.  ! kmd - reset logic for next river boundary segment
         ENDDO
         DEALLOCATE ( DNVELL,DNBVV,DIBTYPE )
         DEALLOCATE ( DNVELLMSG )
         nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
         call memory_dealloc(nbytes)
      ENDIF


c     jgf45.06 end of section written by MEB04/01/04
C
C--Close Global Grid file
C
      CLOSE(14)
      call memory_status()
C
 80   format(a80)
c
      return
C----------------------------------------------------------------------
      end subroutine prepRead14
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C                        S U B R O U T I N E
C            G E T   M E S H   P A R A M E T E R S   X D M F
C----------------------------------------------------------------------
C     jgf51.21.29: If the mesh file is in XDMF format, the adcirc module
C     mesh.F is used to read the data. In order to get the data into the
C     pre_global module, this subroutine is required. Someday, adcprep
C     will no longer use pre_global, and will just rely on the adcirc
C     global module, and we will be able to remove this subroutine.
C----------------------------------------------------------------------
      subroutine getMeshParametersXDMF()
      use memory_usage
      use presizes, only : fluxBoundary, mne, mnp, mnei, mneta, mnope,
     &   skipNeta
      use pre_global, only : x, y, dp, nneg, ibtypee, nvdllmsg, nvdll,
     &   nbdv, nvellmsg, nvell, lbcode, ibtype, nbvv, ibconnr,
     &   bar1, bar2, bar3, weir, weird, flux14_ary, bcts14_ary, agrid,
     &   aperiodic_bc_ts, bndbcriver, exist_bc_ts, exist_flux, iden,
     &   mnbou, mnvel, nbou, nboumsg, nelg, neta, netamsg, nfover,
     &   nnodg, nope, nopemsg, nvel, nvelmsg, nweir, river_above_msl,
     &   alloc_main1
      use mesh, only : getMeshSizesForPrep, getMeshForPrep, nm
      use boundaries, only : getBoundarySizesForPrep,
     &    getBoundariesForPrep
      implicit none
      integer, allocatable :: nneigh(:)
      integer :: i, j, k, nbytes
      !
      !
      call getMeshSizesForPrep(nnodg, nelg)
      call getBoundarySizesForPrep(nbou, nvel, neta, nope,
     &   nweir, exist_flux)
      write(netamsg,'(i0)') neta
      write(nopemsg,'(i0)') nope
      write(nboumsg,'(i0)') nbou
      write(nvelmsg,'(i0)') nvel
      mne = nelg
      mnp = nnodg
      mnope = nope
      mneta = neta
      skipNeta = neta
      mnbou = nbou
      mnvel = 2*nvel
      if (exist_flux.ne.0) then
         fluxBoundary = .true.
      endif
      write(*,'(a,i0)') "neta =", neta
      call alloc_main1()
      !
      ! get the element table
      do i=1,3
         nneg(i,:) = nm(:,i)
      end do
      !
      ! get mesh parameters
      call getMeshForPrep(agrid, x, y, dp)
      call getBoundariesForPrep(nvdll, nvell, ibtypee, nbdv, lbcode, ibtype, nbvv,
     &   ibconnr, bar1, bar2, bar3, weir, weird, iden, nfover,
     &   exist_flux, flux14_ary)
      do k=1,nope
         write(nvdllmsg(k),'(i0)') nvdll(k)
      end do
      do k=1,nbou
         write(nvellmsg(k),'(i0)') nvell(k)
      end do
      !
      ! compute number of neighbors of each node
      allocate (nneigh(mnp))
      nbytes = 4*mnp
      call memory_alloc(nbytes)
      NNeigh(:) = 0
      DO I=1,nelg
        NNeigh(nneg(1,i))=NNeigh(nneg(1,i))+1
        NNeigh(nneg(2,i))=NNeigh(nneg(2,i))+1
        NNeigh(nneg(3,i))=NNeigh(nneg(3,i))+1
      ENDDO
C     determine the maximum NNeigh
      mnei=maxval(nneigh)
      mnei = mnei+1
      deallocate(nneigh)
C----------------------------------------------------------------------
      end subroutine getMeshParametersXDMF
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C     Reads the global ADCIRC Input Data File for data decomposition
C     program ADCPP.  This version compatible with ADCIRC_v34.03
C     vjp  3/28/98
C---------------------------------------------------------------------
      SUBROUTINE READ15()
      USE PRE_GLOBAL
      USE SIZES, ONLY : MNHARF
      USE MESH, ONLY : CYLINDERMAP
      USE GLOBAL, ONLY : parse, a2f, hss, screenUnit, STATNAME,
     &                   STATNAMEV, STATNAMEM, STATNAMEC, ERROR, DEBUG,
     &                   allMessage, setMessageSource, deg2rad,
     &                   unsetMessageSource, INFO, logMessage,
     &                   scratchMessage, screenMessage
      USE HARM, ONLY : NHAINC, ALLOC_HA, NAMEFR, HAFREQ, HAFF, HAFACE,
     &                 THAS, THAF, FMV, NHASE, NHASV, NHAGE, NHAGV,
     &                 IHARIND, ALLOC_MAIN14, CHARMV
      USE PRESIZES, ONLY : SZ
      USE PREP_WEIR
      USE NODALATTRIBUTES, ONLY : outputTau0
      use memory_usage
      IMPLICIT NONE

C    ===================================================================
C      MCF ADDITION FOR READING AND WRITING STATION INFO IN NETCDF
C
C      REAL, EXTERNAL              :: a2f
      CHARACTER(132) STATLINE
      CHARACTER(50) LVAR(3)
C
C    ===================================================================


      integer :: nbytes = 0
      REAL(SZ) RSTIMINC
      INTEGER N1, N2, N3, KMIN, JG, INDX
      INTEGER I,J,K,M
      INTEGER IG1,IG2,IG3,IL1,IL2,IL3
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3
      REAL(8) AE, AEMIN, AREASK, AA
      INTEGER NBV(MNVEL)
      INTEGER ios  ! status of read operations, 0=success
Casey 090825: Added dummy variables.
      INTEGER IDUM
      REAL(8) RDUM
! kmd - add a variable for boundary conditions
      INTEGER ITYPE
C
!   tcm v50.66.02 -- added timebathycontrol namelist related variables
      INTEGER :: ios_nddt
      INTEGER :: ios_metCon
      INTEGER :: ios_tvw
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      INTEGER :: IOS_WC
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      INTEGER :: IOS_OWINC
#endif
C...tcm v51.20.03 Local variables for external station location files
      INTEGER :: NSTAE2,NSTAV2,NSTAC2,NSTAM2
      INTEGER :: IOS_STATIONS
      INTEGER :: ios_wetDry
      namelist /wetDryControl/ outputNodeCode, outputNOFF,
     &          noffActive, StatPartWetFix, How2FixStatPartWet,
C.... DMW 202207 Add options for controlling slope/wind limiting and
C.... "new" formula for computing vel in wetting/drying routine 
     &          slim, windlim, directvelWD, useHF
C.... DMW
      ! jgf52.08.02 : For inundation output files.
      namelist /inundationOutputControl/ inundationOutput, inunThresh

C
C.....tcm v52.30.01 added for Smag_Control
      INTEGER :: IOS_Smag = 0
      ! jgf: Added support for water level correction
      NAMELIST /dynamicWaterLevelCorrectionControl/ dynamicWaterLevelCorrectionFileName,
     &  dynamicWaterLevelCorrectionSkipSnaps, dynamicWaterlevelCorrectionMultiplier,
     &  dynamicWaterLevelCorrectionRampStart, dynamicWaterLevelCorrectionRampEnd,
     &  dynamicWaterLevelCorrectionRampReferenceTime

c.....tcm v51.20.03 initialize external station location specifications
C..... to false, meaning if locations are supplied they are in the fort.15
C..... file.
      USE_ELEV_STAT_FILE = .FALSE.
      USE_VEL_STAT_FILE = .FALSE.
      USE_MET_STAT_FILE = .FALSE.
      USE_CONC_STAT_FILE = .FALSE.
C
      found_tbc_nml = .false.   !flag to determine if the timebathycontrol namelist was present
      found_metCon_nml = .false.   !flag to determine if the metControl namelist was present
      found_tvw_nml = .false.
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      FOUND_WC_NML = .FALSE.
      foundCorrectionControlNamelist = .false.
#ifdef ADCNETCDF
      FOUND_OWINC_NML = .FALSE.
#endif
      call setMessageSource("read15")
#ifdef READ_GLOBAL_TRACE
      call screenMessage(DEBUG,"Enter.")
#endif
C
C--The Run Info file was opened in SIZEUP, performed a rewind, and is ready here
C
C--Run Description and Run Identification
C

c...  tcm v50.66.02 Addtions for Time Varying Bathymetry
c...  read through the fort.15 file for the namelist (TimeBathyControl) for
c...  the time varying bathymetry.  This namelist must be at the bottom of the
c...  fort.15 file. If found, then set the appropriate values (btiminc,bchgtiminc,
c...  and nddt).  If the namelist is not there, then the time varying bathymetry
c...  will not be used.
C...
C...  After this search and read, we will close the file and then reopen it
c...  for further processing the traditional non-namelist components.
C...
      NDDT = 0  !Set default value to be no time varying bathymetry
      READ(UNIT = 15,NML = TimeBathyControl,IOSTAT = IOS_NDDT)

      IF (IOS_NDDT < 0) THEN
c.....   it is possible for the namelist to be present in the file and occuring at the end
c.....   of the file with no line breaks after the ending "\" which causes the iostat to
c.....   return a negative value.  By checking to be sure a namelist variable was set to
c....    a non-default value we can determine this was the case.
         IF (NDDT.NE.0) THEN
!            WRITE(*,*) 'NAMELIST PRESENT, BUT AT THE END OF FILE WITH',
!     &                 ' NO ADVANCING CHARACTER'
            found_tbc_nml = .true.
         ELSE
            call logMessage(INFO,
     &         'The timeBathyControl namelist was not found.')
         ENDIF
      ELSEIF (IOS_NDDT == 0) THEN
!         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
         found_tbc_nml = .true.
      ELSE
         found_tbc_nml = .true.
         WRITE(*,*)
     & 'THERE WAS A PROBLEM PROCESSING THE TimeBathyControl NAMELIST'
         WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCIRC NOW.'
         CALL EXIT(1)   !THERE IS a STOP HERE
      ENDIF

      REWIND(15)  !Return to the top of the fort.15 file (instead of closing then opening)
c..

c... tcm v50.79 added
      READ(UNIT = 15, NML = metControl,IOSTAT = ios_metCon)
      if ( ios_metCon.ge.0) then
         found_metCon_nml = .true.
         call logMessage(INFO,
     &         'The metControl namelist was found.')
      else
         call logMessage(INFO,
     &         'The metControl namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file

Cobell - Added TVW Namelist
      READ(15,NML=tvwControl,IOSTAT=IOS_TVW)
      IF(IOS_TVW.EQ.0)THEN
        FOUND_TVW_NML = .TRUE.
      ELSE
        FOUND_TVW_NML = .FALSE.
      ENDIF
      REWIND(15)

Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      READ(UNIT=15,NML=waveCoupling,IOSTAT=IOS_WC)
      IF(IOS_WC.EQ.0)THEN
         FOUND_WC_NML = .TRUE.
         call logMessage(INFO,
     &         'The waveCoupling namelist was found.')
      else
         call logMessage(INFO,
     &         'The waveCoupling namelist was not found.')
      ENDIF
      REWIND(15)
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      READ(UNIT=15,NML=owiWindNetcdf,IOSTAT=IOS_OWINC)
      IF(IOS_OWINC.EQ.0)THEN
         FOUND_OWINC_NML = .TRUE.
         call logMessage(INFO,
     &         'The owiWindNetcdf namelist was found.')
      else
         call logMessage(INFO,
     &         'The owiWindNetcdf namelist was not found.')
      ENDIF
      REWIND(15)
#endif

      READ(UNIT = 15, NML = wetDryControl,IOSTAT = ios_wetDry)
      if ( ios_wetDry.ge.0) then
         foundWetDryControlNameList = .true.
         call logMessage(INFO,
     &         'The wetDryControl namelist was found.')
      else
         call logMessage(INFO,
     &         'The wetDryControl namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file

      read(unit=15,nml=inundationOutputControl,iostat=ios)
      if ( ios.ge.0) then
         foundInundationOutputControlNamelist = .true.
         call logMessage(INFO,
     &         'The inundationOutput namelist was found.')
      else
         call logMessage(INFO,
     &         'The inundationOutput namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file
      !
      !jgf:   O F F S E T   C O N T R O L   N A M E L I S T
      read(unit=15,nml=dynamicWaterLevelCorrectionControl,iostat=ios)
      select case(ios)
      case(0) ! zero indicates success
         foundCorrectionControlNamelist = .true.
         call logMessage(INFO,
     &         'The dynamicWaterLevelCorrectionControl namelist was found.')
      case(1:)
         write(scratchMessage,'(a,i0,a)')
     &      'Could not read dynamicWaterLevelCorrectionControl'//
     &     ' namelist. The Fortran i/o error code was ',ios,'.'
         call allMessage(ERROR,scratchMessage)
         CALL EXIT(1)
      case default
         call logMessage(INFO,
     &      'End-of-file when searching for dynamicWaterLevelCorrectionControl namelist.')
      end select
      rewind(15)  !Return to the top of the fort.15 file

c...  tcm v52.30.01 Addtions for Smagorinski Limits
c...  read through the fort.15 file for the namelist (Smag_Control) for
c...  the inputs for Smagorinski eddy viscosity limits.  This namelist must be
c...  at the bottom of the fort.15 file. If found, then set the appropriate
c...  values (SMAG_COMP_FLAG,SMAG_UPPER_LIM,SMAG_LOWER_LIM)
c...  If the namelist is not there, then the Smagorinski Limits
c...  will not be used.
C...
C...  After this search and read, we will close the file and then reopen it
c...  for further processing the traditional non-namelist components.
C...
      READ(UNIT = 15,NML = Smag_Control,IOSTAT = IOS_Smag)

      IF (IOS_SMAG < 0) THEN
c.....   it is possible for the namelist to be present in the file and occuring at the end
c.....   of the file with no line breaks after the ending "\" which causes the iostat to
c.....   return a negative value.  By checking to be sure a namelist variable was set to
c....    a non-default value we can determine this was the case.
         IF (SMAG_COMP_FLAG .eqv. .TRUE.) THEN
!            WRITE(*,*) 'NAMELIST PRESENT, BUT AT THE END OF FILE WITH',
!     &                 ' NO ADVANCING CHARACTER'
            found_smag_nml = .true.
         ELSE
            call logMessage(INFO,
     &         'The Smag_Control namelist was not found.')
         ENDIF
      ELSEIF (IOS_SMAG == 0) THEN
!         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
         found_smag_nml = .true.
      ELSE
         found_smag_nml = .true.
         WRITE(*,*)
     & 'THERE WAS A PROBLEM PROCESSING THE Smag_Control NAMELIST'
         WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCPREP NOW.'
         CALL EXIT(1)   !THERE IS a STOP HERE
      ENDIF

      REWIND(15)  !Return to the top of the fort.15 file (instead of closing then opening)
c..
      read(unit=15,nml=WarnElevControl,iostat=ios)
      if ( ios.ge.0) then
         FOUND_WARNELEV_NML = .true.
         call logMessage(INFO,
     &         'The WarnElev namelist was found.')
      else
         call logMessage(INFO,
     &         'The WarnElev namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file
      
      read(unit=15,nml=AliDispersionControl,iostat=ios)
      if ( ios.ge.0) then
         FOUND_ALIDISP_NML = .true.
         call logMessage(INFO,
     &         'The AliDispersion namelist was found.')
      else
         call logMessage(INFO,
     &         'The AliDispersion namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file
      
      read(unit=15,nml=densityControl,iostat=ios)
      if ( ios.ge.0) then
         FOUND_DENSITY_NML = .true.
         call logMessage(INFO,
     &         'The density control namelist was found.')
      else
         call logMessage(INFO,
     &         'The density control namelist was not found.')
      endif
      rewind(15)  !Return to the top of the fort.15 file

      READ(15,80) RUNDES
C
      READ(15,80) RUNID
C
      READ(15,80) OVERMSG
      READ(OVERMSG,*) NFOVER
      IF (NFOVER.EQ.1) THEN
c       print *, "Non-fatal errors will be corrected"
      ELSE
c      print *, "Non-fatal errors will stop execution"
      ENDIF
C
      READ(15,80) ABOUTMSG
      READ(ABOUTMSG,*) NABOUT
C
      READ(15,80) SCREENMSG
      READ(SCREENMSG,*) NSCREEN
      screenUnit = 6 ! always write to screen in adcprep
C
      READ(15,80) HOTMSG
      READ(HOTMSG,*) IHOT

C
      READ(15,80) ICSMSG
      READ(ICSMSG,*) ICS
      ICS = abs(ICS) ! WJP: in case we use rotated coordinates
      SELECT CASE( ICS )
      CASE(0,1,2,20:24)
        ! ICS is set correctly
      CASE DEFAULT
        print *, "ICS set incorrectly"
        CALL EXIT(1)
      END SELECT
C
      READ(15,80) IMMSG
      READ(IMMSG,*) IM
      IF(IM.EQ.2) THEN
        PRINT *, "DSS Model type not presently supported"
        CALL EXIT(1)
      ENDIF
C     jgf46.28 Read IDENMSG if necessary
      IF (CBaroclinic) READ(15,80) IDENMSG
C
      READ(15,80) IBFMSG
      READ(IBFMSG,*) NOLIBF
      IF((NOLIBF.LT.0).OR.(NOLIBF.GT.3)) THEN
         print *, "Value for NOLIBF not allowed"
         CALL EXIT(1)
      ENDIF
C
      READ(15,80) IFAMSG
      READ(IFAMSG,*) NOLIFA
      IF ((NOLIFA.LT.0).OR.(NOLIFA.GT.3)) THEN
         print *, "Value for NOLIFA not allowed"
         CALL EXIT(1)
      ENDIF
C
      READ(15,80) ICAMSG
      READ(ICAMSG,*) NOLICA
      IF ((NOLICA.LT.0).OR.(NOLICA.GT.1)) THEN
         print *, "Value for NOLICA not allowed"
         CALL EXIT(1)
      ENDIF
C
      READ(15,80) ICATMSG
      READ(ICATMSG,*) NOLICAT
      IF ((NOLICAT.LT.0).OR.(NOLICAT.GT.1)) THEN
         print *, "Value for NOLICAT not allowed"
         CALL EXIT(1)
      ENDIF
      !
      ! jgf52.05: Removed this error message because it is alarming
      ! and because it contradicts common practice.
C      IF ((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
C         print *, "NOLIFA and NOLICAT are inconsistent"
C         print *, "May lead to mass balance problems"
C         IF(NFOVER.EQ.1) THEN
C            print *, "Since NFOVER=1, Program will continue"
C        ELSE
C            CALL EXIT(1)
C         ENDIF
C      ENDIF
C
      READ(15,80) NWPMSG
      READ(NWPMSG,*) NWP
      IF (NWP.gt.0) THEN !jgf46.00 read nodal attribute labels
         ALLOCATE(useNodalAttrNames(NWP))
         nbytes = 4*nwp
         call memory_alloc(nbytes)
         DO I=1, NWP
            READ(15,80) useNodalAttrNames(I)
         ENDDO
      ENDIF
C
      READ(15,80) NCORMSG
      READ(NCORMSG,*) NCOR
      IF ((NCOR.NE.0).AND.(NCOR.NE.1)) THEN
         print *, "Value for NCOR not allowed"
         IF (NFOVER.EQ.1) THEN
           NCOR = 0
           print *, "NCOR has been reset to 0"
         ELSE
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      IF ((ICS.EQ.1).AND.(NCOR.EQ.1)) THEN
        print *, "ICS=1 and NCOR=1 may lead to geometric distortions"
         IF(NFOVER.EQ.1) THEN
           print *, "Program will continue with these input values"
           print *, "for large domains it is recommended to use ICS=2"
         ELSE
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      READ(15,80) NTIPMSG
      READ(NTIPMSG,*) NTIP
      IF ((NTIP.LT.0).OR.(NTIP.GT.2)) THEN
         print *, "Value for NTIP not allowed"
         IF(NFOVER.EQ.1) THEN
           NTIP = 0
           print *, "NTIP has been reset to 0"
         ELSE
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      IF ((ICS.EQ.1).AND.(NTIP.GE.1)) THEN
         print *, "ICS=1 & NTIP >= 1 may lead to geometric distortions"
         print *, "for large domains it is recommended to use ICS=2"
         IF (NFOVER.EQ.1) THEN
            print *, "Program will continue with these input values"
         ELSE
            CALL EXIT(1)
         ENDIF
      ENDIF
C
      NRS=0
      NCICE = 0  !tcm v49.64.01 -- added for ice
      READ(15,80) NWSMSG
      READ(NWSMSG,*) NWS

C....tcm v49.46 -- added logic to handle multiple NRS types (100's place)
C.....tcm v49.64.01 Additions for ice
      IF(NWS.EQ.0) THEN
         NWS = 0
         NRS = 0
         NCICE = 0
      ELSE
         NCICE = INT(ABS(NWS)/1000)
         NRS=INT((ABS(NWS) - NCICE*1000)/100)
         NWS=INT((ABS(NWS)- NCICE*1000 - NRS*100))*INT(NWS/ABS(NWS))
      ENDIF


C     jgfdebug46.02 added NWS=45
C     jgf46.02 Added NWS=8.
C     jgf46.16 merged:
C     cf & cm added NWS=9: asymmetric hurricane wind model
C     jie added NWS=20: generalized asymmetric vortex model
C     sb46.28sb01 added NWS=12: OWI format
C     jgf49.0804 Added NWS29 VortexOWI format.
C     jgf: Added NWS30.
Casey 180318: Added NWS=13
        IF((NWS.NE.0).AND.    (NWS.NE.1 ) .AND.(ABS(NWS).NE.2).AND.
     &       (NWS.NE.3).AND.(ABS(NWS).NE.4) .AND.(ABS(NWS).NE.5).AND.
     &       (ABS(NWS).NE.45).AND.(ABS(NWS).NE.6).AND.
     &       (ABS(NWS).NE.7).AND.   ! xyc52.30
     &       (ABS(NWS).NE.8).AND.(ABS(NWS).NE.15).AND. !jgf50.38.05 Added NWS=15
     &       (ABS(NWS).NE.12).AND.(ABS(NWS).NE.13).AND.
     &       (ABS(NWS).NE.14).AND.(ABS(NWS).NE.19).AND.
     &       (ABS(NWS).NE.20) .AND. (ABS(NWS).ne.30) .AND. (NWS.NE.10) .AND.
     &       (ABS(NWS).NE.29) .AND. (ABS(NWS).ne.16) .AND.
     &       (ABS(NWS).NE.17))THEN                             !PV added NWS=17 for NUOPC atm. forcing
           print *,"ERROR: Value for NWS not supported by parallel code"
           CALL EXIT(1)
        ENDIF
C
C... TCM v49.64.02 -- added
C...  BE SURE NWS AND NCICE ARE COMPATABLE
C...  xyc v52.30: note:Not sure if NWS.EQ.7 here should be removed 
      IF((NCICE.GT.0).AND.
     &       ((NWS.EQ.1).OR.(NWS.EQ.2).OR.(NWS.EQ.7))) THEN
         PRINT*,'NCICE = ',NCICE
         PRINT*,'NWS = ', NWS
         PRINT*, "Your selection of NCICE (a UNIT 15 input ",
     &        "parameter) is not compatable"
         PRINT*,"with your NWS value.  NCICE is not allowed for",
     &         " abs(NWS)=1,2, or 7."
         CALL EXIT(1)
      ENDIF

C     jgf46.08 Modified to accomodate fine grained ramp functions.
      READ(15,80) RAMPMSG
      READ(RAMPMSG,*) NRAMP
Corbitt 120322: Includes Zach's Wind Ramping
      IF ((NRAMP.NE.0).AND.(NRAMP.GT.8)) THEN
         print *, "Value for NRAMP not allowed"
         IF (NFOVER.EQ.1) THEN
            print *, "Program will override and use NRAMP = 0"
            NRAMP = 0
         ELSE
            CALL EXIT(1)
         ENDIF
      ENDIF
C
      READ(15,80) GMSG
      READ(GMSG,*) G
      IF ((ICS.EQ.1).AND.(G.NE.9.81d0)) THEN
         IF ((NCOR.EQ.1).OR.(NTIP.EQ.1)) THEN
            print *, "G not consistent with ICS=1"
            print *, "in conjunction with NTIP=1 and/or NCOR=1"
            IF(NFOVER.EQ.1) THEN
               print *, "Program will override and set G=9.81"
               print *, "check to see that all input has SI units"
               G = 9.81d0
            ELSE
               CALL EXIT(1)
            ENDIF
         ENDIF
      ENDIF
C
      IF ((ICS.NE.1).AND.(G.NE.9.81d0)) THEN
         print *, "G not consistent with ICS=2"
         IF(NFOVER.EQ.1) THEN
            print *, "Program will override and set G = 9.81 m/sec*sec"
            print *, "check to see that all input has SI units"
            print *, "execution will continue"
            G = 9.81d0
         ELSE
            CALL EXIT(1)
         ENDIF
      ENDIF
C
C     jgf47.11 Allow the limits of time varying tau0 to be read in from
C     the fort.15 file
      READ(15,80) TAU0MSG
      READ(TAU0MSG,*) TAU0
      if ( abs((int(tau0)-tau0-0.1d0)).lt.1d-6 ) then
         outputTau0 = .true.
      else
         outputTau0 = .false.
      endif
      IF ( (TAU0.LE.-5.d0).AND.(TAU0.GT.-6.d0) ) THEN
         READ(15,80) TAU0LIMMSG
         READ(TAU0LIMMSG,*) Tau0FullDomainMin, Tau0FullDomainMax
      ENDIF

C
      READ(15,80) DTMSG
      READ(DTMSG,*) DT
C
      READ(15,80) STATMSG
      READ(STATMSG,*) STATIM
C
      READ(15,80) REFTMSG
      READ(REFTMSG,*) REFTIM
C
C--If wind stress and surface pressures are applied process this.
C
      IF((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMMSG ! sb46.28sb03
      IF((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMMSG ! sb46.28sb03
C
      IF(NWS.EQ.3) THEN
         READ(15,80) WSMSG1
         READ(15,80) WSMSG2
C....    TCM V49.64.01 ADDITIONS FOR ICE
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC
         ELSEIF((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
         ENDIF
      ENDIF
C
C     jgfdebug46.02 Added NWS=45
C     jgf46.02 Added NWS=8
C     jgf46.16 Merged:
C     cf & cm Added NWS=9: asymmetric hurricane wind model
C     rjw Added NWS=19: asymmetric hurricane wind model v2.0
C     jie Added NWS=20: generalized asymmetric vortex model
C     sb46.28sb01 added NWS=12: OWI format
C
C....TCM V49.64.01 ADDITIONS FOR ICE, broke nws=2 out
      IF(ABS(NWS).EQ.2) THEN
         READ(15,80) WSMSG1  !tcm v51.17 -- added to correct bug
         IF(NRS.EQ.0) READ(WSMSG1,*) WTIMINC
         IF(NRS.GE.1) READ(WSMSG1,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF

      IF((ABS(NWS).EQ.4).OR.(ABS(NWS).EQ.5)
     &     .OR.(ABS(NWS).EQ.45).OR.(ABS(NWS).EQ.8).OR.
     &   (ABS(NWS).EQ.9).OR.(ABS(NWS).EQ.12).OR.(ABS(NWS).EQ.14)
     &     .OR.(NWS.EQ.10) !yf: Added NWS=10 for GFS Met
     &     .OR.(ABS(NWS).EQ.15) !jgf50.38.05: Added NWS=15 for HWind.
     &     .OR.(ABS(NWS).EQ.20) !jie: Added NWS=20 for GAHM
     &     .OR.(ABS(NWS).EQ.13) ! arc added nws13 190110
     &     .OR.(ABS(NWS).EQ.19).OR.(ABS(NWS).EQ.16)) THEN
         READ(15,80) WSMSG1
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) READ(WSMSG1,*) WTIMINC
         IF((NCICE.EQ.0).AND.
     &   ((NRS.EQ.1).OR.(NRS.EQ.2).OR.(NRS.EQ.4))) READ(WSMSG1,*) WTIMINC,RSTIMINC
         IF((NCICE.EQ.0).AND.(NRS.EQ.3)) THEN  !Casey 090825: Fix for NRS = 3.
            IF(ABS(NWS).EQ.8)THEN
               READ(WSMSG1,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
            ELSE
               READ(WSMSG1,*) WTIMINC,RSTIMINC
            ENDIF
         ENDIF
         IF((NCICE.GE.1).AND.(NRS.EQ.0)) READ(WSMSG1,*) WTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.((NRS.EQ.1).OR.(NRS.EQ.2).OR.(NRS.EQ.4)))
     &                       READ(WSMSG1,*) WTIMINC,RSTIMINC,CICE_TIMINC
         IF((NCICE.GE.1).AND.(NRS.EQ.3)) THEN  !Casey 090825: Fix for NRS = 3.
            IF(ABS(NWS).EQ.8)THEN
               READ(WSMSG1,*) IDUM,IDUM,IDUM,IDUM,IDUM,
     &                        RDUM,RSTIMINC,CICE_TIMINC
            ELSE
               READ(WSMSG1,*) WTIMINC,RSTIMINC,CICE_TIMINC
            ENDIF
         ENDIF
      ENDIF
C
C     jgf49.0804 Added NWS29 VortexOWI -- don't think we need to
C     parse out the WTIMINC or other parameters.
      ! jgf: Added NWS30.
      IF (((ABS(NWS)).EQ.29).or.(abs(NWS).eq.30)) THEN
         READ(15,80) WSMSG1
      ENDIF
C
C saeed added  TODO: I should remove this partialy: cap should provide the time of data exchanges
      IF (ABS(NWS).EQ.17) THEN
         READ(15,80) WSMSG1
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
             READ(WSMSG1,*) WTIMINC
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
             READ(WSMSG1,*) WTIMINC , RSTIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
             READ(WSMSG1,*) WTIMINC , CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
             READ(WSMSG1,*) WTIMINC , RSTIMINC, CICE_TIMINC
         END IF
      ENDIF
C
C...  xyc52.30 added NWS=7 to follow NWS=6
      IF((ABS(NWS).EQ.6) .OR. (ABS(NWS).EQ.7)) THEN
         READ(15,80) WSMSG1
C..... tcm v49.64.01 additions for ICE
         IF((NCICE.EQ.0).AND.(NRS.EQ.0)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC
         ELSEIF ((NCICE.EQ.0).AND.(NRS.GE.1)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC
         ELSEIF((NCICE.GE.1).AND.(NRS.EQ.0)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,CICE_TIMINC
         ELSEIF ((NCICE.GE.1).AND.(NRS.GE.1)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
         ENDIF
      ENDIF
C
      READ(15,80) RNDAYMSG
      READ(RNDAYMSG,*) RNDAY
C
      READ(15,80) DRAMPMSG
      READ(DRAMPMSG,*) DRAMP
C
      READ(15,80) COEFMSG
      READ(COEFMSG,*) A00,B00,C00
C
      READ(15,80) H0MSG
      IF (NOLIFA.NE.2) THEN
         READ(H0MSG,*) H0
      ELSE
         READ(H0MSG,*) H0,NODEDRYMIN,NODEWETRMP,VELMIN
      ENDIF
C
      READ(15,80) SLMSG
      READ(SLMSG,*) SLAM0,SFEA0
C
      SL0=DEG2RAD*SLAM0
      SF0=DEG2RAD*SFEA0

      DO I = 1,NNODG
         SLAM(I) = X(I)
         SFEA(I) = Y(I)
      ENDDO

!      WRITE(6,*)"SLAM0,SFEA0", SLAM0,SFEA0
C
C--If ICS = 2,20-24 then apply cylindrical projection
C
      IF (ICS.NE.1) THEN
         DO I = 1,NNODG
            SL1(I) = DEG2RAD*SLAM(I)
            SF1(I) = DEG2RAD*SFEA(I)
            call CYLINDERMAP(slam(i),sfea(i),sl1(i),sf1(i),sl0,sf0,ICS)
         ENDDO
      ELSE
         ! jgf: Need these to be populated for ICS=1 for NetCDF file
         SL1 = SLAM
         SF1 = SFEA
      ENDIF
C
      READ(15,80) TAUMSG
      IF (NOLIBF.EQ.0) THEN
         READ(TAUMSG,*) TAU
      ELSEIF (NOLIBF.EQ.1) THEN
         READ(TAUMSG,*) CF
      ELSEIF (NOLIBF.EQ.2) THEN
         READ(TAUMSG,*) CF,HBREAK,FTHETA,FGAMMA
      ENDIF
C
      READ(15,80) ESLMSG
      IF (IM.EQ.10) THEN
         READ(ESLMSG,*) ESLM,ESLC
      ELSE
         READ(ESLMSG,*) ESLM
      ENDIF
C
      READ(15,80) CORIMSG
      READ(CORIMSG,*) CORI
C
      READ(15,80) NTIFMSG
      READ(NTIFMSG,*) NTIF
      IF (NTIF.GT.MNTIF) THEN
         print *, "NTIF = ",NTIF, " exceeds parameter MNTIF = ",MNTIF
         CALL EXIT(1)
      ENDIF
C
      DO I=1,NTIF
        READ(15,80)  TIPOTAG(I)
        READ(15,80)  TPKMSG(I)
        TIPOTAG(I)=ADJUSTL(TIPOTAG(I))
        READ(TPKMSG(I),*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
      ENDDO
C
      IF (((NTIP.EQ.0).AND.(NTIF.NE.0)).OR.((NTIP.NE.0).AND.
     &    (NTIF.EQ.0))) THEN
        print *, "NTIF and NTIP are not consistent"
        IF(NFOVER.EQ.1) THEN
          print *, "Program will reset NTIP = 0 and continue"
          NTIP = 0
        ELSE
          CALL EXIT(1)
        ENDIF
      ENDIF
C
      READ(15,80) NBFRMSG
      READ(NBFRMSG,*) NBFR
      IF (NBFR.GT.MNBFR) THEN
        print *, "NBFR = ",NBFR, " exceeds parameter MNBFR = ",MNBFR
        CALL EXIT(1)
      ENDIF
C
      DO I=1,NBFR
        READ(15,80) BOUNTAG(I)
        READ(15,80) AMIGMSG(I)
        READ(AMIGMSG(I),*) AMIG(I),FF(I),FACE(I)
      ENDDO
C
      DO I=1,NBFR
        READ(15,80) ALPHA1(I)
        DO J=1,NETA
          READ(15,80) EMOMSG(I,J)
          READ(EMOMSG(I,J),*) EMO(I,J),EFA(I,J)
        ENDDO
      ENDDO
C
      READ(15,80) ANGMSG
      READ(ANGMSG,*) ANGINN
      !
      ! jgf51.52.23: There was code here to detemine if there are
      ! external flux boundaries in the mesh file, but this info
      ! was already determined in read14() above, so this redundant
      ! code was removed.
      ! Also, the number of external flux boundary nodes was already
      ! determined in read14() above, so the loop through the LBCODE
      ! array was not necessary. Also, the loop that was removed was
      ! bounded by NVEL, which is larger than the number of values
      ! in the LBCODE array, causing uninitialized values to be
      ! examined in the loop.
C
Ckmd Added NFFR=-1 so the river information does not have to start at
C    time=0.
      IF (NFLUXF.EQ.1) THEN
         READ(15,80) NFFRMSG
         READ(NFFRMSG,*) NFFR
         IF (NFFR.GT.MNFFR) THEN
            print *, "NFFR = ",NFFR, " exceeds parameter MNFFR = ",MNFFR
            CALL EXIT(1)
         ENDIF
C
         IF ((NFFR.NE.0).AND.(NFFR.NE.-1)) THEN
            DO I=1,NFFR
               READ(15,80) FBOUNTAG(I)
               READ(15,80) FREQMSG(I)
               READ(FREQMSG(I),*) FAMIG(I),FFF(I),FFACE(I)
            ENDDO
            DO I=1,NFFR
               READ(15,80) ALPHA2(I)
               ! jgf51.52.23: Changed the loop extent from NVEL to
               ! numLBCodeValues
               DO J=1,numLBCodeValues
                  IF ((LBCODE(J).EQ.2).OR.(LBCODE(J).EQ.12)
     &                 .OR.(LBCODE(J).EQ.32)
     &                 .OR.(LBCODE(J).EQ.22).OR.(LBCODE(J).EQ.52)) THEN
                     READ(15,80) QNMSG(I,J)
                     READ(QNMSG(I,J),*) QNAM(I,J),QNPH(I,J)
c     dbug             print *, "disc lbnode index = ",J
                  ENDIF
               ENDDO
            ENDDO

Cvjp 5/1/99  added to help localize the flow boundary nodes
C  NBV    =  global node number of all boundary nodes
C  NFLBN  =  number of flow boundary nodes
C  FLBN   =  global node number of flow boundary nodes
C  FLBNX  =  index of flow boundary nodes as per NBV
C
            JG = 1
            DO K = 1,NBOU
               DO I=1, NVELL(K)
                  INDX = NBVV(K,I)
                  NBV(JG) = INDX
                  JG = JG + 1
               ENDDO
            ENDDO
C
            NFLBN = 0
            ! jgf51.52.23: Changed the loop extent from NVEL to
            ! numLBCodeValues.
            DO J=1,numLBCodeValues
               IF ((LBCODE(J).EQ.2).OR.(LBCODE(J).EQ.12)
     &              .OR.(LBCODE(J).EQ.32)
     &              .OR.(LBCODE(J).EQ.22).OR.(LBCODE(J).EQ.52)) THEN
                  NFLBN = NFLBN+1
                  FLBN(NFLBN) = NBV(J)
                  FLBNX(NFLBN) = J
               ENDIF
            ENDDO
         ELSE
            APERIODIC_FLOW_BC = .TRUE. !jgf45.09 need to break up a fort.20
         ENDIF
      ENDIF
C      PRINT *, "APERIODIC_FLOW_BC = ", APERIODIC_FLOW_BC


Cobell...READING TIME VARYING BOUNDARY CONDITIONS
C........CHECK FOR WEIR FLUX BOUNDARIES
      NFLUXB = 0
      NFLUXIB = 0
      NFLUXIBP = 0
      DO I = 1,NVEL
        SELECT CASE(LBCODE(I))
            CASE(3,13,23)
                NFLUXB = 1
            CASE(4,24,64)
                NFLUXIB = 1
            CASE(5,25)
                NFLUXIBP = 1
        END SELECT
      ENDDO
      IF((NFLUXIB.EQ.1) .OR. (NFLUXIBP.EQ.1) .OR. (NFLUXB.EQ.1) )THEN
        IF(USE_TVW)CALL PARSE_TIME_VARYING_WEIR_INFO()
      ENDIF

C
C--Read Elevation Recording Stations
C...
      READ(15,80) STAEMSG
      READ(STAEMSG,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      IF (ABS(NOUTE).GT.5) THEN
        print *, "The value of NOUTE is not allowed"
        CALL EXIT(1)
      ENDIF
C
      READ(15,80) NSTAEMSG
      READ(NSTAEMSG,*) NSTAE
      write(6,*) 'nstae=',nstae ! jgfdebug
      !tcm v51.20.03 -- addition for external elevation station file
      if (NSTAE < 0 ) then
         WRITE(*,*) "External File Used for Elevation Station Locations"
         USE_ELEV_STAT_FILE = .TRUE.
         IOS_STATIONS = 0
         NSTAE = abs(NSTAE)  !reset to positive
         NSTAE2 = 0
         OPEN(unit=151,file='elev_stat.151',status='old',err=7690,iostat=ios_stations)
         READ(151,*) NSTAE2
         IF (ABS(NSTAE2).NE. ABS(NSTAE)) THEN
            NSTAE = ABS(NSTAE2)  !RESET THE VALUE TO WHAT'S IN THE FILE
         ENDIF
7690     IF (IOS_STATIONS .NE. 0) THEN
            WRITE(*,*) "Error in Reading Elevation Station File: elev_stat.151"
            CALL EXIT(1)  ! there is a stop here
         ENDIF
      ELSE
         WRITE(*,*) "Elevation Station Locations contained in fort.15"
      ENDIF
      IF (ABS(NSTAE).GT.MNSTAE) THEN
         print *, "NSTAE = ",abs(NSTAE), " exceeds parameter MNSTAE = ",
     &              MNSTAE
          CALL EXIT(1)
      ENDIF
#ifdef ADCNETCDF
C
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
!      IF(.NOT. ALLOCATED(STATNUMB))ALLOCATE(STATNUMB(abs(NSTAE),SNUMLEN))
      IF(.NOT. ALLOCATED(STATNAME))ALLOCATE(STATNAME(abs(NSTAE)))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
#endif
      DO I=1,ABS(NSTAE)
         IF (USE_ELEV_STAT_FILE ) then  !read from external file
            READ(151,'(A132)') STAELOC(I)
         ELSE
            READ(15,'(A132)') STAELOC(I)
         ENDIF
            IF(ICS.EQ.1) THEN

C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
               IF((abs(NOUTE).NE.3).AND.(abs(NOUTE).NE.5))THEN
                  READ(STAELOC(I),*) XEL(I),YEL(I)
               ELSE
#ifdef ADCNETCDF
!                READ(15,'(A132)') STATLINE
!                call parse(STATLINE, LVAR)
                 call parse(STAELOC(I), LVAR)
                 XEL(I)=a2f(LVAR(1))
                 YEL(I)=a2f(LVAR(2))
                 STATNAME(I)=LVAR(3)
C                WRITE(6,9911) XEL(I), YEL(I), STATNAME(I)
#endif
               ENDIF

            ELSE
               IF((abs(NOUTE).NE.3).AND.(abs(NOUTE).NE.5))THEN

                  READ(STAELOC(I),*) SLEL(I),SFEL(I)
               ELSE
#ifdef ADCNETCDF
!                 READ(15,'(A132)') STATLINE
!                 call parse(STATLINE, LVAR)
                  call parse(STAELOC(I), LVAR)
                  SLEL(I)=a2f(LVAR(1))
                  SFEL(I)=a2f(LVAR(2))
                  STATNAME(I)=LVAR(3)
!                 WRITE(6,*)"parsed stations"
!                 WRITE(6,9911) SLEL(I),SFEL(I), STATNAME(I,SNAMLEN)
#endif
               ENDIF


            SLEL(I)=DEG2RAD*SLEL(I)
            SFEL(I)=DEG2RAD*SFEL(I)
            CALL CYLINDERMAP(XEL(I),YEL(I),SLEL(i),SFEL(i),SL0,SF0,ICS)
         ENDIF
         ENDDO !loop over elevation stations
         if (USE_ELEV_STAT_FILE) CLOSE(151) !close external file if one was used
 9911    FORMAT(F12.3,2X,F12.3, 6X, A50)

C
C--For Each Elevation Station:
C  Find the Global Index of the element it lies in.
C
c...tcm v49.48.01 replace with new call to kdtsearch
!      CALL CoordToEle(XEL,YEL,NNSEG,NSTAE,
!     &     'Elevation recording station   ')
      call setup_kdt_search()
      call kdtsearch(XEL,YEL,NNSEG,abs(NSTAE),
     &     'Elevation recording station   ')
C
C--Read Velocity Recording Stations
C
      READ(15,80) STAVMSG
      READ(STAVMSG,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      IF (ABS(NOUTV).GT.5) THEN
        print *, "Value for NOUTV is not allowable"
        CALL EXIT(1)
      ENDIF
C
      READ(15,80) NSTAVMSG
      READ(NSTAVMSG,*) NSTAV
      !tcm v51.20.03 -- addition for external velocity station location file vel_stat.151
      IF (NSTAV < 0 ) then
         WRITE(*,*) "External File Used for Velocity Station Locations"
         USE_VEL_STAT_FILE = .TRUE.
         IOS_STATIONS = 0
         NSTAV = ABS(NSTAV) ! reset to positive value
         NSTAV2 = 0
         OPEN(unit=151,file='vel_stat.151',status='old',err=7691,iostat=ios_stations)
         READ(151,*) NSTAV2
         IF (ABS(NSTAV2).NE. ABS(NSTAV)) THEN
            NSTAV = ABS(NSTAV2)  !RESET THE VALUE TO WHAT'S IN THE FILE
         ENDIF
7691     IF (IOS_STATIONS .NE. 0) THEN
            WRITE(*,*) "Error in Reading Velocity Station File: vel_stat.151"
            CALL EXIT(1)  ! there is a stop here
         ENDIF
      ELSE
         WRITE(*,*) "Velocity Station Locations Contained in fort.15"
      ENDIF

      IF (ABS(NSTAV).GT.MNSTAV) THEN
        print *, "NSTAV = ",abs(NSTAV), " exceeds parameter MNSTAV = ",
     &             MNSTAV
        CALL EXIT(1)
      ENDIF
#ifdef ADCNETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
!      IF(.NOT. ALLOCATED(STATNUMBV))ALLOCATE(STATNUMBV(NSTAV,SNUMLEN))
      IF(.NOT. ALLOCATED(STATNAMEV))ALLOCATE(STATNAMEV(abs(NSTAV)))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
C
#endif

      DO I=1,ABS(NSTAV)
         if (USE_VEL_STAT_FILE ) then !read from external file
            READ(151,'(A132)') STAVLOC(I)
         else
            READ(15,'(A132)') STAVLOC(I)
         ENDIF
         IF (ICS.eq.1) then
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
            IF((abs(NOUTV).NE.3).AND.(abs(NOUTV).NE.5)) THEN
               READ(STAVLOC(I),*) XEV(I),YEV(I)
            ELSE
#ifdef ADCNETCDF
               call parse(STAVLOC(I), LVAR)
               XEV(I)=a2f(LVAR(1))
               YEV(I)=a2f(LVAR(2))
               STATNAMEV(I)=LVAR(3)
C              WRITE(6,9911) XEV(I), YEV(I), STATNAMEV(I)
#endif
            ENDIF

         ELSE

            IF((abs(NOUTV).NE.3).AND.(abs(NOUTV).NE.5)) THEN
               READ(STAVLOC(I),*) SLEV(I),SFEV(I)

            ELSE
#ifdef ADCNETCDF
               call parse(STAVLOC(I), LVAR)
               SLEV(I)=a2f(LVAR(1))
               SFEV(I)=a2f(LVAR(2))
               STATNAMEV(I)=LVAR(3)
!              WRITE(6,*)"parsed stations"
!              WRITE(6,9911) SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
!              WRITE(6,*)"STATION INFO - AFTER *************",
!     &           SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
#endif
            ENDIF
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================


           SLEV(I)=DEG2RAD*SLEV(I)
           SFEV(I)=DEG2RAD*SFEV(I)
           CALL CYLINDERMAP(XEV(I),YEV(I),SLEV(i),SFEV(i),SL0,SF0,ICS)
         ENDIF  !test for ICS
      ENDDO !loop over velocity stations
      IF (USE_VEL_STAT_FILE) CLOSE(151) !close external file if one was used
C
C--For Each Velocity Station:
C  Find the Global Index of the element it lies in.
C
c...tcm v49.48.01 replace with new call to kdtsearch
!      CALL CoordToEle(XEV,YEV,NNSVG,NSTAV,
!     &     'Velocity recording station    ')
      call kdtsearch(XEV,YEV,NNSVG,abs(NSTAV),
     &     'Velocity recording station    ')

C
C     If Passive Transport is indicated, then read Concentration Station Info
C
      IF (C2D_PTrans.or.C3D_PTrans) THEN
C
        READ(15,80) STACMSG
        READ(STACMSG,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
        IF (ABS(NOUTC).GT.2) THEN
          print *, "Value of NOUTC is not allowable"
          CALL EXIT(1)
        ENDIF
C
        READ(15,80) NSTACMSG
        READ(NSTACMSG,*) NSTAC
       !tcm v51.20.03 -- addition for external Concentration station file
        if (NSTAC < 0 ) then
           WRITE(*,*) "External File Used for Concentration Station Locations"
           USE_CONC_STAT_FILE = .TRUE.
           IOS_STATIONS = 0
           NSTAC = ABS(NSTAC) !RESET TO POSITIVE
           NSTAC2 = 0 !INITIALIZE
           OPEN(unit=151,file='conc_stat.151',status='old',err=7692,iostat=ios_stations)
           READ(151,*) NSTAC2
           IF (ABS(NSTAC2).NE. ABS(NSTAC)) THEN
              NSTAC = ABS(NSTAC2)  !RESET THE VALUE TO WHAT'S IN THE FILE
           ENDIF
7692       IF (IOS_STATIONS .NE. 0) THEN
              WRITE(*,*) "Error in Reading Concentration Station File: conc_stat.151"
              CALL EXIT(1)  ! THERE IS A STOP HERE
           ENDIF
        ELSE
           WRITE(*,*) "Concentration Station Locations Contained in fort.15"
        ENDIF

        IF (ABS(NSTAC).GT.MNSTAC) THEN
          print *, "NSTAC = ",abs(NSTAC), " exceeds parameter MNSTAC = ",
     &              MNSTAC
          CALL EXIT(1)
        ENDIF
C
        DO I=1,ABS(NSTAC)
           IF (USE_CONC_STAT_FILE) THEN !READ FROM EXTERNAL FILE
              READ(151,'(A132)') STACLOC(I)
           ELSE
              READ(15,'(A132)') STACLOC(I)
           ENDIF
           IF (ICS.EQ.1) THEN
             READ(STACLOC(I),*) XEC(I),YEC(I)
           ELSE
             READ(STACLOC(I),*) SLEC(I),SFEC(I)
             SLEC(I)=DEG2RAD*SLEC(I)
             SFEC(I)=DEG2RAD*SFEC(I)
             CALL CYLINDERMAP(XEC(I),YEC(I),SLEC(i),SFEC(i),SL0,SF0,ICS)
           ENDIF !test for ICS
        ENDDO  !loop over stations
        IF (USE_CONC_STAT_FILE) CLOSE(151) !CLOSE EXTERNAL FILE IF ONE WAS USED
C
C--For Each Concentration Recording Station:
C  Find the Global Index of the element it lies in.
C
c...tcm v49.48.01 replace with new call to kdtsearch
!         CALL CoordToEle(XEC,YEC,NNSCG,NSTAC,
!     &        'Concentration station         ')
         call kdtsearch(XEC,YEC,NNSCG,abs(NSTAC),
     &        'Concentration station         ')
C
      ENDIF
C
C--If NWS <> 0 , then read Meteorlogical Station Info
C
      NOUTM = 0
      IF (NWS.NE.0) THEN
C
         READ(15,80) STAMMSG
         READ(STAMMSG,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         IF (ABS(NOUTM).GT.5) THEN
            print *, "Value of NOUTM is not allowable"
            CALL EXIT(1)
         ENDIF
C
         READ(15,80) NSTAMMSG
         READ(NSTAMMSG,*) NSTAM
         !tcm v51.20.03 -- addition for external met station file
         IF (NSTAM < 0 ) THEN
            WRITE(*,*) "External File Used for MET Station Locations"
            USE_MET_STAT_FILE = .TRUE.
            IOS_STATIONS = 0
            NSTAM = ABS(NSTAM) !RESET TO POSITIVE VALUE
            NSTAM2 = 0 !INITIALIZE
            OPEN(unit=151,file='met_stat.151',status='old',
     &           err=7693,iostat=ios_stations)
            READ(151,*) NSTAM2
            IF (ABS(NSTAM2).NE. ABS(NSTAM)) THEN
               NSTAM = ABS(NSTAM2)  !RESET THE VALUE TO WHAT'S IN THE FILE
            ENDIF
7693        IF (IOS_STATIONS .NE. 0) THEN
               WRITE(*,*) "Error in Reading MET Station File: met_stat.151"
               CALL EXIT(1)  ! THERE IS A STOP HERE
            ENDIF
         ELSE
            WRITE(*,*) "MET Station Locations Contained in fort.15"
         ENDIF

         IF (ABS(NSTAM).GT.MNSTAM) THEN
            print *, "NSTAM = ",abs(NSTAM), " exceeds parameter MNSTAM = ",
     &           MNSTAM
            CALL EXIT(1)
         ENDIF
#ifdef ADCNETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C

!         IF(.NOT.ALLOCATED(STATNUMBM))ALLOCATE(STATNUMBM(NSTAM,SNUMLEN))
         IF(.NOT.ALLOCATED(STATNAMEM))ALLOCATE(STATNAMEM(abs(NSTAM)))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
C
#endif
         DO I=1,ABS(NSTAM)
            IF (USE_MET_STAT_FILE) THEN !READ FROM EXTERNAL FILE
               READ(151,'(A132)') STAMLOC(I)
            ELSE
               READ(15,'(A132)') STAMLOC(I)
            ENDIF
C
            IF (ICS.EQ.1) THEN

C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
               IF((ABS(NOUTM).NE.3).AND.(ABS(NOUTM).NE.5)) THEN
                  READ(STAMLOC(I),*) XEM(I),YEM(I)
               ELSE
#ifdef ADCNETCDF
                  call parse(STAMLOC(I), LVAR)
                  XEM(I)=a2f(LVAR(1))
                  YEM(I)=a2f(LVAR(2))
                  STATNAMEM(I)=LVAR(3)
C                 WRITE(6,9911) XEM(I), YEM(I), STATNAMEM(I)
#endif
               ENDIF
            ELSE
               IF((ABS(NOUTM).NE.3).AND.(ABS(NOUTM).NE.5)) THEN

!                 READ(15,80) STAMLOC(I)
                  READ(STAMLOC(I),*) SLEM(I),SFEM(I)
               ELSE
#ifdef ADCNETCDF
                  call parse(STAMLOC(I), LVAR)
                  SLEM(I)=a2f(LVAR(1))
                  SFEM(I)=a2f(LVAR(2))
                  STATNAMEM(I)=LVAR(3)
!                 WRITE(6,9911) SLEM(I), SFEM(I), STATNAMEM(I,SNAMLEN)
!                 WRITE(6,*)"STATION INFO STAMLOC",
!     &           SLEM(I),SFEM(I),STATNAMEM(I,SNAMLEN)
!    &            SLEM(I),SFEM(I),STATNUMBM(I),STATNAMEM(I)
#endif
               ENDIF
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================

               SLEM(I)=DEG2RAD*SLEM(I)
               SFEM(I)=DEG2RAD*SFEM(I)
               CALL CYLINDERMAP(XEM(I),YEM(I),SLEM(i),SFEM(i),SL0,SF0,ICS)
            ENDIF  !test for ICS
         ENDDO  !loop over stations
         IF (USE_MET_STAT_FILE) CLOSE(151) !CLOSE EXTERNAL FILE IF ONE WAS USED

C
C--For Each Meterological Recording Station:
C  Find the Global Index of the element it lies in.
C
c...tcm v49.48.01 replace with new call to kdtsearch
!         CALL CoordToEle(XEM,YEM,NNSMG,NSTAM,
!     &        'Meteorological station        ')
         call kdtsearch(XEM,YEM,NNSMG,abs(NSTAM),
     &         'Meteorological station        ')
C
      ENDIF
C
C
C--Read Global Elevation Data Output
C
      READ(15,80) OUTGEMSG
      READ(OUTGEMSG,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      IF (ABS(NOUTGE).GT.7) THEN
         print *, "NOUTGE does not have an allowable value"
         CALL EXIT(1)
      ENDIF
C
C--Read Global Velocity Data Output
C
      READ(15,80) OUTGVMSG
      READ(OUTGVMSG,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      IF (ABS(NOUTGV).GT.7) THEN
        print *, "NOUTGV does not have an allowable value"
        CALL EXIT(1)
      ENDIF
C
C     If Passive Transport is indicated, read Global Concentration Data Output
C
      IF (C2D_PTrans.or.C3D_PTrans) THEN
        READ(15,80) OUTGCMSG
        READ(OUTGCMSG,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
        IF (ABS(NOUTGC).GT.2) THEN
          print *, "NOUTGC does not have an allowable value"
          CALL EXIT(1)
        ENDIF
      ENDIF
C
      IF (NWS.NE.0) THEN
        READ(15,80) OUTGWMSG
        READ(OUTGWMSG,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
Casey 090302: Changed the following line for ABS(NOUTGW).EQ.4
        IF (ABS(NOUTGW).GT.7) THEN
          print *, "NOUTGW does not have an allowable value"
          CALL EXIT(1)
        ENDIF
      ENDIF
C
C--Read Harmonic Analysis Data
C
      READ(15,80) HARFRMSG
      READ(HARFRMSG,*) NHARFR
      MNHARF = NHARFR
      IF (NHARFR.EQ.0) MNHARF = 1
      CALL ALLOC_HA()
      CALL ALLOC_MAIN14()
      ALLOCATE(HAFREMSG(MNHARF))
C
      DO I=1,NHARFR
          READ(15,'(A10)') NAMEFR(I)
          READ(15,80) HAFREMSG(I)
          READ(HAFREMSG(I),*) HAFREQ(I),HAFF(I),HAFACE(I)
      ENDDO
C
      READ(15,80) HARPARMSG
      READ(HARPARMSG,*) THAS,THAF,NHAINC,FMV
      READ(15,80) OUTHARMSG
      READ(OUTHARMSG,*) NHASE,NHASV,NHAGE,NHAGV
      IF ((NHASE.LT.0).OR.(NHASE.GT.5)) THEN
        print *, "NHASE does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHASE=0 "
          NHASE = 0
        ELSE
          CALL EXIT(1)
        ENDIF
      ENDIF
C
      IF ((NHASV.LT.0).OR.(NHASV.GT.5)) THEN
         print *, "NHASV does not have an allowable value"
         IF (NFOVER.EQ.1) THEN
           print *, "Program will override an reset NHASV=0 "
           NHASV = 0
         ELSE
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      IF ((NHAGE.LT.0).OR.(NHAGE.GT.5)) THEN
        print *, "NHAGE does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHAGE=0 "
          NHAGE = 0
        ELSE
          CALL EXIT(1)
        ENDIF
      ENDIF
C
      IF ((NHAGV.LT.0).OR.(NHAGV.GT.5)) THEN
        print *, "NHAGV does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHAGV=0 "
          NHAGV = 0
        ELSE
          CALL EXIT(1)
        ENDIF
      ENDIF
      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1
      IF ((FMV.GT.0.).AND.(NHARFR.GT.0).AND.(C2DDI)) CHARMV = .TRUE.
C
C--Read Hot Start Data
C
C     jgf45.07 added undocumented option to allow ADCIRC to stop after writing
C     hot start file. This is used to test hot starting capabilities.
      READ(15,80) HSTARMSG
      READ(HSTARMSG,*) NHSTAR,NHSINC
C     jgf49.41: Added new options for hotstart file format indicator: 2, 3,
C     367, 368. Not sure we need to check the validity of NHSTAR values in
C     adcprep though.      !wjp allows for single output hotstart (negative value that is not -1) 
      SELECT CASE(abs(NHSTAR))  !tcm v51.26 mods for binary time-stamped hotstart file generation NHSTAR=-1
      CASE(0,1,2,3,5,67,68,367,368,567,568)
         ! do nothing, these are allowable values
      CASE DEFAULT
        print *, "WARNING: NHSTAR=",NHSTAR," is not valid."
        IF (NFOVER.EQ.1) THEN
           print *, "WARNING: Program will override and reset NHSTAR=0."
           NHSTAR = 0
        ELSE
           print *,"ERROR: Halting adcprep due to invalid NHSTAR value."
           CALL EXIT(1)
        ENDIF
      END SELECT
C
C--Read Solver Data
C
      READ(15,80) SOLVMSG
      READ(SOLVMSG,*) ITITER,ISLDIA,CONVCR,ITMAX
      ! jgf49.82: Removed value check on ITITER ... this parameter
      ! is deprecated and its value is now arbitrary.
C
C
C--Read in 3D data
C
      IF(C3DVS) THEN
        CALL READ15_3DVS()
c     ELSEIF(C3DDSS) THEN
c       CALL READ15_3DDSS()
      ENDIF

      IF((useNetCDF.eqv..true.).and.(NETCDF_AVAIL.eqv..false.)) THEN
         CALL allMessage(ERROR,
     &      "NetCDF was selected in the fort.15 file.")
         CALL allMessage(ERROR,
     &      "NetCDF was not compiled into this executable program.")
         CALL EXIT(1)
      ENDIF
#ifdef ADCNETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION for NETCDF I/O - 6/20/07
C
!     --------------------------------------------------------------------
!     IF output is in netCDF format read global attributes for netCDF file
!     --------------------------------------------------------------------

      IF(useNetCDF.eqv..true.) THEN
         IF (NOUTGC.eq.3.or.                               ! fort.93
     &       NOUTC.eq.3) then                                 ! fort.91
            CALL allMessage(ERROR,
     &         "NetCDF output not available for concentration output.")
            CALL EXIT(1)
         ENDIF
C
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) title
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) institution
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) source
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) history
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) references
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) comments
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) host
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) convention
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) contact
         READ(15,'(A80)',ERR=999,END=999,IOSTAT=ios) base_date

999      IF (ios.ne.0) THEN
            CALL allMessage(ERROR,
     &         "Could not read metadata for NetCDF at end of fort.15.")
            CALL EXIT(1)
         ENDIF

      ENDIF
#endif
C
C--Close Global Run Info file
C
      CLOSE(15)
      call memory_status()
C
#ifdef READ_GLOBAL_TRACE
      call screenMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C
      RETURN
 80   FORMAT(A80)
      END SUBROUTINE READ15


C----------------------------------------------------------------------
C     Reads the 3DVS portion of the global ADCIRC Input Data File for
C     data decomposition program ADCPP.  This version compatible with
C     ADCIRC_v41.11a
C     tjc  6/24/02
C---------------------------------------------------------------------
      SUBROUTINE READ15_3DVS()
      USE PRE_GLOBAL
      USE GLOBAL, ONLY : parse, a2f, SNAMLEN, deg2rad, day2sec
      USE GLOBAL_3DVS, ONLY : STATNAMED, STATNAMEV3D, STATNAMET, SIGMA,
     &    FEGRIDS
      USE MESH, only : CYLINDERMAP
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,K,M,N                                        ! loop counters
      CHARACTER(50) LVAR(3)
      REAL(8) :: STATIME
      REAL(sz) :: HH1 ! domain averaged depth
      INTEGER :: NH   ! horizontal node counter
C
C     jgf49.43.19: Moved memory allocation for number of 3D stations
C     on each subdomain to here and made it unconditional.
C     The NNSTA3D*P arrays will still be needed later b/c they are
C     written to the fort.18 message passing file.
      ALLOCATE(NNSTA3DDP(MNPROC))
      nbytes = 4*mnproc
      call memory_alloc(nbytes)
      ALLOCATE(NNSTA3DVP(MNPROC))
      nbytes = 4*mnproc
      call memory_alloc(nbytes)
      ALLOCATE(NNSTA3DTP(MNPROC))
      nbytes = 4*mnproc
      call memory_alloc(nbytes)
      ! initialize
      NNSTA3DDP = 0
      NNSTA3DVP = 0
      NNSTA3DTP = 0
C
C
 350  FORMAT(//,2X,'***** INVALID INPUT IN THE PRIMARY VERTICAL INPUT',
     &          ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

C     jgf45.10 removed IDIAG
C... SPECIFY WHETHER A BAROTROPIC OR BAROCLINIC RUN
      READ(15,80) IDENMSG
      READ(IDENMSG,*) IDEN
      IF((IDEN.GT.4).OR.(IDEN.LT.-4)) THEN
         WRITE(*,*) "ERROR: IDEN=",IDEN
         WRITE(*,423)
         WRITE(*,350)
 423     FORMAT(/,2X,' IDEN must be an integer between -4 and 4.')
         CALL EXIT(1)
      ENDIF

C     jgf45.12 Set flag to read in unit 11 file (initial density) if
C     necessary.
      IF ( IDEN .ne. 0 ) THEN
         CBaroclinic = .true.
      ENDIF
C     jgf45.12 Set flag to read in additional parameters from unit 15
C     file that will be used when solving for transport of salinity,
C     temperature, etc in prognostic baroclinic simulation.
      IF ( IDEN .gt. 0 ) THEN
         C3D_BTrans = .true.
      ENDIF

C... READ IN THE TYPE OF BOTTOM BOUNDARY CONDITION AND THE SLIP COEFFICIENTS

      READ(15,80) SLIPMSG
      READ(SLIPMSG,*) ISLIP,KP
      IF((ISLIP.LT.0).OR.(ISLIP.GT.3)) THEN
         WRITE(6,350)
         WRITE(6,360)
 360     FORMAT(/,2X,'    THE SLIP CODE MUST = 0,1,2,OR 3.')
         CALL EXIT(1)
      ENDIF

C... READ IN THE SURFACE AND BOTTOM ROUGHNESSES

      READ(15,80) Z0MSG
      READ(Z0MSG,*) Z0S, Z0B

C... READ IN THE TIME STEPPING COEFFICIENTS

      READ(15,80) ALPMSG
      READ(ALPMSG,*) ALP1,ALP2,ALP3

C... READ IN IGC & NFEN: F.E. GRID CODE & # NODES IN F.E. GRID

      READ(15,80) FEMSG
      READ(FEMSG,*) IGC,NFEN

C     jgf45.12 Add code to read in nondimensional thicknesses of
C     vertical layers.
      ALLOCATE ( Sigma(NFEN) )
      nbytes = 8*nfen
      call memory_alloc(nbytes)
      IF(IGC.EQ.0) THEN
         DO N=1,NFEN
            READ(15,*) Sigma(N)
         ENDDO
      ELSE  ! jgf50.60.10 Populate sigma() for netcdf output files
         HH1=0.d0
         DO NH=1,MNP
            HH1=HH1+DP(NH)
         ENDDO
         HH1=HH1/MNP                        !domain averaged depth
         CALL FEGRIDS(IGC,HH1)
      ENDIF



C... SPECIFY TYPE OF EDDY VISCOSITY PROFILE
      READ(15,80) EVCMSG
      READ(EVCMSG,*) IEVC,EVMIN,EVCON
      IF((IEVC.NE.0 ).and.(IEVC.NE.1 ).and.
     &   (IEVC.NE.10).and.(IEVC.NE.11).and.
     &   (IEVC.NE.20).and.(IEVC.NE.21).and.
     &   (IEVC.NE.22).and.(IEVC.NE.23).and.
     &   (IEVC.NE.30).and.(IEVC.NE.31).and.(IEVC.NE.32).and.
     &   (IEVC.NE.33).and.
     &   (IEVC.NE.40).and.(IEVC.NE.41).and.(IEVC.NE.42).and.
     &   (IEVC.NE.43).and.
     &   (IEVC.NE.50).and.(IEVC.NE.51)) THEN
         WRITE(*,350)
         WRITE(*,411)
 411     FORMAT(/,2X,'    IEVC MUST BE 0,1,10,11,20,21,22,23,',
     &              '30,31,32,33,40,41,42,43,50,51')
         CALL EXIT(1)
      ENDIF
      IF((IEVC.EQ.50).OR.(IEVC.EQ.51)) THEN
         READ(15,80) THETAMSG
         READ(THETAMSG,*) THETA1,THETA2
      ENDIF
C     jgf45.12 Add code to read in vertical eddy viscosity profile.
      IF(IEVC.EQ.0) THEN
         ALLOCATE ( EVTot(NFEN) )
         nbytes = 8*nfen
         call memory_alloc(nbytes)
         DO N=1,NFEN
            READ(15,*) EVTot(N)
         ENDDO
      ENDIF
C     -----------------------------------------------------
C.... Station 3D Density, Temperature, Salinity output
C     jgf45.11 switched from node numbers to coordinates

      READ(15,80) DSDMSG
      READ(DSDMSG,*) I3DSD,TO3DSDS,TO3DSDF,NSPO3DSD
      READ(15,80) NSTA3DDMSG
      READ(NSTA3DDMSG,*) NSTA3DD
      ALLOCATE(STA3DDLOC(NSTA3DD))
      nbytes = 8*nsta3dd
      call memory_alloc(nbytes)
      ALLOCATE(X3DD(NSTA3DD),Y3DD(NSTA3DD))
      nbytes = 16*nsta3dd
      call memory_alloc(nbytes)
      ALLOCATE(SL3DD(NSTA3DD),SF3DD(NSTA3DD))
      nbytes = 16*nsta3dd
      call memory_alloc(nbytes)
      ALLOCATE(NNS3DDG(NSTA3DD))
      nbytes = 4*nsta3dd
      call memory_alloc(nbytes)
      ALLOCATE(IMAP_STA3DD_LG(NSTA3DD,MNPROC))
      nbytes = 4*nsta3dd*mnproc
      call memory_alloc(nbytes)
      IF(.NOT. ALLOCATED(STATNAMED))ALLOCATE(STATNAMED(NSTA3DD))
C     ALLOCATE(STATNAMED(NSTA3DD))
      nbytes = 4*SNAMLEN*NSTA3DD
      call memory_alloc(nbytes)
C  kmd48.33bc changed to -2 and I3DSD not equal to 0
      IF(ABS(I3DSD).GT.5) THEN
         WRITE(*,350)
         WRITE(*,511)
 511     FORMAT(/,2X,'    I3DSD MUST BE -5,-3,-2,-1,0,1,2, 3, or 5')
         CALL EXIT(1)
      ENDIF
      IF(NSTA3DD.NE.0) THEN  ! kmd : changed to match 2D option
         DO I=1,NSTA3DD
            READ(15,80) STA3DDLOC(I)
            call parse(STA3DDLOC(I), LVAR)
            STATNAMED(I)=LVAR(3)
            IF (ICS.EQ.1) THEN
               X3DD(I)=a2f(LVAR(1))
               Y3DD(I)=a2f(LVAR(2))
            ELSE
               SL3DD(I)=a2f(LVAR(1))
               SF3DD(I)=a2f(LVAR(2))
               SL3DD(I)=DEG2RAD*SL3DD(I)
               SF3DD(I)=DEG2RAD*SF3DD(I)
               call CYLINDERMAP(x3dd(i),y3dd(i),sl3dd(i),sf3dd(i),
     &                          sl0,sf0,ICS)
            ENDIF
         ENDDO
C
C     jgf45.11 For each 3D density station, find the full-domain index of
C     the corresponding element.
c...tcm v49.48.01 replace with new call to kdtsearch
!         CALL CoordToEle(X3DD,Y3DD,NNS3DDG,NSTA3DD,
!     &        '3D density recording station  ')
          call kdtsearch(X3DD,Y3DD,NNS3DDG,NSTA3DD,
     &           '3D density recording station  ')

      ENDIF

C     -----------------------------------------------------
C.... Station 3D Velocity output
C     jgf45.11 switched from node numbers to coordinates
      READ(15,80) DSVMSG
      READ(DSVMSG,*) I3DSV,TO3DSVS,TO3DSVF,NSPO3DSV
      READ(15,80) NSTA3DVMSG
      READ(NSTA3DVMSG,*) NSTA3DV
      ALLOCATE(IMAP_STA3DV_LG(NSTA3DV,MNPROC))
      nbytes = 4*nsta3dv*mnproc
      call memory_alloc(nbytes)
      ALLOCATE(STA3DVLOC(NSTA3DV))
      nbytes = 8*nsta3dv
      call memory_alloc(nbytes)
      ALLOCATE(X3DV(NSTA3DV),Y3DV(NSTA3DV))
      nbytes = 16*nsta3dv
      call memory_alloc(nbytes)
      ALLOCATE(NNS3DVG(NSTA3DV))
      nbytes = 4*nsta3dv
      call memory_alloc(nbytes)
      IF(.NOT. ALLOCATED(STATNAMEV3D))ALLOCATE(STATNAMEV3D(NSTA3DV))
C      ALLOCATE(STATNAMEV3D(NSTA3DV))
      nbytes = 4*SNAMLEN*NSTA3DV
      call memory_alloc(nbytes)
       ALLOCATE(SL3DV(NSTA3DV),SF3DV(NSTA3DV))
        nbytes = 16*nsta3dv
        call memory_alloc(nbytes)

C  kmd48.33 changed to -2
      IF(ABS(I3DSV).GT.5) THEN
         WRITE(*,350)
         WRITE(*,512)
 512     FORMAT(/,2X,'    I3DSV MUST BE -5,-3,-2,-1,0,1,2,3 or 5')
         CALL EXIT(1)
      ENDIF
      IF(NSTA3DV.NE.0) THEN ! kmd : Changed to match 2D option
         DO I=1,NSTA3DV
            READ(15,80) STA3DVLOC(I)
            call parse(STA3DVLOC(I), LVAR)
            STATNAMEV3D(I)=LVAR(3)
            IF (ICS.EQ.1) THEN
               X3DV(I)=a2f(LVAR(1))
               Y3DV(I)=a2f(LVAR(2))
            ELSE
               SL3DV(I)=a2f(LVAR(1))
               SF3DV(I)=a2f(LVAR(2))
               SL3DV(I)=DEG2RAD*SL3DV(I)
               SF3DV(I)=DEG2RAD*SF3DV(I)
               call CYLINDERMAP(x3dv(i),y3dv(i),sl3dv(i),sf3dv(i),
     &                          sl0,sf0,ics)
            ENDIF
         ENDDO
C     jgf45.11 For each 3D velocity station, find the full-domain index of
C     the corresponding element.
c...tcm v49.48.01 replace with new call to kdtsearch
!         CALL CoordToEle(X3DV,Y3DV,NNS3DVG,NSTA3DV,
!     &        '3D velocity recording station ')
          call kdtsearch(X3DV,Y3DV,NNS3DVG,NSTA3DV,
     &           '3D velocity recording station ')
      ENDIF
C     -----------------------------------------------------
C.... Station 3D turbulence output
C     jgf45.11 switched from node numbers to coordinates
      READ(15,80) DSTMSG
      READ(DSTMSG,*) I3DST,TO3DSTS,TO3DSTF,NSPO3DST
      READ(15,80) NSTA3DTMSG
      READ(NSTA3DTMSG,*) NSTA3DT
      ALLOCATE(IMAP_STA3DT_LG(NSTA3DT,MNPROC))
      nbytes = 4*nsta3dt*mnproc
      call memory_alloc(nbytes)
      ALLOCATE(STA3DTLOC(NSTA3DT))
      nbytes = 4*nsta3dt
      call memory_alloc(nbytes)
      ALLOCATE(X3DT(NSTA3DT),Y3DT(NSTA3DT))
      nbytes = 16*nsta3dt
      call memory_alloc(nbytes)
      IF(.NOT. ALLOCATED(STATNAMET))ALLOCATE(STATNAMET(NSTA3DT))
C      ALLOCATE(STATNAMET(NSTA3DT))
      nbytes = 4*SNAMLEN*NSTA3DT
      call memory_alloc(nbytes)
      ALLOCATE(SL3DT(NSTA3DT),SF3DT(NSTA3DT))
      nbytes = 16*nsta3dv
      call memory_alloc(nbytes)
C  kmd48.33 changed to -2
      IF(ABS(I3DST).GT.5) THEN
         WRITE(*,350)
         WRITE(*,513)
 513     FORMAT(/,2X,'    I3DST MUST BE -5,-3,-2,-1,0,1,2,3 or 5')
         CALL EXIT(1)
      ENDIF
      IF(NSTA3DT.NE.0) THEN ! kmd - Changed to match 2D option
         DO I=1,NSTA3DT
            READ(15,80) STA3DTLOC(I)
            call parse(STA3DTLOC(I), LVAR)
            STATNAMET(I)=LVAR(3)
            IF (ICS.EQ.1) THEN
               X3DT(I)=a2f(LVAR(1))
               Y3DT(I)=a2f(LVAR(2))
            ELSE
               SL3DT(I)=a2f(LVAR(1))
               SF3DT(I)=a2f(LVAR(2))
               SL3DT(I)=DEG2RAD*SL3DT(I)
               SF3DT(I)=DEG2RAD*SF3DT(I)
               call CYLINDERMAP(x3dt(i),y3dt(i),sl3dt(i),sf3dt(i),
     &                          sl0,sf0,ICS)
            ENDIF
         ENDDO

C     jgf45.11 For each 3D turbulence station, find the full-domain index of
C     the corresponding element.

         ALLOCATE(NNS3DTG(NSTA3DT))
         nbytes = 4*nsta3dt
         call memory_alloc(nbytes)
c...tcm v49.48.01 replace with new call to kdtsearch
!         CALL CoordToEle(X3DT,Y3DT,NNS3DTG,NSTA3DT,
!     &        '3D turbulence rec. station    ')
          CALL KDTSEARCH(X3DT,Y3DT,NNS3DTG,NSTA3DT,
     &        '3D turbulence rec. station    ')

      ENDIF
C....  GLOBAL 3D DENSITY, TEMPERATURE, SALINITY OUTPUT

      READ(15,80) DGDMSG
      READ(DGDMSG,*) I3DGD,TO3DGDS,TO3DGDF,NSPO3DGD
C   kmd48.33bc changed to -2
      IF(ABS(I3DGD).GT.5) THEN  !tcm 20110511 should this be .GT.3 instead of .GT.2 ?
         WRITE(*,350)
         WRITE(*,514)
 514     FORMAT(/,2X,'    I3DGD MUST BE -3,-2,-1,0,1,2, or 3')
         CALL EXIT(1)
      ENDIF

C....  GLOBAL 3D VELOCITY OUTPUT

      READ(15,80) DGVMSG
      READ(DGVMSG,*) I3DGV,TO3DGVS,TO3DGVF,NSPO3DGV
C   kmd48.33bc changed to -2
      IF(ABS(I3DGV).GT.5) THEN
         WRITE(*,350)
         WRITE(*,515)
 515     FORMAT(/,2X,'    I3DGV MUST BE -3,-2,-1,0,1,2 or 3')
         CALL EXIT(1)
      ENDIF

C....  GLOBAL 3D TURBULENCE OUTPUT

      READ(15,80) DGTMSG
      READ(DGTMSG,*) I3DGT,TO3DGTS,TO3DGTF,NSPO3DGT
C   kmd48.33bc changed to -2
      IF(ABS(I3DGT).GT.5) THEN
         WRITE(*,350)
         WRITE(16,516)
 516     FORMAT(/,2X,'    I3DGT MUST BE -3,-2,-1,0,1,2, or 3')
         CALL EXIT(1)
      ENDIF

C   kmd48.33bc add in information for 3D boundary condition
      IF (CBAROCLINIC) THEN
        READ(15,80) RESBCFLAGMSG
        READ(RESBCFLAGMSG,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP
        IF (RES_BC_FLAG.NE.0) THEN
           IF (NOPE.GT.0) THEN
             READ(15,80) BCTIMEMSG
             READ(15,80) BCSTATMSG
           END IF
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(15,80) TBCTIMEMSG
           END IF
        END IF
        IF ((RES_BC_FLAG.LT.0).OR.(RES_BC_FLAG.EQ.1))  THEN
           IF (NOPE.GT.0) THEN
             READ(BCTIMEMSG,*) RBCTIMEINC
             READ(BCSTATMSG,*) BCSTATIM
           END IF
        ELSE IF (RES_BC_FLAG.EQ.2) THEN
           IF (NOPE.GT.0) THEN
             READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC
             READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM
           END IF
        ELSE IF (RES_BC_FLAG.EQ.3) THEN
           IF (NOPE.GT.0) THEN
             READ(BCTIMEMSG,*) RBCTIMEINC, TBCTIMEINC
             READ(BCSTATMSG,*) BCSTATIM, TBCSTATIM
           END IF
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
           END IF
        ELSE IF (RES_BC_FLAG.EQ.4) THEN
           IF (NOPE.GT.0) THEN
             READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC, TBCTIMEINC
             READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM, TBCSTATIM
           END IF
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
           END IF
        END IF
      END IF

      IF (CBAROCLINIC) THEN
        READ(15,80) SPONGEDISTMSG
        READ(15,80) EqnstateMSG
      END IF

C
C     jgf45.12: Read in the parameters for the transport equation, if
C     necessary.
      IF (C3D_BTrans) THEN

C     Lateral and vertical diffusion coefficients.
         READ(15,*) NLSD, NVSD
         READ(15,*) NLTD, NVTD
C     Time stepping coefficient for the transport equation terms.
         READ(15,*) ALP4
C     Temperature boundary condition file type, if necessary
!         if ( IDEN .eq. 3 .or. IDEN .eq. 4 ) then
!            READ(15,*) NTF
!         endif
         endif

      call memory_status()
      RETURN
C
 80   FORMAT(A80)
C----------------------------------------------------------------------
      END SUBROUTINE READ15_3DVS
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C         S U B R O U T I N E    S E T U P K D T S E A R C H
C----------------------------------------------------------------------
C
C  tcm v49.48.01 Subroutine that creates a global list of element
C             center's and a global list of the size of an
C             element's radius
C  Written by:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39056
C
C----------------------------------------------------------------------
      subroutine setup_kdt_search
      use pre_global, only : SLAM,SFEA,NNEG,RMAX,BCXY,NELG

      implicit none
      integer :: i
      integer :: ielm(3)
      real(8) :: xelm(3),yelm(3),shplocal(7)

      do i=1,NELG
         ielm(:) = NNEG((/1,2,3/),i)  !element's node numbers
         xelm(:) = slam(ielm(:))      !element's vertex x-values
         yelm(:) = sfea(ielm(:))      !element's vertex y-values
         shplocal(1)= yelm(2)-yelm(3)
         shplocal(2)= yelm(3)-yelm(1)
         shplocal(3)= yelm(1)-yelm(2)
         shplocal(4)= xelm(3)-xelm(2)
         shplocal(5)= xelm(1)-xelm(3)
         shplocal(6)= xelm(2)-xelm(1)
         shplocal(7)= shplocal(1)*shplocal(5)
     &                   - shplocal(2)*shplocal(4)

!        compute the radius of the circle that circumscribes
!        the triangle then scale it by 50% larger to allow for
!        a buffer later on
         rmax(i)=1.50D0*( sqrt(shplocal(6)*shplocal(6)+
     &                         shplocal(3)*shplocal(3))*
     &                    sqrt(shplocal(4)*shplocal(4)+
     &                         shplocal(1)*shplocal(1))*
     &                    sqrt(shplocal(5)*shplocal(5)+
     &                         shplocal(2)*shplocal(2))/
     &                     (2.d0*shplocal(7)) )
!        Compute the barycenter of the element
          bcxy(1,i) = sum(xelm)/3.d0
          bcxy(2,i) = sum(yelm)/3.d0
       enddo !i

      return
C----------------------------------------------------------------------
      end subroutine setup_kdt_search
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C         S U B R O U T I N E    K D T S E A R C H
C----------------------------------------------------------------------
C
C  tcm v49.48.01 Subroutine that uses the KDTREE2 algorithm for
C             finding what element a point lies in.
C  Written by:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39056
C
C----------------------------------------------------------------------
      subroutine kdtsearch(XCoords,YCoords,FDEle,NumofStations,Desc)
      use PRE_Global
      use kdtree2_module

      implicit NONE
      INTEGER, intent(in) :: NumOfStations                     ! total
      REAL(8), intent(in), dimension(NumOfStations) :: XCoords ! cartesian
      REAL(8), intent(in), dimension(NumOfStations) :: YCoords ! cartesian
      INTEGER, intent(out), dimension(NumOfStations):: FDEle   ! FullDomain
      CHARACTER(len=30), intent(in) :: Desc                    ! description

      INTEGER :: I,itc,iek
      INTEGER :: ielm(3)
      real(8) :: Xsta,Ysta,dist
      real(8) :: x1,x2,x3,y1,y2,y3,A1,A2,A3,AA,AREASK,AE
      real(8) :: elmmin(2),xelm(3),yelm(3)
      TYPE(KDTREE2), POINTER :: TREE
      TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)
      LOGICAL ElementFound  ! .true. when a corresponding element is found
      INTEGER ClosestElement  ! element with closest match
      REAL(sz), PARAMETER :: Tolerance = 1.0d-5     ! area difference for match

      ElementFound = .false.

!    Be sure the maximum search depth is not larger than
!    the number of elements being kept
      if (nelg.lt.srchdp) srchdp = nelg

!     Create the search tree
      tree => kdtree2_create(bcxy,rearrange=.true.,sort=.true.)

!     allocate space for the search results from the tree
      allocate(KDRESULTS(srchdp))


      do I = 1, NumOfStations
         Xsta = XCoords(I)
         Ysta = YCoords(I)
         ElementFound = .false.

        ! Find the srchdp nearest elements to this point
         call kdtree2_n_nearest(tp=tree,qv=(/Xsta,Ysta/),
     &                  nn=srchdp,results=KDRESULTS)
         itc = 1
         ClosestElement = KDRESULTS(itc)%idx
!       Check to see if the point lies with rmax of any of these elements
         elmmin = minval(sqrt(KDRESULTS(1:srchdp)%dis)
     &                              - rmax(KDRESULTS(1:srchdp)%idx) )
         if(elmmin(1).le.0.0D0) then  ! Point lies within search radius of an element
!           loop through the elements in the search list
            do while ((ElementFound.eqv..false.).and.(itc.le.srchdp))
               iek = KDRESULTS(itc)%idx  !Current search element number
!              Get the distance from this point to the barycenter of the
!              current element
               dist = sqrt(KDRESULTS(itc)%dis)
!              If the distance is less than or equal to rmax (rmax=1.5*element radius)
!              Then the point is near the element and might be in it
!              Proceed with the weights test
               if(dist-rmax(iek).le.0.0d0) then
                  !get the shape function for this element
                  ielm(:) = NNEG((/1,2,3/),iek)  !element's node numbers
                  xelm(:) = slam(ielm(:))      !element's vertex x-values
                  yelm(:) = sfea(ielm(:))      !element's vertex y-values
                  X1=xelm(1)
                  X2=xelm(2)
                  X3=xelm(3)
                  Y1=yelm(1)
                  Y2=yelm(2)
                  Y3=yelm(3)
                  A1=(Xsta-X3)*(Y2-Y3)+(X2-X3)*(Y3-Ysta)
                  A2=(Xsta-X1)*(Y3-Y1)-(Ysta-Y1)*(X3-X1)
                  A3=(Ysta-Y1)*(X2-X1)-(Xsta-X1)*(Y2-Y1)
                  AA=ABS(A1)+ABS(A2)+ABS(A3)
                  AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
                  AE=ABS(AA-AREASK)/AREASK
               IF (AE.LT.Tolerance) THEN
                  ElementFound = .true.
                  ClosestElement = iek
                  FDEle(I) = ClosestElement
                  else !not in this element keep looking
                    itc = itc + 1
                  endif !End area ratio test
               else !
!                point is too far away from the barycenter of the
!                element to possibly be in the element, so move to
!                the next element
                 itc = itc + 1
               endif !end Radius test
            enddo !end the while loop
         endif
         if(.not.ElementFound) then !We did not find the element
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER.EQ.1) THEN
               print *, "The program will estimate nearest element."
               print *, "WARNING. Distance to nearest element is ",
     &                   sqrt(KDRESULTS(1)%dis)
               print *, " "
               FDEle(I) = ClosestElement
            ELSE
               print *, "ERROR. Distance to nearest element is ",
     &                   sqrt(KDRESULTS(1)%dis)
               CALL EXIT(1)
            ENDIF
         ENDIF
      enddo !loop over station points

!     Deallocate the tree
      call kdtree2_destroy(tp=tree)

 1234 format(A30,1x,I4,1x,'does not lie in the grid.')

      return

      end subroutine kdtsearch


C-----------------------------------------------------------------------
C          S U B R O U T I N E   C O O R D  T O  E L E
C-----------------------------------------------------------------------
C
C     jgf45.11 Subroutine to take a list of X and Y cartesian
C     coordinates and find the corresponding elements.
C
C-----------------------------------------------------------------------
      SUBROUTINE CoordToEle(XCoords,YCoords,FDEle,NumOfStations,Desc)
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: NumOfStations                     ! total
      REAL(8), intent(in), dimension(NumOfStations) :: XCoords ! cartesian
      REAL(8), intent(in), dimension(NumOfStations) :: YCoords ! cartesian
      INTEGER, intent(out), dimension(NumOfStations):: FDEle   ! FullDomain
      CHARACTER(len=30), intent(in) :: Desc                    ! description

      INTEGER I,J,K,M,N                                        ! loop counters
      INTEGER N1, N2, N3, KMIN
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3       ! geometry
      REAL(8) AE, AEMIN, AREASK, AA                            ! area
      LOGICAL ElementFound  ! .true. when a corresponding element is found

      DO I = 1, NumOfStations
         ElementFound = .false.
         AEMIN=1.0E+25
         KMIN=0
         DO K=1,NELG
            N1=NNEG(1,K)
            N2=NNEG(2,K)
            N3=NNEG(3,K)
            X1=SLAM(N1)
            X2=SLAM(N2)
            X3=SLAM(N3)
            X4=XCoords(I)
            Y1=SFEA(N1)
            Y2=SFEA(N2)
            Y3=SFEA(N3)
            Y4=YCoords(I)
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
            AE=ABS(AA-AREASK)/AREASK
            IF (AE.LT.AEMIN) THEN
               AEMIN=AE
               KMIN=K
            ENDIF
            IF (AE.LT.1.0E-5) THEN
               ElementFound = .true.
               FDEle(I)=K
            ENDIF
         ENDDO
         IF ( ElementFound .eqv. .false.) THEN
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER.EQ.1) THEN
               print *, "The program will estimate nearest element."
               print *, "WARNING. Proximity index is ",AEMIN
               print *, " "
               FDEle(I) = KMIN
            ELSE
               print *, "ERROR. Proximity index is ",AEMIN
               CALL EXIT(1)
            ENDIF
         ENDIF
      ENDDO

 1234 format(A30,1x,I4,1x,'does not lie in the grid.')

      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE CoordToEle
C-----------------------------------------------------------------------
