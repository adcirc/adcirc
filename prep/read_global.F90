!----------------------------------------------------------------------------

!                           MODULE READ_GLOBAL

!----------------------------------------------------------------------------

!                  For use with ADCPREP Version 2.1 ( 03/17/2006 )

!                     current for ADCIRC v43.03   5/20/2003
!                 jgf Updated for ADCIRC v45.07  11/08/2005
!                 jgf Updated for ADCIRC v45.09  01/09/2006
!                 jgf Updated for ADCIRC v45.11  01/09/2006
!                 jgf Updated for ADCIRC v45.12  03/17/2006
!----------------------------------------------------------------------------


!  version 1.2 vjp 12/7/99
!  version 1.3 ral 10/10/01 revisions
!  version 1.6 meb 3/03 & ral 5/21/03
!  version 1.7 meb 4/04 added by jgf


!----------------------------------------------------------------------
!         S U B R O U T I N E     P R E P   R E A D   1 4
!----------------------------------------------------------------------
!     Reads the global ADCIRC grid data file for data decomposition
!     program ADCPP.  This version compatible with ADCIRC_v34.03
!     vjp  2/28/98
!---------------------------------------------------------------------
    subroutine prepRead14()
    use pre_global, only : x, y, dp, nneg, ibtypee, nvdllmsg, nvdll, &
    nbdv, nvellmsg, nvell, lbcode, ibtype, nbvv, ibconnr, &
    bar1, bar2, bar3, bar4, bar5, bar6, weir, weird, flux14_ary, &
    bcts14_ary, agrid, aperiodic_bc_ts, bndbcriver, exist_bc_ts, &
    exist_flux, iden, mnbou, mnvel, nbou, nboumsg, nelg, neta, &
    netamsg, nfover, nnodg, nope, nopemsg, nvel, nvelmsg, nweir, &
    river_above_msl, numLBCodeValues, nfluxf

    use memory_usage
    implicit none

    integer :: nbytes = 0
    INTEGER :: I,J,JW,K,ITEMP,ITYPE
    INTEGER :: DISC,BBN,IBP,I1
    INTEGER :: MAXNEIGH
    INTEGER :: NUMLINES

! Dummy variables

    INTEGER :: DNBOU,DNVEL,DITYPE
    INTEGER,ALLOCATABLE :: DNVELL(:),DNBVV(:,:)
    INTEGER,ALLOCATABLE :: DIBTYPE(:)
    CHARACTER(80) DNBOUMSG,DNVELMSG
    CHARACTER(80),ALLOCATABLE ::  DNVELLMSG(:)
    CHARACTER(80) DUM
! kmd - variable for river in baroclinic simulation
    LOGICAL :: BC_Rivers

    ALLOCATE ( DNVELL(MNBOU),DNBVV(MNBOU,0:MNVEL),DIBTYPE(MNBOU) )
    ALLOCATE ( DNVELLMSG(MNBOU+1) )
    nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
    call memory_alloc(nbytes)

!--The Grid file was opened in SIZEUP performed a rewind, and is ready here

!--Read Grid Title

    EXIST_FLUX = 0
    READ(14,80) AGRID

!--Read Total Number of Elements and Nodes

    READ(14,*) NELG,NNODG
!     READ(14,80) NSIZES
!     READ(NSIZES,*) NELG,NNODG
!     CALL GETMSG(NSIZES,SIZEMSG)

!--Read Nodal Coordinates and Bathymetry
!  If ICS=2 Will Convert later in read15

    NUMLINES=2
    DO I = 1,NNODG
        READ(14,*) J,X(I),Y(I),DP(I)
        IF (J /= I) THEN
            print *, I,J
            STOP 'Node Numbering not in Sequential Order'
        ENDIF
    ENDDO

!--Read Element Connectivity Table

    DO I = 1,NELG
        READ(14,*) J,ITEMP,NNEG(1,I),NNEG(2,I),NNEG(3,I)
        IF (J /= I) THEN
            print *, I,J
            STOP 'Element Numbering not Sequential'
        ENDIF
    ENDDO
    NUMLINES=NUMLINES+NNODG+NELG
!...
!...Read Total Number of Open Boundary Segments
!...
    READ(14,80) NOPEMSG
    READ(NOPEMSG,*) NOPE
    NUMLINES=NUMLINES+1
!...
!...Read Total Number of Open Boundary Forcing Nodes
!...
    READ(14,80) NETAMSG
    READ(NETAMSG,*) NETA
    NUMLINES=NUMLINES+1
!...
!...Read Number of Nodes on Open Boundary Segment
!...and Segment Nodes Numbers
!...
    J=0
    IBTYPEE(:) = 0 !jgf49.32 initialize all elements of the array to zero
    DO K=1,NOPE
        READ(14,80) NVDLLMSG(K) !jgf49.32: TODO: IBTYPEE should be read here
        NUMLINES=NUMLINES+1
        READ(NVDLLMSG(K),*) NVDLL(K)
        DO I=1,NVDLL(K)
            READ(14,*) NBDV(K,I)
            NUMLINES=NUMLINES+1
        ENDDO
        J=J+NVDLL(K)
    ENDDO

    IF (NETA /= J) THEN
        print *, "Total Number of Boundary Nodes = ",J
        print *, "This exceeds NETA = ",NETA
        IF (NFOVER == 1) THEN  !jgf51.21.27: Has NFOVER been read at this point?
            NETA = J
            print *, "ADCPP corrected this error"
        ELSE
            stop
        ENDIF
    ENDIF

!--Read Total Number of Land Boundary Segments

    READ(14,80) NBOUMSG
    READ(NBOUMSG,*) NBOU

!--Read Total of Land Boundary Nodes

    READ(14,80) NVELMSG
    READ(NVELMSG,*) NVEL

!--Read Number of Nodes in the Land Boundary Segment and Boundary Type
!  and construct LBCODE array for read15 routine

    J=0
    NWEIR = 0
! jgf51.52.23: The LBCODE array is dimensioned by MNVEL which is
! larger than the number of values in the LBCODE array. So we
! need to initialize to a value that will hopefully cause obvious
! isses if an initialized value is accidentally used.
    LBCODE(:) = -99
    DO K = 1,NBOU
    
        READ(14,80) NVELLMSG(K)
        READ(NVELLMSG(K),*) NVELL(K),IBTYPE(K)
        ITYPE = IBTYPE(K)
    ! kmd - added for rivers in a baroclinic simulation
        IF (ABS(ITYPE/100) == 1) THEN
            IF (IDEN <= 0) THEN
                ITYPE=(ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
            ELSE IF (IDEN > 0) THEN
                APERIODIC_BC_TS= .TRUE. 
                BndBCRiver= .TRUE. 
                ITYPE=(ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
            END IF
        END IF
    
        DO I=1, NVELL(K)
            J = J+1
            LBCODE(J) = ITYPE
        ENDDO
    
    ! jgf51.21: Cleaned up reading of boundary arrays by replacing
    ! a complex set of if/thens with a simple select case.
        select case(ITYPE)
        case(0,10,20,30,40,1,11,21,41)
        DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I)
            IBCONNR(K,I) = 0
        ENDDO
        case(3,13,23)
        DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I),BAR1(K,I),BAR2(K,I)
            IBCONNR(K,I) = 0
        ENDDO
        case(4,24)
        DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I),IBCONNR(K,I), &
            BAR1(K,I),BAR2(K,I),BAR3(K,I)
        !--Construct List of WEIR nodes and their duals
            NWEIR = NWEIR+1
            WEIR(NWEIR) = NBVV(K,I)
            WEIRD(NWEIR) = IBCONNR(K,I)
        ENDDO
        case(5,25)
        DO I =1,NVELL(K)
            READ(14,*) NBVV(K,I),IBCONNR(K,I), &
            BAR1(K,I),BAR2(K,I),BAR3(K,I), &
            BAR4(K,I),BAR5(K,I),BAR6(K,I)
        !--Construct List of WEIR nodes and their duals
            NWEIR = NWEIR+1
            WEIR(NWEIR) = NBVV(K,I)
            WEIRD(NWEIR) = IBCONNR(K,I)
        ENDDO
        case(2,12,22,32,52)
        NFLUXF=1 !jgf51.52.23
        DO I=1,NVELL(K)
            EXIST_FLUX = EXIST_FLUX + 1
            READ(14,*) NBVV(K,I)
            IBCONNR(K,I) = 0
        ! md - added river boundary conditions for baroclinic
            IF (BndBCRiver) THEN
                EXIST_BC_TS=EXIST_BC_TS+1
            END IF
        ENDDO
        BndBCRiver= .FALSE. 
        case default
        print *, "IBTYPE not set correctly for segment ",K
        stop
        end select
    ENDDO
! jgf51.52.23: Record the actual number of values stored in the
! BCODE array.
    numLBCodeValues = J


!     jgf45.06 added the following section
! MEB --------------------------------------------------
! MEB This entire section written 04/01/04 to accomodate
! MEB division of fort.20 files between PE directories.
! MEB This must be modified if fort.14 format changes.
! MEB --------------------------------------------------
    IF (EXIST_FLUX > 0) THEN
        REWIND(14)
        ALLOCATE ( FLUX14_ARY(EXIST_FLUX) )
    ! md - added river boundary conditions for baroclinic
        ALLOCATE ( BCTS14_ARY(EXIST_BC_TS) )
        EXIST_FLUX=0
        EXIST_BC_TS=0
    
    !  SKIP OVER UNNEEDED LINES
    
        DO I=1,NUMLINES
            READ(14,*) DUM
        ENDDO
    
        READ(14,80) DNBOUMSG
        READ(NBOUMSG,*) DNBOU
    
    !--Read Total of Land Boundary Nodes
    
        READ(14,80) DNVELMSG
        READ(NVELMSG,*) DNVEL
    
        DO K=1,DNBOU
            READ(14,80) DNVELLMSG(K)
            READ(DNVELLMSG(K),*) DNVELL(K), DIBTYPE(K)
            DITYPE=DIBTYPE(K)
        ! kmd - added in boundary conditions for rivers in a baroclinic simulation
            IF (ABS(DITYPE/100) == 1) THEN
                IF (IDEN <= 0) THEN
                    DITYPE=(ABS(DITYPE)-100)*(DITYPE/ABS(DITYPE))
                ELSE IF (IDEN > 0) THEN
                    BC_Rivers= .TRUE. 
                    DITYPE=(ABS(DITYPE)-100)*(DITYPE/ABS(DITYPE))
                END IF
            END IF
        ! jgf50.21 Use select case and add IBTYPE 32.
            select case(DITYPE)
            case(2,12,22,32,52)
            DO I=1, DNVELL(K)
                EXIST_FLUX = EXIST_FLUX + 1
                READ(14,*) DNBVV(K,I)
                FLUX14_ARY(EXIST_FLUX) = DNBVV(K,I)
            ! kmd - added in boundary conditions for rivers in a baroclinic simulation
                IF (BC_Rivers) THEN
                    EXIST_BC_TS=EXIST_BC_TS+1
                    BCTS14_ARY(EXIST_BC_TS)=DNBVV(K,I)
                END IF
                IF (DP(DNBVV(K,I)) < 0.d0) THEN
                    River_above_MSL= .TRUE.  ! kmd - Evan's changes for rivers above MSL
                END IF
            ENDDO
            case default
            DO I=1, DNVELL(K)
                READ(14,*) DNBVV(K,I)
            ENDDO
            end select
            BC_Rivers = .FALSE.  ! kmd - reset logic for next river boundary segment
        ENDDO
        DEALLOCATE ( DNVELL,DNBVV,DIBTYPE )
        DEALLOCATE ( DNVELLMSG )
        nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
        call memory_dealloc(nbytes)
    ENDIF


!     jgf45.06 end of section written by MEB04/01/04

!--Close Global Grid file

    CLOSE(14)
    call memory_status()

    80 format(a80)

    return
!----------------------------------------------------------------------
    end subroutine prepRead14
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!                        S U B R O U T I N E
!            G E T   M E S H   P A R A M E T E R S   X D M F
!----------------------------------------------------------------------
!     jgf51.21.29: If the mesh file is in XDMF format, the adcirc module
!     mesh.F is used to read the data. In order to get the data into the
!     pre_global module, this subroutine is required. Someday, adcprep
!     will no longer use pre_global, and will just rely on the adcirc
!     global module, and we will be able to remove this subroutine.
!----------------------------------------------------------------------
    subroutine getMeshParametersXDMF()
    use memory_usage
    use presizes, only : fluxBoundary, mne, mnp, mnei, mneta, mnope, &
    skipNeta
    use pre_global, only : x, y, dp, nneg, ibtypee, nvdllmsg, nvdll, &
    nbdv, nvellmsg, nvell, lbcode, ibtype, nbvv, ibconnr, &
    bar1, bar2, bar3, weir, weird, flux14_ary, bcts14_ary, agrid, &
    aperiodic_bc_ts, bndbcriver, exist_bc_ts, exist_flux, iden, &
    mnbou, mnvel, nbou, nboumsg, nelg, neta, netamsg, nfover, &
    nnodg, nope, nopemsg, nvel, nvelmsg, nweir, river_above_msl, &
    alloc_main1
    use mesh, only : getMeshSizesForPrep, getMeshForPrep, nm
    use boundaries, only : getBoundarySizesForPrep, &
    getBoundariesForPrep
    implicit none
    integer, allocatable :: nneigh(:)
    integer :: i, j, k, nbytes


    call getMeshSizesForPrep(nnodg, nelg)
    call getBoundarySizesForPrep(nbou, nvel, neta, nope, &
    nweir, exist_flux)
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
    if (exist_flux /= 0) then
        fluxBoundary = .TRUE. 
    endif
    write(*,'(a,i0)') "neta =", neta
    call alloc_main1()

! get the element table
    do i=1,3
        nneg(i,:) = nm(:,i)
    end do

! get mesh parameters
    call getMeshForPrep(agrid, x, y, dp)
    call getBoundariesForPrep(nvdll, nvell, ibtypee, nbdv, lbcode, ibtype, nbvv, &
    ibconnr, bar1, bar2, bar3, weir, weird, iden, nfover, &
    exist_flux, flux14_ary)
    do k=1,nope
        write(nvdllmsg(k),'(i0)') nvdll(k)
    end do
    do k=1,nbou
        write(nvellmsg(k),'(i0)') nvell(k)
    end do

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
!     determine the maximum NNeigh
    mnei=maxval(nneigh)
    mnei = mnei+1
    deallocate(nneigh)
!----------------------------------------------------------------------
    end subroutine getMeshParametersXDMF
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!     Reads the global ADCIRC Input Data File for data decomposition
!     program ADCPP.  This version compatible with ADCIRC_v34.03
!     vjp  3/28/98
!---------------------------------------------------------------------
    SUBROUTINE READ15()
    USE PRE_GLOBAL
    USE SIZES, ONLY : MNHARF
    USE MESH, ONLY : CPP
    USE GLOBAL, ONLY : parse, a2f, hss, screenUnit, STATNAME, &
    STATNAMEV, STATNAMEM, STATNAMEC, ERROR, DEBUG, &
    allMessage, setMessageSource, deg2rad, &
    unsetMessageSource, INFO, logMessage
    USE HARM, ONLY : NHAINC, ALLOC_HA, NAMEFR, HAFREQ, HAFF, HAFACE, &
    THAS, THAF, FMV, NHASE, NHASV, NHAGE, NHAGV, &
    IHARIND, ALLOC_MAIN14, CHARMV
    USE PRESIZES, ONLY : SZ
    USE PREP_WEIR
    USE NODALATTRIBUTES, ONLY : outputTau0
    use memory_usage
    IMPLICIT NONE

!    ===================================================================
!      MCF ADDITION FOR READING AND WRITING STATION INFO IN NETCDF

!      REAL, EXTERNAL              :: a2f
    CHARACTER(132) STATLINE
    CHARACTER(50) LVAR(3)

!    ===================================================================


    integer :: nbytes = 0
    REAL(SZ) RSTIMINC
    INTEGER :: N1, N2, N3, KMIN, JG, INDX
    INTEGER :: I,J,K,M
    INTEGER :: IG1,IG2,IG3,IL1,IL2,IL3
    REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3
    REAL(8) AE, AEMIN, AREASK, AA
    INTEGER :: NBV(MNVEL)
    INTEGER :: ios  ! status of read operations, 0=success
! Casey 090825: Added dummy variables.
    INTEGER :: IDUM
    REAL(8) RDUM
! kmd - add a variable for boundary conditions
    INTEGER :: ITYPE

!   tcm v50.66.02 -- added timebathycontrol namelist related variables
    INTEGER :: ios_nddt
    INTEGER :: ios_metCon
    INTEGER :: ios_tvw
! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    INTEGER :: IOS_WC
!...tcm v51.20.03 Local variables for external station location files
    INTEGER :: NSTAE2,NSTAV2,NSTAC2,NSTAM2
    INTEGER :: IOS_STATIONS
    integer :: ios_wetDry
    namelist /wetDryControl/ outputNodeCode, outputNOFF, noffActive
! jgf52.08.02 : For inundation output files.
    namelist /inundationOutputControl/ inundationOutput, inunThresh

!.....tcm v51.20.03 initialize external station location specifications
!..... to false, meaning if locations are supplied they are in the fort.15
!..... file.
    USE_ELEV_STAT_FILE = .FALSE. 
    USE_VEL_STAT_FILE = .FALSE. 
    USE_MET_STAT_FILE = .FALSE. 
    USE_CONC_STAT_FILE = .FALSE. 

    found_tbc_nml = .FALSE.   !flag to determine if the timebathycontrol namelist was present
    found_metCon_nml = .FALSE.   !flag to determine if the metControl namelist was present
    found_tvw_nml = .FALSE. 
! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    FOUND_WC_NML = .FALSE. 


    call setMessageSource("read15")
#ifdef READ_GLOBAL_TRACE
    call screen(DEBUG,"Enter.")
#endif

!--The Run Info file was opened in SIZEUP, performed a rewind, and is ready here

!--Run Description and Run Identification


!...  tcm v50.66.02 Addtions for Time Varying Bathymetry
!...  read through the fort.15 file for the namelist (TimeBathyControl) for
!...  the time varying bathymetry.  This namelist must be at the bottom of the
!...  fort.15 file. If found, then set the appropriate values (btiminc,bchgtiminc,
!...  and nddt).  If the namelist is not there, then the time varying bathymetry
!...  will not be used.
!...
!...  After this search and read, we will close the file and then reopen it
!...  for further processing the traditional non-namelist components.
!...
    NDDT = 0  !Set default value to be no time varying bathymetry
    READ(UNIT = 15,NML = TimeBathyControl,IOSTAT = IOS_NDDT)

    IF (IOS_NDDT < 0) THEN
    !.....   it is possible for the namelist to be present in the file and occuring at the end
    !.....   of the file with no line breaks after the ending "\" which causes the iostat to
    !.....   return a negative value.  By checking to be sure a namelist variable was set to
    !....    a non-default value we can determine this was the case.
        IF (NDDT /= 0) THEN
        !            WRITE(*,*) 'NAMELIST PRESENT, BUT AT THE END OF FILE WITH',
        !     &                 ' NO ADVANCING CHARACTER'
            found_tbc_nml = .TRUE. 
        ELSE
            call logMessage(INFO, &
            'The timeBathyControl namelist was not found.')
        ENDIF
    ELSEIF (IOS_NDDT == 0) THEN
    !         WRITE(*,*) 'NAME LIST PRESENT AND CORRECT'
        found_tbc_nml = .TRUE. 
    ELSE
        found_tbc_nml = .TRUE. 
        WRITE(*,*) &
        'THERE WAS A PROBLEM PROCESSING THE TimeBathyControl NAMELIST'
        WRITE(*,*) 'IN THE FORT.15 FILE.  SHUTTING DOWN ADCIRC NOW.'
        STOP   !THERE IS a STOP HERE
    ENDIF

    REWIND(15)  !Return to the top of the fort.15 file (instead of closing then opening)
!..

!... tcm v50.79 added
    READ(UNIT = 15, NML = metControl,IOSTAT = ios_metCon)
    if ( ios_metCon >= 0) then
        found_metCon_nml = .TRUE. 
        call logMessage(INFO, &
        'The metControl namelist was found.')
    else
        call logMessage(INFO, &
        'The metControl namelist was not found.')
    endif
    rewind(15)  !Return to the top of the fort.15 file

! bell - Added TVW Namelist
    READ(15,NML=tvwControl,IOSTAT=IOS_TVW)
    IF(IOS_TVW == 0)THEN
        FOUND_TVW_NML = .TRUE. 
    ELSE
        FOUND_TVW_NML = .FALSE. 
    ENDIF
    REWIND(15)

! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
    READ(UNIT=15,NML=waveCoupling,IOSTAT=IOS_WC)
    IF(IOS_WC == 0)THEN
        FOUND_WC_NML = .TRUE. 
        call logMessage(INFO, &
        'The waveCoupling namelist was found.')
    else
        call logMessage(INFO, &
        'The waveCoupling namelist was not found.')
    ENDIF
    REWIND(15)
     
    READ(UNIT = 15, NML = wetDryControl,IOSTAT = ios_wetDry)
    if ( ios_wetDry >= 0) then
        foundWetDryControlNameList = .TRUE. 
        call logMessage(INFO, &
        'The wetDryControl namelist was found.')
    else
        call logMessage(INFO, &
        'The wetDryControl namelist was not found.')
    endif
    rewind(15)  !Return to the top of the fort.15 file

    read(unit=15,nml=inundationOutputControl,iostat=ios)
    if ( ios >= 0) then
        foundInundationOutputControlNamelist = .TRUE. 
        call logMessage(INFO, &
        'The inundationOutput namelist was found.')
    else
        call logMessage(INFO, &
        'The inundationOutput namelist was not found.')
    endif
               
    rewind(15)  !Return to the top of the fort.15 file


    READ(15,80) RUNDES

    READ(15,80) RUNID

    READ(15,80) OVERMSG
    READ(OVERMSG,*) NFOVER
    IF (NFOVER == 1) THEN
    !       print *, "Non-fatal errors will be corrected"
    ELSE
    !      print *, "Non-fatal errors will stop execution"
    ENDIF

    READ(15,80) ABOUTMSG
    READ(ABOUTMSG,*) NABOUT

    READ(15,80) SCREENMSG
    READ(SCREENMSG,*) NSCREEN
    screenUnit = 6 ! always write to screen in adcprep

    READ(15,80) HOTMSG
    READ(HOTMSG,*) IHOT


    READ(15,80) ICSMSG
    READ(ICSMSG,*) ICS
    IF ((ICS /= 1) .AND. (ICS /= 2)) THEN
        print *, "ICS set incorrectly"
        STOP
    ENDIF

    READ(15,80) IMMSG
    READ(IMMSG,*) IM
    IF(IM == 2) THEN
        PRINT *, "DSS Model type not presently supported"
        STOP
    ENDIF
!     jgf46.28 Read IDEN if necessary
    IF (CBaroclinic) READ(15,*) IDEN

    READ(15,80) IBFMSG
    READ(IBFMSG,*) NOLIBF
    IF((NOLIBF < 0) .OR. (NOLIBF > 2)) THEN
        print *, "Value for NOLIBF not allowed"
        stop
    ENDIF

    READ(15,80) IFAMSG
    READ(IFAMSG,*) NOLIFA
    IF ((NOLIFA < 0) .OR. (NOLIFA > 3)) THEN
        print *, "Value for NOLIFA not allowed"
        stop
    ENDIF

    READ(15,80) ICAMSG
    READ(ICAMSG,*) NOLICA
    IF ((NOLICA < 0) .OR. (NOLICA > 1)) THEN
        print *, "Value for NOLICA not allowed"
        stop
    ENDIF

    READ(15,80) ICATMSG
    READ(ICATMSG,*) NOLICAT
    IF ((NOLICAT < 0) .OR. (NOLICAT > 1)) THEN
        print *, "Value for NOLICAT not allowed"
        stop
    ENDIF

! jgf52.05: Removed this error message because it is alarming
! and because it contradicts common practice.
!      IF ((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
!         print *, "NOLIFA and NOLICAT are inconsistent"
!         print *, "May lead to mass balance problems"
!         IF(NFOVER.EQ.1) THEN
!            print *, "Since NFOVER=1, Program will continue"
!        ELSE
!            stop
!         ENDIF
!      ENDIF

    READ(15,80) NWPMSG
    READ(NWPMSG,*) NWP
    IF (NWP > 0) THEN !jgf46.00 read nodal attribute labels
        ALLOCATE(useNodalAttrNames(NWP))
        nbytes = 4*nwp
        call memory_alloc(nbytes)
        DO I=1, NWP
            READ(15,80) useNodalAttrNames(I)
        ENDDO
    ENDIF

    READ(15,80) NCORMSG
    READ(NCORMSG,*) NCOR
    IF ((NCOR /= 0) .AND. (NCOR /= 1)) THEN
        print *, "Value for NCOR not allowed"
        IF (NFOVER == 1) THEN
            NCOR = 0
            print *, "NCOR has been reset to 0"
        ELSE
            stop
        ENDIF
    ENDIF

    IF ((ICS == 1) .AND. (NCOR == 1)) THEN
        print *, "ICS=1 and NCOR=1 may lead to geometric distortions"
        IF(NFOVER == 1) THEN
            print *, "Program will continue with these input values"
            print *, "for large domains it is recommended to use ICS=2"
        ELSE
            stop
        ENDIF
    ENDIF

    READ(15,80) NTIPMSG
    READ(NTIPMSG,*) NTIP
    IF ((NTIP < 0) .OR. (NTIP > 2)) THEN
        print *, "Value for NTIP not allowed"
        IF(NFOVER == 1) THEN
            NTIP = 0
            print *, "NTIP has been reset to 0"
        ELSE
            stop
        ENDIF
    ENDIF

    IF ((ICS == 1) .AND. (NTIP >= 1)) THEN
        print *, "ICS=1 & NTIP >= 1 may lead to geometric distortions"
        print *, "for large domains it is recommended to use ICS=2"
        IF (NFOVER == 1) THEN
            print *, "Program will continue with these input values"
        ELSE
            stop
        ENDIF
    ENDIF

    NRS=0
    NCICE = 0  !tcm v49.64.01 -- added for ice
    READ(15,80) NWSMSG
    READ(NWSMSG,*) NWS

!....tcm v49.46 -- added logic to handle multiple NRS types (100's place)
!.....tcm v49.64.01 Additions for ice
    IF(NWS == 0) THEN
        NWS = 0
        NRS = 0
        NCICE = 0
    ELSE
        NCICE = INT(ABS(NWS)/1000)
        NRS=INT((ABS(NWS) - NCICE*1000)/100)
        NWS=INT((ABS(NWS)- NCICE*1000 - NRS*100))*INT(NWS/ABS(NWS))
    ENDIF

!     jgfdebug46.02 added NWS=45
!     jgf46.02 Added NWS=8.
!     jgf46.16 merged:
!     cf & cm added NWS=9: asymmetric hurricane wind model
!     jie added NWS=20: generalized asymmetric vortex model
!     sb46.28sb01 added NWS=12: OWI format
!     jgf49.0804 Added NWS29 VortexOWI format.
    IF((NWS /= 0) .AND.    (NWS /= 1 ) .AND. (ABS(NWS) /= 2) .AND. &
    (NWS /= 3) .AND. (ABS(NWS) /= 4) .AND. (ABS(NWS) /= 5) .AND. &
    (ABS(NWS) /= 45) .AND. (ABS(NWS) /= 6) .AND. &
    (ABS(NWS) /= 8) .AND. (ABS(NWS) /= 15) .AND.  & !jgf50.38.05 Added NWS=15
    (ABS(NWS) /= 12) .AND. (ABS(NWS) /= 19) .AND. &
    (ABS(NWS) /= 20) .AND. &
    (ABS(NWS) /= 29) .AND. (abs(NWS) /= 16))  THEN
    print *,"ERROR: Value for NWS not supported by parallel code"
    stop
ENDIF

!... TCM v49.64.02 -- added
!...  BE SURE NWS AND NCICE ARE COMPATABLE
    IF((NCICE > 0) .AND. &
    ((NWS == 1) .OR. (NWS == 2) .OR. (NWS == 7))) THEN
        PRINT*,'NCICE = ',NCICE
        PRINT*,'NWS = ', NWS
        PRINT*, "Your selection of NCICE (a UNIT 15 input ", &
        "parameter) is not compatable"
        PRINT*,"with your NWS value.  NCICE is not allowed for", &
        " abs(NWS)=1,2, or 7."
        STOP
    ENDIF

!     jgf46.08 Modified to accomodate fine grained ramp functions.
    READ(15,80) RAMPMSG
    READ(RAMPMSG,*) NRAMP
! Corbitt 120322: Includes Zach's Wind Ramping
    IF ((NRAMP /= 0) .AND. (NRAMP > 8)) THEN
        print *, "Value for NRAMP not allowed"
        IF (NFOVER == 1) THEN
            print *, "Program will override and use NRAMP = 0"
            NRAMP = 0
        ELSE
            stop
        ENDIF
    ENDIF

    READ(15,80) GMSG
    READ(GMSG,*) G
    IF ((ICS == 1) .AND. (G /= 9.81d0)) THEN
        IF ((NCOR == 1) .OR. (NTIP == 1)) THEN
            print *, "G not consistent with ICS=1"
            print *, "in conjunction with NTIP=1 and/or NCOR=1"
            IF(NFOVER == 1) THEN
                print *, "Program will override and set G=9.81"
                print *, "check to see that all input has SI units"
                G = 9.81d0
            ELSE
                stop
            ENDIF
        ENDIF
    ENDIF

    IF ((ICS == 2) .AND. (G /= 9.81d0)) THEN
        print *, "G not consistent with ICS=2"
        IF(NFOVER == 1) THEN
            print *, "Program will override and set G = 9.81 m/sec*sec"
            print *, "check to see that all input has SI units"
            print *, "execution will continue"
            G = 9.81d0
        ELSE
            stop
        ENDIF
    ENDIF

!     jgf47.11 Allow the limits of time varying tau0 to be read in from
!     the fort.15 file
    READ(15,80) TAU0MSG
    READ(TAU0MSG,*) TAU0
    if ( abs((int(tau0)-tau0-0.1d0)) < 1d-6 ) then
        outputTau0 = .TRUE. 
    else
        outputTau0 = .FALSE. 
    endif
    IF ( (TAU0 <= -5.d0) .AND. (TAU0 > -6.d0) ) THEN
        READ(15,80) TAU0LIMMSG
        READ(TAU0LIMMSG,*) Tau0FullDomainMin, Tau0FullDomainMax
    ENDIF


    READ(15,80) DTMSG
    READ(DTMSG,*) DT

    READ(15,80) STATMSG
    READ(STATMSG,*) STATIM

    READ(15,80) REFTMSG
    READ(REFTMSG,*) REFTIM

!--If wind stress and surface pressures are applied process this.

    IF((NWS == 0) .AND. (NRS >= 1)) READ(15,*) RSTIMMSG ! sb46.28sb03
    IF((NWS == 1) .AND. (NRS >= 1)) READ(15,*) RSTIMMSG ! sb46.28sb03

    IF(NWS == 3) THEN
        READ(15,80) WSMSG1
        READ(15,80) WSMSG2
    !....    TCM V49.64.01 ADDITIONS FOR ICE
        IF((NCICE == 0) .AND. (NRS == 0)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC
        ELSEIF ((NCICE == 0) .AND. (NRS >= 1)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,RSTIMINC
        ELSEIF((NCICE >= 1) .AND. (NRS == 0)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,CICE_TIMINC
        ELSEIF ((NCICE >= 1) .AND. (NRS >= 1)) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
        ENDIF
    ENDIF

!     jgfdebug46.02 Added NWS=45
!     jgf46.02 Added NWS=8
!     jgf46.16 Merged:
!     cf & cm Added NWS=9: asymmetric hurricane wind model
!     rjw Added NWS=19: asymmetric hurricane wind model v2.0
!     jie Added NWS=20: generalized asymmetric vortex model
!     sb46.28sb01 added NWS=12: OWI format

!....TCM V49.64.01 ADDITIONS FOR ICE, broke nws=2 out
    IF(ABS(NWS) == 2) THEN
        READ(15,80) WSMSG1  !tcm v51.17 -- added to correct bug
        IF(NRS == 0) READ(WSMSG1,*) WTIMINC
        IF(NRS >= 1) READ(WSMSG1,*) WTIMINC,RSTIMINC ! sb46.28sb03
    ENDIF
         
    IF((ABS(NWS) == 4) .OR. (ABS(NWS) == 5) &
     .OR. (ABS(NWS) == 45) .OR. (ABS(NWS) == 8) .OR. &
    (ABS(NWS) == 9) .OR. (ABS(NWS) == 12) &
     .OR. (ABS(NWS) == 15)  & !jgf50.38.05: Added NWS=15 for HWind.
     .OR. (ABS(NWS) == 20)  & !jie: Added NWS=20 for GAHM
     .OR. (ABS(NWS) == 19) .OR. (ABS(NWS) == 16)) THEN
    READ(15,80) WSMSG1
    IF((NCICE == 0) .AND. (NRS == 0)) READ(WSMSG1,*) WTIMINC
    IF((NCICE == 0) .AND. &
    ((NRS == 1) .OR. (NRS == 2) .OR. (NRS == 4))) READ(WSMSG1,*) WTIMINC,RSTIMINC
    IF((NCICE == 0) .AND. (NRS == 3)) THEN  !Casey 090825: Fix for NRS = 3.
        IF(ABS(NWS) == 8)THEN
            READ(WSMSG1,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
        ELSE
            READ(WSMSG1,*) WTIMINC,RSTIMINC
        ENDIF
    ENDIF
    IF((NCICE >= 1) .AND. (NRS == 0)) READ(WSMSG1,*) WTIMINC,CICE_TIMINC
    IF((NCICE >= 1) .AND. ((NRS == 1) .OR. (NRS == 2) .OR. (NRS == 4))) &
    READ(WSMSG1,*) WTIMINC,RSTIMINC,CICE_TIMINC
    IF((NCICE >= 1) .AND. (NRS == 3)) THEN  !Casey 090825: Fix for NRS = 3.
        IF(ABS(NWS) == 8)THEN
            READ(WSMSG1,*) IDUM,IDUM,IDUM,IDUM,IDUM, &
            RDUM,RSTIMINC,CICE_TIMINC
        ELSE
            READ(WSMSG1,*) WTIMINC,RSTIMINC,CICE_TIMINC
        ENDIF
    ENDIF
ENDIF

!     jgf49.0804 Added NWS29 VortexOWI -- don't think we need to
!     parse out the WTIMINC or other parameters.
    IF ((ABS(NWS)) == 29) THEN
        READ(15,80) WSMSG1
    ENDIF

    IF(ABS(NWS) == 6) THEN
        READ(15,80) WSMSG1
    !..... tcm v49.64.01 additions for ICE
        IF((NCICE == 0) .AND. (NRS == 0)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC
        ELSEIF ((NCICE == 0) .AND. (NRS >= 1)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,RSTIMINC
        ELSEIF((NCICE >= 1) .AND. (NRS == 0)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,CICE_TIMINC
        ELSEIF ((NCICE >= 1) .AND. (NRS >= 1)) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, &
            WLONINC,WTIMINC,RSTIMINC,CICE_TIMINC
        ENDIF
    ENDIF

    READ(15,80) RNDAYMSG
    READ(RNDAYMSG,*) RNDAY

    READ(15,80) DRAMPMSG
    READ(DRAMPMSG,*) DRAMP

    READ(15,80) COEFMSG
    READ(COEFMSG,*) A00,B00,C00

    READ(15,80) H0MSG
    IF (NOLIFA /= 2) THEN
        READ(H0MSG,*) H0
    ELSE
        READ(H0MSG,*) H0,NODEDRYMIN,NODEWETRMP,VELMIN
    ENDIF

    READ(15,80) SLMSG
    READ(SLMSG,*) SLAM0,SFEA0

    SL0=DEG2RAD*SLAM0
    SF0=DEG2RAD*SFEA0

    DO I = 1,NNODG
        SLAM(I) = X(I)
        SFEA(I) = Y(I)
    ENDDO

!      WRITE(6,*)"SLAM0,SFEA0", SLAM0,SFEA0

!--If ICS = 2 then apply CPP projection

    IF (ICS == 2) THEN
        DO I = 1,NNODG
            SL1(I) = DEG2RAD*SLAM(I)
            SF1(I) = DEG2RAD*SFEA(I)
            call cpp(slam(i),sfea(i),sl1(i),sf1(i),sl0,sf0)
        ENDDO
    ELSE
    ! jgf: Need these to be populated for ICS=1 for NetCDF file
        SL1 = SLAM
        SF1 = SFEA
    ENDIF

    READ(15,80) TAUMSG
    IF (NOLIBF == 0) THEN
        READ(TAUMSG,*) TAU
    ELSEIF (NOLIBF == 1) THEN
        READ(TAUMSG,*) CF
    ELSEIF (NOLIBF == 2) THEN
        READ(TAUMSG,*) CF,HBREAK,FTHETA,FGAMMA
    ENDIF

    READ(15,80) ESLMSG
    IF (IM == 10) THEN
        READ(ESLMSG,*) ESLM,ESLC
    ELSE
        READ(ESLMSG,*) ESLM
    ENDIF

    READ(15,80) CORIMSG
    READ(CORIMSG,*) CORI

    READ(15,80) NTIFMSG
    READ(NTIFMSG,*) NTIF
    IF (NTIF > MNTIF) THEN
        print *, "NTIF = ",NTIF, " exceeds parameter MNTIF = ",MNTIF
        stop
    ENDIF

    DO I=1,NTIF
        READ(15,80)  TIPOTAG(I)
        READ(15,80)  TPKMSG(I)
        READ(TPKMSG(I),*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
    ENDDO

    IF (((NTIP == 0) .AND. (NTIF /= 0)) .OR. ((NTIP /= 0) .AND. &
    (NTIF == 0))) THEN
        print *, "NTIF and NTIP are not consistent"
        IF(NFOVER == 1) THEN
            print *, "Program will reset NTIP = 0 and continue"
            NTIP = 0
        ELSE
            stop
        ENDIF
    ENDIF

    READ(15,80) NBFRMSG
    READ(NBFRMSG,*) NBFR
    IF (NBFR > MNBFR) THEN
        print *, "NBFR = ",NBFR, " exceeds parameter MNBFR = ",MNBFR
        stop
    ENDIF

    DO I=1,NBFR
        READ(15,80) BOUNTAG(I)
        READ(15,80) AMIGMSG(I)
        READ(AMIGMSG(I),*) AMIG(I),FF(I),FACE(I)
    ENDDO

    DO I=1,NBFR
        READ(15,80) ALPHA1(I)
        DO J=1,NETA
            READ(15,80) EMOMSG(I,J)
            READ(EMOMSG(I,J),*) EMO(I,J),EFA(I,J)
        ENDDO
    ENDDO

    READ(15,80) ANGMSG
    READ(ANGMSG,*) ANGINN

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

! md Added NFFR=-1 so the river information does not have to start at
!    time=0.
    IF (NFLUXF == 1) THEN
        READ(15,80) NFFRMSG
        READ(NFFRMSG,*) NFFR
        IF (NFFR > MNFFR) THEN
            print *, "NFFR = ",NFFR, " exceeds parameter MNFFR = ",MNFFR
            stop
        ENDIF
    
        IF ((NFFR /= 0) .AND. (NFFR /= -1)) THEN
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
                    IF ((LBCODE(J) == 2) .OR. (LBCODE(J) == 12) &
                     .OR. (LBCODE(J) == 22) .OR. (LBCODE(J) == 52)) THEN
                        READ(15,80) QNMSG(I,J)
                        READ(QNMSG(I,J),*) QNAM(I,J),QNPH(I,J)
                    !     dbug             print *, "disc lbnode index = ",J
                    ENDIF
                ENDDO
            ENDDO

        ! jp 5/1/99  added to help localize the flow boundary nodes
        !  NBV    =  global node number of all boundary nodes
        !  NFLBN  =  number of flow boundary nodes
        !  FLBN   =  global node number of flow boundary nodes
        !  FLBNX  =  index of flow boundary nodes as per NBV
        
            JG = 1
            DO K = 1,NBOU
                DO I=1, NVELL(K)
                    INDX = NBVV(K,I)
                    NBV(JG) = INDX
                    JG = JG + 1
                ENDDO
            ENDDO
        
            NFLBN = 0
        ! jgf51.52.23: Changed the loop extent from NVEL to
        ! numLBCodeValues.
            DO J=1,numLBCodeValues
                IF ((LBCODE(J) == 2) .OR. (LBCODE(J) == 12) &
                 .OR. (LBCODE(J) == 22) .OR. (LBCODE(J) == 52)) THEN
                    NFLBN = NFLBN+1
                    FLBN(NFLBN) = NBV(J)
                    FLBNX(NFLBN) = J
                ENDIF
            ENDDO
        ELSE
            APERIODIC_FLOW_BC = .TRUE. !jgf45.09 need to break up a fort.20
        ENDIF
    ENDIF
!      PRINT *, "APERIODIC_FLOW_BC = ", APERIODIC_FLOW_BC


! bell...READING TIME VARYING BOUNDARY CONDITIONS
!........CHECK FOR WEIR FLUX BOUNDARIES
    NFLUXB = 0
    NFLUXIB = 0
    NFLUXIBP = 0
    DO I = 1,NVEL
        SELECT CASE(LBCODE(I))
        CASE(3,13,23)
        NFLUXB = 1
        CASE(4,24)
        NFLUXIB = 1
        CASE(5,25)
        NFLUXIBP = 1
        END SELECT
    ENDDO
    IF((NFLUXIB == 1) .OR. (NFLUXIBP == 1) .OR. (NFLUXB == 1) )THEN
        IF(USE_TVW)CALL PARSE_TIME_VARYING_WEIR_INFO()
    ENDIF


!--Read Elevation Recording Stations
!...
    READ(15,80) STAEMSG
    READ(STAEMSG,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
    IF (ABS(NOUTE) > 5) THEN
        print *, "The value of NOUTE is not allowed"
        stop
    ENDIF

    READ(15,80) NSTAEMSG
    READ(NSTAEMSG,*) NSTAE
! cm v51.20.03 -- addition for external elevation station file
    if (NSTAE < 0 ) then
        WRITE(*,*) "External File Used for Elevation Station Locations"
        USE_ELEV_STAT_FILE = .TRUE. 
        IOS_STATIONS = 0
        NSTAE = abs(NSTAE)  !reset to positive
        NSTAE2 = 0
        OPEN(unit=151,file='elev_stat.151',status='old',err=7690,iostat=ios_stations)
        READ(151,*) NSTAE2
        IF (ABS(NSTAE2) /= ABS(NSTAE)) THEN
            NSTAE = ABS(NSTAE2)  !RESET THE VALUE TO WHAT'S IN THE FILE
        ENDIF
        7690 IF (IOS_STATIONS /= 0) THEN
            WRITE(*,*) "Error in Reading Elevation Station File: elev_stat.151"
            stop  ! there is a stop here
        ENDIF
    ELSE
        WRITE(*,*) "Elevation Station Locations contained in fort.15"
    ENDIF
    IF (ABS(NSTAE) > MNSTAE) THEN
        print *, "NSTAE = ",abs(NSTAE), " exceeds parameter MNSTAE = ", &
        MNSTAE
        stop
    ENDIF
#ifdef ADCNETCDF

!     ==================================================================
!     START ADDITION BY MCF FOR NETCDF OPTION

!      IF(.NOT. ALLOCATED(STATNUMB))ALLOCATE(STATNUMB(abs(NSTAE),SNUMLEN))
    IF( .NOT. ALLOCATED(STATNAME))ALLOCATE(STATNAME(abs(NSTAE)))

!     END ADDITION BY MCF FOR NETCDF OPTION
!     ==================================================================
#endif
    DO I=1,ABS(NSTAE)
        IF (USE_ELEV_STAT_FILE ) then  !read from external file
            READ(151,'(A132)') STAELOC(I)
        ELSE
            READ(15,'(A132)') STAELOC(I)
        ENDIF
        IF(ICS == 1) THEN

        !     ==================================================================
        !     START ADDITION BY MCF FOR NETCDF OPTION
        
            IF((abs(NOUTE) /= 3) .AND. (abs(NOUTE) /= 5))THEN
                READ(STAELOC(I),*) XEL(I),YEL(I)
            ELSE
#ifdef ADCNETCDF
            !                READ(15,'(A132)') STATLINE
            !                call parse(STATLINE, LVAR)
                call parse(STAELOC(I), LVAR)
                XEL(I)=a2f(LVAR(1))
                YEL(I)=a2f(LVAR(2))
                STATNAME(I)=LVAR(3)
            !                WRITE(6,9911) XEL(I), YEL(I), STATNAME(I)
#endif
            ENDIF

        ELSE
            IF((abs(NOUTE) /= 3) .AND. (abs(NOUTE) /= 5))THEN

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
            CALL CPP(XEL(I),YEL(I),SLEL(i),SFEL(i),SL0,SF0)
        ENDIF
    ENDDO !loop over elevation stations
    if (USE_ELEV_STAT_FILE) CLOSE(151) !close external file if one was used
    9911 FORMAT(F12.3,2X,F12.3, 6X, A50)


!--For Each Elevation Station:
!  Find the Global Index of the element it lies in.

!...tcm v49.48.01 replace with new call to kdtsearch
!      CALL CoordToEle(XEL,YEL,NNSEG,NSTAE,
!     &     'Elevation recording station   ')
    call setup_kdt_search()
    call kdtsearch(XEL,YEL,NNSEG,abs(NSTAE), &
    'Elevation recording station   ')

!--Read Velocity Recording Stations

    READ(15,80) STAVMSG
    READ(STAVMSG,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
    IF (ABS(NOUTV) > 5) THEN
        print *, "Value for NOUTV is not allowable"
        stop
    ENDIF

    READ(15,80) NSTAVMSG
    READ(NSTAVMSG,*) NSTAV
! cm v51.20.03 -- addition for external velocity station location file vel_stat.151
    IF (NSTAV < 0 ) then
        WRITE(*,*) "External File Used for Velocity Station Locations"
        USE_VEL_STAT_FILE = .TRUE. 
        IOS_STATIONS = 0
        NSTAV = ABS(NSTAV) ! reset to positive value
        NSTAV2 = 0
        OPEN(unit=151,file='vel_stat.151',status='old',err=7691,iostat=ios_stations)
        READ(151,*) NSTAV2
        IF (ABS(NSTAV2) /= ABS(NSTAV)) THEN
            NSTAV = ABS(NSTAV2)  !RESET THE VALUE TO WHAT'S IN THE FILE
        ENDIF
        7691 IF (IOS_STATIONS /= 0) THEN
            WRITE(*,*) "Error in Reading Velocity Station File: vel_stat.151"
            stop  ! there is a stop here
        ENDIF
    ELSE
        WRITE(*,*) "Velocity Station Locations Contained in fort.15"
    ENDIF
          
    IF (ABS(NSTAV) > MNSTAV) THEN
        print *, "NSTAV = ",abs(NSTAV), " exceeds parameter MNSTAV = ", &
        MNSTAV
        stop
    ENDIF
#ifdef ADCNETCDF
!     ==================================================================
!     START ADDITION BY MCF FOR NETCDF OPTION

!      IF(.NOT. ALLOCATED(STATNUMBV))ALLOCATE(STATNUMBV(NSTAV,SNUMLEN))
    IF( .NOT. ALLOCATED(STATNAMEV))ALLOCATE(STATNAMEV(abs(NSTAV)))

!     END ADDITION BY MCF FOR NETCDF OPTION
!     ==================================================================

#endif

    DO I=1,ABS(NSTAV)
        if (USE_VEL_STAT_FILE ) then !read from external file
            READ(151,'(A132)') STAVLOC(I)
        else
            READ(15,'(A132)') STAVLOC(I)
        ENDIF
        IF (ICS == 1) then
        !     ==================================================================
        !     START ADDITION BY MCF FOR NETCDF OPTION
        
            IF((abs(NOUTV) /= 3) .AND. (abs(NOUTV) /= 5)) THEN
                READ(STAVLOC(I),*) XEV(I),YEV(I)
            ELSE
#ifdef ADCNETCDF
                call parse(STAVLOC(I), LVAR)
                XEV(I)=a2f(LVAR(1))
                YEV(I)=a2f(LVAR(2))
                STATNAMEV(I)=LVAR(3)
            !              WRITE(6,9911) XEV(I), YEV(I), STATNAMEV(I)
#endif
            ENDIF

        ELSE
                     
            IF((abs(NOUTV) /= 3) .AND. (abs(NOUTV) /= 5)) THEN
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
        
        !     END ADDITION BY MCF FOR NETCDF OPTION
        !     ==================================================================


            SLEV(I)=DEG2RAD*SLEV(I)
            SFEV(I)=DEG2RAD*SFEV(I)
            CALL CPP(XEV(I),YEV(I),SLEV(i),SFEV(i),SL0,SF0)
        ENDIF  !test for ICS
    ENDDO !loop over velocity stations
    IF (USE_VEL_STAT_FILE) CLOSE(151) !close external file if one was used

!--For Each Velocity Station:
!  Find the Global Index of the element it lies in.

!...tcm v49.48.01 replace with new call to kdtsearch
!      CALL CoordToEle(XEV,YEV,NNSVG,NSTAV,
!     &     'Velocity recording station    ')
    call kdtsearch(XEV,YEV,NNSVG,abs(NSTAV), &
    'Velocity recording station    ')


!     If Passive Transport is indicated, then read Concentration Station Info

    IF (C2D_PTrans .OR. C3D_PTrans) THEN
    
        READ(15,80) STACMSG
        READ(STACMSG,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
        IF (ABS(NOUTC) > 2) THEN
            print *, "Value of NOUTC is not allowable"
            stop
        ENDIF
    
        READ(15,80) NSTACMSG
        READ(NSTACMSG,*) NSTAC
    ! cm v51.20.03 -- addition for external Concentration station file
        if (NSTAC < 0 ) then
            WRITE(*,*) "External File Used for Concentration Station Locations"
            USE_CONC_STAT_FILE = .TRUE. 
            IOS_STATIONS = 0
            NSTAC = ABS(NSTAC) !RESET TO POSITIVE
            NSTAC2 = 0 !INITIALIZE
            OPEN(unit=151,file='conc_stat.151',status='old',err=7692,iostat=ios_stations)
            READ(151,*) NSTAC2
            IF (ABS(NSTAC2) /= ABS(NSTAC)) THEN
                NSTAC = ABS(NSTAC2)  !RESET THE VALUE TO WHAT'S IN THE FILE
            ENDIF
            7692 IF (IOS_STATIONS /= 0) THEN
                WRITE(*,*) "Error in Reading Concentration Station File: conc_stat.151"
                STOP  ! THERE IS A STOP HERE
            ENDIF
        ELSE
            WRITE(*,*) "Concentration Station Locations Contained in fort.15"
        ENDIF
                
        IF (ABS(NSTAC) > MNSTAC) THEN
            print *, "NSTAC = ",abs(NSTAC), " exceeds parameter MNSTAC = ", &
            MNSTAC
            stop
        ENDIF
    
        DO I=1,ABS(NSTAC)
            IF (USE_CONC_STAT_FILE) THEN !READ FROM EXTERNAL FILE
                READ(151,'(A132)') STACLOC(I)
            ELSE
                READ(15,'(A132)') STACLOC(I)
            ENDIF
            IF (ICS == 1) THEN
                READ(STACLOC(I),*) XEC(I),YEC(I)
            ELSE
                READ(STACLOC(I),*) SLEC(I),SFEC(I)
                SLEC(I)=DEG2RAD*SLEC(I)
                SFEC(I)=DEG2RAD*SFEC(I)
                CALL CPP(XEC(I),YEC(I),SLEC(i),SFEC(i),SL0,SF0)
            ENDIF !test for ICS
        ENDDO  !loop over stations
        IF (USE_CONC_STAT_FILE) CLOSE(151) !CLOSE EXTERNAL FILE IF ONE WAS USED
    
    !--For Each Concentration Recording Station:
    !  Find the Global Index of the element it lies in.
    
    !...tcm v49.48.01 replace with new call to kdtsearch
    !         CALL CoordToEle(XEC,YEC,NNSCG,NSTAC,
    !     &        'Concentration station         ')
        call kdtsearch(XEC,YEC,NNSCG,abs(NSTAC), &
        'Concentration station         ')
    
    ENDIF

!--If NWS <> 0 , then read Meteorlogical Station Info

    NOUTM = 0
    IF (NWS /= 0) THEN
    
        READ(15,80) STAMMSG
        READ(STAMMSG,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
        IF (ABS(NOUTM) > 5) THEN
            print *, "Value of NOUTM is not allowable"
            stop
        ENDIF
    
        READ(15,80) NSTAMMSG
        READ(NSTAMMSG,*) NSTAM
    ! cm v51.20.03 -- addition for external met station file
        IF (NSTAM < 0 ) THEN
            WRITE(*,*) "External File Used for MET Station Locations"
            USE_MET_STAT_FILE = .TRUE. 
            IOS_STATIONS = 0
            NSTAM = ABS(NSTAM) !RESET TO POSITIVE VALUE
            NSTAM2 = 0 !INITIALIZE
            OPEN(unit=151,file='met_stat.151',status='old', &
            err=7693,iostat=ios_stations)
            READ(151,*) NSTAM2
            IF (ABS(NSTAM2) /= ABS(NSTAM)) THEN
                NSTAM = ABS(NSTAM2)  !RESET THE VALUE TO WHAT'S IN THE FILE
            ENDIF
            7693 IF (IOS_STATIONS /= 0) THEN
                WRITE(*,*) "Error in Reading MET Station File: met_stat.151"
                STOP  ! THERE IS A STOP HERE
            ENDIF
        ELSE
            WRITE(*,*) "MET Station Locations Contained in fort.15"
        ENDIF

        IF (ABS(NSTAM) > MNSTAM) THEN
            print *, "NSTAM = ",abs(NSTAM), " exceeds parameter MNSTAM = ", &
            MNSTAM
            stop
        ENDIF
#ifdef ADCNETCDF
    !     ==================================================================
    !     START ADDITION BY MCF FOR NETCDF OPTION
    

    !         IF(.NOT.ALLOCATED(STATNUMBM))ALLOCATE(STATNUMBM(NSTAM,SNUMLEN))
        IF( .NOT. ALLOCATED(STATNAMEM))ALLOCATE(STATNAMEM(abs(NSTAM)))
    
    !     END ADDITION BY MCF FOR NETCDF OPTION
    !     ==================================================================
    
#endif
        DO I=1,ABS(NSTAM)
            IF (USE_MET_STAT_FILE) THEN !READ FROM EXTERNAL FILE
                READ(151,'(A132)') STAMLOC(I)
            ELSE
                READ(15,'(A132)') STAMLOC(I)
            ENDIF
        
            IF (ICS == 1) THEN

            !     ==================================================================
            !     START ADDITION BY MCF FOR NETCDF OPTION
            
                IF((ABS(NOUTM) /= 3) .AND. (ABS(NOUTM) /= 5)) THEN
                    READ(STAMLOC(I),*) XEM(I),YEM(I)
                ELSE
#ifdef ADCNETCDF
                    call parse(STAMLOC(I), LVAR)
                    XEM(I)=a2f(LVAR(1))
                    YEM(I)=a2f(LVAR(2))
                    STATNAMEM(I)=LVAR(3)
                !                 WRITE(6,9911) XEM(I), YEM(I), STATNAMEM(I)
#endif
                ENDIF
            ELSE
                IF((ABS(NOUTM) /= 3) .AND. (ABS(NOUTM) /= 5)) THEN

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
            
            !     END ADDITION BY MCF FOR NETCDF OPTION
            !     ==================================================================

                SLEM(I)=DEG2RAD*SLEM(I)
                SFEM(I)=DEG2RAD*SFEM(I)
                CALL CPP(XEM(I),YEM(I),SLEM(i),SFEM(i),SL0,SF0)
            ENDIF  !test for ICS
        ENDDO  !loop over stations
        IF (USE_MET_STAT_FILE) CLOSE(151) !CLOSE EXTERNAL FILE IF ONE WAS USED

    
    !--For Each Meterological Recording Station:
    !  Find the Global Index of the element it lies in.
    
    !...tcm v49.48.01 replace with new call to kdtsearch
    !         CALL CoordToEle(XEM,YEM,NNSMG,NSTAM,
    !     &        'Meteorological station        ')
        call kdtsearch(XEM,YEM,NNSMG,abs(NSTAM), &
        'Meteorological station        ')
    
    ENDIF


!--Read Global Elevation Data Output

    READ(15,80) OUTGEMSG
    READ(OUTGEMSG,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
    IF (ABS(NOUTGE) > 7) THEN
        print *, "NOUTGE does not have an allowable value"
        stop
    ENDIF

!--Read Global Velocity Data Output

    READ(15,80) OUTGVMSG
    READ(OUTGVMSG,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
    IF (ABS(NOUTGV) > 7) THEN
        print *, "NOUTGV does not have an allowable value"
        stop
    ENDIF

!     If Passive Transport is indicated, read Global Concentration Data Output

    IF (C2D_PTrans .OR. C3D_PTrans) THEN
        READ(15,80) OUTGCMSG
        READ(OUTGCMSG,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
        IF (ABS(NOUTGC) > 2) THEN
            print *, "NOUTGC does not have an allowable value"
            stop
        ENDIF
    ENDIF

    IF (NWS /= 0) THEN
        READ(15,80) OUTGWMSG
        READ(OUTGWMSG,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
    ! Casey 090302: Changed the following line for ABS(NOUTGW).EQ.4
        IF (ABS(NOUTGW) > 7) THEN
            print *, "NOUTGW does not have an allowable value"
            stop
        ENDIF
    ENDIF

!--Read Harmonic Analysis Data

    READ(15,80) HARFRMSG
    READ(HARFRMSG,*) NHARFR
    MNHARF = NHARFR
    IF (NHARFR == 0) MNHARF = 1
    CALL ALLOC_HA()
    CALL ALLOC_MAIN14()
    ALLOCATE(HAFREMSG(MNHARF))

    DO I=1,NHARFR
        READ(15,'(A10)') NAMEFR(I)
        READ(15,80) HAFREMSG(I)
        READ(HAFREMSG(I),*) HAFREQ(I),HAFF(I),HAFACE(I)
    ENDDO

    READ(15,80) HARPARMSG
    READ(HARPARMSG,*) THAS,THAF,NHAINC,FMV
    READ(15,80) OUTHARMSG
    READ(OUTHARMSG,*) NHASE,NHASV,NHAGE,NHAGV
    IF ((NHASE < 0) .OR. (NHASE > 1)) THEN
        print *, "NHASE does not have an allowable value"
        IF (NFOVER == 1) THEN
            print *, "Program will override an reset NHASE=0 "
            NHASE = 0
        ELSE
            stop
        ENDIF
    ENDIF

    IF ((NHASV < 0) .OR. (NHASV > 1)) THEN
        print *, "NHASV does not have an allowable value"
        IF (NFOVER == 1) THEN
            print *, "Program will override an reset NHASV=0 "
            NHASV = 0
        ELSE
            stop
        ENDIF
    ENDIF

    IF ((NHAGE < 0) .OR. (NHAGE > 1)) THEN
        print *, "NHAGE does not have an allowable value"
        IF (NFOVER == 1) THEN
            print *, "Program will override an reset NHAGE=0 "
            NHAGE = 0
        ELSE
            stop
        ENDIF
    ENDIF

    IF ((NHAGV < 0) .OR. (NHAGV > 1)) THEN
        print *, "NHAGV does not have an allowable value"
        IF (NFOVER == 1) THEN
            print *, "Program will override an reset NHAGV=0 "
            NHAGV = 0
        ELSE
            stop
        ENDIF
    ENDIF
    IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
    IF(IHARIND > 0) IHARIND=1
    IF ((FMV > 0.) .AND. (NHARFR > 0) .AND. (C2DDI)) CHARMV = .TRUE. 

!--Read Hot Start Data

!     jgf45.07 added undocumented option to allow ADCIRC to stop after writing
!     hot start file. This is used to test hot starting capabilities.
    READ(15,80) HSTARMSG
    READ(HSTARMSG,*) NHSTAR,NHSINC
!     jgf49.41: Added new options for hotstart file format indicator: 2, 3,
!     367, 368. Not sure we need to check the validity of NHSTAR values in
!     adcprep though.
    SELECT CASE(NHSTAR)  !tcm v51.26 mods for binary time-stamped hotstart file generation NHSTAR=-1
    CASE(0,-1,1,2,3,5,67,68,367,368,567,568)
! do nothing, these are allowable values
    CASE DEFAULT
    print *, "WARNING: NHSTAR=",NHSTAR," is not valid."
    IF (NFOVER == 1) THEN
        print *, "WARNING: Program will override and reset NHSTAR=0."
        NHSTAR = 0
    ELSE
        print *,"ERROR: Halting adcprep due to invalid NHSTAR value."
        stop
    ENDIF
    END SELECT

!--Read Solver Data

    READ(15,80) SOLVMSG
    READ(SOLVMSG,*) ITITER,ISLDIA,CONVCR,ITMAX
! jgf49.82: Removed value check on ITITER ... this parameter
! is deprecated and its value is now arbitrary.


!--Read in 3D data

    IF(C3DVS) THEN
        CALL READ15_3DVS()
    !     ELSEIF(C3DDSS) THEN
    !       CALL READ15_3DDSS()
    ENDIF

    IF((useNetCDF.eqv. .TRUE. ) .AND. (NETCDF_AVAIL.eqv. .FALSE. )) THEN
        CALL allMessage(ERROR, &
        "NetCDF was selected in the fort.15 file.")
        CALL allMessage(ERROR, &
        "NetCDF was not compiled into this executable program.")
        stop
    ENDIF
#ifdef ADCNETCDF
!     ==================================================================
!     START ADDITION BY MCF FOR NETCDF OPTION for NETCDF I/O - 6/20/07

!     --------------------------------------------------------------------
!     IF output is in netCDF format read global attributes for netCDF file
!     --------------------------------------------------------------------

    IF(useNetCDF.eqv. .TRUE. ) THEN
        IF (NOUTGC == 3 .OR.                                & ! fort.93
        NOUTC == 3) then                                 ! fort.91
        CALL allMessage(ERROR, &
        "NetCDF output not available for concentration output.")
        STOP
    ENDIF

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

    999 IF (ios /= 0) THEN
        CALL allMessage(ERROR, &
        "Could not read metadata for NetCDF at end of fort.15.")
        STOP
    ENDIF

ENDIF
#endif

!--Close Global Run Info file

    CLOSE(15)
    call memory_status()

#ifdef READ_GLOBAL_TRACE
    call screen(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    RETURN
    80 FORMAT(A80)
    END SUBROUTINE READ15


!----------------------------------------------------------------------
!     Reads the 3DVS portion of the global ADCIRC Input Data File for
!     data decomposition program ADCPP.  This version compatible with
!     ADCIRC_v41.11a
!     tjc  6/24/02
!---------------------------------------------------------------------
    SUBROUTINE READ15_3DVS()
    USE PRE_GLOBAL
    USE GLOBAL, ONLY : parse, a2f, SNAMLEN, deg2rad, day2sec
    USE GLOBAL_3DVS, ONLY : STATNAMED, STATNAMEV3D, STATNAMET, SIGMA, &
    FEGRIDS
    USE MESH, only : cpp
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,K,M,N                                        ! loop counters
    CHARACTER(50) LVAR(3)
    REAL(8) :: STATIME
    REAL(sz) :: HH1 ! domain averaged depth
    INTEGER :: NH   ! horizontal node counter

!     jgf49.43.19: Moved memory allocation for number of 3D stations
!     on each subdomain to here and made it unconditional.
!     The NNSTA3D*P arrays will still be needed later b/c they are
!     written to the fort.18 message passing file.
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


    350 FORMAT(//,2X,'***** INVALID INPUT IN THE PRIMARY VERTICAL INPUT', &
    ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

!     jgf45.10 removed IDIAG
!... SPECIFY WHETHER A BAROTROPIC OR BAROCLINIC RUN
    READ(15,80) IDENMSG
    READ(IDENMSG,*) IDEN
    IF((IDEN > 4) .OR. (IDEN < -4)) THEN
        WRITE(*,*) "ERROR: IDEN=",IDEN
        WRITE(*,423)
        WRITE(*,350)
        423 FORMAT(/,2X,' IDEN must be an integer between -4 and 4.')
        STOP
    ENDIF

!     jgf45.12 Set flag to read in unit 11 file (initial density) if
!     necessary.
    IF ( IDEN /= 0 ) THEN
        CBaroclinic = .TRUE. 
    ENDIF
!     jgf45.12 Set flag to read in additional parameters from unit 15
!     file that will be used when solving for transport of salinity,
!     temperature, etc in prognostic baroclinic simulation.
    IF ( IDEN > 0 ) THEN
        C3D_BTrans = .TRUE. 
    ENDIF

!... READ IN THE TYPE OF BOTTOM BOUNDARY CONDITION AND THE SLIP COEFFICIENTS

    READ(15,80) SLIPMSG
    READ(SLIPMSG,*) ISLIP,KP
    IF((ISLIP < 0) .OR. (ISLIP > 3)) THEN
        WRITE(6,350)
        WRITE(6,360)
        360 FORMAT(/,2X,'    THE SLIP CODE MUST = 0,1,2,OR 3.')
        STOP
    ENDIF

!... READ IN THE SURFACE AND BOTTOM ROUGHNESSES

    READ(15,80) Z0MSG
    READ(Z0MSG,*) Z0S, Z0B

!... READ IN THE TIME STEPPING COEFFICIENTS

    READ(15,80) ALPMSG
    READ(ALPMSG,*) ALP1,ALP2,ALP3

!... READ IN IGC & NFEN: F.E. GRID CODE &# NODES IN F.E. GRID

    READ(15,80) FEMSG
    READ(FEMSG,*) IGC,NFEN

!     jgf45.12 Add code to read in nondimensional thicknesses of
!     vertical layers.
    ALLOCATE ( Sigma(NFEN) )
    nbytes = 8*nfen
    call memory_alloc(nbytes)
    IF(IGC == 0) THEN
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



!... SPECIFY TYPE OF EDDY VISCOSITY PROFILE
    READ(15,80) EVCMSG
    READ(EVCMSG,*) IEVC,EVMIN,EVCON
    IF((IEVC /= 0 ) .AND. (IEVC /= 1 ) .AND. &
    (IEVC /= 10) .AND. (IEVC /= 11) .AND. &
    (IEVC /= 20) .AND. (IEVC /= 21) .AND. &
    (IEVC /= 22) .AND. (IEVC /= 23) .AND. &
    (IEVC /= 30) .AND. (IEVC /= 31) .AND. (IEVC /= 32) .AND. &
    (IEVC /= 33) .AND. &
    (IEVC /= 40) .AND. (IEVC /= 41) .AND. (IEVC /= 42) .AND. &
    (IEVC /= 43) .AND. &
    (IEVC /= 50) .AND. (IEVC /= 51)) THEN
        WRITE(*,350)
        WRITE(*,411)
        411 FORMAT(/,2X,'    IEVC MUST BE 0,1,10,11,20,21,22,23,', &
        '30,31,32,33,40,41,42,43,50,51')
        STOP
    ENDIF
    IF((IEVC == 50) .OR. (IEVC == 51)) THEN
        READ(15,80) THETAMSG
        READ(THETAMSG,*) THETA1,THETA2
    ENDIF
!     jgf45.12 Add code to read in vertical eddy viscosity profile.
    IF(IEVC == 0) THEN
        ALLOCATE ( EVTot(NFEN) )
        nbytes = 8*nfen
        call memory_alloc(nbytes)
        DO N=1,NFEN
            READ(15,*) EVTot(N)
        ENDDO
    ENDIF
!     -----------------------------------------------------
!.... Station 3D Density, Temperature, Salinity output
!     jgf45.11 switched from node numbers to coordinates

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
    IF( .NOT. ALLOCATED(STATNAMED))ALLOCATE(STATNAMED(NSTA3DD))
!     ALLOCATE(STATNAMED(NSTA3DD))
    nbytes = 4*SNAMLEN*NSTA3DD
    call memory_alloc(nbytes)
!  kmd48.33bc changed to -2 and I3DSD not equal to 0
    IF(ABS(I3DSD) > 5) THEN
        WRITE(*,350)
        WRITE(*,511)
        511 FORMAT(/,2X,'    I3DSD MUST BE -5,-3,-2,-1,0,1,2, 3, or 5')
        STOP
    ENDIF
    IF(NSTA3DD /= 0) THEN  ! kmd : changed to match 2D option
        DO I=1,NSTA3DD
            READ(15,80) STA3DDLOC(I)
            call parse(STA3DDLOC(I), LVAR)
            STATNAMED(I)=LVAR(3)
            IF (ICS == 1) THEN
                X3DD(I)=a2f(LVAR(1))
                Y3DD(I)=a2f(LVAR(2))
            ELSE
                SL3DD(I)=a2f(LVAR(1))
                SF3DD(I)=a2f(LVAR(2))
                SL3DD(I)=DEG2RAD*SL3DD(I)
                SF3DD(I)=DEG2RAD*SF3DD(I)
                call cpp(x3dd(i),y3dd(i),sl3dd(i),sf3dd(i),sl0,sf0)
            ENDIF
        ENDDO
    
    !     jgf45.11 For each 3D density station, find the full-domain index of
    !     the corresponding element.
    !...tcm v49.48.01 replace with new call to kdtsearch
    !         CALL CoordToEle(X3DD,Y3DD,NNS3DDG,NSTA3DD,
    !     &        '3D density recording station  ')
        call kdtsearch(X3DD,Y3DD,NNS3DDG,NSTA3DD, &
        '3D density recording station  ')

    ENDIF

!     -----------------------------------------------------
!.... Station 3D Velocity output
!     jgf45.11 switched from node numbers to coordinates
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
    IF( .NOT. ALLOCATED(STATNAMEV3D))ALLOCATE(STATNAMEV3D(NSTA3DV))
!      ALLOCATE(STATNAMEV3D(NSTA3DV))
    nbytes = 4*SNAMLEN*NSTA3DV
    call memory_alloc(nbytes)
    ALLOCATE(SL3DV(NSTA3DV),SF3DV(NSTA3DV))
    nbytes = 16*nsta3dv
    call memory_alloc(nbytes)

!  kmd48.33 changed to -2
    IF(ABS(I3DSV) > 5) THEN
        WRITE(*,350)
        WRITE(*,512)
        512 FORMAT(/,2X,'    I3DSV MUST BE -5,-3,-2,-1,0,1,2,3 or 5')
        STOP
    ENDIF
    IF(NSTA3DV /= 0) THEN ! kmd : Changed to match 2D option
        DO I=1,NSTA3DV
            READ(15,80) STA3DVLOC(I)
            call parse(STA3DVLOC(I), LVAR)
            STATNAMEV3D(I)=LVAR(3)
            IF (ICS == 1) THEN
                X3DV(I)=a2f(LVAR(1))
                Y3DV(I)=a2f(LVAR(2))
            ELSE
                SL3DV(I)=a2f(LVAR(1))
                SF3DV(I)=a2f(LVAR(2))
                SL3DV(I)=DEG2RAD*SL3DV(I)
                SF3DV(I)=DEG2RAD*SF3DV(I)
                call cpp(x3dv(i),y3dv(i),sl3dv(i),sf3dv(i),sl0,sf0)
            ENDIF
        ENDDO
    !     jgf45.11 For each 3D velocity station, find the full-domain index of
    !     the corresponding element.
    !...tcm v49.48.01 replace with new call to kdtsearch
    !         CALL CoordToEle(X3DV,Y3DV,NNS3DVG,NSTA3DV,
    !     &        '3D velocity recording station ')
        call kdtsearch(X3DV,Y3DV,NNS3DVG,NSTA3DV, &
        '3D velocity recording station ')
    ENDIF
!     -----------------------------------------------------
!.... Station 3D turbulence output
!     jgf45.11 switched from node numbers to coordinates
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
    IF( .NOT. ALLOCATED(STATNAMET))ALLOCATE(STATNAMET(NSTA3DT))
!      ALLOCATE(STATNAMET(NSTA3DT))
    nbytes = 4*SNAMLEN*NSTA3DT
    call memory_alloc(nbytes)
    ALLOCATE(SL3DT(NSTA3DT),SF3DT(NSTA3DT))
    nbytes = 16*nsta3dv
    call memory_alloc(nbytes)
!  kmd48.33 changed to -2
    IF(ABS(I3DST) > 5) THEN
        WRITE(*,350)
        WRITE(*,513)
        513 FORMAT(/,2X,'    I3DST MUST BE -5,-3,-2,-1,0,1,2,3 or 5')
        STOP
    ENDIF
    IF(NSTA3DT /= 0) THEN ! kmd - Changed to match 2D option
        DO I=1,NSTA3DT
            READ(15,80) STA3DTLOC(I)
            call parse(STA3DTLOC(I), LVAR)
            STATNAMET(I)=LVAR(3)
            IF (ICS == 1) THEN
                X3DT(I)=a2f(LVAR(1))
                Y3DT(I)=a2f(LVAR(2))
            ELSE
                SL3DT(I)=a2f(LVAR(1))
                SF3DT(I)=a2f(LVAR(2))
                SL3DT(I)=DEG2RAD*SL3DT(I)
                SF3DT(I)=DEG2RAD*SF3DT(I)
                call cpp(x3dt(i),y3dt(i),sl3dt(i),sf3dt(i),sl0,sf0)
            ENDIF
        ENDDO

    !     jgf45.11 For each 3D turbulence station, find the full-domain index of
    !     the corresponding element.

        ALLOCATE(NNS3DTG(NSTA3DT))
        nbytes = 4*nsta3dt
        call memory_alloc(nbytes)
    !...tcm v49.48.01 replace with new call to kdtsearch
    !         CALL CoordToEle(X3DT,Y3DT,NNS3DTG,NSTA3DT,
    !     &        '3D turbulence rec. station    ')
        CALL KDTSEARCH(X3DT,Y3DT,NNS3DTG,NSTA3DT, &
        '3D turbulence rec. station    ')

    ENDIF
!....  GLOBAL 3D DENSITY, TEMPERATURE, SALINITY OUTPUT

    READ(15,80) DGDMSG
    READ(DGDMSG,*) I3DGD,TO3DGDS,TO3DGDF,NSPO3DGD
!   kmd48.33bc changed to -2
    IF(ABS(I3DGD) > 5) THEN  !tcm 20110511 should this be > 3 instead of > 2 ?
        WRITE(*,350)
        WRITE(*,514)
        514 FORMAT(/,2X,'    I3DGD MUST BE -3,-2,-1,0,1,2, or 3')
        STOP
    ENDIF

!....  GLOBAL 3D VELOCITY OUTPUT

    READ(15,80) DGVMSG
    READ(DGVMSG,*) I3DGV,TO3DGVS,TO3DGVF,NSPO3DGV
!   kmd48.33bc changed to -2
    IF(ABS(I3DGV) > 5) THEN
        WRITE(*,350)
        WRITE(*,515)
        515 FORMAT(/,2X,'    I3DGV MUST BE -3,-2,-1,0,1,2 or 3')
        STOP
    ENDIF

!....  GLOBAL 3D TURBULENCE OUTPUT

    READ(15,80) DGTMSG
    READ(DGTMSG,*) I3DGT,TO3DGTS,TO3DGTF,NSPO3DGT
!   kmd48.33bc changed to -2
    IF(ABS(I3DGT) > 5) THEN
        WRITE(*,350)
        WRITE(16,516)
        516 FORMAT(/,2X,'    I3DGT MUST BE -3,-2,-1,0,1,2, or 3')
        STOP
    ENDIF

!   kmd48.33bc add in information for 3D boundary condition
    IF (CBAROCLINIC) THEN
        READ(15,80) RESBCFLAGMSG
        READ(RESBCFLAGMSG,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP
        IF (RES_BC_FLAG /= 0) THEN
            IF (NOPE > 0) THEN
                READ(15,80) BCTIMEMSG
                READ(15,80) BCSTATMSG
            END IF
            IF (BCFLAG_TEMP /= 0) THEN
                READ(15,80) TBCTIMEMSG
            END IF
        END IF
        IF ((RES_BC_FLAG < 0) .OR. (RES_BC_FLAG == 1))  THEN
            IF (NOPE > 0) THEN
                READ(BCTIMEMSG,*) RBCTIMEINC
                READ(BCSTATMSG,*) BCSTATIM
            END IF
        ELSE IF (RES_BC_FLAG == 2) THEN
            IF (NOPE > 0) THEN
                READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC
                READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM
            END IF
        ELSE IF (RES_BC_FLAG == 3) THEN
            IF (NOPE > 0) THEN
                READ(BCTIMEMSG,*) RBCTIMEINC, TBCTIMEINC
                READ(BCSTATMSG,*) BCSTATIM, TBCSTATIM
            END IF
            IF (BCFLAG_TEMP /= 0) THEN
                READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
            END IF
        ELSE IF (RES_BC_FLAG == 4) THEN
            IF (NOPE > 0) THEN
                READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC, TBCTIMEINC
                READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM, TBCSTATIM
            END IF
            IF (BCFLAG_TEMP /= 0) THEN
                READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
            END IF
        END IF
    END IF

    IF (CBAROCLINIC) THEN
        READ(15,80) SPONGEDISTMSG
        READ(15,80) EqnstateMSG
    END IF


!     jgf45.12: Read in the parameters for the transport equation, if
!     necessary.
    IF (C3D_BTrans) THEN

    !     Lateral and vertical diffusion coefficients.
        READ(15,*) NLSD, NVSD
        READ(15,*) NLTD, NVTD
    !     Time stepping coefficient for the transport equation terms.
        READ(15,*) ALP4
    !     Temperature boundary condition file type, if necessary
    !         if ( IDEN .eq. 3 .or. IDEN .eq. 4 ) then
    !            READ(15,*) NTF
    !         endif
    endif

    call memory_status()
    RETURN

    80 FORMAT(A80)
!----------------------------------------------------------------------
    END SUBROUTINE READ15_3DVS
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!         S U B R O U T I N E    S E T U P K D T S E A R C H
!----------------------------------------------------------------------

!  tcm v49.48.01 Subroutine that creates a global list of element
!             center's and a global list of the size of an
!             element's radius
!  Written by:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39056

!----------------------------------------------------------------------
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
        shplocal(7)= shplocal(1)*shplocal(5) &
        - shplocal(2)*shplocal(4)

    !        compute the radius of the circle that circumscribes
    !        the triangle then scale it by 50% larger to allow for
    !        a buffer later on
        rmax(i)=1.50D0*( sqrt(shplocal(6)*shplocal(6)+ &
        shplocal(3)*shplocal(3))* &
        sqrt(shplocal(4)*shplocal(4)+ &
        shplocal(1)*shplocal(1))* &
        sqrt(shplocal(5)*shplocal(5)+ &
        shplocal(2)*shplocal(2))/ &
        (2.d0*shplocal(7)) )
    !        Compute the barycenter of the element
        bcxy(1,i) = sum(xelm)/3.d0
        bcxy(2,i) = sum(yelm)/3.d0
    enddo !i

    return
!----------------------------------------------------------------------
    end subroutine setup_kdt_search
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!         S U B R O U T I N E    K D T S E A R C H
!----------------------------------------------------------------------

!  tcm v49.48.01 Subroutine that uses the KDTREE2 algorithm for
!             finding what element a point lies in.
!  Written by:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39056

!----------------------------------------------------------------------
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
    LOGICAL :: ElementFound  ! .TRUE. when a corresponding element is found
    INTEGER :: ClosestElement  ! element with closest match
    REAL(sz), PARAMETER :: Tolerance = 1.0d-5     ! area difference for match

    ElementFound = .FALSE. 

!    Be sure the maximum search depth is not larger than
!    the number of elements being kept
    if (nelg < srchdp) srchdp = nelg

!     Create the search tree
    tree => kdtree2_create(bcxy,rearrange= .TRUE. ,sort= .TRUE. )

!     allocate space for the search results from the tree
    allocate(KDRESULTS(srchdp))


    do I = 1, NumOfStations
        Xsta = XCoords(I)
        Ysta = YCoords(I)
        ElementFound = .FALSE. 

    ! Find the srchdp nearest elements to this point
        call kdtree2_n_nearest(tp=tree,qv=(/Xsta,Ysta/), &
        nn=srchdp,results=KDRESULTS)
        itc = 1
        ClosestElement = KDRESULTS(itc)%idx
    !       Check to see if the point lies with rmax of any of these elements
        elmmin = minval(sqrt(KDRESULTS(1:srchdp)%dis) &
        - rmax(KDRESULTS(1:srchdp)%idx) )
        if(elmmin(1) <= 0.0D0) then  ! Point lies within search radius of an element
        !           loop through the elements in the search list
            do while ((ElementFound.eqv. .FALSE. ) .AND. (itc <= srchdp))
                iek = KDRESULTS(itc)%idx  !Current search element number
            !              Get the distance from this point to the barycenter of the
            !              current element
                dist = sqrt(KDRESULTS(itc)%dis)
            !              If the distance is less than or equal to rmax (rmax=1.5*element radius)
            !              Then the point is near the element and might be in it
            !              Proceed with the weights test
                if(dist-rmax(iek) <= 0.0d0) then
                ! et the shape function for this element
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
                    IF (AE < Tolerance) THEN
                        ElementFound = .TRUE. 
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
        if( .NOT. ElementFound) then !We did not find the element
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER == 1) THEN
                print *, "The program will estimate nearest element."
                print *, "WARNING. Distance to nearest element is ", &
                sqrt(KDRESULTS(1)%dis)
                print *, " "
                FDEle(I) = ClosestElement
            ELSE
                print *, "ERROR. Distance to nearest element is ", &
                sqrt(KDRESULTS(1)%dis)
                stop
            ENDIF
        ENDIF
    enddo !loop over station points

!     Deallocate the tree
    call kdtree2_destroy(tp=tree)

    1234 format(A30,1x,I4,1x,'does not lie in the grid.')

    return

    end subroutine kdtsearch


!-----------------------------------------------------------------------
!          S U B R O U T I N E   C O O R D  T O  E L E
!-----------------------------------------------------------------------

!     jgf45.11 Subroutine to take a list of X and Y cartesian
!     coordinates and find the corresponding elements.

!-----------------------------------------------------------------------
    SUBROUTINE CoordToEle(XCoords,YCoords,FDEle,NumOfStations,Desc)
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: NumOfStations                     ! total
    REAL(8), intent(in), dimension(NumOfStations) :: XCoords ! cartesian
    REAL(8), intent(in), dimension(NumOfStations) :: YCoords ! cartesian
    INTEGER, intent(out), dimension(NumOfStations):: FDEle   ! FullDomain
    CHARACTER(len=30), intent(in) :: Desc                    ! description

    INTEGER :: I,J,K,M,N                                        ! loop counters
    INTEGER :: N1, N2, N3, KMIN
    REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3       ! geometry
    REAL(8) AE, AEMIN, AREASK, AA                            ! area
    LOGICAL :: ElementFound  ! .TRUE. when a corresponding element is found

    DO I = 1, NumOfStations
        ElementFound = .FALSE. 
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
            IF (AE < AEMIN) THEN
                AEMIN=AE
                KMIN=K
            ENDIF
            IF (AE < 1.0E-5) THEN
                ElementFound = .TRUE. 
                FDEle(I)=K
            ENDIF
        ENDDO
        IF ( ElementFound .eqv. .FALSE. ) THEN
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER == 1) THEN
                print *, "The program will estimate nearest element."
                print *, "WARNING. Proximity index is ",AEMIN
                print *, " "
                FDEle(I) = KMIN
            ELSE
                print *, "ERROR. Proximity index is ",AEMIN
                stop
            ENDIF
        ENDIF
    ENDDO

    1234 format(A30,1x,I4,1x,'does not lie in the grid.')

    RETURN
!-----------------------------------------------------------------------
    END SUBROUTINE CoordToEle
!-----------------------------------------------------------------------



