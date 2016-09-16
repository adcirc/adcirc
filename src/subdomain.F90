
!                      ADCIRC - SUBDOMAIN MODELING MODULE

!    ========================================================================
!    |                                                                      |
!    |   This file contains the subroutines required by Subdomain Modeling, |
!    |   an approach to reduce the total runtime of a series of hurricane   |
!    |   storm surge simulations on a smaller grid, a subdomain grid,       |
!    |   extracted from the original grid.                                  |
!    |                                                                      |
!    |   Written by Alper Altuntas, aaltunt@ncsu.edu                        |
!    |   North Carolina State University,                                   |
!    |   2013                                                               |
!    |                                                                      |
!    ========================================================================

    module subdomain
    use sizes, only : sz, globaldir
    use global, only : DEBUG, ECHO, INFO, WARNING, ERROR, logMessage, &
    allMessage, scratchMessage, setMessageSource, &
    unsetMessageSource
    use adcirc_mod, only : adcirc_terminate

    logical :: subdomainOn
    integer :: NOUTGS, NSPOOLGS, enforceBN
    integer :: ncbnr, nobnr, nibnr, localncbnr, localnobnr, localnibnr
    integer,allocatable :: cbnr(:), obnr(:), ibnr(:)
    integer,allocatable :: localcbnr(:), localobnr(:), localibnr(:)
    integer :: sbtiminc, nobn, nibn, ncbn
    real(sz),allocatable :: eobn1(:),uobn1(:),vobn1(:),eibn1(:)
    real(sz),allocatable :: eobn2(:),uobn2(:),vobn2(:),eibn2(:)
    real(sz),allocatable :: ecbn1(:),ucbn1(:),vcbn1(:)
    real(sz),allocatable :: ecbn2(:),ucbn2(:),vcbn2(:)
    real(sz),allocatable :: setEob(:), setUob(:), setEcb(:), setUcb(:)
    real(sz),allocatable :: setVob(:), setVcb(:), setEib(:)
    integer,allocatable :: setWDob(:), setWDcb(:)
    integer,allocatable :: wdobn1(:),wdobn2(:),wdcbn1(:),wdcbn2(:)
    integer,allocatable :: cbn(:),obn(:),ibn(:)
    integer :: nlines, bchange

    contains

    SUBROUTINE readFort015()
!    ========================================================================
!      This subroutine reads in additional modeling parameters from fort.015
!      and initiates the Subdomain Modeling

!       - NOUTGS: Type of the run:
!     NOUTGS=0 => full run
!     NOUTGS=1 => subdomain run (old file formatting)
!     NOUTGS=2 => subdomain run (new file formatting)
!       - NSPOOLGS: The number of timesteps at which information is
!                   written to output files fort.06*
!       - enforceBN: Boundary enforcing flag.
!     enforceBN = 0 => no forcing (full domain)
!     enforceBN = 1 => forcing (subdomain, old file formatting)
!     enforceBN = 2 => forcing (subdomain, new file formatting)
!       - nobnr: The number of outer boundary nodes of subdomain grids to be
!                recorded to fort.06* during a full run.
!       - obnr(i): i.th outer boundary node to be recorded to fort.065 and
!                  fort.066
!       - NIBNR: The number of inner boundary nodes of subdomain grids to be
!                recorded to fort.065 and fort.066 during a full run.
!       - IBNR(i): i.th inner boundary node to be recorded to fort.065 and
!                  fort.066
!    ========================================================================

#ifndef CMPI
    Use GLOBAL, only: myproc
#else 
    USE GLOBAL, only: myproc,nodes_lg
    use mesh, only: np
    Use MESSENGER
#endif
    implicit none
    integer :: i,j
    logical :: fileFound

    call setMessageSource("readFort015")
#if defined(SUBDOMAIN_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif
            
    fileFound = .FALSE. 
    inquire(file=trim(globaldir)//'/'//'fort.015', Exist=fileFound)
    if (fileFound.eqv. .FALSE. ) then
        call allMessage(ERROR, &
        'Subdomain modeling was activated but the fort.015 file ' &
        // ' was not found in "' // trim(globaldir) // '".')
        call adcirc_terminate()
    endif

    if (myproc == 0) print *, "Subdomain Active"
    open(1015, file=trim(globaldir)//'/'//'fort.015')
    read(1015,*) NOUTGS
    read(1015,*) NSPOOLGS
    read(1015,*) enforceBN
    select case(noutgs)
    case(0)
! subdomain run (no b.c. recording)
    case(1)
! open ocean boundary nodes:
    read(1015,*) ncbnr
    allocate(cbnr(ncbnr))
    do i=1, ncbnr
        read(1015,*) cbnr(i)
    enddo
    case(2)
! outer boundary nodes:
    read(1015,*) nobnr
    allocate(obnr(nobnr))
    do i=1, nobnr
        read(1015,*) obnr(i)
    enddo
! inner boundary nodes:
    read(1015,*) nibnr
    allocate(ibnr(nibnr))
    do i=1, nibnr
        read(1015,*) ibnr(i)
    enddo
    case default
    call allMessage(ERROR,'The NOUTGS value is invalid.')
    call adcirc_terminate()
    end select
#ifdef CMPI
! localize global record-node arrays:
    select case(noutgs)
    case(0)
! subdomain run (no b.c. recording)
    case(1)
! open ocean boundary nodes:
    localncbnr = 0
    do i=1,np
        if(any(cbnr == nodes_lg(i))) then
            localncbnr = localncbnr+1
        endif
    enddo
    allocate( localcbnr(localncbnr))
    j=1
    do i=1,np
        if(any(cbnr == nodes_lg(i))) then
            localcbnr(j) = i
            j=j+1
        endif
    enddo
    case(2)
! outer boundary nodes:
    localnobnr = 0
    do i=1,np
        if(any(obnr == nodes_lg(i))) then
            localnobnr = localnobnr+1
        endif
    enddo
    allocate( localobnr(localnobnr))
    j=1
    do i=1,np
        if(any(obnr == nodes_lg(i))) then
            localobnr(j) = i
            j=j+1
        endif
    enddo
! inner boundary nodes:
    localnibnr = 0
    do i=1,np
        if(any(ibnr == nodes_lg(i))) then
            localnibnr = localnibnr+1
        endif
    enddo
    allocate( localibnr(localnibnr))
    j=1
    do i=1,np
        if(any(ibnr == nodes_lg(i))) then
            localibnr(j) = i
            j=j+1
        endif
    enddo
    case default
    call allMessage(ERROR,'The NOUTGS value is invalid.')
    call adcirc_terminate()
    end select
#endif

#if defined(SUBPREP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()
          
    END SUBROUTINE readFort015


    SUBROUTINE writeFort065(it)

!     | This subroutine writes ETA2, x/y velocity and wet/dry states of
!     | open ocean boundary nodes of subdomain grids to fort.065 file,
!     | during a full run. This subroutine is called within timestep.F

#ifndef CMPI
    Use GLOBAL, only : myproc, dt, rnday, eta2, uu2, vv2, nodecode
#else
    Use GLOBAL, only : myproc, dt, rnday, eta2, uu2, vv2, nodecode, nodes_lg
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer,intent(in) :: it
    integer :: i,n,gn

! open fort.065 at the first timestep
    if (it == 1) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1065,file=procLoc//'/'//'fort.065', status='REPLACE')
        write(1065,*) 'Subdomain Modeling CBN output'
        write(1065,*) nspoolgs,localncbnr,int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLGS,lncbnr,tsteps"
#else
        open(1065,file='fort.065',status='REPLACE')
        write(1065,*) 'Subdomain Modeling CBN output'
        write(1065,*) nspoolgs,ncbnr,int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLGS,ncbnr,tsteps"
#endif
        close(1065)
    endif

    if(mod(it,NSPOOLGS) == 0) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1065,file=procLoc//'/'//'fort.065', &
        access='SEQUENTIAL', position='APPEND')
        write(1065,*) it, "! Timestep"
        do i=1,localncbnr
            n = localcbnr(i)
            gn = nodes_lg(n)
            write(1065,*) gn, eta2(n), uu2(n)
            write(1065,*) vv2(n), nodecode(n)
        enddo
#else
        open(1065,file='fort.065',access='SEQUENTIAL', &
        position='APPEND')
        write(1065,*) it, "! Timestep"
        do i=1,ncbnr
            n = cbnr(i)
            write(1065,*) n, eta2(n), uu2(n)
            write(1065,*) vv2(n), nodecode(n)
        enddo
#endif
        close(1065)
    endif

    END SUBROUTINE writeFort065






    SUBROUTINE writeFort066(it)

!     | This subroutine writes ETAS, x/y velocity and wet/dry states of
!     | outer boundary nodes of subdomain grids to fort.066 file,
!     | during a full run. This subroutine is called within timestep.F

#ifndef CMPI
    Use GLOBAL, only : myproc, dt, rnday, etas, uu2, vv2, nodecode
#else 
    Use GLOBAL, only : myproc, dt, rnday, etas, uu2, vv2, nodecode,&
                       nodes_lg
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer,intent(in) :: it
    integer :: i,n,gn

! open fort.066 at the first timestep
    if (it == 1) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1066,file=procLoc//'/'//'fort.066', status='REPLACE')
        write(1066,*) 'Subdomain Modeling OBN output'
        write(1066,*) nspoolgs,localnobnr,int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLSG,lnobnr,tsteps"
#else
        open(1066,file='fort.066',status='REPLACE')
        write(1066,*) 'Subdomain Modeling OBN output'
        write(1066,*) nspoolgs, nobnr,int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLGS,nobnr,tsteps"
#endif
        close(1066)
    endif

    if(mod(it,NSPOOLGS) == 0) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1066,file=procLoc//'/'//'fort.066', &
        access='SEQUENTIAL', position='APPEND')
        write(1066,*) it, "! Timestep"
        do i=1,localnobnr
            n = localobnr(i)
            gn = nodes_lg(n)
            write(1066,*) gn, etas(n), uu2(n)
            write(1066,*) vv2(n), nodecode(n)
        enddo
#else
        open(1066,file='fort.066',access='SEQUENTIAL', &
        position='APPEND')
        write(1066,*) it, "! Timestep"
        do i=1,nobnr
            n = obnr(i)
            write(1066,*) n, etas(n), uu2(n)
            write(1066,*) vv2(n), nodecode(n)
        enddo
#endif
        close(1066)
    endif

    END SUBROUTINE writeFort066






    SUBROUTINE writeFort067(it)

!     | This subroutine writes ETAS of inner boundary nodes of subdomain
!     | grids to fort.067 file during a full run. This subroutine is called
!     | within timestep.F

#ifndef CMPI
    Use GLOBAL, only : myproc, dt, rnday, etas
#else
    Use GLOBAL, only : myproc, dt, rnday, etas, nodes_lg
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer,intent(in) :: it
    integer :: i,n,gn

! open fort.067 at the first timestep
    if (it == 1) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1067,file=procLoc//'/'//'fort.067', status='REPLACE')
        write(1067,*) 'Subdomain Modeling IBN output'
        write(1067,*) nspoolgs,localnibnr,int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLGS,lnibnr,tsteps"
#else
        open(1067,file='fort.067',status='REPLACE')
        write(1067,*) 'Subdomain Modeling IBN output'
        write(1067,*) nspoolgs, nibnr, int(rnday*86400/(dt*nspoolgs)), &
        " ! NSPOOLGS,nibnr,tsteps"
#endif
        close(1067)
    endif

    if(mod(it,NSPOOLGS) == 0) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1067,file=procLoc//'/'//'fort.067', &
        access='SEQUENTIAL', position='APPEND')
        write(1067,*) it, "! Timestep"
        do i=1,localnibnr
            n = localibnr(i)
            gn = nodes_lg(n)
            write(1067,*) gn, etas(n)
        enddo
#else
        open(1067,file='fort.067',access='SEQUENTIAL', &
        position='APPEND')
        write(1067,*) it, "! Timestep"
        do i=1,nibnr
            n = ibnr(i)
            write(1067,*) n, etas(n)
        enddo
#endif
        close(1067)
    endif

    END SUBROUTINE writeFort067



    SUBROUTINE openFort019H(TimeLoc)
!     | This subroutine opens fort.019 for a hot-started model and
!     | finds the proper place in the b.c. file.

    Use GLOBAL, only : myproc, iths
    USE BOUNDARIES, ONLY : neta
#ifdef CMPI
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer :: i,j,n,it,itread
    REAL(8), intent(in) :: TimeLoc
#ifdef CMPI
    write(procLoc,'(A2,I4.4)') "PE",myproc
    open(1019, file=procLoc//'/'//'fort.019')
#else
    open(1019, file=trim(globaldir)//'/'//'fort.019')
#endif
    read(1019,*) ! Header
    read(1019,*) sbtiminc, ncbn
    allocate(cbn(ncbn))
    do i=1,ncbn
        read(1019,*) cbn(i)
    enddo
    if (abs(ncbn-neta) > 1) then
        print *, "WARNING: ncbn=!neta"
        print *,  ncbn,neta,myproc
        CALL ADCIRC_Terminate()
    endif
    allocate(ecbn1(ncbn),ucbn1(ncbn),vcbn1(ncbn),wdcbn1(ncbn))
    allocate(ecbn2(ncbn),ucbn2(ncbn),vcbn2(ncbn),wdcbn2(ncbn))
    allocate(setEcb(ncbn),setUcb(ncbn))
    allocate(setVcb(ncbn),setWDcb(ncbn))
    do i=1,ncbn
        ecbn1(i)=0.0
        ucbn1(i)=0.0
        vcbn1(i)=0.0
        wdcbn1(i)=0
    enddo
! Read the first set of boundary conditions:
    read (1019,*) n
    do i=1,ncbn
        read (1019,*) n,ecbn2(i),ucbn2(i)
        read (1019,*) vcbn2(i),wdcbn2(i)
        wdcbn1(i)=wdcbn2(i) !ensure that the wet nodes begin wet
    enddo

    END SUBROUTINE openFort019H




    SUBROUTINE openFort019C()
!     | This subroutine opens fort.019 for a cold-started model

    Use GLOBAL, only : myproc
    USE BOUNDARIES, ONLY : neta
    USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
#ifdef CMPI
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer :: i,n

#ifdef CMPI
    write(procLoc,'(A2,I4.4)') "PE",myproc
    open(1019, file=procLoc//'/'//'fort.019')
#else
    open(1019, file=trim(globaldir)//'/'//'fort.019')
#endif
    read(1019,*) ! Header
    read(1019,*) sbtiminc, ncbn
    allocate(cbn(ncbn))
    do i=1,ncbn
        read(1019,*) cbn(i)
    enddo
    if (abs(ncbn-neta) > 1) then
        print *, "WARNING: ncbn=!neta"
        print *,  ncbn,neta,myproc
        CALL ADCIRC_Terminate()
    endif
    allocate(ecbn1(ncbn),ucbn1(ncbn),vcbn1(ncbn),wdcbn1(ncbn))
    allocate(ecbn2(ncbn),ucbn2(ncbn),vcbn2(ncbn),wdcbn2(ncbn))
    allocate(setEcb(ncbn),setUcb(ncbn))
    allocate(setVcb(ncbn),setWDcb(ncbn))
    do i=1,ncbn
        ecbn1(i)=0.0
        ucbn1(i)=0.0
        vcbn1(i)=0.0
        wdcbn1(i)=0
    enddo

! Read the first set of boundary conditions:
    read (1019,*) n
    do i=1,ncbn
        read (1019,*) n,ecbn2(i),ucbn2(i)
        read (1019,*) vcbn2(i),wdcbn2(i)
        wdcbn1(i)=wdcbn2(i) !ensure that the wet nodes begin wet
    enddo


    END SUBROUTINE openFort019C



    SUBROUTINE readFort019(it)
!     | This subroutine reads in ETA2, x/y velocities and wet/dry states of
!     | outer boundary nodes of the subdomain grid from fort.019 file during
!     | a subdomain run. This subroutine is called within timestep.F.

    Use GLOBAL, only : ihot
    USE ADCIRC_MOD, ONLY : ADCIRC_Terminate
#ifdef CMPI
    Use MESSENGER
#endif
    implicit none
    character(6) :: procLoc
    integer,intent(in) :: it
    double precision :: rateTS
    integer :: i,n

    if (it == 1 .AND. ihot == 0) then  ! coldstart
        call openFort019C()
    endif

    if (mod(it,sbtiminc) == 0) then
        read (1019,*) n
    !             print *, "Reading fort.019 at timestep",n
        do i=1,ncbn
            ecbn1(i) = ecbn2(i)
            ucbn1(i) = ucbn2(i)
            vcbn1(i) = vcbn2(i)
            wdcbn1(i) = wdcbn2(i)
            read (1019,*) n,ecbn2(i),ucbn2(i)
            read (1019,*) vcbn2(i),wdcbn2(i)
        enddo
    endif

! Iteration
    rateTS = mod(it,sbtiminc)/dble(sbtiminc)
    do i=1,ncbn
        setEcb(i) = ecbn1(i) + (ecbn2(i)-ecbn1(i))*ratets
        setUcb(i) = ucbn1(i) + (ucbn2(i)-ucbn1(i))*ratets
        setVcb(i) = vcbn1(i) + (vcbn2(i)-vcbn1(i))*ratets
        setWDcb(i) = wdcbn1(i)
    enddo

    END SUBROUTINE readFort019





    SUBROUTINE readFort020(it)

!     | This subroutine reads in ETAS, x/y velocities and wet/dry states of
!     | outer boundary nodes of the subdomain grid from fort.020 file during
!     | a subdomain run. This subroutine is called within timestep.F.

    Use GLOBAL, only : myproc
#ifdef CMPI
    Use MESSENGER
#endif
    implicit none
    integer,intent(in) :: it
    character(6) :: procLoc
    integer :: i,n
            
    if (it == 1) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1020, file=procLoc//'/'//'fort.020')
#else
        open(1020, file=trim(globaldir)//'/'//'fort.020')
#endif
        read(1020,*) ! Header
        read(1020,*) sbtiminc, nobn
        allocate(obn(nobn))
        do i=1,nobn
            read(1020,*) obn(i)
        enddo
        allocate(eobn2(nobn),uobn2(nobn),vobn2(nobn))
        allocate(wdobn1(nobn),wdobn2(nobn))
        allocate(setEob(nobn),setUob(nobn))
        allocate(setVob(nobn),setWDob(nobn))
    endif

! f (mod(it,sbtiminc).eq.0) then
    read (1020,*) n
    do i=1,nobn
        wdobn1(i) = wdobn2(i)
        read (1020,*) n,eobn2(i),uobn2(i)
        read (1020,*) vobn2(i),wdobn2(i)
        setEob(i) = eobn2(i)
        setUob(i) = uobn2(i)
        setVob(i) = vobn2(i)
        setWDob(i) = wdobn2(i)
    enddo
! ndif

    END SUBROUTINE readFort020

          


    SUBROUTINE readFort021(it)

!     | This subroutine reads in ETAS of inner boundary nodes of the
!     | subdomain grid from fort.021 file during a subdomain run.
!     | This subroutine is called within timestep.F.

    Use GLOBAL, only : myproc
#ifdef CMPI
    Use MESSENGER
#endif
    implicit none
    integer,intent(in) :: it
    character(6) :: procLoc
    integer :: i,n
            
    if (it == 1) then
#ifdef CMPI
        write(procLoc,'(A2,I4.4)') "PE",myproc
        open(1021, file=procLoc//'/'//'fort.021')
#else
        open(1021, file=trim(globaldir)//'/'//'fort.021')
#endif
        read(1021,*) ! Header
        read(1021,*) sbtiminc, nibn
        allocate(ibn(nibn))
        do i=1,nibn
            read(1021,*) ibn(i)
        enddo
        allocate(eibn2(nibn),setEib(nibn))
    endif

!         if (mod(it,sbtiminc).eq.0) then
    read (1021,*) n
!             print *, "Reading fort.021 at timestep",n
    do i=1,nibn
        read (1021,*) n,eibn2(i)
        setEib(i) = eibn2(i)
    enddo
!         endif

    END SUBROUTINE readFort021



    SUBROUTINE enforceEcb()

!     | This subroutine enforces ETA2 of opean ocean boundary nodes
!     | using setEcb array. Called within timestep.F

    Use GLOBAL, only : eta2
    implicit none
    integer :: i,n

    do i=1,ncbn
        n = cbn(i)
        ETA2(n) = setEcb(i)
    enddo

    END SUBROUTINE enforceEcb




    SUBROUTINE enforceEob()

!     | This subroutine enforces ETAS of outer boundary nodes
!     | using setEob array. Called within timestep.F

    Use GLOBAL, only : etas
    implicit none
    integer :: i,n
           
    do i=1,nobn
        n = obn(i)
        ETAS(n) = setEob(i)
    enddo

    END SUBROUTINE enforceEob




    SUBROUTINE enforceUVcb()

!     | This subroutine enforces U and V velocities of open ocean boundary
!     | nodes using setUcb and setVcb arrays. Called within timestep.F

    Use GLOBAL, only : uu2, vv2
    implicit none
    integer :: i,n

    do i=1,ncbn
        n = cbn(i)
        UU2(n) = setUcb(i)
        VV2(n) = setVcb(i)
    enddo

    END SUBROUTINE enforceUVcb




    SUBROUTINE enforceUVob()

!     | This subroutine enforces U and V velocities of outer boundary
!     | nodes using setUob and setVob arrays. Called within timestep.F

    Use GLOBAL, only : uu2, vv2
    implicit none
    integer :: i,n
           
    do i=1,nobn
        n = obn(i)
        UU2(n) = setUob(i)
        VV2(n) = setVob(i)
    enddo

    END SUBROUTINE enforceUVob



    SUBROUTINE enforceWDcb()

!     | This subroutine enforces wet/dry flags of outer boundary
!     | nodes using setWDcb. Called within timestep.F

    Use GLOBAL, only : nnodecode
    implicit none
    integer :: i,n

    do i=1,ncbn
        n = cbn(i)
        NNODECODE(n) = setWDcb(i)
    enddo

    END SUBROUTINE enforceWDcb





    SUBROUTINE enforceWDob()

!     | This subroutine enforces wet/dry flags of outer boundary
!     | nodes using setWDob. Called within timestep.F

    Use GLOBAL, only : nnodecode
    implicit none
    integer :: i,n
           
    do i=1,nobn
        n = obn(i)
        NNODECODE(n) = setWDob(i)
    enddo

    END SUBROUTINE enforceWDob






    SUBROUTINE enforceEib()

!     | This subroutine enforces ETAS of inner boundary nodes
!     | using setEib array. Called within timestep.F

    Use GLOBAL, only : etas
    implicit none
    integer :: i,n

    do i=1,nibn
        n = ibn(i)
        ETAS(n) = setEib(i)
    enddo

    END SUBROUTINE enforceEib




    SUBROUTINE enforceGWCELVob()

!     | This subroutine changes the value of the GWCE_LV vector
!     | elements of the outer boundary nodes of subdomains so
!     | that jcg solver leads to the enforced values of ETA for
!     | the outer boundary nodes.

    use sizes, only : mnei
    use global, only : coef, etas, gwce_lv
    use mesh, only : neitab
    implicit none
    integer :: i,j,n,neighbor
    real(8) :: newGWCElv

    do i=1,nobn
        n = obn(i)
        newGWCElv = 0.
        do j=1,mnei
            neighbor = neitab(n,j)
            if (neighbor /= 0) then
                newGWCElv = newGWCElv + COEF(n,j)*ETAS(neighbor)
            endif
        enddo
        GWCE_LV(n) = newGWCElv
    enddo

    END SUBROUTINE enforceGWCELVob




    SUBROUTINE checkChange()

!     | This subroutine checks if w/d status of any o.b.n. has changed.
!     | If a change has detected, ncchange is set to 1 at the next timestep
!     | to recalculate COEF of outer-inner boundary nodes

    Use GLOBAL, only : ilump
#ifdef CMPI
    Use MESSENGER
#endif
    integer :: i

    if (bchange == 1) then
        ncchange = 1
        bchange = 0
    endif
    do i=1,nobn
        if (wdobn2(i) /= wdobn1(i)) then
            bchange = 1
        endif
    enddo
#ifdef CMPI
    IF ( ILump == 0 ) THEN
        call WetDrySum(NCCHANGE)
    ENDIF
#endif

    END SUBROUTINE checkChange

    end module subdomain
