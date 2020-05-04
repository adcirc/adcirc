C
C                      ADCIRC - SUBDOMAIN MODELING MODULE
C       
C    ========================================================================
C    |                                                                      |
C    |   This file contains the subroutines required by Subdomain Modeling, |
C    |   an approach to reduce the total runtime of a series of hurricane   |
C    |   storm surge simulations on a smaller grid, a subdomain grid,       |
C    |   extracted from the original grid.                                  |
C    |                                                                      |
C    |   Written by Alper Altuntas, aaltunt@ncsu.edu                        |
C    |   North Carolina State University,                                   |
C    |   2013                                                               |
C    |                                                                      |
C    ========================================================================

      module subdomain
      use sizes, only : sz, globaldir
      use global, only : DEBUG, ECHO, INFO, WARNING, ERROR, logMessage,
     &                    allMessage, scratchMessage, setMessageSource,
     &                    unsetMessageSource 

      logical subdomainOn                                         
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
C    ========================================================================
C      This subroutine reads in additional modeling parameters from fort.015
C      and initiates the Subdomain Modeling
C     
C       - NOUTGS: Type of the run:
C     NOUTGS=0 => full run
C     NOUTGS=1 => subdomain run (old file formatting)
C     NOUTGS=2 => subdomain run (new file formatting)
C       - NSPOOLGS: The number of timesteps at which information is
C                   written to output files fort.06*
C       - enforceBN: Boundary enforcing flag.
C     enforceBN = 0 => no forcing (full domain)
C     enforceBN = 1 => forcing (subdomain, old file formatting)
C     enforceBN = 2 => forcing (subdomain, new file formatting)
C       - nobnr: The number of outer boundary nodes of subdomain grids to be
C                recorded to fort.06* during a full run.
C       - obnr(i): i.th outer boundary node to be recorded to fort.065 and
C                  fort.066 
C       - NIBNR: The number of inner boundary nodes of subdomain grids to be
C                recorded to fort.065 and fort.066 during a full run.
C       - IBNR(i): i.th inner boundary node to be recorded to fort.065 and
C                  fort.066 
C    ========================================================================


         Use sizes, only : myproc
#ifdef CMPI
         use global, only : nodes_lg
         use mesh, only : np
         Use MESSENGER 
#endif
         implicit none
         integer :: i,j
         logical :: fileFound

      call setMessageSource("readFort015")
#if defined(SUBDOMAIN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
        
         fileFound = .false.
         inquire(file=trim(globaldir)//'/'//'fort.015', Exist=fileFound)
         if (fileFound.eqv..false.) then
            call allMessage(ERROR,
     &       'Subdomain modeling was activated but the fort.015 file '
     &       // ' was not found in "' // trim(globaldir) // '".')
            call subdomainTerminate()
         endif

         if (myproc.eq.0) print *, "Subdomain Active"
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
            call subdomainTerminate()
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
               if(any(cbnr.eq.nodes_lg(i))) then
                  localncbnr = localncbnr+1
               endif
            enddo
            allocate( localcbnr(localncbnr))
            j=1
            do i=1,np
               if(any(cbnr.eq.nodes_lg(i))) then
                  localcbnr(j) = i
                  j=j+1
               endif
            enddo
         case(2)
            ! outer boundary nodes:
            localnobnr = 0 
            do i=1,np
               if(any(obnr.eq.nodes_lg(i))) then
                  localnobnr = localnobnr+1
               endif
            enddo
            allocate( localobnr(localnobnr))
            j=1
            do i=1,np
               if(any(obnr.eq.nodes_lg(i))) then
                  localobnr(j) = i
                  j=j+1
               endif
            enddo
            ! inner boundary nodes:
            localnibnr = 0 
            do i=1,np
               if(any(ibnr.eq.nodes_lg(i))) then
                  localnibnr = localnibnr+1
               endif
            enddo
            allocate( localibnr(localnibnr))
            j=1
            do i=1,np
               if(any(ibnr.eq.nodes_lg(i))) then
                  localibnr(j) = i
                  j=j+1
               endif
            enddo
         case default
            call allMessage(ERROR,'The NOUTGS value is invalid.')
            call subdomainTerminate()
         end select
#endif

#if defined(SUBPREP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      
      END SUBROUTINE readFort015


      SUBROUTINE writeFort065(it)

C     | This subroutine writes ETA2, x/y velocity and wet/dry states of
C     | open ocean boundary nodes of subdomain grids to fort.065 file,
C     | during a full run. This subroutine is called within timestep.F

         Use SIZES, only : myproc
         Use GLOBAL, only : dt, rnday, eta2, uu2, vv2, nodecode
#ifdef CMPI
     &                     , nodes_lg
         Use MESSENGER
#endif
         implicit none
         character(6) :: procLoc
         integer,intent(in) :: it
         integer :: i,n,gn

         ! open fort.065 at the first timestep
         if (it.eq.1) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1065,file=procLoc//'/'//'fort.065', status='REPLACE')
             write(1065,*) 'Subdomain Modeling CBN output'
             write(1065,*) nspoolgs,localncbnr,int(rnday*86400/(dt*nspoolgs)),
     &                      " ! NSPOOLGS,lncbnr,tsteps"
#else
             open(1065,file='fort.065',status='REPLACE')
             write(1065,*) 'Subdomain Modeling CBN output'
             write(1065,*) nspoolgs,ncbnr,int(rnday*86400/(dt*nspoolgs)), 
     &                     " ! NSPOOLGS,ncbnr,tsteps"
#endif 
             close(1065)
         endif

         if(mod(it,NSPOOLGS).eq.0) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1065,file=procLoc//'/'//'fort.065',
     &            access='SEQUENTIAL', position='APPEND')
             write(1065,*) it, "! Timestep"
             do i=1,localncbnr
                 n = localcbnr(i)
                 gn = nodes_lg(n)
                 write(1065,*) gn, eta2(n), uu2(n)
                 write(1065,*) vv2(n), nodecode(n)
             enddo
#else
             open(1065,file='fort.065',access='SEQUENTIAL',
     &            position='APPEND')
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

C     | This subroutine writes ETAS, x/y velocity and wet/dry states of
C     | outer boundary nodes of subdomain grids to fort.066 file,
C     | during a full run. This subroutine is called within timestep.F

         Use SIZES, only : myproc
         Use GLOBAL, only : dt, rnday, etas, uu2, vv2, nodecode
#ifdef CMPI
     &                      , nodes_lg
         Use MESSENGER 
#endif
         implicit none
         character(6) :: procLoc
         integer,intent(in) :: it
         integer :: i,n,gn

         ! open fort.066 at the first timestep
         if (it.eq.1) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1066,file=procLoc//'/'//'fort.066', status='REPLACE')
             write(1066,*) 'Subdomain Modeling OBN output'
             write(1066,*) nspoolgs,localnobnr,int(rnday*86400/(dt*nspoolgs)),
     &                     " ! NSPOOLSG,lnobnr,tsteps"
#else
             open(1066,file='fort.066',status='REPLACE')
             write(1066,*) 'Subdomain Modeling OBN output'
             write(1066,*) nspoolgs, nobnr,int(rnday*86400/(dt*nspoolgs)),
     &                     " ! NSPOOLGS,nobnr,tsteps"
#endif 
             close(1066)
         endif

         if(mod(it,NSPOOLGS).eq.0) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1066,file=procLoc//'/'//'fort.066',
     &            access='SEQUENTIAL', position='APPEND')
             write(1066,*) it, "! Timestep"
             do i=1,localnobnr
                 n = localobnr(i)
                 gn = nodes_lg(n)
                 write(1066,*) gn, etas(n), uu2(n)
                 write(1066,*) vv2(n), nodecode(n)
             enddo
#else
             open(1066,file='fort.066',access='SEQUENTIAL',
     &            position='APPEND')
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

C     | This subroutine writes ETAS of inner boundary nodes of subdomain  
C     | grids to fort.067 file during a full run. This subroutine is called
C     | within timestep.F

         Use sizes, only : myproc
         Use global, only : dt, rnday, etas
#ifdef CMPI
     &                      , nodes_lg
         Use MESSENGER 
#endif
         implicit none
         character(6) :: procLoc
         integer,intent(in) :: it
         integer :: i,n,gn

         ! open fort.067 at the first timestep
         if (it.eq.1) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1067,file=procLoc//'/'//'fort.067', status='REPLACE')
             write(1067,*) 'Subdomain Modeling IBN output'
             write(1067,*) nspoolgs,localnibnr,int(rnday*86400/(dt*nspoolgs)),
     &                     " ! NSPOOLGS,lnibnr,tsteps"
#else
             open(1067,file='fort.067',status='REPLACE')
             write(1067,*) 'Subdomain Modeling IBN output'
             write(1067,*) nspoolgs, nibnr, int(rnday*86400/(dt*nspoolgs)),
     &                     " ! NSPOOLGS,nibnr,tsteps"
#endif 
             close(1067)
         endif

         if(mod(it,NSPOOLGS).eq.0) then
#ifdef CMPI
             write(procLoc,'(A2,I4.4)') "PE",myproc
             open(1067,file=procLoc//'/'//'fort.067',
     &            access='SEQUENTIAL', position='APPEND')
             write(1067,*) it, "! Timestep"
             do i=1,localnibnr
                 n = localibnr(i)
                 gn = nodes_lg(n)
                 write(1067,*) gn, etas(n)
             enddo
#else
             open(1067,file='fort.067',access='SEQUENTIAL',
     &            position='APPEND')
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
C     | This subroutine opens fort.019 for a hot-started model and 
C     | finds the proper place in the b.c. file.

         Use SIZES, only : myproc
         Use GLOBAL, only : iths
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
         if (abs(ncbn-neta).gt.1) then
             print *, "WARNING: ncbn=!neta"
             print *,  ncbn,neta,myproc
             CALL subdomainTerminate()
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
C     | This subroutine opens fort.019 for a cold-started model

         Use SIZES, only : myproc
         USE BOUNDARIES, ONLY : neta
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
         if (abs(ncbn-neta).gt.1) then
            print *, "WARNING: ncbn=!neta"
            print *,  ncbn,neta,myproc
            CALL subdomainTerminate()
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
C     | This subroutine reads in ETA2, x/y velocities and wet/dry states of
C     | outer boundary nodes of the subdomain grid from fort.019 file during
C     | a subdomain run. This subroutine is called within timestep.F.

         Use GLOBAL, only : ihot 
#ifdef CMPI
         Use MESSENGER
#endif
         implicit none
         character(6) :: procLoc
         integer,intent(in) :: it
         double precision rateTS
         integer :: i,n

         if (it.eq.1.and.ihot.eq.0) then  ! coldstart
                 call openFort019C() 
         endif

         if (mod(it,sbtiminc).eq.0) then
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

C     | This subroutine reads in ETAS, x/y velocities and wet/dry states of
C     | outer boundary nodes of the subdomain grid from fort.020 file during
C     | a subdomain run. This subroutine is called within timestep.F.

         Use SIZES, only : myproc
#ifdef CMPI
         Use MESSENGER 
#endif
         implicit none
         integer,intent(in) :: it
         character(6) :: procLoc
         integer :: i,n
        
         if (it.eq.1) then
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

         !if (mod(it,sbtiminc).eq.0) then
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
         !endif

      END SUBROUTINE readFort020

      


      SUBROUTINE readFort021(it)

C     | This subroutine reads in ETAS of inner boundary nodes of the
C     | subdomain grid from fort.021 file during a subdomain run.
C     | This subroutine is called within timestep.F.

         Use SIZES, only : myproc
#ifdef CMPI
         Use MESSENGER 
#endif
         implicit none
         integer,intent(in) :: it
         character(6) :: procLoc
         integer :: i,n
        
         if (it.eq.1) then
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

C     | This subroutine enforces ETA2 of opean ocean boundary nodes 
C     | using setEcb array. Called within timestep.F

         Use GLOBAL, only : eta2
         implicit none
         integer :: i,n

         do i=1,ncbn
             n = cbn(i)
             ETA2(n) = setEcb(i)
         enddo

      END SUBROUTINE enforceEcb




      SUBROUTINE enforceEob()

C     | This subroutine enforces ETAS of outer boundary nodes 
C     | using setEob array. Called within timestep.F

         Use GLOBAL, only : etas
         implicit none
         integer :: i,n
       
         do i=1,nobn
             n = obn(i)
             ETAS(n) = setEob(i)
         enddo

      END SUBROUTINE enforceEob




      SUBROUTINE enforceUVcb()

C     | This subroutine enforces U and V velocities of open ocean boundary
C     | nodes using setUcb and setVcb arrays. Called within timestep.F

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

C     | This subroutine enforces U and V velocities of outer boundary
C     | nodes using setUob and setVob arrays. Called within timestep.F

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

C     | This subroutine enforces wet/dry flags of outer boundary 
C     | nodes using setWDcb. Called within timestep.F

         Use GLOBAL, only : nnodecode
         implicit none
         integer :: i,n

         do i=1,ncbn
             n = cbn(i)
             NNODECODE(n) = setWDcb(i)
         enddo

      END SUBROUTINE enforceWDcb





      SUBROUTINE enforceWDob()

C     | This subroutine enforces wet/dry flags of outer boundary 
C     | nodes using setWDob. Called within timestep.F

         Use GLOBAL, only : nnodecode
         implicit none
         integer :: i,n
       
         do i=1,nobn
             n = obn(i)
             NNODECODE(n) = setWDob(i)
         enddo

      END SUBROUTINE enforceWDob






      SUBROUTINE enforceEib()

C     | This subroutine enforces ETAS of inner boundary nodes 
C     | using setEib array. Called within timestep.F

         Use GLOBAL, only : etas
         implicit none
         integer :: i,n

         do i=1,nibn
             n = ibn(i)
             ETAS(n) = setEib(i)
         enddo

      END SUBROUTINE enforceEib




      SUBROUTINE enforceGWCELVob()

C     | This subroutine changes the value of the GWCE_LV vector
C     | elements of the outer boundary nodes of subdomains so 
C     | that jcg solver leads to the enforced values of ETA for 
C     | the outer boundary nodes.

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
                 if (neighbor.ne.0) then
                     newGWCElv = newGWCElv + COEF(n,j)*ETAS(neighbor)
                 endif
             enddo
             GWCE_LV(n) = newGWCElv
         enddo

      END SUBROUTINE enforceGWCELVob




      SUBROUTINE checkChange()

C     | This subroutine checks if w/d status of any o.b.n. has changed.
C     | If a change has detected, ncchange is set to 1 at the next timestep
C     | to recalculate COEF of outer-inner boundary nodes

         Use GLOBAL, only : ilump
#ifdef CMPI
         Use MESSENGER         
#endif
         integer :: i

         if (bchange.eq.1) then
             ncchange = 1
             bchange = 0
         endif
         do i=1,nobn
             if (wdobn2(i).ne.wdobn1(i)) then
                 bchange = 1
             endif
         enddo
#ifdef CMPI
         IF ( ILump.eq.0 ) THEN
            call WetDrySum(NCCHANGE)
         ENDIF
#endif

      END SUBROUTINE checkChange


      !jgf: Added this subroutine to eliminate dependence on adcirc_mod
      ! module, which should be built last to simplify the build system.
      SUBROUTINE subdomainTerminate(NO_MPI_FINALIZE)
#ifdef CMPI
         Use MESSENGER         
#endif
      implicit none
c
      logical, optional :: no_mpi_finalize
C
      call setMessageSource("subdomainTerminate")
#if defined(SUBDOMAIN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1)
C
#if defined(SUBDOMAIN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      END SUBROUTINE subdomainTerminate



      end module subdomain
