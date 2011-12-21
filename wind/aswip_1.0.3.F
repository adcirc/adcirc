C***********************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                     *
C last changes in this file VERSION 45.05                              *
C                                                                      *
C***********************************************************************
C
      PROGRAM ASWIP
      ! Import custom precision types for cross-platform portability
      USE precision
      ! Import global constants
      USE constants
      ! Import asymmetric hurricane vortex class
      USE vortex
      IMPLICIT NONE

      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====              PROGRAM ASWIP                =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      INTEGER  :: np, ics,ient,icyc,itpc,method
      REAL(sz) :: time, time1,time2, wtratio,timest
      INTEGER,DIMENSION(:),ALLOCATABLE::iFcstInc ! hours between forecasts
      INTEGER,DIMENSION(:),ALLOCATABLE :: iyear, imth, iday, ihr
      INTEGER,DIMENSION(:),ALLOCATABLE::ilat,ilon
      REAL(sz),DIMENSION(:),ALLOCATABLE::rlat,rlon,CycleTime
      INTEGER,DIMENSION(:),ALLOCATABLE::ispd,icpress,atcfRMW
      REAL(sz), SAVE ::  dt
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE:: castType! hindcast/nowcast or forecast?
      REAL(sz),DIMENSION(:),ALLOCATABLE :: HollB
      INTEGER ,DIMENSION(:),ALLOCATABLE :: advr
      REAL(sz) :: Vmax !,speedmax,Rmax,Rmaxmax
      REAL(sz) :: Vr,speedmax2
      REAL(sz), DIMENSION(4) :: r
      REAL(sz) :: Pn
      REAL(sz) :: Pc
      REAL(sz) :: cLat
      REAL(sz) :: cLon
      REAL(sz) :: offset
      REAL(sz) :: dx,dy

      INTEGER ,DIMENSION(:),ALLOCATABLE :: ivr, dir, speed
      INTEGER , DIMENSION(:,:), ALLOCATABLE :: ir   ! given isotach radii
      INTEGER , DIMENSION(:,:), ALLOCATABLE :: irad ! working isotach radii
      INTEGER , DIMENSION(:,:), ALLOCATABLE :: quadflag
      REAL(sz) , DIMENSION(:,:), ALLOCATABLE :: rmaxw
      INTEGER ,DIMENSION(:),ALLOCATABLE :: count_newcycle
      INTEGER ,DIMENSION(:),ALLOCATABLE :: cycle_num
      INTEGER ,DIMENSION(:),ALLOCATABLE :: isotachs_per_cycle
      INTEGER ,DIMENSION(:),ALLOCATABLE :: ipn
      CHARACTER (LEN = 10),DIMENSION(:),ALLOCATABLE :: stormname
      REAL(sz), PARAMETER :: outFreq = thirtySixHundred
      INTEGER :: i, nwi,screenunit,num_cycles,numdiffcyc,num_entry

      REAL(sz) :: clatOld,clonOld, clatNew,clonNew,
     &                    timeOld,timeNew
      REAL(sz) ,DIMENSION(:),ALLOCATABLE :: uTrans, vTrans
      REAL(sz) ,DIMENSION(:),ALLOCATABLE :: uTrans2, vTrans2
      REAL(sz) ,DIMENSION(:),ALLOCATABLE :: h_speed, h_dir
      CHARACTER*1 ew,E,W,ns,N,S
      REAL(sz) :: VmaxBL ! max sustained wind at top of atm. b.l.
      REAL(sz) :: VrBL   ! wind vel at top of atm. b.l. at wind radius 
      INTEGER, DIMENSION(0:5) :: lookup_radii ! periodic interpolation 
      INTEGER :: radii_sum ! record radius values for filling in missing vals
      INTEGER :: num_nonzero ! number of nonzero isotach radii 
      INTEGER :: first_entry    ! first entry in the cycle
      INTEGER :: last_entry     ! last entry in the cycle
      LOGICAL :: higher_isotach ! true if a valid Rmax exists at higher isotach
C
C     jgf49.30: New variables to support command line options
      INTEGER :: IARGC   ! function to return number of command line arguments
      INTEGER :: ARGCOUNT ! number of command line arguments
      CHARACTER(2048) :: CMDLINEARG
      CHARACTER(2) :: CMDLINEOPTION ! the option flag itself
      LOGICAL :: fileFound ! .true. if the input file is there
      REAL(sz) :: percent  ! desired percentage change in Rmax relative to the
                           ! value calculated from isotach radii data (0 to 100)
      REAL(sz) :: rMax     ! constant Rmax value for the duration of the 
                           ! forecast (nautical miles)
      LOGICAL :: constantRmax ! .true. if rMax was specified as something .gt. 
                              ! zero ... takes precedence over percent
                              ! change in Rmax
      LOGICAL :: percentSpecified ! .true. if percent was specified and is
                                  ! something other than 100 percent
C   
C     jgf49.30 Initialize to reasonable defaults.
      screenunit=6
      rMax = 25.d0         ! default to a reasonable number, just in case
      constantRmax=.false. ! rMax is normally calculated from isotach data
      percentSpecified=.false. 
      method=2             ! calc rMaxes based on largest available isotach
                           ! radius for each quadrant, a.k.a. the fancy method
      !
      ! jgf49.30: Process command line options, if any.
      ARGCOUNT = IARGC()
      IF (ARGCOUNT.ne.0) THEN
         I=0
         DO WHILE (I.lt.ARGCOUNT)
            I=I+1
            CALL GETARG(I,CMDLINEARG) 
            CMDLINEOPTION = CMDLINEARG(1:2)
            SELECT CASE(CMDLINEOPTION)
            CASE("-M","-m")  ! method of selecting quadrant Rmax values
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'." 
               READ(CMDLINEARG,*) method
            CASE("-P","-p")  ! percent change of rmax relative to calculated val
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'." 
               READ(CMDLINEARG,*) percent
               IF ( percent.ne.100 ) THEN
                  percentSpecified = .true.
               ENDIF            
            CASE("-R","-r")  ! constant rMax value for forecast
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               READ(CMDLINEARG,*) rMax
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'." 
               IF ( rMax.gt.0 ) THEN
                  constantRmax=.true.
                  write(screenunit,*) 
     &               "INFO: All quadrant Rmaxes will be set to ",rMax,
     &               " throughout the forecast." 
               ENDIF
            CASE DEFAULT
               write(screenunit,*) "WARNING: The command line option '",
     &             CMDLINEARG(1:2),"' is not valid and will be ignored."
            END SELECT
         END DO
      ENDIF
C
      ! jgf49.30 Direct Rmax specification takes precedence over percentage
      ! variation of Rmax.
      IF ((constantRmax.eqv..true.).and.
     &           (percentSpecified.eqv..true.)) THEN 
         write(screenunit,*) 
     &     "WARNING: Rmax variation was specified as a percent (-P ",
     &     percent,") as well as a constant (-R ",rMax,")."
         write(screenunit,*) 
     & "WARNING: The percentage variation for Rmax will not be used."
      ENDIF        
C
      E='E'
      W='W'
      N='N'
      S='S'
      !------------
      ! open output file
      !------------
      OPEN(222,file="NWS_19_fort.22",status="unknown")
C
C     jgf49.30 Check to see if the fort.22 file is there and if not, emit error 
C     message and terminate.
      INQUIRE(file="fort.22",exist=fileFound)
      IF (fileFound.eqv..true.) THEN
         OPEN(22,file="fort.22",status="old")
      ELSE
         write(screenunit,*) 
     &      "ERROR: aswip: The fort.22 file was not found!"
         write(screenunit,*) "ERROR: aswip: Execution terminated."
         STOP 1
      ENDIF
      i=0
      DO 
         READ(22,*,END=991) 
         i=i+1
      ENDDO    
991   continue
      num_entry=i
      close(22)
      OPEN(22,file="fort.22",status="old")

      i=num_entry 
      allocate(advr(i),iyear(i),imth(i),iday(i),ihr(i),castType(i))
      allocate(iFcstInc(i),ilat(i),ilon(i), ispd(i))
      allocate(rlat(i),rlon(i))
      allocate(irad(i,4))
      allocate(icpress(i), ivr(i),ir(i,4), rmaxw(i,4), ipn(i))
      allocate(dir(i),speed(i),count_newcycle(i),cycle_num(i))
      allocate(isotachs_per_cycle(i),CycleTime(i))
      allocate(quadflag(i,4),stormname(i))
      allocate(uTrans(i), vTrans(i),HollB(i),atcfRMW(i) )
      allocate(uTrans2(i), vTrans2(i)) !,HollB(i),atcfRMW(i) )
      allocate(h_speed(i), h_dir(i) )
      
      ient=1
      READ(22,22) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient),atcfRMW(ient),
     &      dir(ient),speed(ient), stormname(ient)

      num_cycles = 1
      count_newcycle(1) = 1
      cycle_num(ient)=count_newcycle(1)
           
!     The logic for the code is set for W and N so if E or S
!     multiply lats and lons by -1
      IF(ew.EQ.E)THEN
         ilon(ient)=(-1)*ilon(ient)
      ENDIF
      IF(ns.EQ.S)THEN
         ilat(ient)=(-1)*ilat(ient)
      ENDIF
      rlon(ient) = (-1.0d0) * ilon(ient) * oneTenth
      rlat(ient) = ilat(ient) * oneTenth

      CycleTime(ient) =  iFcstInc(ient) * hour2sec
      !------------
      ! Set ipn to background pressure of 1013
      !------------
      ipn(ient)=1013
      if (icpress(ient) == 0  ) icpress(ient) = ipn(ient)

      DO ient=2,num_entry

         READ(22,22) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &      dir(ient),speed(ient), stormname(ient)


!        The logic for the code is set for W and N so if E or S
!        multiply lats and lons by -1
         IF(ew.EQ.E)THEN
            ilon(ient)=(-1)*ilon(ient)
         ENDIF
         IF(ns.EQ.S)THEN
            ilat(ient)=(-1)*ilat(ient)
         ENDIF
         rlon(ient) = (-1.0d0) * ilon(ient) * oneTenth
         rlat(ient) = ilat(ient) * oneTenth
C
C        jgf49.30 Commented this line out so that the code that handles
C        the "no isotach" case (near the end of this program) will be active.
C        if (atcfRMW(ient) .eq. 0.0 .AND. atcfRMW(ient-1) .ne. 0.0 )
C     &                                    atcfRMW(ient)=atcfRMW(ient-1)
      !------------
      ! Set ipn to background pressure of 1013
      !---------===
         ipn(ient)=1013
c        if (ipn(ient) == 0 )then
c            ipn(ient)=ipn(ient-1)
c         endif
         if (icpress(ient) == 0  ) icpress(ient) = icpress(ient-1)
C
         if ( (iday(ient).eq.iday(ient-1)).AND.
     &                     (ihr(ient).eq.ihr(ient-1))) then
            if ( iFcstInc(ient) .ne. iFcstInc(ient-1)) then
               num_cycles = num_cycles + 1
               count_newcycle(num_cycles)=ient
            endif
         else
            num_cycles = num_cycles + 1
            count_newcycle(num_cycles)=ient              
         endif
         cycle_num(ient)=num_cycles
         CycleTime(ient) =  iFcstInc(ient) * hour2sec

      ENDDO

        !-------------------------------
         ! Calculate translation velocity
         !-------------------------------
      DO icyc = 1,num_cycles
        if (icyc .lt. num_cycles ) then
           isotachs_per_cycle(icyc) = 
     &                       count_newcycle(icyc+1)-count_newcycle(icyc)
          DO ient =  count_newcycle(icyc), count_newcycle(icyc+1)-1
              
             if (icyc==1 ) then
      !----------------------------------------------------------
                dx = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(cos(rLAT(count_newcycle(icyc+1))*DEG2RAD)
     & * cos(rLAT(count_newcycle(icyc))*DEG2RAD)
     & * sin((rLON(count_newcycle(icyc+1))*DEG2RAD
     & - rLON(count_newcycle(icyc))*DEG2RAD)/2.0d0)**2.0d0))))
     & ,(rLon(count_newcycle(icyc+1))-rLon(count_newcycle(icyc)))  )  
      ! get correct sign

                dy = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(sin((rLAT(count_newcycle(icyc+1))*DEG2RAD
     & -rLAT(count_newcycle(icyc))*DEG2RAD)/2.0d0)**2.0d0))) )
     & ,(rLat(count_newcycle(icyc+1))-rLat(count_newcycle(icyc)))  ) 
      ! get correct sign

                dt = (CycleTime(count_newcycle(icyc+1)) -
     &                             CycleTime(count_newcycle(icyc)))
                if(dt == 0.0d0) dt=6.0d0* hour2sec
      !----------------------------------------------------------
             else
                dx = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(cos(rLAT(count_newcycle(icyc))*DEG2RAD)
     & * cos(rLAT(count_newcycle(icyc-1))*DEG2RAD)
     & * sin((rLON(count_newcycle(icyc))*DEG2RAD
     & - rLON(count_newcycle(icyc-1))*DEG2RAD)/2.0d0)**2.0d0))))
     & ,(rLon(count_newcycle(icyc))-rLon(count_newcycle(icyc-1)))  )  
       ! get correct sign

                dy = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(sin((rLAT(count_newcycle(icyc))*DEG2RAD
     & -rLAT(count_newcycle(icyc-1))*DEG2RAD)/2.0d0)**2.0d0))) )
     & ,(rLat(count_newcycle(icyc))-rLat(count_newcycle(icyc-1)))  )  
      ! get correct sign

                dt = (CycleTime(count_newcycle(icyc)) -
     &                             CycleTime(count_newcycle(icyc-1)))
                if(dt == 0.0d0) dt=6.0d0* hour2sec

             endif
             uTrans(ient) = dx/dt * ms2kt ! now in kts
             vTrans(ient) = dy/dt * ms2kt ! now in kts
          ENDDO ! ient

       ELSE  ! icyc = num_cycles
          isotachs_per_cycle(icyc) =
     &            num_entry - count_newcycle(icyc) + 1
            
          DO ient = count_newcycle(icyc), num_entry
             ! get correct sign
             dx = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(cos(rLAT(count_newcycle(icyc))*DEG2RAD)
     & * cos(rLAT(count_newcycle(icyc-1))*DEG2RAD)
     & * sin((rLON(count_newcycle(icyc))*DEG2RAD
     & - rLON(count_newcycle(icyc-1))*DEG2RAD)/2.0d0)**2.0d0))))
     & ,(rLon(count_newcycle(icyc))-rLon(count_newcycle(icyc-1)))  )

             ! get correct sign
             dy = SIGN( (Rearth*(2.0d0*
     & ASIN(sqrt(sin((rLAT(count_newcycle(icyc))*DEG2RAD
     & -rLAT(count_newcycle(icyc-1))*DEG2RAD)/2.0d0)**2.0d0))) )
     & ,(rLat(count_newcycle(icyc))-rLat(count_newcycle(icyc-1)))  )  

             dt = (CycleTime(count_newcycle(icyc)) -
     &                             CycleTime(count_newcycle(icyc-1)))
             if(dt == 0.0d0) dt=6.0d0* hour2sec

             uTrans(ient) = dx/dt  * ms2kt ! now in kts
             vTrans(ient) = dy/dt  * ms2kt ! now in kts

          ENDDO ! ient
       ENDIF  ! icyc = num_cycles
      ENDDO ! numcycles

      !-------------------------------------------   
      ! convert Utrans and Vtrans to speed and direction  
      ! direction is in compass coordinates 0 == North
      ! increasing clockwise
      !-------------------------------------------   

      DO ient=1,num_entry

         h_speed(ient) = SQRT(uTrans(ient) * uTrans(ient)+ 
     &                          vTrans(ient) * vTrans(ient))           ! in kts


         IF (h_speed(ient) .lt. 1.0d0 ) then
            !jgf50.29: The vortex module can't handle speed and direction 
            ! being zero; it will return NaNs as a result. Persist the
            ! direction from the previous cycle, and make the storm translation
            ! speed small but nonzero.
            h_speed(ient) = 1.0d0
            if (ient.gt.1) then 
               h_dir(ient) = h_dir(ient-1)
            else 
               h_dir(ient)=0.0 
            endif
         ELSE 
            h_dir(ient)= asin(abs(uTrans(ient))/h_speed(ient))*rad2deg  ! in deg
            if ( uTrans(ient) .gt. 0.0d0 ) then
               ! I or IV quadrant
               if( vTrans(ient) .lt. 0.0d0 ) then
                  !  IV quadrant
                  h_dir(ient)=180.0d0 - h_dir(ient) 
               endif
            else !uTrans  < 0
               ! II or III quadrant
               if( vTrans(ient) .gt. 0.0d0 ) then
                  ! II  quadrant
                  h_dir(ient)=-h_dir(ient)
               else
                  h_dir(ient)=h_dir(ient)-180.0d0
               endif
            endif
            if (h_dir(ient) .lt. 0.0d0 ) then 
               h_dir(ient)=h_dir(ient) + 360.0d0
            endif
         ENDIF ! hspeed == 0
      ENDDO

      !-----------------------------------------
      !  Now using the calculated translational velocities 
      !  call the vortex module and compute the Rmax's
      !  to be used in the new input file        
      !----------------------------------------------
      irad(:,:) = ir(:,:)      
      DO ient=1,num_entry

         !--------------------------------------------------------------
         ! transform variables from integers
         ! to real numbers for hurricane vortex calcualtions.
         !--------------------------------------------------------------
         Vmax =  one * ( ispd(ient)  )
         Pn   =  one * ( ipn(ient)   )
         Pc   =  one * ( icpress(ient)  )
         cLat = rlat(ient) 
         cLon = rlon(ient)  

         !---------------------------------------------
         !  need to get some logic incase Vr is  zero
         !  if so we will also be setting ir(:) to Rmax 
         !---------------------------------------------
         if ( ivr(ient) == 0 ) then
            Vr = Vmax
         else
            Vr = one* ivr(ient)
         endif
         
         ! jgf50.36: fill missing values based on number of nonzero radii.
         lookup_radii(0)=ir(ient,4)
         lookup_radii(5)=ir(ient,1)
         radii_sum = 0
         num_nonzero = 0
         do i=1,4
            lookup_radii(i) = ir(ient,i)
            radii_sum = radii_sum + ir(ient,i)
            if (ir(ient,i).gt.0) then
               num_nonzero = num_nonzero + 1
               quadflag(ient,i) = 1   ! use the Rmax resulting from this
            else
               quadflag(ient,i) = 0   ! don't use Rmax resulting from this 
            endif
         end do
         ! fill missing values based on how many are missing
         select case(num_nonzero)
         case(0) ! no isotachs reported, use overall Rmax; set isotach to Vmax
            quadflag(ient,:) = 1
            if (atcfRMW(ient).ne.0) then
               irad(ient,:) = atcfRMW(ient)
            else 
               irad(ient,:) = 40 ! need a nonzero value for Rmax calcs,
                                 ! this result will be thrown away later
            endif   
            Vr=Vmax 
         case(1) ! set all radii equal to the nonzero radius
            where (ir(ient,:).eq.0) irad(ient,:) = radii_sum
         case(2) ! avg the 2 radii that are given to get the 2 that are missing
            where (ir(ient,:).eq.0) irad(ient,:) = 0.5 * radii_sum 
         case(3) ! average the radii on either side of the missing radius
            do i=1,4
               if (ir(ient,i).eq.0) then
                  irad(ient,i) = 0.5*
     &               (lookup_radii(i+1)+lookup_radii(i-1))
               endif
            enddo 
         case(4) 
            ! use all these radii as-is
         case default
            ! the following error message should be unreachable
            write(16,*) "ERROR: aswip: Number of nonzero radii on ",
     &      "line ",ient," not in range 0 to 4."
         end select

        DO i=1,4
           r(i) = one * irad(ient,i)
        END DO
        HollB(ient)=one
        !-------------------------------------------------------
        ! Create a new asymmetric hurricane vortex.
        !
        ! Note: Subtract translational speed from Vmax, then
        ! scale (Vmax - Vt) and Vr up to the top of the surface,
        ! where the cylcostrophic wind balance is valid.
        !-------------------------------------------------------
        VmaxBL = (Vmax-SQRT(uTrans(ient)*uTrans(ient)
     &                     +vTrans(ient)*vTrans(ient)))/windReduction
        VrBL = Vr/windReduction
        ! jgfdebug
        do i=1,4
           if (r(i).lt.1.d0) then
              print *,"ERROR: ient=",ient,"r(",i,")=",r(i)
           endif
        enddo
        CALL newVortex_rmax_b_OUT(stormname(ient),
     &   advr(ient),iyear(ient),imth(ient),iday(ient),ihr(ient),
     &   iFcstInc(ient),Pn,Pc,cLat,cLon,VmaxBL,VrBL,r,
     &   rmaxw(ient,:),HollB(ient))

C    reset rmax to zero if there was a zero radius to the isotach for all
C    isotachs EXCEPT the 34 kt isotach.  in that case leave the radius that
C    has been substituted.
C       jgf49.30: The isotach wind speed can sometimes be zero in cases
C       where all radii are zero (this has been observed in the BEST
C       track file for IGOR2010). Including this possibility in the if
C       statement, so that we can avoid setting the quadrant Rmax to zero
C       if ivr was zero.
        DO i=1,4
           if ( (ivr(ient).ne.34).AND.(ivr(ient).ne.0).AND.
     &          (ir(ient,i).eq.0) ) then
              rmaxw(ient,i) = 0.0
           endif
        END DO

      ENDDO ! num_entry (main do loop)

C     !-------------------------------------
      ! Now indicate which isotach quadrant radius 
      ! that the user desires ADCIRC to read in 
      ! for the final calculation of RMX in the 
      ! Asymmetric Holland wind calculations
      !
      ! 34... - 0 0 0 0 ...
      ! 50... - 0 0 1 1 ...
      ! 64... - 1 1 0 0 ...
      !  
      ! would indicate -
      ! use NO radii from the 34 kt isotach
      ! use the 3 & 4 radii form the 50 kt isotach
      ! use the 1 & 2 radii form the 64 kt isotach  
      
      ! users can then modify the input file
      ! to indicate which set of radii to use
      ! for each cycle
      !
      !  Loop through each cycle and choose 
      !  the isotach radii to use
      ! 
      !  method 1
      !  use the 34kt isotach only (like original NWS=9)
      ! 
      !  method 2 
      !  use the fancy way of taking the highest 
      !  isotach Rmax that exists
      ! 
      !------------------------------------
      
      !------------------------------------
      ! SET METHOD HERE!
      !------------------------------------
C     method = 2  ! jgf49.30: this is now set at the top of the file 
                  ! and can be optionally selected via a command line option
      select case(method)
      case(1) ! just use the Rmaxes from the 34kt isotach
         do ient=1,num_entry
            if ((ivr(ient).eq.34).or.(ivr(ient).eq.0)) then
               quadflag(ient,:)=1   
            else
               quadflag(ient,:)=0
            endif
         enddo
      case(2) ! use the Rmax from the highest isotach in each quadrant
         do icyc=1,num_cycles
            last_entry = sum(isotachs_per_cycle(1:icyc))
            first_entry = last_entry-(isotachs_per_cycle(icyc)-1)
            if (isotachs_per_cycle(icyc).eq.1) then
               quadflag(last_entry,:)=1
            else 
               ! loop over quadrants
               do i=1,4
                  num_nonzero = 
     &               count(quadflag(first_entry:last_entry,i).ne.0)
                  select case(num_nonzero)
                     case(0,1) ! none, or only 34kt isotach has a radius value
                        quadflag(first_entry,i) = 1
                     case(2) ! the 34kt and 50kt isotachs have radius value
                        quadflag(first_entry,i) = 0
                     case(3) ! the 34kt, 50kt, and 64kt isotachs have values
                        quadflag(first_entry:first_entry+1,i) = 0
                     case default ! zero isotachs have been flagged
                        write(*,*) 
     &                  "ERROR: aswip: ",num_nonzero," isotachs were ",
     &                      "nonzero."
                  end select
               enddo
            endif
         enddo
      case default      
         write(*,*) "ERROR: method=",method,
     &      " is not valid for setting rmax in quadrants."
         write(*,*) "ERROR: Execution terminated."
      end select

      ! persist last good 34kt Rmax values if all radii are missing
      do icyc=1,num_cycles
         if (isotachs_per_cycle(icyc).eq.1) then
            ient = sum(isotachs_per_cycle(1:icyc))
            if ((all(ir(ient,:).eq.0)).and.(atcfRMW(ient).eq.0)) then
               if ((icyc-1).ge.1) then
                  rmaxw(ient,:) = 
     &                     rmaxw(ient-isotachs_per_cycle(icyc-1),:)
               else
                  rmaxw(ient,:) = 25 ! default value when all else fails
               endif
            endif
         endif
      enddo

C!-------------------------------------------------------------------

         !------------------------
         ! WRITE OUTPUT
         !------------------------
     
      DO ient=1,num_entry
         !-----------------------------------
         ! NHC advisory best-track i/o format
         !-----------------------------------
22       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 1x, i5,a1, 2x,
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10)
         !-----------------------------------
         ! NHC advisory best-track i/o format
         ! WITH nws=19 ADDED DATA
         !-----------------------------------
26       FORMAT("AL,",1x,i2.2,",",i5,3i2.2,",",3x,",",1x,a4,",",
     &         i4,",",i4,a1,",",i5,a1,",",i4,",",i5,",",3x,",",i4,",",
     &       " NEQ,",5(i5,","),5x,",",i4,", ",5(4x,","),i3,", ",i3,",",
     &         a10," ",1x,",",i4,",",i5,",",1x, 4(i1,",",1x),2x, 
     &         4(f6.1,",",1x), 2x, f8.4)!,2(2x,f10.4))


         ! jgf49.30: Modify the Rmax as specified on the command line.
         IF (constantRmax.eqv..true.) THEN
            DO i=1,4
               rmaxw(ient,i) = rMax
            ENDDO
         ELSEIF (percentSpecified.eqv..true.) THEN
            DO i=1,4
               rmaxw(ient,i) = rmaxw(ient,i)*(percent/100.d0)
            ENDDO
         ENDIF
C
         WRITE(222,26) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &      nint(h_dir(ient)),nint(h_speed(ient)), stormname(ient), 
     &      cycle_num(ient),isotachs_per_cycle(cycle_num(ient)),
     &     (quadflag(ient,i),i=1,4),(rmaxw(ient,i),i=1,4),HollB(ient)

      ENDDO

      close(222)

      END PROGRAM
