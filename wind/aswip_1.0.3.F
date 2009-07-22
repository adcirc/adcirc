C***********************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                     *
C last changes in this file VERSION 45.05                              *
C                                                                      *
C***********************************************************************
C
       PROGRAM COMPUTE_RMAX_B

         !-------------------------------------------------------------
         ! Import custom precision types for cross-platform portability
         !-------------------------------------------------------------
         USE precision

         !------------------------
         ! Import global constants
         !------------------------
         USE constants

         !-----------------------------------------
         ! Import asymmetric hurricane vortex class
         !-----------------------------------------
         USE vortex

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         IMPLICIT NONE

      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====      PROGRAM ASYMMETRIC_RMAX_B            =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

         INTEGER  :: np, ics,ient,icyc,itpc,itpc2,method
         REAL(sz) :: time, time1,time2, wtratio,timest

         INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: iyear, imth, iday, ihr
         INTEGER , SAVE :: iyear1, imth1, iday1, ihr1
         INTEGER , SAVE :: iyear2, imth2, iday2, ihr2
         INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ilat,ilon
         INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ispd,icpress,atcfRMW
         INTEGER , SAVE :: ilat1,ilon1, ispd1, icpress1!, irmw1
         INTEGER , SAVE :: ilat2,ilon2, ispd2, icpress2!, irmw2
         REAL(sz), SAVE ::  lat1, lon1,  spd1,  cpress1!,  rmw1
         REAL(sz), SAVE ::  lat2, lon2,  spd2,  cpress2!,  rmw2
         REAL(sz), SAVE ::  dt, cLat1, cLon1
         CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE:: type! hindcast/nowcast or forecast?
         CHARACTER(LEN=4) :: type1              ! hindcast/nowcast or forecast?
         CHARACTER(LEN=4) :: type2              ! hindcast/nowcast or forecast?
         INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::iFcstInc ! hours between forecasts
         INTEGER,SAVE::iFcstInc1, iFcstInc2 ! hours between forecasts
         INTEGER , SAVE :: firstCall = 0
         LOGICAL , SAVE :: firstTime = .TRUE.

         REAL(sz),DIMENSION(:),ALLOCATABLE,SAVE :: HollB
         REAL(sz), SAVE :: HollB1
         INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: advr
         INTEGER , SAVE :: adv
         REAL(sz), SAVE :: Vmax,speedmax,Rmax,Rmaxmax
         REAL(sz), SAVE :: Vr,speedmax2
         REAL(sz), SAVE, DIMENSION(4) :: r
         REAL(sz), SAVE :: Pn
         REAL(sz), SAVE :: Pc
         REAL(sz), SAVE :: cLat
         REAL(sz), SAVE :: cLon
         REAL(sz), SAVE :: offset
         REAL(sz), SAVE :: dx,dy

         INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: ivr, dir, speed
         INTEGER , SAVE :: ivr1, dir1, speed1
         INTEGER , SAVE :: ivr2, dir2, speed2
         INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE :: ir
         INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE :: quadflag
         INTEGER , DIMENSION(4), SAVE :: ir1,ir2
         REAL(sz) , DIMENSION(:,:), ALLOCATABLE,SAVE :: rmaxw
         REAL(sz) , DIMENSION(4), SAVE :: rmaxw1
         INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: count_newcycle
         INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: cycle_num
         INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle
         INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: ipn
         INTEGER , SAVE :: ipn1,ipn2
         CHARACTER (LEN = 10),DIMENSION(:),ALLOCATABLE, SAVE :: name
         REAL(sz), PARAMETER :: outFreq = thirtySixHundred
!        REAL(sz), PARAMETER :: outFreq = one
         INTEGER :: i, nwi,screenunit,num_cycles,numdiffcyc,num_entry

         REAL(sz),SAVE :: clatOld,clonOld, clatNew,clonNew,
     &                    timeOld,timeNew
         REAL(sz) :: uTrans1 = zero
         REAL(sz) :: vTrans1 = zero
         REAL(sz) ,DIMENSION(:),ALLOCATABLE, SAVE :: uTrans, vTrans

         INTEGER , PARAMETER :: node2print = 6205
         CHARACTER*1 ew,E,W,ns,N,S
   
          screenunit=6

         E='E'
         W='W'
         N='N'
         S='S'
        !------------
        ! open output file
        !------------
         OPEN(222,file="NWS_19_fort.22",status="unknown")

         OPEN(22,file="fort.22",status="old")
            i=0
         DO 
            READ(22,*,END=991) 
            i=i+1
         ENDDO    
991    continue
           num_entry=i
        close(22)
         OPEN(22,file="fort.22",status="old")

              i=num_entry 
         allocate(advr(i),iyear(i),imth(i),iday(i),ihr(i),type(i))
         allocate(iFcstInc(i),ilat(i),ilon(i), ispd(i))
         allocate(icpress(i), ivr(i),ir(i,4), rmaxw(i,4), ipn(i))
         allocate(dir(i),speed(i),count_newcycle(i),cycle_num(i))
         allocate(isotachs_per_cycle(i))
         allocate(quadflag(i,4),name(i))
         allocate(uTrans(i), vTrans(i),HollB(i),atcfRMW(i) )
      
          ient=1
            READ(22,22) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), type(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient),atcfRMW(ient),
     &      dir(ient),speed(ient), name(ient)

          num_cycles = 1
          count_newcycle(1) = 1
          cycle_num(ient)=count_newcycle(1)
         DO ient=2,num_entry

            READ(22,22) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), type(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &      dir(ient),speed(ient), name(ient)

  
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
             

       ENDDO


         !-------------------------------
         ! Calculate translation velocity
         !-------------------------------
          DO icyc = 1,num_cycles
        if (icyc .lt. num_cycles ) then
           isotachs_per_cycle(icyc) = 
     &                       count_newcycle(icyc+1)-count_newcycle(icyc)
          DO ient =  count_newcycle(icyc), count_newcycle(icyc+1)-1
              
           write(*,*) icyc, ient
        if (icyc==1 ) then
            clatOld =  ilat(count_newcycle(icyc+1)-1) * oneTenth
            clatNew =  ilat(count_newcycle(icyc+1)) * oneTenth
            clonOld = -ilon(count_newcycle(icyc+1)-1) * oneTenth
            clonNew = -ilon(count_newcycle(icyc+1)) * oneTenth

            time1 =  iFcstInc(count_newcycle(icyc+1)-1) * hour2sec
            time2 =  iFcstInc(count_newcycle(icyc+1)) * hour2sec

               dx = deg2rad * Rearth * (clonNew - clonOld)
     &            * COS(deg2rad*half * (clatOld + clatNew))
               dy = deg2rad * Rearth * (clatNew - clatOld)
               dt = (time2 - time1)
               if(dt == 0.0d0) dt=6.0d0* hour2sec
         else
            clatOld =  ilat(count_newcycle(icyc)-1) * oneTenth
            clatNew =  ilat(count_newcycle(icyc+1)) * oneTenth
            clonOld = -ilon(count_newcycle(icyc)-1) * oneTenth
            clonNew = -ilon(count_newcycle(icyc+1)) * oneTenth

            time1 =  iFcstInc(count_newcycle(icyc)-1) * hour2sec
            time2 =  iFcstInc(count_newcycle(icyc+1)) * hour2sec

               dx = deg2rad * Rearth * (clonNew - clonOld)
     &            * COS(deg2rad*half * (clatOld + clatNew))
               dy = deg2rad * Rearth * (clatNew - clatOld)
               dt = (time2 - time1)
               if(dt == 0.0d0) dt=12.0d0* hour2sec
         endif
               uTrans(ient) = dx/dt
               vTrans(ient) = dy/dt
          ENDDO

        ELSE  ! icyc = num_cycles
        isotachs_per_cycle(icyc) =
     &            num_entry - count_newcycle(icyc) + 1
            
           DO ient = count_newcycle(icyc), num_entry
            clatOld =  ilat(count_newcycle(icyc)-1) * oneTenth
            clatNew =  ilat(count_newcycle(icyc)) * oneTenth
            clonOld = -ilon(count_newcycle(icyc)-1) * oneTenth
            clonNew = -ilon(count_newcycle(icyc)) * oneTenth

            time1 =  iFcstInc(count_newcycle(icyc)-1) * hour2sec
            time2 =  iFcstInc(count_newcycle(icyc)) * hour2sec

               dx = deg2rad * Rearth * (clonNew - clonOld)
     &            * COS(deg2rad*half * (clatOld + clatNew))
               dy = deg2rad * Rearth * (clatNew - clatOld)
               dt = (time2 - time1)
               if(dt == 0.0d0) dt=6.0d0* hour2sec
               uTrans(ient) = dx/dt
               vTrans(ient) = dy/dt
            ENDDO ! ient
          ENDIF  ! icyc = num_cycles
          ENDDO ! numcycles

          DO ient=1,num_entry
          ENDDO
            DO  icyc = 1,num_cycles
            ENDDO

         !-----------------------------------------
         !  Now using the calculated translational velocities 
         !  call the vortex module and compute the Rmax's
         !  to be used in the new input file        
         !----------------------------------------------
            
         DO ient=1,num_entry
      !
      !  need to get some logic incase Vr is  zero
      !
         if (ipn(ient) == 0 )then
             ipn(ient)=ipn(ient-1)
          endif

            adv    = advr(ient)
            iyear1 = iyear(ient)
            imth1  = imth(ient)
            iday1  = iday(ient)
            ihr1   = ihr(ient)
            type1  = type(ient)
            iFcstInc1 = iFcstInc(ient)
            ilat1  = ilat(ient)
            ilon1  = ilon(ient)
            ispd1  = ispd(ient)
            icpress1 = icpress(ient)
            ivr1   = ivr(ient)
            ir1(:) = ir(ient,:)
            ipn1   = ipn(ient)
            dir1   = dir(ient)
            speed1 = speed(ient)
            uTrans1 =uTrans(ient)
            vTrans1 =vTrans(ient)
! 
!         The logic for the code is set for W and N so if E or S
!         multiply lats and lons by -1
!
            IF(ew.EQ.E)THEN
               ilon1=(-1)*ilon1 
               ilon2=(-1)*ilon2
            ENDIF
            IF(ns.EQ.S)THEN
               ilat1=(-1)*ilat1
               ilat2=(-1)*ilat2
            ENDIF
 
           !--------------------------------------------------------------
           ! transform variables from integers
           ! to real numbers for hurricane vortex calcualtions.
           !--------------------------------------------------------------
          Vmax =  one * ( ispd1  )
          Pn   =  one * ( ipn1   )
          Pc   =  one * ( icpress1  )
          cLat =  ilat1 * oneTenth 
          cLon = -ilon1 * oneTenth 

      !---------------------------------------------
      !  need to get some logic incase Vr is  zero
      !  if so we will also be setting ir(:) to Rmax 
      !---------------------------------------------
          if ( ivr1 == 0 ) ivr1 = Vmax

          Vr = one* ivr1

C         GO TO 3333
         
      !-------------------------------------------------------
      !  need to get some logic incase any of radii are zero
      !  in the end we will keep the adjustments to the 34kt
      !  isotach...and reset the other isotachs to their initial
      !  zero value...
      !  We will then indicate which quadrant radii from which
      !  isotach we will have ADCIRC use in calculating the 
      !  Asymmetric Holland wind
      !-------------------------------------------------------
          if (ir1(1) .eq. 0 ) then                !rad1 = 0
             quadflag(ient,1)=0 
      ! comput a radius from nonzero entries
              if (ir1(2) .ne. 0 ) then            !rad1 = 0,rad2 != 0
                quadflag(ient,2)=1 
                    if (ir1(3) .ne. 0 ) then      !rad1 = 0,rad2 != 0,rad3 != 0
                     quadflag(ient,3)=1 
                        ir1(1) = 0.5 * (ir1(2)+ir1(3))
                    elseif (ir1(4) .ne. 0 ) then  !rad1 = 0,rad2 != 0,rad3 = 0,rad4 != 0
                     quadflag(ient,3)=0 
                     quadflag(ient,4)=1 
                        ir1(1) = 0.5 * (ir1(2)+ir1(4))
                        ir1(3) = 0.5 * (ir1(2)+ir1(4))
                    else                          !rad1 = 0,rad2 != 0,rad3 = 0,rad4 = 0
                     quadflag(ient,3)=0 
                     quadflag(ient,4)=0 
                        ir1(1) = (ir1(2))
                        ir1(3) = (ir1(2))
                        ir1(4) = (ir1(2))
                    endif                      
              elseif (ir1(3) .ne. 0 ) then    !rad1 = 0,rad2 = 0,rad3 != 0 
                  quadflag(ient,2)=0
                  quadflag(ient,3)=1
                    if (ir1(4) .ne. 0 ) then  !rad1 = 0,rad2 = 0,rad3 != 0,rad4 != 0
                      quadflag(ient,4)=1
                       ir1(1) = 0.5 * (ir1(3)+ir1(4))  
                       ir1(2) = 0.5 * (ir1(3)+ir1(4))
                    else                      !rad1 = 0,rad2 = 0,rad3 != 0,rad4 = 0
                      quadflag(ient,4)=0
                        ir1(1) = (ir1(3))
                        ir1(2) = (ir1(3))
                        ir1(4) = (ir1(3))
                    endif     
              elseif (ir1(4) .ne. 0 ) then     !rad1 = 0,rad2 = 0,rad3 = 0,rad4 != 0
                      quadflag(ient,2)=0
                      quadflag(ient,3)=0
                      quadflag(ient,4)=1
                      ir1(1) = (ir1(4))
                      ir1(2) = (ir1(4))
                      ir1(3) = (ir1(4))
              else                             !rad1 = 0,rad2 = 0,rad3 = 0,rad4 = 0
              ! will only happen if no isotachs are reported in that caseuse RMW and Vmax
                      quadflag(ient,1)=1
                      quadflag(ient,2)=1
                      quadflag(ient,3)=1
                      quadflag(ient,4)=1
                      ir1(1) = (ir1(4))
                     write(*,*) "NO NON ZERO QUAD RADII"
                     write(*,*) "USE RMW from ATCF"
                       do i=1,4
                        ir1(i) = atcfRMW(ient)
                       enddo
              endif
          
          elseif (ir1(2) .eq. 0 ) then       !rad1 != 0,rad2 = 0
                      quadflag(ient,1)=1
                      quadflag(ient,2)=0
               if (ir1(3) .ne. 0 ) then      !rad1 != 0,rad2 = 0,rad3 != 0
                      quadflag(ient,3)=1
                     ir1(2) =  0.5 * (ir1(1)+ir1(3))
               elseif( ir1(4) .ne. 0 ) then  !rad1 != 0,rad2 = 0,rad3 = 0 ,rad4 != 0
                      quadflag(ient,3)=0
                      quadflag(ient,4)=1
                      ir1(2) =  0.5 * (ir1(1)+ir1(4))
                      ir1(3) =  0.5 * (ir1(1)+ir1(4))
               else                          !rad1 != 0,rad2 = 0,rad3 = 0 ,rad4 = 0
                      quadflag(ient,3)=0
                      quadflag(ient,4)=0
                      ir1(2) = (ir1(1))
                      ir1(3) = (ir1(1))
                      ir1(4) = (ir1(1))
               endif

          elseif (ir1(3) .eq. 0 ) then       !rad1 != 0,rad2 != 0,rad3 = 0
                      quadflag(ient,1)=1
                      quadflag(ient,2)=1
                      quadflag(ient,3)=0
              if (ir1(4) .ne. 0 ) then       !rad1 != 0,rad2 != 0,rad3 = 0 ,rad4 != 0
                      quadflag(ient,4)=1
                    ir1(3) =  0.5 * (ir1(2)+ir1(4))
              else                           !rad1 != 0,rad2 != 0,rad3 = 0 ,rad4 = 0
                      quadflag(ient,4)=0
                    ir1(3) =  0.5 * (ir1(1)+ir1(2))
                    ir1(4) =  0.5 * (ir1(1)+ir1(2))
              endif                  

           elseif (ir1(4) .eq. 0 ) then      !rad1 != 0,rad2 != 0,rad3 = !0 ,rad4 = 0
                      quadflag(ient,1)=1
                      quadflag(ient,2)=1
                      quadflag(ient,3)=1
                      quadflag(ient,4)=0
                    ir1(4) =  0.5 * (ir1(2)+ir1(3))

           else 
                     
                      quadflag(ient,1)=1
                      quadflag(ient,2)=1
                      quadflag(ient,3)=1
                      quadflag(ient,4)=1
                  write(*,*) " all 4 nonzero radii ! "

           endif
C3333    continue
              
 
            DO i=1,4
               r(i) = one * ir1(i)
            END DO
            DO i=1,4
               rmaxw1(i) = atcfRMW(ient)
            END DO
                HollB1=one
         !-------------------------------------------------------
         ! Create a new asymmetric hurricane vortex.
         !
         ! Note: Subtract translational speed from Vmax, then
         ! scale (Vmax - Vt) and Vr up to the top of the surface,
         ! where the cylcostrophic wind balance is valid.
         !-------------------------------------------------------
       CALL newVortex_rmax_b_OUT(
     &   name(ient/2),adv,iyear1,imth1,iday1,ihr1,iFcstInc1,
     &   Pn,Pc,cLat,cLon,((Vmax-ms2kt*
     &   SQRT(uTrans1*uTrans1+vTrans1*vTrans1))/windReduction),
     &   (Vr/windReduction),r,rmaxw1,HollB1)
         !---------------------------
         ! Print the hurricane vortex
         !---------------------------

                  HollB(ient)=HollB1
                   write(*,'(f10.4)')HollB(ient) 
                  rmaxw(ient,:)=rmaxw1(:)
                   write(*,'(f10.4)') rmaxw(ient,:)
C    reset rmax to zero if there was a zero radius to the isotach for all isotachs
C    EXCEPT the 34 kt isotach.  in that case leave the radius that has been substituted.

            DO i=1,4
              if ( ivr(ient) .ne. 34 .AND. ir(ient,i) .eq.  0 ) then
               rmaxw(ient,i) = 0.0
              endif

            END DO

          ENDDO ! num_entry (main do loop)



C       !-------------------------------------
        ! Now indicate which isotach quadrant radius 
        ! that the user desires ADCIRC to read in 
        ! for teh final calculation of RMX in the 
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
        !
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
           method = 1

          
         ient = 0
       DO icyc = 1,num_cycles
       IF (method == 1) then 
               ient = ient+isotachs_per_cycle(icyc)
          IF ( isotachs_per_cycle(icyc) .gt. 1) then
            DO i = 1,4  ! quadrants
               DO itpc=1,isotachs_per_cycle(icyc)-1 ! isotachs
                    quadflag((ient-(itpc-1)),i) = 0
               ENDDO
            ENDDO
          ELSE
            DO i = 1,4  ! quadrants
               quadflag((ient),i) = 1
            enddo
          ENDIF  ! isotachs per cycle > 1

       ELSEIF (method == 2) then
             ient = ient+isotachs_per_cycle(icyc)
            !-----------------------
            ! only check isotachs greater than 34 kt
            !----------------------
          IF ( isotachs_per_cycle(icyc) .gt. 1) then 
                  
            DO i = 1,4  ! quadrants

               DO itpc=1,isotachs_per_cycle(icyc) ! isotachs
 
                if ( quadflag((ient-(itpc-1)),i) == 1 ) then ! the highest vel isotach has a value here so zero out the others

                     DO itpc2=itpc,isotachs_per_cycle(icyc)-1

                    quadflag((ient-itpc2),i) = 0

                     ENDDO
                 GOTO 9992
                 
                 endif

              ENDDO  ! isotachs
9992      CONTINUE
              ENDDO ! quardants
          ELSE
            
            DO i = 1,4  ! quadrants
               quadflag((ient),i) = 1
            enddo          
          ENDIF  ! isotachs per cycle > 1
       ELSE ! method type if statement 
       write(*,*)"not a recognized method for setting rmax in quadrants"
       ENDIF ! method type if statement
       ENDDO    ! num_cycles        
C!-------------------------------------------------------------------
          !------------------------
          ! WRITE OUTPUT
          !------------------------
     
        DO ient=1,num_entry
         !-----------------------------------
         ! NHC advisory best-track i/o format
         !-----------------------------------
22      FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 3x, i3,a1, 2x,
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10)

26      FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 3x, i3,a1, 2x,
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10, 2x, i4, 2x, i4,
     &         2x, 4(i1,2x),2x, 4(f6.1,2x), 2x, f8.4)

            WRITE(222,26) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), type(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &      dir(ient),speed(ient), name(ient), 
     &      cycle_num(ient),isotachs_per_cycle(cycle_num(ient)),
     &     (quadflag(ient,i),i=1,4),(rmaxw(ient,i),i=1,4),HollB(ient)

         ENDDO

          close(222)
 


      END PROGRAM

