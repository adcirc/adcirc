!-----------------------------------------------------------------------
! rainfall.F
!   Written by Zachary Cobell, 2018/08/20
!              ARCADIS, INC.
!              Zachary.Cobell@Arcadis.com
!
!   Contains the routine to apply rainfall to the ADCIRC mesh
!   by addition to the ETA2 variable
!
!-----------------------------------------------------------------------
        MODULE RAINFALL
            USE WIND,ONLY:EyeLat,EyeLon,FoundEye,FoundRMax,
     &                    rMaxLon,rMaxLat,PRN2,SZ
            USE GLOBAL,ONLY: DTDP,ERROR,INFO,WARNING,scratchMessage,
     &                       allMessage,logmessage,openFileForRead,
     &                       setMessageSource,unsetMessageSource,DEBUG,
     &                       nnodecode,nodecode
            USE nodalattributes, ONLY : STARTDRY
            USE MESH,ONLY:NP,NE,ICS,DP
            IMPLICIT NONE

            !...public variables (set via namelist)
            LOGICAL              :: use_rain = .false.
            INTEGER              :: rainfallType = 0 !...Default to ipet
            REAL(SZ)             :: rainfallMultiplier = 1.0
            INTEGER              :: rainfallTimeIncrement = 1
            REAL(SZ)             :: rainfallThreshold = 0.0254D0    ! 0.0254 m = 1 inch
            REAL(SZ)             :: rainfallDepth = 0.1D0
            INTEGER              :: rainfallBlankSnaps = 0
            CHARACTER(200)       :: rainfallFile = ""

            NAMELIST /rainfallControl/ use_rain,rainfallMultiplier,
     &                    rainfallFile,rainfallType,rainfallTimeIncrement,
     &                    rainfallThreshold,rainfallDepth,rainfallBlankSnaps

            !...private variables
            INTEGER,PARAMETER    :: ipetRainfall     = 0
            INTEGER,PARAMETER    :: owiRainfall      = 1
            INTEGER,PARAMETER    :: combinedRainfall = 2

            INTEGER,PARAMETER    :: fid = 622
            INTEGER              :: iSnap,snapCounter
            LOGICAL              :: updateMapping

            INTEGER              :: iLatR,iLongR
            INTEGER              :: WDCHECK
            INTEGER              :: RAINCHECK
            INTEGER,ALLOCATABLE  :: rainCode(:)
            REAL(SZ)             :: dxR,dyR,swlatR,swlongR
            REAL(SZ)             :: RTIME1,RTIME2
            REAL(SZ),ALLOCATABLE :: rain1(:),rain2(:)
            REAL(SZ),ALLOCATABLE :: raingrid(:,:)
            REAL(SZ),ALLOCATABLE :: latR(:),longR(:)
            REAL(SZ),ALLOCATABLE :: swpointsR(:,:),wr(:,:)

            CHARACTER(200)       :: errorVar

            PRIVATE ipetRainfall,owiRainfall,fid,iSnap,snapCounter
            PRIVATE updateMapping,iLatR,iLongR,dxR,dyR,swLatR,swLongR
            PRIVATE rain1,rain2,raingrid,latr,longr,swpointsr,wr,errorvar
            PRIVATE computeOwiWeighting,computeOwiRainfall,computeIpetRainfall
            PRIVATE owiRainfallInit,radius_earth,haversine,owiRainfallGet
            PRIVATE readOwiRainfallHeader,owiRainfallMapping,check_err
            PRIVATE WDCHECK,RAINCHECK,rainCode

            CONTAINS

            
            SUBROUTINE rainfallInit(time)
                USE GLOBAL,ONLY:NWS
                USE WIND,ONLY:DragLawString
#ifdef CMPI
                USE MESSENGER, ONLY : MSG_FINI
#endif
                IMPLICIT NONE
                REAL(SZ),INTENT(IN) :: time

                call setMessageSource("rainfallInit")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                ALLOCATE(rainCode(1:NP))
                rainCode(:) = 0

                IF(rainfallType.EQ.ipetRainfall.OR.rainfallType.EQ.combinedRainfall)THEN
                    IF(NWS.EQ.0)THEN
                        call allMessage(ERROR,"Must use NWS with IPET rainfall model.")
#ifdef CMPI
                        call msg_fini()
#endif                        
                        STOP
                    ENDIF
                    IF(TRIM(DragLawString).NE."Powell".OR.
     &                   TRIM(DragLawString).EQ."powell".OR.
     &                   TRIM(DragLawString).EQ."POWELL")THEN
                        call allMessage(ERROR,"Must use Powell drag law string with IPET rainfall")
#ifdef CMPI
                        call msg_fini()
#endif
                        STOP
                    ENDIF
                ENDIF

                IF(rainfallType.EQ.owiRainfall.OR.
     &             rainfallType.EQ.combinedRainfall)THEN
                    CALL owiRainfallInit(time)
                    CALL owiRainfallGet()
                    CALL owiRainfallGet()
                ENDIF

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

            END SUBROUTINE rainfallInit


            SUBROUTINE applyRainfall(time)
                USE GLOBAL,ONLY:DTDP,ETA2,rainfallRate,rainfallTotal,
     &                          rainfallAccounted,rainfallUnaccounted,RAD2DEG,NOFF
                USE MESH,ONLY: X,Y,SFEA,SLAM,DP,NM
                IMPLICIT NONE           
                REAL(SZ),INTENT(IN) :: time
                INTEGER             :: i,j
                REAL(SZ)            :: radiusToMaxWinds
                REAL(SZ)            :: r1,r2,rain,weight

                call setMessageSource("applyRainfall")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                r1 = 0.0D0
                r2 = 0.0D0
                rain = 0.0D0

                IF(rainfallType.EQ.ipetRainfall.OR.
     &             rainfallType.EQ.combinedRainfall)THEN
                    radiusToMaxWinds = 
     &                  HAVERSINE(EyeLon(3),EyeLat(3),rMaxLon,rMaxLat)/1000D0
                ENDIF

                IF(rainfallType.EQ.owiRainfall.OR.
     &             rainfallType.EQ.combinedRainfall)THEN
                    call computeOwiWeighting(time,weight)
                ENDIF

                DO I = 1,NP
                    IF(rainfallType.EQ.ipetRainfall.OR.
     &                      rainfallType.EQ.combinedRainfall)THEN
                        IF(.NOT.foundEye.OR..NOT.foundRMax)THEN
                            r1 = 0D0
                        ELSE
                            IF(ICS.EQ.1)THEN
                                r1 = computeIpetRainfall(X(I),
     &                                  Y(I),PRN2(I),radiusToMaxWinds)
                            ELSE
                                r1 = computeIpetRainfall(SLAM(I)*RAD2DEG,
     &                                  SFEA(I)*RAD2DEG,PRN2(I),radiusToMaxWinds)
                            ENDIF
                        ENDIF
                    ENDIF

                    IF(rainfallType.EQ.owiRainfall.OR.
     &                     rainfallType.EQ.combinedRainfall)THEN
                        r2 = computeOWIRainfall(i,weight)
                    ENDIF

                    IF(rainfallType.EQ.ipetRainfall)THEN
                        rain = r1
                    ELSEIF(rainfallType.EQ.owiRainfall)THEN
                        rain = r2
                    ELSEIF(rainfallType.EQ.combinedRainfall)THEN
                        rain = MAX(r1,r2)
                    ENDIF
                    
                    rainfallTotal(I) = rainfallTotal(I)+RAIN                   !..m
                    rainfallRate(I)  = RAIN*3600D0*1000D0/DTDP                 !...mm/h

                    IF(STARTDRY(I).EQ.0)THEN    !...Don't collect water in leveed regions

                        IF(NNODECODE(I).EQ.1)THEN  !...Collect water only on "wet" nodes
                            ETA2(I)              = ETA2(I) + RAIN
                            rainfallAccounted(I) = rainfallAccounted(I) + RAIN     !..m
                            rainCode(I) = 0
                        ELSEIF(RAIN.GT.0.0D0)THEN
                            IF(rainfallUnaccounted(I).GT.rainfallThreshold)THEN
                                rainCode(I) = 1
                            ELSE
                                rainCode(I) = 0
                            ENDIF
                            rainfallUnaccounted(I) = rainfallUnaccounted(I) + RAIN !..m
                        ELSE
                            rainCode(I) = 0
                            rainfallUnaccounted(I) = rainfallUnaccounted(I) + RAIN !..m
                        ENDIF

                    ELSE

                        rainfallUnaccounted(I) = rainfallUnaccounted(I) + RAIN !..m

                    ENDIF
                        
                ENDDO

                DO I = 1,NE
                    IF(NOFF(I).EQ.0)THEN
                        WDCHECK = NNODECODE(NM(I,1)) + NNODECODE(NM(I,2)) + NNODECODE(NM(I,3))
                        IF(WDCHECK.NE.3)THEN
                            RAINCHECK = rainCode(NM(I,1)) + rainCode(NM(I,2)) + rainCode(NM(I,3))
                            IF(RAINCHECK.EQ.3)THEN
                                DO J = 1, 3
                                    ETA2(NM(I,J)) = ABS(DP(NM(I,J))) + rainfallDepth
                                    NNODECODE(NM(I,J)) = 1
                                ENDDO
                            ENDIF
                        ENDIF
                    ENDIF
                ENDDO

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

            END SUBROUTINE applyRainfall


            SUBROUTINE computeOwiWeighting(time,weight)
                IMPLICIT NONE
                REAL(SZ),INTENT(IN)  :: time
                REAL(SZ),INTENT(OUT) :: weight

                call setMessageSource("computeOwiWeighting")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                IF(Time.GT.RTIME2)THEN
                    CALL owiRainfallGet()
                ENDIF
                weight = (Time-RTime1)/rainfallTimeIncrement
                
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

            END SUBROUTINE computeOwiWeighting


            REAL(SZ) FUNCTION computeOwiRainfall(node,weight)
                IMPLICIT NONE
                INTEGER,INTENT(IN)  :: node
                REAL(SZ),INTENT(IN) :: weight
                call setMessageSource("computeOwiRainfall")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif
                computeOwiRainfall = (rain1(node)+weight*(rain2(node)-rain1(node))) *
     &            rainfallMultiplier*DTDP*(1D0/(3600D0*1000D0)) !... mm/hr -> m/s -> m/dt
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()
            END FUNCTION computeOwiRainfall


            REAL(8) FUNCTION computeIpetRainfall(X,Y,P,radiusToMaxWinds)
                !...Rainfall is computed in m/s
                !   This function uses the methodology
                !   from the IPET study in 
                !   Louisiana. The methodology approximates
                !   the rainfall rate based upon radius to maximum
                !   winds, pressure differential from 1013mb, distance
                !   to the storm center, and azimuth. It is only
                !   a very rough approximation to historical data and
                !   likely does not perform well in real world scenarios
                USE GLOBAL,ONLY:RAD2DEG
                USE WIND,ONLY:StormDirection
                IMPLICIT NONE
                REAL(8),INTENT(IN) :: X
                REAL(8),INTENT(IN) :: Y
                REAL(8),INTENT(IN) :: P
                REAL(8),INTENT(IN) :: radiusToMaxWinds
                REAL(8)            :: dp
                REAL(8)            :: r
                REAL(8)            :: multiplier
                REAL(8)            :: dTheta
                REAL(8)            :: NodeDirection

                call setMessageSource("computeIpetRainfall")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                !...Calculate the node direction
                NodeDirection = 90D0-(ATAN2(Y-EyeLat(3),X-EyeLon(3)))*RAD2DEG
                IF(NodeDirection.LT.0.D0)THEN
                   NodeDirection = NodeDirection + 360.D0
                ENDIF
                
                !...Calculate which side of the storm track
                !   this node falls on
                dTheta = StormDirection - NodeDirection
                IF(dTheta.LT.0D0)THEN
                    dTheta = dTheta + 360D0
                ENDIF

                !...Select the wind multiplier. 
                !   -> Right side of track => 1.5
                !   -> Left side of track  => 1.0
                IF(dTheta.GT.180D0)THEN
                    multiplier=1.5D0
                ELSE
                    multiplier=1.0D0
                ENDIF
                
                !...Pressure difference in millibars
                dp = 1013D0 - (p*98.0638D0)

                !...km distance to storm center
                r  = haversine(eyeLon(3),eyeLat(3),X,Y)/1000D0

                !...Check against radius to maximum winds
                !   and then compute rainfall rate in mm/h
                IF(r.LT.radiusToMaxWinds)THEN
                     computeIpetRainfall = 1.14D0+0.12D0*dp*multiplier
                ELSE
                    computeIpetRainfall = (1.14D0+0.12D0*dp) * 
     &                  EXP(-0.3D0*((r-radiusTomaxWinds)/radiusToMaxWinds))*multiplier
                ENDIF

                !...Convert mm/h to m/s
                computeIpetRainfall = computeIpetRainfall * (1D0/(3600D0*1000D0))

                !...Compute rainfall using multiplier and adcirc time step
                computeIpetRainfall = computeIpetRainfall*DTDP*rainfallMultiplier 

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                RETURN

            END FUNCTION computeIpetRainfall


            REAL(8) FUNCTION RADIUS_EARTH(LATITUDE)
                IMPLICIT NONE 

                REAL(8),PARAMETER  :: ER = 6378137.0D0
                REAL(8),PARAMETER  :: PR = 6356752.3D0
                REAL(8),INTENT(IN) :: LATITUDE
                REAL(8)            :: LAT
                REAL(8)            :: DEG2RAD
                REAL(8)            :: PI

                call setMessageSource("radius_earth")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                PI      = 4D0*ATAN2(1D0,1D0)
                DEG2RAD = PI / 180D0
                LAT     = LATITUDE * DEG2RAD

                RADIUS_EARTH = SQRT(                      
     &              (ER**4D0*COS(LAT)*COS(LAT) +        
     &               PR**4D0*SIN(LAT)*SIN(LAT))          
     &              /                                    
     &              (ER**2D0*COS(LAT)*COS(LAT) +         
     &               PR**2D0*SIN(LAT)*SIN(LAT) ) )

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                RETURN
            
            END FUNCTION RADIUS_EARTH
            
            REAL(8) FUNCTION HAVERSINE(LON1d,LAT1d,LON2d,LAT2d)
                IMPLICIT NONE
                REAL(8),INTENT(IN) :: LON1d
                REAL(8),INTENT(IN) :: LAT1d
                REAL(8),INTENT(IN) :: LON2d
                REAL(8),INTENT(IN) :: LAT2d
                REAL(8)            :: LON1
                REAL(8)            :: LAT1
                REAL(8)            :: LON2
                REAL(8)            :: LAT2
                REAL(8)            :: RADIUS
                REAL(8)            :: PI
                REAL(8)            :: DEG2RAD

                call setMessageSource("haversine")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                PI       = 4D0*ATAN2(1D0,1D0)
                DEG2RAD  = PI / 180D0
                LAT1     = LAT1d*DEG2RAD 
                LAT2     = LAT2d*DEG2RAD
                LON1     = LON1d*DEG2RAD
                LON2     = LON2d*DEG2RAD
                RADIUS   = RADIUS_EARTH((LAT1+LAT2)/2D0)

                HAVERSINE = 2D0*RADIUS*ASIN(SQRT(                      
     &               SIN( (LAT2-LAT1)/2D0 )**2D0 +                      
     &               COS(LAT1)*COS(LAT2)*SIN( (LON2-LON1)/2D0)**2D0 ) )

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                RETURN

            END FUNCTION HAVERSINE
            
            
            SUBROUTINE owiRainfallInit(time)
                USE SIZES,ONLY: gblInputDir
                IMPLICIT NONE
                REAL(SZ),INTENT(IN) :: time
                CHARACTER(200)      :: filename
                INTEGER             :: errorIO
                INTEGER             :: I
                INTEGER             :: numSkipSnaps

                call setMessageSource("owiRainfallInit")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                filename = TRIM(GBLINPUTDIR)//"/"//TRIM(rainfallFile)
                call openFileForRead(fid,filename,errorIO)
                errorVar = "opening rainfall field"
                call check_err(errorIO)
                call readOwiRainfallHeader()

                if(allocated(swpointsR))deallocate(swpointsR)
                if(allocated(wr))deallocate(wr)

                allocate(swpointsR(np,2))
                allocate(wr(np,4))
                allocate(rain1(np))
                allocate(rain2(np))

                isnap = 0
                updateMapping = .TRUE.
                snapCounter = 0
                rtime2 = time 

                IF(rainfallBlankSnaps.LT.0)THEN
                    numSkipSnaps = ABS(rainfallBlankSnaps)
                    rainfallBlankSnaps = 0
                ELSE
                    numSkipSnaps = 0
                ENDIF

                DO I=1,numSkipSnaps
                    write(scratchMessage,41) i
                    call allMessage(DEBUG,scratchMessage)
                    CALL owiRainfallGet()
                ENDDO

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

41              format("Skipping snap '",I6,"' in rainfall data.")                    

            END SUBROUTINE owiRainfallInit


            SUBROUTINE owiRainfallGet()
                IMPLICIT NONE
                
                INTEGER       :: errorIO
                INTEGER       :: iCYMDH
                INTEGER       :: iMin
                INTEGER       :: i,j
                INTEGER       :: xi,yi
                INTEGER       :: iLat,iLong
                REAL(SZ)      :: dx,dy,swlat,swlong
                REAL(SZ)      :: rr,w1,w2,w3,w4

                call setMessageSource("owiRainfallGet")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                snapCounter = snapCounter + 1
                IF(snapCounter.GT.1)THEN
                    RTIME1=RTIME2
                    RTIME2=RTIME2+rainfallTimeIncrement
                ELSE
                    RTIME1=0.0D0
                    RTIME2=0.0D0
                ENDIF
                RAIN1(:)=RAIN2(:)
                
                IF(snapCounter.LE.rainfallBlankSnaps)THEN
                    rain2(:) = 0D0
                    write(scratchMessage,10) snapCounter
10                  format('Inserting a blank rainfall snap, count = ',I4)
                    call allMessage(INFO,trim(scratchMessage))
                    call unsetMessageSource()
                    return
                ENDIF

                iSnap = iSnap + 1
      
                errorVar = "grid specifications/date in rainfall file"
                read (fid,11,end=10000,err=9999,iostat=errorIO)
     &             iLat,iLong,dx,dy,swlat,swlong,iCYMDH,iMin
                call check_err(errorIO)

                IF(iLat.NE.iLatR.OR.iLong.NE.iLongR.OR.dx.NE.dxR.OR.
     &             dy.NE.dyR.OR.swlat.NE.swlatR.OR.swlong.NE.swlongR)THEN
                    updateMapping = .TRUE.
                    IF(iSnap.GT.1)THEN
                        CALL allMessage(INFO,"Rainfall grid has changed.")
                    ENDIF
                    iLatR = iLat
                    iLongR = iLong
                    dxR = dx
                    dyR = dy
                    swlatR = swlat
                    swlongR = swlong
                ENDIF

                IF(updateMapping)call owiRainfallMapping()

                write(scratchMessage,15) iCYMDH,iMin
15              format("Processing rainfall data ",I12," ",I2,".")
                call allMessage(INFO,scratchMessage)

                errorVar="rainfall snapshot"
                read(fid,22,end=10000,err=9999,iostat=errorIO)
     &              ((raingrid(i,j),i=1,iLongR),j=1,iLatR)
                call check_err(errorIO)

                do i =1,np
                    if(swpointsR(i,1).gt.0)then
                        xi=swpointsR(i,1)
                        yi=swpointsR(i,2)
                        w1=wR(i,1)
                        w2=wR(i,2)
                        w3=wR(i,3)
                        w4=wR(i,4)
                        rr=w1*raingrid(xi,yi)+w2*raingrid(xi+1,yi)+
     &                      w3*raingrid(xi+1,yi+1)+w4*raingrid(xi,yi+1)
                    else
                        rr = 0.0
                    endif

                    rain2(i)=rr

                enddo
                
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                return

11              format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &                 t44,f8.0,t58,f8.0,t69,i10,i2)
22              format(8f10.0)

9999            call check_err(1) ! ERR during read jumps to here

                ! END during read of required data jumps to here
9997            call allMessage(ERROR,"Unexpected end-of-file reached.")
                call check_err(1)

10000           continue  ! END during read of u, v, p data jumps to here
                write(scratchMessage,'(A)') "Unexpected end-of-file while reading rainfall"//
     &                                      " Precipitation set to zero."
                call allMessage(WARNING,trim(scratchMessage))

                rain2 = 0.0D0
                


            END SUBROUTINE owiRainfallGet


            SUBROUTINE owiRainfallMapping()
                USE GLOBAL,ONLY:RAD2DEG
                USE MESH,ONLY: slam,sfea,x,y
                IMPLICIT NONE

                INTEGER  :: i,j,k
                INTEGER  :: xi,yi
                REAL(SZ) :: adcLong,adcLat
                REAL(SZ) :: w,w1,w2,w3,w4

                call setMessageSource("owiRainfallMapping")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif

                if(allocated(raingrid))deallocate(raingrid)
                if(allocated(longR))deallocate(longR)
                if(allocated(latR))deallocate(latR)
                
                allocate(raingrid(iLongR,iLatR))
                allocate(latR(iLatR))
                allocate(longR(iLongR))

                do i=1,iLatR
                    latR(i)=swLatR+(i-1)*dyR
                enddo

                do i=1,iLongR
                    longR(i)=swLongR+(i-1)*dxR
                enddo

                do i=1,np
                    IF(ICS.EQ.2)THEN
                        adcLat  = RAD2DEG*SFEA(i)
                        adcLong = RAD2DEG*SLAM(i)
                    ELSE
                        adcLat  = y(i)
                        adcLong = x(i)
                    ENDIF
                    
                    if(adcLong.GE.longR(1).and.adcLong.LT.longR(iLongR).and.
     &                 adcLat.GE.latR(1).and.adcLat.LT.latR(ilatR))then
                        
                        do j =1,iLongR-1
                            if(adcLong.GE.longR(j).and.adcLong.LT.longR(j+1))then
                                xi=j
                                exit
                            endif
                        enddo
                        
                        do k=1,iLatR-1
                            if(adcLat.GE.latR(k).and.adcLat.LT.latR(k+1))then
                                yi=k
                                exit
                            endif
                        enddo

                        swpointsR(i,1)=xi
                        swpointsR(i,2)=yi
                        w=(longR(xi+1)-longR(xi))*(latR(yi+1)-latR(yi))
                        w1=(longR(xi+1)-adcLong)*(latR(yi+1)-adcLat)
                        w2=(adcLong-longR(xi))*(latR(yi+1)-adcLat)
                        w3=(adcLong-longR(xi))*(adcLat-latR(yi))
                        w4=(longR(xi+1)-adcLong)*(adcLat-latR(yi))

                        wR(i,1)=w1/w
                        wR(i,2)=w2/w
                        wR(i,3)=w3/w
                        wR(i,4)=w4/w

                    else
                        swpointsR(i,1)=0
                        swpointsR(i,2)=0
                    endif
                enddo

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                return

            END SUBROUTINE owiRainfallMapping


            SUBROUTINE readOwiRainfallHeader()
                IMPLICIT NONE
                CHARACTER(80) :: owiheader
                INTEGER       :: date1,date2
                INTEGER       :: errorIO

                call setMessageSource("readOwiRainfallHeader")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter")
#endif
                
                owiHeader(:) = " "

                errorVar="rainfall header"
                read(fid, fmt='(a80)',
     &              end=99998,err=99999,iostat=errorIO) owiheader
                call check_err(errorIO)
                read(owiheader(56:65),'(I10)',
     &              end=99998,err=99999,iostat=errorIO) date1
                call check_err(errorIO)
                read(owiheader(71:80),'(I10)',
     &              end=99998,err=99999,iostat=errorIO) date2
                call check_err(errorIO)

                write(scratchMessage,'(A,I10,A,I10,A)') 
     &               "Rainfall data found for dates: ",date1," to ",date2,"."
                call allMessage(INFO,TRIM(scratchMessage))

#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return")
#endif
                call unsetMessageSource()

                RETURN

99998           call allMessage(ERROR,"Unexpectedly reached end-of-file.") 
99999           call check_err(1) 

            END SUBROUTINE readOwiRainfallHeader

            subroutine check_err(iret)
#ifdef CMPI
                USE MESSENGER, ONLY : MSG_FINI
#endif
                IMPLICIT NONE
                INTEGER, intent(in) :: iret
                call setMessageSource("check_err")
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Enter.")
#endif
                if (iret.ne.0) then
                    if (errorVar.ne."") then
                        write(scratchMessage,888)
     &                  trim(errorVar), trim(rainfallFile), fid
888                     format("Failed to read '",A,"' from '",A,
     &                      "' (unit number ",I3,").")
                        call allMessage(ERROR,trim(scratchMessage))
                    endif
                    call allMessage(ERROR,"ADCIRC execution terminated.")
#ifdef CMPI
                    call msg_fini()
#endif
                    stop
                endif
#if defined(ALL_TRACE)
                call allMessage(DEBUG,"Return.")
#endif
                call unsetMessageSource()
                RETURN
            end subroutine check_err

        END MODULE RAINFALL
