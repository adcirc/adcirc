C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      MODULE Couple2BC3D
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      USE SIZES,  ONLY: SZ, MNEI
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage, NB, NM, MLD, SigTS,
     &          VIDBCPDXOH, VIDBCPDYOH, VIDISPDXOH, VIDISPDYOH, IDEN
      USE CONSTANTS, ONLY: SigT0
#ifdef DATETIME
      use datetime_module, only: datetime, timedelta, strptime
#endif
      use nodalattributes, only: LoadIT_Fric, NOLIBF, HBREAK 
      implicit none

!     WJP notes on variables and interpolation.
!     McDougall & Parker (2014) (TEOS-10 creators) note that 
!     interpolation of nonconservative variables such as in-situ temp
!     (T), in-situ density (S), potential density, etc, is not good 
!     oceanographic practice. Instead we should interpolate Conservative
!     Temperature (theta) and Absolute Salinity (SA) before then 
!     calculating N (buoyancy frequency) and rho (in-situ density) on 
!     our grid. Thus, we read in T and S from the third-party 3D
!     baroclinic model, interpolate these to our grid and in the
!     temporal dimension. To save space we only read in two timesteps at
!     a time and interpolate between these. SA and theta are then
!     converted to N and rho (or sigma_T) as required. 

      CHARACTER(LEN=100) :: densityFileName ! Name of the BC3D NETCDF file (fort.11.nc by default)
      INTEGER  :: NX, NY, NYc, NYs, NXc, NYY, NYYY, NXo, NT, SL ! dimensions of 3D BC model required
      INTEGER  :: XI, YI, XIc, YIc, YIs  ! the 'start' locations in the NETCDF file
      REAL(SZ) :: BC3D_DT    ! Delta T for the 3D Baroclinic inputs
      INTEGER  :: BC3D_IT    ! 3D BC temporal iteration number 
      INTEGER  :: NP11       ! Number of nodes when reading from fort.11
      INTEGER  :: densityTimeIterator  ! BC3D time step iteraton skipper (usually 1)
      real(sz),allocatable,dimension(:) :: BC3D_Lon, BC3D_Lat, BXE, BYE,
     &                         DXE, DYE, BC3D_Lats, BC3D_Lonc, BC3D_Latc
      real(sz),allocatable,dimension(:,:) :: BC3D_NM, BC3D_NB, BC3D_BX,
     &          BC3D_SigTS, BC3D_BY, BC3D_MLD, BC3D_CD, BC3D_DX, BC3D_DY
      logical,allocatable,dimension(:,:) :: NSkip
      real(sz),parameter :: LatUL = 89d9, DFV = 1d0 
#ifdef DATETIME
      type(datetime) :: CurDT
#endif
!     Variables for temporal integration
      real(sz) :: DTIME1, DTIME2 
      real(sz),allocatable,dimension(:,:) :: VIDBCPDXOH_F, VIDBCPDYOH_F,
     &      SigTS_F, NB_F, NM_F, MLD_F, CD_F, VIDISPDXOH_F, VIDISPDYOH_F
      integer,allocatable  :: indm(:,:), indx(:,:), indy(:,:), inds(:,:)
      real(sz),allocatable :: weightsm(:,:), weightsx(:,:),
     &                        weightsy(:,:), weightss(:,:)

C-----------------------------------------------------------------------
      CONTAINS

C-----------------------------------------------------------------------
C...  Allocate space for Arrays dimensioned by MNP
C... 
      SUBROUTINE ALLOC_BC3D_to_2D()
      use sizes, only : mnp, mne
      IMPLICIT NONE
C  
      allocate ( VIDBCPDXOH_F(2,MNP),VIDBCPDYOH_F(2,MNP))
      VIDBCPDXOH_F = 0.0d0 ; VIDBCPDYOH_F = 0.0d0
      allocate (SigTS_F(2,MNP))
      SigTS_F = SigT0 
      allocate ( MLD_F(2,MNP), CD_F(2,MNP))
      MLD_F = DFV 
      CD_F = 0.0d0
      allocate ( VIDISPDXOH_F(2,MNE),VIDISPDYOH_F(2,MNE))
      VIDISPDXOH_F = 0.0d0 ; VIDISPDYOH_F = 0.0d0
      allocate (NB_F(2,MNP),NM_F(2,MNP))
      NB_F = 0.0d0; NM_F = 0.0d0 
      IF (IDEN.LT.8) THEN
         allocate(indm(4,MNP),weightsm(4,MNP))
         allocate(indx(4,MNP),indy(4,MNP),inds(4,MNP))
         allocate(weightsx(4,MNP),weightsy(4,MNP),weightss(4,MNP))
      ENDIF
      allocate(BXE(MNP),BYE(MNP),DXE(MNP),DYE(MNP))
      BXE = 0.0d0; BYE = 0.0d0; DXE = 0.0d0; DYE = 0.0d0;
C-----------------------------------------------------------------------
      END SUBROUTINE ALLOC_BC3D_to_2D
C-----------------------------------------------------------------------
C
#ifdef ADCNETCDF
C-----------------------------------------------------------------------
C     S U B R O U T I N E  R E A D _ B C 3 D _ N E T C D F
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Read_BC3D_NetCDF()
      use netcdf, only: nf90_open, nf90_close, nf90_nowrite 
      implicit none
      integer :: nc_id
 
      call setMessageSource("Read_BC3D_NetCDF")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      WRITE(16,*) 'inside read_bc3d_netcdf'
      ! Open NETCDF 
      call Check_err(nf90_open(densityFileName,nf90_nowrite,nc_id))
    
      ! On first call: 
      ! Get the dimensions and time and spatial arrays of the data,
      ! and calculate the interpolant coefficients
      ! On > first call: check the datetimes and update BC3D_DT
      call Get_LonLatDepthTime(nc_id)
      IF (abs(IDEN).GE.8) THEN
         WRITE(16,*) 'IDEN >= 8'
         CALL Read_BC3D_NetCDF_on_ADCIRC_Grid(nc_id)
         ! Close NETCDF
         call Check_err(nf90_close(nc_id))
#ifdef ALL_TRACE
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      ENDIF
      if (abs(iden).ne.7) then     
         ! Read all the necessary temporal varying data 
         ! Depth-averaged east-west baroclinic pressure graidient
         BC3D_BX = read_nc_var(NC_ID,'BPGX',BC3D_IT,XIc,NXc,YI,NY)
         ! Depth-averaged north-south baroclinic pressure graidient
         BC3D_BY = read_nc_var(NC_ID,'BPGY',BC3D_IT,XI,NX,YIc,NYc)
         ! Free surface sigmat density
         BC3D_SigTS = read_nc_var(NC_ID,'SigTS',BC3D_IT,XI,NX,YI,NY)
      endif
      !
      if (LoadIT_Fric) then
         ! Buoyancy frequency at the seabed
         BC3D_NB = read_nc_var(NC_ID,'NB',BC3D_IT,XI,NX,YI,NY)
         ! Depth-averaged nuoyancy frequency
         BC3D_NM = read_nc_var(NC_ID,'NM',BC3D_IT,XI,NX,YI,NY)
      endif
      !
      if (NOLIBF.eq.3) then
         if (HBREAK.gt.0d0) then
            ! Mixed-layer Depth
            BC3D_MLD = read_nc_var(NC_ID,'MLD',BC3D_IT,XI,NX,YI,NY)
         else
            ! CD due to dispersion
            BC3D_CD = read_nc_var(NC_ID,'CDisp',BC3D_IT,XI,NX,YIs,NYs)
            where(BC3D_CD > 1.d0) BC3D_CD = 1.d0
         endif
      endif
      if (abs(iden).eq.6) then
         ! Depth-averaged east-west baroclinic pressure graidient
         BC3D_DX = read_nc_var(NC_ID,'DispX',BC3D_IT,XI,NX,YIs,NYs)
         ! Depth-averaged north-south baroclinic pressure graidient
         BC3D_DY = read_nc_var(NC_ID,'DispY',BC3D_IT,XI,NX,YIs,NYs)
      endif
      ! Close NETCDF
      call Check_err(nf90_close(nc_id))

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      
C-----------------------------------------------------------------------
      END SUBROUTINE Read_BC3D_NetCDF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Get_LonLatDepthTime(NC_ID)
      use netcdf
      use mesh, ONLY: slam, sfea, binarysearch, bl_interp, np, ne, nm,
     &                conserve_interp, neitabele, nodele !indmax
      use global, ONLY: rad2deg, np_g, ne_g
      implicit none
      integer,intent(in) :: NC_ID
      integer  :: Temp_ID, i, j, ic, jc, js, indt(4), numfound !ind(2,indmax), cc
      real(sz) :: Lon_s, Lat_s, Lat_e, Lon_e, xx, yy, wt(4), 
     &            Lonn_s, Lonn_e, x1, y1, xl, yl, xu, yu !, x2, y2, x3, y3 !, ww(indmax)
      real(sz),allocatable :: Lonc(:), Lon(:), Latc(:), Lat(:), Lats(:)
#ifdef DATETIME       
      type(timedelta) :: DT
      type(datetime)  :: TNDT, TNDT1
      character(len=16) :: TSS
#endif        

      if (BC3D_IT > 1) then
         WRITE(16,*) 'BC3D_IT = ',BC3D_IT
         ! Check times if not first iteration
#ifdef DATETIME
         call Check_err(NF90_INQ_VARID(NC_ID,'time',Temp_ID))
         ! Get current and next date times and update BC3D_DT
         call check_err(nf90_get_var(nc_id, Temp_ID,        
     &                  TSS,start=[1, BC3D_IT],count=[SL, 1]))
         TNDT = strptime(trim(TSS),"%Y-%m-%d %H:%M")
         WRITE(16,*) 'TNDT = ',TSS
         call check_err(nf90_get_var(nc_id, Temp_ID, TSS,
     &      start=[1, BC3D_IT+densityTimeIterator],count=[SL, 1]))
         TNDT1 = strptime(trim(TSS),"%Y-%m-%d %H:%M")
         WRITE(16,*) 'TNDT1 = ',TSS
         DT = TNDT1 - TNDT 
         if (abs(DT%total_seconds() - BC3D_DT) > 1d-6) then
            BC3D_DT = DT%total_seconds()
            write(16,*) 'BC3D_DT is changed to ',BC3D_DT
         endif
#endif
         return
      endif
      ! First call
      call check_err(nf90_inq_dimid(nc_id,'time', Temp_ID))
      call Check_err(nf90_inquire_dimension(nc_id,Temp_ID,len=NT))
      call check_err(nf90_inq_dimid(nc_id,'strlen', Temp_ID))
      call Check_err(nf90_inquire_dimension(nc_id,Temp_ID,len=SL))
      ! CPB 9/2022: moved the timestep check here to limit how much
      ! rewriting of this subroutine I need to do.
      ! Get the timestep and starting time
      call Check_err(NF90_INQ_VARID(NC_ID,'time',Temp_ID))
      BC3D_DT = 0d0
      if (iden.gt.0) then
#ifdef DATETIME
      ! Let's get the start time index
      do BC3D_IT = 1,NT
         call check_err(nf90_get_var(nc_id, Temp_ID,        
     &                  TSS,start=[1, BC3D_IT],count=[SL, 1]))
         TNDT = strptime(trim(TSS),"%Y-%m-%d %H:%M")
         DT   = CurDT - TNDT
         if (DT%total_seconds() <= 0) exit
      enddo
      call check_err(nf90_get_var(nc_id, Temp_ID, TSS,           
     &     start=[1, BC3D_IT+densityTimeIterator],count=[SL, 1]))
      TNDT1 = strptime(trim(TSS),"%Y-%m-%d %H:%M")
      DT = TNDT1 - TNDT 
      BC3D_DT = DT%total_seconds()
#else
      ! Just set at 3 hrs
      BC3D_DT = 3600d0*3d0
#endif
      endif
      ! CPB 2/2022: in case where baroclinic forces are calculated on
      ! ADCIRC grid (IDEN = 8) we don't need to check the nx, ny, etc
      ! but we do want to verify that the input file has the correct
      ! number of nodes and such.
      IF (IDEN.GE.8) THEN
         ! check if the number of nodes is correct
         call Check_err(NF90_INQ_DIMID(NC_ID,'node',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
         IF (NY.NE.NP_G) THEN
            WRITE(16,*) 'ERROR: Baroclinic forcing input (fort.11.nc) '
     &                //'does not have the correct number of nodes for '
     &                //'IDEN = 8. Check your input files to make sure '
     &                //'you have BC3D forcing for each node in the '
     &                //'model or change which baroclinic forcing '
     &                //'method you are using.'
            call BC3DTerminate()
         ENDIF
         ! check if the number of elements is correct
         call Check_err(NF90_INQ_DIMID(NC_ID,'nele',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
         IF (NY.NE.NE_G) THEN
            WRITE(16,*) 'ERROR: Baroclinic forcing input (fort.11.nc) '
     &                //'does not have the correct number of elements'
     &                //' for IDEN = 8. Check your input files to make '
     &                //'sure you have BC3D forcing for each node in '
     &                //'the model or change which baroclinic forcing '
     &                //'method you are using.'
            call BC3DTerminate()
         ENDIF
         return
      ENDIF
      call Check_err(NF90_INQ_DIMID(NC_ID,'NY',Temp_ID))
      call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
      call Check_err(NF90_INQ_DIMID(NC_ID,'NX',Temp_ID))
      call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))
      if (abs(IDEN).ne.7) then
         call Check_err(NF90_INQ_DIMID(NC_ID,'NYY',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NYY))
      endif
      if (abs(IDEN).eq.6) then 
         call Check_err(NF90_INQ_DIMID(NC_ID,'NYYY',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NYYY))
      endif   
      ! Allocate the lat and lon  arrays first 
      allocate(BC3D_Lat(NY),BC3D_Lon(NX),
     &         BC3D_Latc(NYY),BC3D_Lonc(NX),BC3D_Lats(NYYY))

      ! Read the latitude, longitude variables 
      call Check_err(NF90_INQ_VARID(NC_ID,'lat',Temp_ID))
      call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lat))
      call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
      call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lon))
      if (abs(IDEN).ne.7) then
         call Check_err(NF90_INQ_VARID(NC_ID,'latc',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Latc))
         call Check_err(NF90_INQ_VARID(NC_ID,'lonc',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lonc))
      endif
      if (abs(IDEN).eq.6) then 
         call Check_err(NF90_INQ_VARID(NC_ID,'lats',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lats))
      endif     
 
      ! Get bounds of the mesh
      Lat_s = max(-79.99d0,rad2deg*minval(sfea))
      Lat_e = rad2deg*maxval(sfea)
      Lon_s = rad2deg*minval(slam)
      Lon_e = min(179.91d0,rad2deg*maxval(slam))
      ! Update the longitude if it passes the 180/-180
      if (Lon_s < 0d0 .and. Lon_e >= 0d0) then
         Lonn_e = rad2deg*minval(slam,slam >= 0d0)
         Lonn_s = rad2deg*maxval(slam,slam <= 0d0)
         ! If the difference is smaller than the original orientation
         if (360d0 - (Lonn_e - Lonn_s) < Lon_e - Lon_s) then 
            Lon_s = Lonn_e; Lon_e = Lonn_s;
         endif
      endif

      ! Correct mesh bounds if BC3D is 0 to 360 
      Lonn_s = Lon_s; Lonn_e = Lon_e;
      if (BC3D_Lon(NX).gt.180d0) then
         if (Lon_s < 0d0) Lon_s = Lon_s + 360d0;
         if (Lon_e < 0d0) Lon_e = Lon_e + 360d0;
         if (Lon_s > Lon_e) then
            Lonn_s = Lon_e; Lonn_e = Lon_s;
         else
            Lonn_s = Lon_s; Lonn_e = Lon_e;
         endif
      endif
      
      ! Test the bounds 
      if (BC3D_Lon(1) > Lonn_s.or.BC3D_Lat(1) > Lat_s.or.
     &    Lonn_e > BC3D_Lon(NX).or.Lat_e > BC3D_Lat(NY)) then
          call allMessage(Error,'Part of 3DBC domain is not '
     &                   //'contained within this subdomain. '
     &                   //'Try download again')
          write(16,*) BC3D_Lon(1), Lon_s
          write(16,*) BC3D_Lon(NX), Lon_e
          write(16,*) BC3D_Lat(1), Lat_s
          write(16,*) BC3D_Lat(NY), Lat_e
          call BC3DTerminate()
      else
          ! Now find only the portions of lat lon required for this
          ! sub-domain and update NX, NY
          ! For longitude (with buffer)
          NXo = NX
          XI  = max(1,binarysearch(NXo, BC3D_Lon, Lon_s) - 1)
          i   = min(NXo,binarysearch(NXo, BC3D_Lon, Lon_e) + 2) 
          NX  = i - XI + 1; if (NX <  0) NX = NX + NXo;
          
          ! For latitude (with buffer)
          YI  = max(1,binarysearch(NY, BC3D_Lat, Lat_s) - 1)
          j   = min(NY,binarysearch(NY, BC3D_Lat, Lat_e) + 2)
          NY  = j - YI + 1

          write(16,*) 'BC3D_IT = ',BC3D_IT, 'BC3D_DT = ',BC3D_DT
          write(16,*) 'Lon: ',BC3D_Lon(XI), BC3D_Lon(i)
          write(16,*) 'Lat: ',BC3D_Lat(YI), BC3D_Lat(j)
          write(16,*) 'XI = ',XI, 'i = ', i, 'NX = ',NX
          write(16,*) 'YI = ',YI, 'j = ', j, 'NY = ',NY
          allocate(Lon(NX),Lat(NY))
          Lat  = BC3D_Lat(YI:j)
          if (i > XI) then
             Lon  = BC3D_Lon(XI:i)
          else
             Lon  = [BC3D_Lon(XI:NXo)  , BC3D_Lon(1:i)   + 360d0 ]
          endif
          
          if (abs(IDEN).ne.7) then
             XIc = max(1,binarysearch(NXo, BC3D_Lonc, Lon_s) - 1)
             ic  = min(NXo,binarysearch(NXo, BC3D_Lonc, Lon_e) + 2) 
             NXc = ic - XIc + 1; if (NXc <  0) NXc = NXc + NXo;

             YIc = max(1,binarysearch(NYY, BC3D_Latc, Lat_s) - 1)
             jc  = min(NYY,binarysearch(NYY, BC3D_Latc, Lat_e) + 2)
             NYc = jc - YIc + 1
             write(16,*) 'Lonc: ',BC3D_Lonc(XIc), BC3D_Lonc(ic)
             write(16,*) 'Latc: ',BC3D_Latc(YIc), BC3D_Latc(jc)
             write(16,*) 'XIc = ',XIc, 'ic = ', ic, 'NXc = ',NXc
             write(16,*) 'YIc = ',YIc, 'jc = ', jc, 'NYc = ',NYc
             allocate(Lonc(NXc),Latc(NYc))
             allocate(BC3D_BX(NXc,NY),BC3D_BY(NX,NYc),BC3D_SigTS(NX,NY))
             Latc = BC3D_Latc(YIc:jc)
             if (ic > XIc) then 
                Lonc = BC3D_Lonc(XIc:ic)
             else
                Lonc = [BC3D_Lonc(XIc:NXo), BC3D_Lonc(1:ic) + 360d0 ]
             endif
          endif
          if (abs(IDEN).eq.6) then 
             YIs = max(1,binarysearch(NYYY, BC3D_Lats, Lat_s) - 1)
             js  = min(NYYY,binarysearch(NYYY, BC3D_Lats, Lat_e) + 2)
             NYs = js - YIs + 1
             allocate(Lats(NYs))
             allocate(BC3D_DX(NX,NYs),BC3D_DY(NX,NYs))
             Lats = BC3D_Lats(YIs:js)
          endif
  
          if (LoadIT_Fric) allocate(BC3D_NB(NX,NY),BC3D_NM(NX,NY))
          if (NOLIBF.eq.3) then
             if (HBREAK.gt.0d0) then
                allocate(BC3D_MLD(NX,NY))
             else
                allocate(BC3D_CD(NX,NYs)) 
             endif
          endif
          
          ! Compute the interpolants
          do i = 1,NP  
             xx = rad2deg*slam(i)
             ! Convert our numbers if grids are 0 to 360
             if (minval(Lon) >= 0d0.and.xx < 0d0) xx = xx + 360d0
             yy = rad2deg*sfea(i)
             ! For non-gradient variables just use linear interp
             call bl_interp(NX,Lon,NY,Lat,xx,yy,indt,wt)
             indm(:,i) = indt
             weightsm(:,i) = wt
             ! For gradient variables lets use conserve interp
             if (abs(iden).ne.7) then
                j = 0; numfound = 0
                xl = xx; xu = xx
                yl = yy; yu = yy
                do j = 1,MNEI
                   if (NEITABELE(i,j).ne.0) then
                   numfound = numfound + 1
                   ! Get neighboring element mid-point
                   x1 = (rad2deg*slam(nm(NEITABELE(i,j),1))
     &                + rad2deg*slam(nm(NEITABELE(i,j),2))
     &                + rad2deg*slam(nm(NEITABELE(i,j),3)))/3d0
                   if (minval(Lon) >= 0d0) then
                      if (x1 < 0d0) x1 = x1 + 360d0
                   endif
                   y1 = (rad2deg*sfea(nm(NEITABELE(i,j),1))
     &                + rad2deg*sfea(nm(NEITABELE(i,j),2))
     &                + rad2deg*sfea(nm(NEITABELE(i,j),3)))/3d0
                   ! Set upper/low bounds
                   xl = min(x1,xl); xu = max(x1,xu)
                   yl = min(y1,yl); yu = max(y1,yu)
                   endif
                enddo  
                ! For x-gradient variable (BC3D_BX)
                ! CPB 12/08/2021 switch from conserve_interp to bl_interp for proper 
                ! wrapping around the prime meridian
                call bl_interp(NXc,Lonc,NY,Lat,xx,yy,indt,wt)
                indx(:,i) = indt
                weightsx(:,i) = wt
                ! For y-gradient variable (BC3D_BY)
                ! CPB 12/08/2021 switch from conserve_interp to bl_interp for proper 
                ! wrapping around the prime meridian
                call bl_interp(NX,Lon,NYc,Latc,xx,yy,indt,wt) 
                indy(:,i) = indt
                weightsy(:,i) = wt
                if (abs(iden).eq.6) then
                   ! For other momentum dispersion variables (BC3D_DX, BC3D_DY)
                   call conserve_interp(NX,Lon,NYs,Lats,
     &                                  xx,yy,xl,yl,xu,yu,indt,wt) 
                   inds(:,i) = indt
                   weightss(:,i) = wt
                endif
             endif
          enddo
          deallocate(Lon,Lat)
          if (abs(IDEN).ne.7) deallocate(Lonc,Latc)
          if (abs(IDEN).eq.6) deallocate(Lats)
      endif
      END SUBROUTINE Get_LonLatDepthTime
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Read_BC3D_NetCDF_on_ADCIRC_Grid(NC_ID)
      use netcdf
      use mesh, only : np,sfea,slam 
      use sizes, only : mnproc,myproc, localdir
      use global, only : nodes_lg, np_g
#ifdef CMPI
      USE GL2LOC_MAPPING, ONLY : MAPTOLOCAL_REAL
#endif      
      use global, ONLY: rad2deg
      implicit none
      INTEGER,INTENT(IN) :: NC_ID
      ! global arrays
      REAL(SZ),ALLOCATABLE,DIMENSION(:) :: bpgx_g, bpgy_g, sigts_g, 
     &                                     nb_g, nm_g, mld_g
      ! local array
      REAL(SZ),ALLOCATABLE,DIMENSION(:) :: tmp
      ! temporary looping variable to write out values to see if this is
      ! working
      INTEGER :: ii
      ! on processor 0 we have to read in the values
      IF (MYPROC.EQ.0) THEN
         WRITE(16,*) 'Proc 0 reading in vals'
         ALLOCATE( bpgx_g(np_g), bpgy_g(np_g), sigts_g(np_g) )
         bpgx_g  = read_nc_var_adc(NC_ID,'BPGX',BC3D_IT,np_g)
         bpgy_g  = read_nc_var_adc(NC_ID,'BPGY',BC3D_IT,np_g)
         sigts_g = read_nc_var_adc(NC_ID,'SigTS',BC3D_IT,np_g)
         IF (LoadIT_Fric) THEN
            ALLOCATE( nb_g(np_g), nm_g(np_g) )
            nb_g = read_nc_var_adc(NC_ID,'NB',BC3D_IT,np_g)
            nm_g = read_nc_var_adc(NC_ID,'NM',BC3D_IT,np_g)
         ENDIF
         IF (NOLIBF.EQ.3) THEN
            ALLOCATE( mld_g(np_g) )
            mld_g = read_nc_var_adc(NC_ID,'MLD',BC3D_IT,np_g)
         ENDIF
      ENDIF
#ifdef CMPI
      ! use mapToSubdomainRealMPI to distribute to the processors
      ! allocate dummy array
      ALLOCATE( tmp(NP) )
      tmp = 0d0
      IF (ABS(IDEN).NE.9) THEN
         ! BPGX
         CALL MAPTOLOCAL_REAL(bpgx_g,tmp)
         VIDBCPDXOH_F(2,:) = tmp
         ! BPGY
         CALL MAPTOLOCAL_REAL(bpgy_g,tmp)
         VIDBCPDYOH_F(2,:) = tmp
         ! SigTS
         CALL MAPTOLOCAL_REAL(sigts_g,tmp)
         SigTS_F(2,:)      = tmp
      ENDIF
      IF (LoadIT_Fric) THEN
         ! NB
         CALL MAPTOLOCAL_REAL(nb_g,tmp)
         NB_F(2,:) = tmp
         ! NM
         CALL MAPTOLOCAL_REAL(nb_g,tmp)
         NM_F(2,:) = tmp
      ENDIF
      IF (NOLIBF.EQ.3) THEN
         ! MLD
         CALL MAPTOLOCAL_REAL(mld_g,tmp)
         MLD_F(2,:) = tmp
      ENDIF
      DEALLOCATE( tmp )
#endif
      IF (MYPROC.EQ.0) THEN
         DEALLOCATE( bpgx_g, bpgy_g, sigts_g )
         IF (LoadIT_Fric) THEN
            DEALLOCATE( nb_g, nm_g )
         ENDIF
         IF (NOLIBF.EQ.3) THEN
            DEALLOCATE( mld_g )
         ENDIF
      ENDIF
      END SUBROUTINE Read_BC3D_NetCDF_on_ADCIRC_Grid
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION read_nc_var(NC_ID,varname,BCIT_IN,XIin,NXin,YIin,NYin)
     & result (Var)
      use netcdf
         
      implicit none

      integer,intent(in) :: NC_ID, BCIT_IN, XIin, NXin, YIin, NYin
      character(*),intent(in) :: varname
      integer :: Temp_ID, Xe 
      real(sz) :: Var(NXin,NYin)

      call Check_err(NF90_INQ_VARID(NC_ID,trim(varname),Temp_ID))
      Xe = NXo - XIin + 1;
      if (NXin > Xe) then
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var(1:Xe,:),
     &                  start=[XIin,YIin,BCIT_IN],count=[Xe,NYin,1]))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var(Xe+1:NXin,:),
     &                  start=[1,YIin,BCIT_IN],count=[NXin-Xe,NYin,1]))
      else
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var,
     &                  start=[XIin,YIin,BCIT_IN],count=[NXin,NYin,1]))
      endif
C
C-----------------------------------------------------------------------
      END FUNCTION read_nc_var
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION read_nc_var_adc(NC_ID,varname,BCIT_IN,NumNodes)
     & result (Var)
      use netcdf
         
      implicit none

      integer,intent(in) :: NC_ID, BCIT_IN, NumNodes
      character(*),intent(in) :: varname
      integer :: Temp_ID
      real(sz) :: Var(NumNodes)

      call Check_err(NF90_INQ_VARID(NC_ID,trim(varname),Temp_ID))
      call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var,
     &               start=[1,BCIT_IN],count=[NumNodes,1]))
C
C-----------------------------------------------------------------------
      END FUNCTION read_nc_var_adc
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
C     jgf49.17.02 Checks the return value from netCDF calls; if there
C     was an error, it writes the error message to the screen and to the
C     fort.16 file.
C-----------------------------------------------------------------------
      subroutine check_err(iret)
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      USE NETCDF
      IMPLICIT NONE

      INTEGER, intent(in) :: iret

      !call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call BC3DTerminate()
      endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     S U B R O U T I N E  I N I T I A L _ B C 3 D _ N E T C D F
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Initial_BC3D_NetCDF(TimeLoc)
#ifdef DATETIME
      use global, only: basedatetime
#endif
      implicit none
      real*8,intent(in) :: TimeLoc
C
      !call setMessageSource("Initial_BC3D_NetCDF")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! Read first value in time 
      BC3D_IT = 1
      ! Times are rounded down to nearest minute
      DTIME1  = floor(TimeLoc/60d0)*60d0
#ifdef DATETIME
      CurDT   = basedatetime + timedelta(minutes=floor(DTIME1/60d0))
#endif           
      WRITE(16,*) 'About to read_bc3d_netcdf in initial_...' 
      call Read_BC3D_Netcdf()
      ! Put the first value on ADCIRC grid
      ! CPB 9/2022: added logical for IDEN case 8
      IF (abs(IDEN).LT.8) THEN
         call Put_BC3D_on_ADCIRC_Grid() 
      ENDIF
      call logMessage(ECHO,'Read and processed 1st BC3D '
     &                         //'timesnap'
#ifdef DATETIME
     &                         //': '//CurDT%isoformat(' ')
#endif               
     &                         ) 
      ! WJP Get the start time and the next data time
      ! Put the first values in first position
      if (abs(IDEN).ne.7.AND.abs(IDEN).NE.9) then
         VIDBCPDXOH_F(1,:) = VIDBCPDXOH_F(2,:) 
         VIDBCPDYOH_F(1,:) = VIDBCPDYOH_F(2,:) 
         SigTS_F(1,:)      = SigTS_F(2,:)
      endif
      if (LoadIT_Fric) then
         NB_F(1,:)      = NB_F(2,:)
         NM_F(1,:)      = NM_F(2,:)
      endif
      if (NOLIBF.eq.3) then
         if (HBREAK.gt.0d0) then
            MLD_F(1,:)  = MLD_F(2,:)
         else
            CD_F(1,:)   = CD_F(2,:)
         endif
      endif
      if (abs(IDEN).eq.6) then
         VIDISPDXOH_F(1,:) = VIDISPDXOH_F(2,:) 
         VIDISPDYOH_F(1,:) = VIDISPDYOH_F(2,:) 
      endif
      ! Return if diagnostic
      if (IDEN.lt.0) return
      ! Read next value in time
      BC3D_IT = BC3D_IT + densityTimeIterator
#ifdef DATETIME
      CurDT = CurDT + timedelta(minutes=nint(BC3D_DT/60d0))
#endif
      CALL Read_BC3D_Netcdf()
      ! BC3D_DT may have been updated in Read_BC3D..
      DTIME2 = DTIME1 + BC3D_DT
      ! Put the second value on ADCIRC grid
      IF (abs(IDEN).LT.8) THEN
         call Put_BC3D_on_ADCIRC_Grid() 
      ENDIF
      call logMessage(ECHO,'Read and processed 2nd BC3D '
     &                         //'timesnap'
#ifdef DATETIME
     &                         //': '//CurDT%isoformat(' ')
#endif               
     &                         ) 
C
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE Initial_BC3D_NetCDF
C
#endif
C-----------------------------------------------------------------------
C     S U B R O U T I N E  U P D A T E _ B C 3 D _ I N F O
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Update_BC3D_Info(flag,TimeLoc)

      use mesh, only: np, ne, uvectmp, vvectmp, DRVMAP2DSPVEC
      use global, only: IFSPROTS
      implicit none
      
      integer, intent(in) :: flag
      real(8),intent(in)  :: TimeLoc 
      real(sz)            :: DTRATIO
      integer             :: J,JKI 
C
      !call setMessageSource("Update_BC3D_Info")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! Check if prognostic calculation
      IF (BC3D_DT > 0.d0) THEN
         ! Get new times if time is right
         IF (TimeLoc.GT.DTIME2) THEN
            ! Updating new time and arrays
            DTIME1 = DTIME2
            if (IDEN.ne.7.AND.abs(IDEN).NE.9) then
               SigTS_F(1,:)      = SigTS_F(2,:)
               VIDBCPDXOH_F(1,:) = VIDBCPDXOH_F(2,:) 
               VIDBCPDYOH_F(1,:) = VIDBCPDYOH_F(2,:) 
            endif
            if (LoadIT_Fric) then
               NB_F(1,:)      = NB_F(2,:)
               NM_F(1,:)      = NM_F(2,:)
            endif
            if (NOLIBF.eq.3) then
               if (HBREAK.gt.0d0) then
                  MLD_F(1,:)  = MLD_F(2,:)
               else
                  CD_F(1,:)   = CD_F(2,:)
               endif
            endif
            if (IDEN.eq.6) then
               VIDISPDXOH_F(1,:) = VIDISPDXOH_F(2,:) 
               VIDISPDYOH_F(1,:) = VIDISPDYOH_F(2,:) 
            endif
            IF (flag > 0) THEN
#ifdef ADCNETCDF
               BC3D_IT = BC3D_IT + densityTimeIterator
#ifdef DATETIME
               CurDT = CurDT + timedelta(minutes=nint(BC3D_DT/60d0))
#endif
               ! Read from netcdf 
               CALL Read_BC3D_NetCDF()
               
               ! BC3D_DT may have been changed above..
               DTIME2 = DTIME1 + BC3D_DT
               IF (IDEN.LT.8) THEN
                  ! Calculate the BPG, free surface density, buoyancy
                  ! frequency, and dispersion terms on the ADCIRC grid
                  CALL Put_BC3D_on_ADCIRC_Grid()
               ENDIF
               WRITE(16,*) 'Read and processed new BC3D timesnap: '
     &                     //CurDT%isoformat(' ')
!               call logMessage(ECHO,'Read and processed new BC3D '
!     &                         //'timesnap'
!#ifdef DATETIME
!     &                         //': '//CurDT%isoformat(' ')
!#endif               
!     &                         ) 
#endif
            ELSEIF (flag < 0) THEN
               ! Read next values from fort.11
               DO J = 1,NP11
                  READ(11,*) JKI,VIDBCPDXOH_F(2,JKI),
     &                       VIDBCPDYOH_F(2,JKI),SigTS_F(2,JKI)
               END DO
               DTIME2 = DTIME1 + BC3D_DT
            ENDIF
         ENDIF
         ! Update to the new DT RATIO
         DTRATIO = (TimeLoc-DTIME1)/BC3D_DT
      ELSE
         DTRATIO = 0.0d0
      ENDIF

      ! Linear interpolate all the values in time
      if (abs(IDEN).ne.7.AND.abs(IDEN).NE.9) then
         DO J = 1,NP
           SigTS(J) = SigTS_F(1,J) + DTRATIO*(SigTS_F(2,J)-SigTS_F(1,J))
           ! This is interpolated onto element and put back onto node in
           ! another subroutine (FBPG_Disp_from_BC3D)
           BXE(J) = VIDBCPDXOH_F(1,J) + DTRATIO * 
     &              (VIDBCPDXOH_F(2,J) - VIDBCPDXOH_F(1,J))
           BYE(J) = VIDBCPDYOH_F(1,J) + DTRATIO * 
     &              (VIDBCPDYOH_F(2,J) - VIDBCPDYOH_F(1,J))
         ENDDO
C        DW/WJP: rotate BCP if coordinate transform
         IF ( IFSPROTS .eq. 1 ) THEN
            UVECTMP(1:NP) = BXE(1:NP); 
            VVECTMP(1:NP) = BYE(1:NP);
            CALL DRVMAP2DSPVEC( BXE, BYE,
     &             UVECTMP, VVECTMP, NP, FWD = .TRUE. ) ;
         ENDIF
      endif
      IF (LoadIT_Fric) then
         DO J = 1,NP
            ! This is used for the internal tide energy conversion
            NB(J) = NB_F(1,J) + DTRATIO*(NB_F(2,J)-NB_F(1,J))
            NM(J) = NM_F(1,J) + DTRATIO*(NM_F(2,J)-NM_F(1,J))
         ENDDO
      ENDIF
      IF (NOLIBF.eq.3) then
         IF (HBREAK.gt.0d0) then
            DO J = 1,NP
               ! This is the mixed layer depth
               MLD(J) = MLD_F(1,J) + DTRATIO*(MLD_F(2,J)-MLD_F(1,J))
            ENDDO
         ELSE
            DO J = 1,NP
               ! This is the bottom friction due to momentum dispersion
               MLD(J) = CD_F(1,J) + DTRATIO*(CD_F(2,J)-CD_F(1,J))
            ENDDO
         ENDIF
      ENDIF
      IF (abs(IDEN).eq.6) THEN
         ! This is interpolated onto element and put back onto node in
         ! another subroutine (FBPG_Disp_from_BC3D)
         DO J = 1,NP 
            DXE(J) = VIDISPDXOH_F(1,J) + DTRATIO * 
     &               (VIDISPDXOH_F(2,J) - VIDISPDXOH_F(1,J))
            DYE(J) = VIDISPDYOH_F(1,J) + DTRATIO * 
     &               (VIDISPDYOH_F(2,J) - VIDISPDYOH_F(1,J))
         ENDDO
C        DW/WJP: rotate dispersion gradients if coordinate transform
         IF ( IFSPROTS .eq. 1 ) THEN
            UVECTMP(1:NP) = DXE(1:NP); 
            VVECTMP(1:NP) = DYE(1:NP);
            CALL DRVMAP2DSPVEC( DXE, DYE,
     &            UVECTMP, VVECTMP, NP, FWD = .TRUE. ) ;
         ENDIF
      ENDIF
C
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE Update_BC3D_Info
C
C-----------------------------------------------------------------------
C     S U B R O U T I N E  P U T _ B C 3 D _ O N _ A D C I R C _ G R I D
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Put_BC3D_on_ADCIRC_Grid()

      use mesh, only : np, ne
      implicit none
      
      integer :: i, j, ii, jj, icount 
C
      !call setMessageSource("Put_3DBC_on_ADCIRC_Grid")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
  
      ! Interpolate to our mesh
      IDEN7: if (abs(iden).ne.7) then
      do i = 1,NP  
         if (BC3D_SigTS(indm(1,i),indm(2,i)) < -1d2.or.
     &       BC3D_SigTS(indm(3,i),indm(2,i)) < -1d2.or.
     &       BC3D_SigTS(indm(1,i),indm(4,i)) < -1d2.or.
     &       BC3D_SigTS(indm(3,i),indm(4,i)) < -1d2) then
             SigTS_F(2,i) = -999d0
         else
             SigTS_F(2,i) = 
     &                  BC3D_SigTS(indm(1,i),indm(2,i))*weightsm(1,i) + 
     &                  BC3D_SigTS(indm(3,i),indm(2,i))*weightsm(2,i) +
     &                  BC3D_SigTS(indm(1,i),indm(4,i))*weightsm(3,i) +
     &                  BC3D_SigTS(indm(3,i),indm(4,i))*weightsm(4,i)
         endif
      enddo
      do i = 1,NP
         if (weightsx(1,i).le.1d0) then
            ! Linear interpolation 
            VIDBCPDXOH_F(2,i) = BC3D_BX(indx(1,i),indx(2,i))*weightsx(1,i) + 
     &                  BC3D_BX(indx(3,i),indx(2,i))*weightsx(2,i) +
     &                  BC3D_BX(indx(1,i),indx(4,i))*weightsx(3,i) +
     &                  BC3D_BX(indx(3,i),indx(4,i))*weightsx(4,i)
         else
            ! The grid-averaging "conservative" (approximately) interpolation 
            VIDBCPDXOH_F(2,i) = 0.0d0
            icount = 0
            do ii = indx(1,i),indx(3,i)
               do jj = indx(2,i),indx(4,i)
                  icount = icount + 1
                  VIDBCPDXOH_F(2,i) = VIDBCPDXOH_F(2,i) + BC3D_BX(ii,jj)
               enddo
            enddo
            VIDBCPDXOH_F(2,i) = VIDBCPDXOH_F(2,i)/icount
         endif
         if (weightsy(1,i).le.1d0) then
            ! Linear interpolation 
            VIDBCPDYOH_F(2,i) = BC3D_BY(indy(1,i),indy(2,i))*weightsy(1,i) + 
     &                  BC3D_BY(indy(3,i),indy(2,i))*weightsy(2,i) +
     &                  BC3D_BY(indy(1,i),indy(4,i))*weightsy(3,i) +
     &                  BC3D_BY(indy(3,i),indy(4,i))*weightsy(4,i)
         else
            VIDBCPDYOH_F(2,i) = 0.0d0
            icount = 0
            do ii = indy(1,i),indy(3,i)
               do jj = indy(2,i),indy(4,i)
                  icount = icount + 1
                  VIDBCPDYOH_F(2,i) = VIDBCPDYOH_F(2,i) + BC3D_BY(ii,jj)
               enddo
            enddo
            VIDBCPDYOH_F(2,i) = VIDBCPDYOH_F(2,i)/icount
         endif
      enddo
      endif IDEN7
      if (LoadIT_Fric) then
         do i = 1,NP  
            NB_F(2,i) = BC3D_NB(indm(1,i),indm(2,i))*weightsm(1,i) + 
     &                  BC3D_NB(indm(3,i),indm(2,i))*weightsm(2,i) +
     &                  BC3D_NB(indm(1,i),indm(4,i))*weightsm(3,i) +
     &                  BC3D_NB(indm(3,i),indm(4,i))*weightsm(4,i)
            NM_F(2,i) = BC3D_NM(indm(1,i),indm(2,i))*weightsm(1,i) + 
     &                  BC3D_NM(indm(3,i),indm(2,i))*weightsm(2,i) +
     &                  BC3D_NM(indm(1,i),indm(4,i))*weightsm(3,i) +
     &                  BC3D_NM(indm(3,i),indm(4,i))*weightsm(4,i)
         enddo
      endif
      if (NOLIBF.eq.3) then
         if (HBREAK.gt.0d0) then
            do i = 1,NP
               if (BC3D_MLD(indm(1,i),indm(2,i)).eq.DFV.or.
     &             BC3D_MLD(indm(3,i),indm(2,i)).eq.DFV.or.
     &             BC3D_MLD(indm(1,i),indm(4,i)).eq.DFV.or.
     &             BC3D_MLD(indm(3,i),indm(4,i)).eq.DFV) then
                   MLD_F(2,i) = DFV
              else
                   MLD_F(2,i) = 
     &                  BC3D_MLD(indm(1,i),indm(2,i))*weightsm(1,i) + 
     &                  BC3D_MLD(indm(3,i),indm(2,i))*weightsm(2,i) +
     &                  BC3D_MLD(indm(1,i),indm(4,i))*weightsm(3,i) +
     &                  BC3D_MLD(indm(3,i),indm(4,i))*weightsm(4,i)
              endif
            enddo
         else
            do i = 1,NP
               if (weightss(1,i).le.1d0) then
                  CD_F(2,i) = 
     &                  BC3D_CD(inds(1,i),inds(2,i))*weightss(1,i) + 
     &                  BC3D_CD(inds(3,i),inds(2,i))*weightss(2,i) +
     &                  BC3D_CD(inds(1,i),inds(4,i))*weightss(3,i) +
     &                  BC3D_CD(inds(3,i),inds(4,i))*weightss(4,i)
               else
                  ! The grid-averaging "conservative" (approximately) interpolation 
                  CD_F(2,i) = 0.0d0;
                  icount = 0
                  do ii = inds(1,i),inds(3,i)
                     do jj = inds(2,i),inds(4,i)
                        icount = icount + 1
                        CD_F(2,i) = CD_F(2,i) + BC3D_CD(ii,jj)
                     enddo
                  enddo
                  CD_F(2,i) = CD_F(2,i)/icount
               endif
            enddo
         endif
      endif
      if (abs(iden).eq.6) then
         do i = 1,NP
            if (weightss(1,i).le.1d0) then
               ! Linear interpolation 
               VIDISPDXOH_F(2,i) = BC3D_DX(inds(1,i),inds(2,i))*weightss(1,i) + 
     &                  BC3D_DX(inds(3,i),inds(2,i))*weightss(2,i) +
     &                  BC3D_DX(inds(1,i),inds(4,i))*weightss(3,i) +
     &                  BC3D_DX(inds(3,i),inds(4,i))*weightss(4,i)
               VIDISPDYOH_F(2,i) = BC3D_DY(inds(1,i),inds(2,i))*weightss(1,i) + 
     &                  BC3D_DY(inds(3,i),inds(2,i))*weightss(2,i) +
     &                  BC3D_DY(inds(1,i),inds(4,i))*weightss(3,i) +
     &                  BC3D_DY(inds(3,i),inds(4,i))*weightss(4,i)
            else
               ! The grid-averaging "conservative" (approximately) interpolation 
               VIDISPDXOH_F(2,i) = 0.0d0; VIDISPDYOH_F(2,i) = 0.0d0
               icount = 0
               do ii = inds(1,i),inds(3,i)
                  do jj = inds(2,i),inds(4,i)
                     icount = icount + 1
                     VIDISPDXOH_F(2,i) = VIDISPDXOH_F(2,i) + BC3D_DX(ii,jj)
                     VIDISPDYOH_F(2,i) = VIDISPDYOH_F(2,i) + BC3D_DY(ii,jj)
                  enddo
               enddo
               VIDISPDXOH_F(2,i) = VIDISPDXOH_F(2,i)/icount
               VIDISPDYOH_F(2,i) = VIDISPDYOH_F(2,i)/icount
            endif
         enddo
      endif 
C
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE Put_BC3D_on_ADCIRC_Grid
C
C-----------------------------------------------------------------------
C     S U B R O U T I N E  F B P G _ D i s p _ f r o m _ B C 3 D
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FBPG_Disp_from_BC3D

      USE GLOBAL, ONLY: G, IFNLFA, IFNLCT, RAMP, NODECODE, NOFF, 
     &                  ETA2, RhoWat0, IFSFM 
      USE SPONGELAYER, ONLY: NO_BPG_IN_SPONGE 
      USE NodalAttributes, ONLY: absorblayer_sigma_mnx, 
     &                           absorblayer_sigma_mny
      USE MESH, ONLY: NE, NM, NP, DP, AREAS, NODELE, FDXE, FDYE,
     &      NEITABELE, NEITAB, NNEIGH, SFAC, TOTALAREA, SFacEle,
     &      SFMXEle, SFMYEle
      USE BOUNDARIES, ONLY: NVELME, ME2GW, NBV
 
      IMPLICIT NONE
      
      INTEGER :: IE             !element loop counter
      INTEGER :: J, K, NBDI     !node loop counter
      INTEGER :: NCELE          !element code
      INTEGER :: NEle           !local value of NetTabEle
      INTEGER :: NM1,NM2,NM3    !local node numbers used to compute gradients
      INTEGER :: NC1,NC2,NC3    !local node codes

      REAL(SZ) :: SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac
      REAL(SZ) :: AreaIE2, AreaEle
      REAL(SZ) :: FDX1, FDX2, FDX3
      REAL(SZ) :: FDY1, FDY2, FDY3
      REAL(SZ) :: FDX1O2A, FDX2O2A, FDX3O2A
      REAL(SZ) :: FDY1O2A, FDY2O2A, FDY3O2A
      REAL(SZ) :: EtaN1, EtaN2, EtaN3
      REAL(SZ) :: DEta2DX, DEta2DY
      REAL(SZ) :: DRhoDX, DRhoDY
      REAL(SZ) :: DARhoMRho0N1,DARhoMRho0N2,DARhoMRho0N3
      REAL(SZ) :: VIDBCPDXOHN1,VIDBCPDXOHN2,VIDBCPDXOHN3
      REAL(SZ) :: VIDBCPDYOHN1,VIDBCPDYOHN2,VIDBCPDYOHN3
      REAL(SZ) :: VIDBCPDXOHAvgArea,VIDBCPDYOHAvgArea

C
      !call setMessageSource("FBPG_Disp_from_3DBC")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      
      ! Initialize the nodal values
      VIDISPDXOH = 0.0d0; VIDISPDYOH = 0.0d0;
      VIDBCPDXOH = 0.0d0; VIDBCPDYOH = 0.0d0;
#ifndef NOFSBPG
      ! Loop over each element and get the free surface BCP gradient
      ! and add on the already known elemental gradients
      DO IE = 1,NE
         NM1 = NM(IE,1)
         NM2 = NM(IE,2)
         NM3 = NM(IE,3)
         NC1 = NODECODE(NM1)
         NC2 = NODECODE(NM2)
         NC3 = NODECODE(NM3)
         NCEle = NC1*NC2*NC3*NOFF(IE)
         EtaN1 = IFNLFA*Eta2(NM1)
         EtaN2 = IFNLFA*Eta2(NM2)
         EtaN3 = IFNLFA*Eta2(NM3)
         SFacAvg = SFacEle(IE)
         ! CPB 5/25/2022: correct rotated coordinates
         SFmxAvg = SFMXEle(IE) ;
         SFmyAvg = SFMYEle(IE) ;
         sfdxfac = (1 - IFSFM)*SFacAvg + IFSFM*SFmxAvg ; 
         sfdyfac = (1 - IFSFM)*1.0_SZ  + IFSFM*SFmyAvg ; 

         FDX1 = FDXE(1,IE)*sfdxfac ; !c FDX1=(Y(NM2)-Y(NM3))*SFacAvg !b1
         FDX2 = FDXE(2,IE)*sfdxfac ; !c FDX2=(Y(NM3)-Y(NM1))*SFacAvg !b2
         FDX3 = FDXE(3,IE)*sfdxfac ; !c FDX3=(Y(NM1)-Y(NM2))*SFacAvg !b3
         FDY1 = FDYE(1,IE)*sfdyfac ; !c  FDY1=X(NM3)-X(NM2) !a1
         FDY2 = FDYE(2,IE)*sfdyfac ; !c  FDY2=X(NM1)-X(NM3) !a2
         FDY3 = FDYE(3,IE)*sfdyfac ; !c  FDY3=X(NM2)-X(NM1) !a3

         AreaIE2 = Areas(IE)
         AreaEle = 0.5d0*NCEle*AreaIE2
        
         FDX1O2A = FDX1/AreaIE2 !dphi1/dx
         FDY1O2A = FDY1/AreaIE2 !dphi1/dy
         FDX2O2A = FDX2/AreaIE2 !dphi2/dx
         FDY2O2A = FDY2/AreaIE2 !dphi2/dy
         FDX3O2A = FDX3/AreaIE2 !dphi3/dx
         FDY3O2A = FDY3/AreaIE2 !dphi3/dy 
         !Use free surface density only if all non-zero
         IF (SigTS(NM1) < -1d2.or.SigTS(NM2) < -1d2.or.
     &       SigTS(NM3) < -1d2) THEN
            DARhoMRho0N1 = 0d0
            DARhoMRho0N2 = 0d0
            DARhoMRho0N3 = 0d0
         ELSE
            DARhoMRho0N1 = SigTS(NM1)
            DARhoMRho0N2 = SigTS(NM2) 
            DARhoMRho0N3 = SigTS(NM3)
         ENDIF

         ! Add on free surface baroclinic pressure on the element
         VIDBCPDXOHN1 = G/RhoWat0*(DARhoMRho0N1*EtaN1*FDX1O2A 
     &                + DARhoMRho0N2*EtaN2*FDX2O2A
     &                + DARhoMRho0N3*EtaN3*FDX3O2A)
 
         VIDBCPDYOHN1 = G/RhoWat0*(DARhoMRho0N1*EtaN1*FDY1O2A 
     &                + DARhoMRho0N2*EtaN2*FDY2O2A
     &                + DARhoMRho0N3*EtaN3*FDY3O2A)

         ! Get the element-area weighted contributions
         VIDBCPDXOHAvgArea = AreaEle*VIDBCPDXOHN1
         VIDBCPDYOHAvgArea = AreaEle*VIDBCPDYOHN1

         ! Add element-area weighted contributions to nodal contrib.
         VIDBCPDXOH(NM1) = VIDBCPDXOH(NM1) + VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM2) = VIDBCPDXOH(NM2) + VIDBCPDXOHAvgArea
         VIDBCPDXOH(NM3) = VIDBCPDXOH(NM3) + VIDBCPDXOHAvgArea
         VIDBCPDYOH(NM1) = VIDBCPDYOH(NM1) + VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM2) = VIDBCPDYOH(NM2) + VIDBCPDYOHAvgArea
         VIDBCPDYOH(NM3) = VIDBCPDYOH(NM3) + VIDBCPDYOHAvgArea
      ENDDO
#endif
      ! Loop over each node, divide by the total area and scale by
      ! ramp, gravity and the reference density. Also add the
      ! internal BPG calculated previously outside of this routine
      DO J = 1,NP
         IF (TotalArea(J).NE.0.d0) THEN
            ! Divide BPG by element area
            VIDBCPDXOH(J) = Ramp*(VIDBCPDXOH(J)/TotalArea(J) + BXE(J))
            VIDBCPDYOH(J) = Ramp*(VIDBCPDYOH(J)/TotalArea(J) + BYE(J))
            VIDISPDXOH(J) = Ramp*DXE(J)
            VIDISPDYOH(J) = Ramp*DYE(J)
         ENDIF
         ! If need to omit from sponge..
         IF ( NO_BPG_IN_SPONGE ) THEN
            IF ( absorblayer_sigma_mnx(J,1) > 1.0e-9 ) THEN
               VIDBCPDXOH(J) = 0.0d0; VIDISPDXOH(J) = 0.0d0; 
            ENDIF
            IF ( absorblayer_sigma_mny(J,1) > 1.0e-9 ) THEN
               VIDBCPDYOH(J) = 0.0d0; VIDISPDYOH(J) = 0.0d0; 
            ENDIF
         ENDIF
      ENDDO
      ! Zeroing out on boundary nodes
      !DO J = 1,NVELME
      !   K = ME2GW(J)
      !   !IF (LBCODEI(K).GT.) cycle
      !   NBDI = NBV(K)
      !   VIDBCPDXOH(NBDI) = 0; VIDBCPDYOH(NBDI) = 0;
      !   VIDISPDXOH(NBDI) = 0; VIDISPDYOH(NBDI) = 0;
      !ENDDO
C
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE FBPG_Disp_from_BC3D
C-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E   B C 3 D  T E R M I N A T E
!-----------------------------------------------------------------------
      SUBROUTINE BC3DTerminate(NO_MPI_FINALIZE)
      USE SIZES, ONLY : myproc
#ifdef CMPI 
      USE MESSENGER
#endif
      IMPLICIT NONE 
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
#ifdef ALL_TRACE
      REAL, ALLOCATABLE :: dummy(:)
#endif
C
      !call setMessageSource("BC3DTerminate")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef ALL_TRACE
         ! intentionally create a segmentation fault so that we can get
         ! a stack trace to determine the line number of the netcdf call
         ! that went bad ... this assumes that the code was compiled with
         ! debugging symbols, bounds checking, and stack trace turned on.
      write(6,*) 'myproc = ',myProc      
      dummy(1) = 99.9d0
#endif

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
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      !call unsetMessageSource()
!-----------------------------------------------------------------------
      END SUBROUTINE BC3DTerminate
!-----------------------------------------------------------------------
      END MODULE Couple2BC3D
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
