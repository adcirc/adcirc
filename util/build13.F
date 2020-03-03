      program build13
      implicit none
      logical has_dry, has_fric
!$omp threadprivate(has_dry,has_fric)
      character(80) agrid, agrid2, agrid3
!$omp threadprivate(agrid,agrid2,agrid3)
      integer ne, np, ne2, np2, idumy, i, j, ntot
!$omp threadprivate(ne,np,ne2,np2,idumy,i,j,ntot)
      real(8) xdumy, ydumy
!$omp threadprivate(xdumy,ydumy)
      integer,allocatable :: iz0(:)
!$omp threadprivate(iz0)
      real(8),allocatable :: dp(:), startdry(:), vcanopy(:), z0land(:,:)
!$omp threadprivate(dp,startdry,vcanopy,z0land)
      real(8),allocatable :: tau0var(:), fric(:)
!$omp threadprivate(tau0var,fric)

c-------------------------------------------------------------------------
c  This utility constructs a fort.13 file for adcirc v46
c  vjp 9/1/2006
c-------------------------------------------------------------------------

      open(14+CK_LUN,file='fort.14')
      read(14+CK_LUN,'(A)') agrid
      read(14+CK_LUN,*) ne,np
      allocate( dp(np) )
      do i=1,np
         read(14+CK_LUN,*) idumy,xdumy,ydumy,dp(i)
      enddo
      close(14+CK_LUN)

      inquire(file='fort.12', exist = has_dry)
      if (has_dry)  then
        open(12+CK_LUN,file='fort.12')
        read(12+CK_LUN,'(A)') agrid2
        read(12+CK_LUN,*) ne2,np2
        allocate( startdry(np), vcanopy(np), z0land(np,12), iz0(np) )
        do i=1,np
           read(12+CK_LUN,*) idumy,xdumy,ydumy,startdry(i),
     &       (z0land(i,j),j=1, 12),vcanopy(i)
        enddo  
        close(12+CK_LUN)

        allocate( tau0var(np) )
        do i=1,np
           if(dp(i).le.10.0d0) tau0var(i)=0.02d0
           if(dp(i).gt.10.0d0) tau0var(i)=0.005d0
           if(startdry(i).lt.-60000.0d0) tau0var(i)=0.03d0
        enddo
      endif

c Read local fort.21 file

      inquire(file='fort.21', exist = has_fric)
      if (has_fric) then
        OPEN(21+CK_LUN,FILE='fort.21')
        allocate( fric(np) )
        read(21+CK_LUN,'(A)') agrid3
        do i=1, np
           read(21+CK_LUN,*) idumy, fric(i)
        enddo
        close(21+CK_LUN)
      endif

c-----------------------------------------------------------------------
c   write variable descriptors
c-----------------------------------------------------------------------

      open(13+CK_LUN,file='fort.13')
      write(13+CK_LUN,'(A)') trim(agrid2)
      write(13+CK_LUN,'(I7)')  np
      if (has_dry .and. has_fric) then
        write(13+CK_LUN,'(I2)')  6
      elseif (has_dry .and. .not.has_fric) then
        write(13+CK_LUN,'(I2)')  5
      elseif (has_fric .and. .not.has_dry) then
        write(13+CK_LUN,'(I2)')  1
      else
         print *, "cannot find either a fort.12 or a for.21 file"
         stop 
      endif


      if (has_dry) then
        write(13+CK_LUN,'(A)')  "sea_surface_height_above_geoid"
        write(13+CK_LUN,'(A2)')  " m"
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(A7)')  "0.36576"

        write(13+CK_LUN,'(A)')  "surface_submergence_state"
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(A3)')  "0.0"

        write(13+CK_LUN,'(A)') "surface_canopy_coefficient"
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(A4)')  " 1.0"

        write(13+CK_LUN,'(A)')  "primitive_weighting_in_continuity_equation"
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(A4)')  "0.0"

        write(13+CK_LUN,'(A)') "surface_directional_effective_roughness_length"
        write(13+CK_LUN,'(A2)')  " m"
        write(13+CK_LUN,'(I3)')  12
        write(13+CK_LUN,'(A)') 
     &  "0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0  0.0"
      endif

      if (has_fric) then
        write(13+CK_LUN,'(A)')  "mannings_n_at_sea_floor"
        write(13+CK_LUN,'(A2)')  " m"
        write(13+CK_LUN,'(I2)')  1
        write(13+CK_LUN,'(A4)')  "0.0"
      endif

c-----------------------------------------------------------------------
c  write variable values
c-----------------------------------------------------------------------
        
      if (has_dry) then

! VARABLE 1
      write(13+CK_LUN,'(A)')  "sea_surface_height_above_geoid"
      write(13+CK_LUN,'(I1)')  0

! VARABLE 2
      write(13+CK_LUN,'(A)')  "surface_submergence_state"
      ntot = 0
      do i=1, np
         if(startdry(i).eq.-88888.0d0) ntot = ntot+1
      enddo
      write(13+CK_LUN,'(I7)') ntot
      do i=1, np
         if(startdry(i).eq.-88888.0d0) write(13+CK_LUN,'(I7,2X,F3.1)')  i,1.0
      enddo

! VARABLE 3
      write(13+CK_LUN,'(A)') "surface_canopy_coefficient"
      ntot = 0
      do i=1, np
         if(vcanopy(i).ne.1.0d0) ntot = ntot+1
      enddo
      write(13+CK_LUN,'(I8)') ntot
      do i=1, np
         if(vcanopy(i).ne.1.0d0) write(13+CK_LUN,'(I7,F4.1)') i,0.0
      enddo

! VARABLE 4
      write(13+CK_LUN,'(A)')  "primitive_weighting_in_continuity_equation"
      write(13+CK_LUN,'(I8)') np
      do i=1,np
         write(13+CK_LUN,'(I8,F8.3)') i,tau0var(i)
      enddo

! VARABLE 5
      write(13+CK_LUN,'(A)') "surface_directional_effective_roughness_length"
      ntot = 0
      do i=1, np
         iz0(i) = 0
         do j=1, 12
           if(z0land(i,j).ne.0.0d0) then
             iz0(i) = 1
             ntot = ntot+1
             exit
           endif
         enddo
      enddo
      write(13+CK_LUN,'(I8)') ntot
      do i=1, np
         if (iz0(i).eq.1) then
           write(13+CK_LUN,'(I7,12E16.8)') i,(z0land(i,j),j=1, 12)
         endif
      enddo

      endif

      if (has_fric) then

! VARABLE 6
      write(13+CK_LUN,'(A)') "mannings_n_at_sea_floor"
      write(13+CK_LUN,'(I8)') np
      do i=1, np
         write(13+CK_LUN,'(I7,E16.8)') i, fric(i)
      enddo

      endif

      close(13+CK_LUN)
      stop
      end
