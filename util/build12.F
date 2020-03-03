      program build12
c-------------------------------------------------------------------------
c  This utility constructs fort.12 and fort.21 from fort.13
c  vjp 8/1/2006
c-------------------------------------------------------------------------
      implicit none
      integer ne, np, ne2, np2, idumy, i, j, ntot, nattr, ncount, inode
!$omp threadprivate(ne,np,ne2,np2,idumy,i,j,ntot,nattr,ncount,inode)
      real(8) xdumy, ydumy, value
!$omp threadprivate(xdumy,ydumy,value)
      real(8) sea_surface_hag_def, surface_substate_def
!$omp threadprivate(sea_surface_hag_def,surface_substate_def)
      real(8) surface_canopy_def, prim_weight_def
!$omp threadprivate(surface_canopy_def,prim_weight_def)
      real(8) mannings_def, surface_derl_def(12)
!$omp threadprivate(mannings_def,surface_derl_def)
      character(24) agrid, agrid2, agrid3
!$omp threadprivate(agrid,agrid2,agrid3)
      character(80) attr, inline
!$omp threadprivate(attr,inline)
      character(42) :: prim_weight = "primitive_weighting_in_continuity_equation"
!$omp threadprivate(prim_weight)
      character(23) :: mannings    = "mannings_n_at_sea_floor"
!$omp threadprivate(mannings)
      character(25) :: surface_substate = "surface_submergence_state"
!$omp threadprivate(surface_substate)
      character(26) :: surface_canopy = "surface_canopy_coefficient" 
!$omp threadprivate(surface_canopy)
      character(46) :: surface_derl = "surface_directional_effective_roughness_length"
!$omp threadprivate(surface_derl)
      character(30) :: sea_surface_hag = "sea_surface_height_above_geoid"
!$omp threadprivate(sea_surface_hag)
      real(8),allocatable :: x(:), y(:), dp(:)
!$omp threadprivate(x,y,dp)
      real(8),allocatable :: startdry(:), vcanopy(:), z0land(:,:)
!$omp threadprivate(startdry,vcanopy,z0land)
      real(8),allocatable :: fric(:)
!$omp threadprivate(fric)

      open(14+CK_LUN,file='fort.14')
      read(14+CK_LUN,'(A)') agrid
      read(14+CK_LUN,*) ne,np

      allocate( x(np), y(np), dp(np) )

      do i=1,np
         read(14+CK_LUN,*) idumy,x(i),y(i),dp(i)
      enddo
      close(14+CK_LUN)

      allocate( startdry(np), vcanopy(np), z0land(np,12), fric(np) )

c-----------------------------------------------------------------------
c   read variable descriptors
c-----------------------------------------------------------------------

      open(13+CK_LUN,file='fort.13')
      read(13+CK_LUN,'(A)') agrid2
      print *, agrid2
      read(13+CK_LUN,*)  np
      print *, "np = ",np
      read(13+CK_LUN,*)  nattr
      print *, "nattr = ", nattr
      do i=1, nattr
         read(13+CK_LUN,'(A)') attr
         print *, "found attribute = ",trim(attr)

         read(13+CK_LUN,'(A)') inline
         print *, "units = ", trim(inline)

         read(13+CK_LUN,*) idumy   
         print *, "values per node = ", idumy      

         if (trim(attr).eq. prim_weight) then
           read(13+CK_LUN,*) prim_weight_def 
           print *, "prim_weight default = ", prim_weight_def

         elseif (trim(attr).eq. mannings) then
           read(13+CK_LUN,*) mannings_def 
           print *, "mannings default = ", mannings_def

         elseif (trim(attr).eq. surface_substate) then
           read(13+CK_LUN,*) surface_substate_def 
           print *, "surface_substate default = ", surface_substate_def

         elseif (trim(attr).eq. surface_canopy) then
           read(13+CK_LUN,*) surface_canopy_def 
           print *, "surface_canopy default = ", surface_canopy_def

         elseif (trim(attr).eq. surface_derl) then
           read(13+CK_LUN,*) surface_derl_def(1:12) 
           print *, "surface_derl default = "
           print *, surface_derl(1:12)

         elseif (trim(attr).eq. sea_surface_hag) then
           read(13+CK_LUN,*) sea_surface_hag_def 
           print *, "surface_surface_hag default = ", sea_surface_hag_def
         endif
      enddo

c-----------------------------------------------------------------------
c  read variable values
c-----------------------------------------------------------------------
        
      do i=1, nattr
         read(13+CK_LUN,'(A)') attr
         print *, "processing: ",trim(attr)
         if (trim(attr).eq. prim_weight) then
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(13+CK_LUN,'(A)') inline
           enddo
         elseif (trim(attr).eq. mannings) then
           do j=1, np
              fric(j) = mannings_def
           enddo 
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(13+CK_LUN,*) inode, value
              fric(inode) = value
           enddo
         elseif (trim(attr).eq. surface_substate) then
           do j=1, np
              startdry(j) = surface_substate_def 
           enddo 
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(13+CK_LUN,*) inode, value
              startdry(inode) = 1.0  
           enddo
         elseif (trim(attr).eq. surface_canopy) then
           do j=1, np
              vcanopy(j) = surface_canopy_def   
           enddo 
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(13+CK_LUN,*) inode, value
              vcanopy(inode) = value
           enddo
         elseif (trim(attr).eq. surface_derl) then
           do j=1, np
              z0land(j,1:12) = surface_derl_def(1:12)
           enddo 
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(13+CK_LUN,*) inode, z0land(inode,1:12)
           enddo
         elseif (trim(attr).eq. sea_surface_hag) then
           read(13+CK_LUN,*) ncount           
           print *, "ncount = ", ncount
           if (ncount .gt. 0) then
             do j=1, ncount
              read(13+CK_LUN,'(A)') inline
             enddo
           endif
         endif
      enddo
      close(13+CK_LUN)

c Write fort.12 file

      print *, "writing fort.12 file"
      open(12+CK_LUN,file='fort.12')
      write(12+CK_LUN,'(A)') agrid2
      write(12+CK_LUN,*) ne,np
      do i=1,np
         write(12+CK_LUN,'(I8,16E17.9)') i,x(i),y(i),startdry(i),
     &     (z0land(i,j),j=1, 12),vcanopy(i)
      enddo  
      close(12+CK_LUN)

c Write fort.21 file

      print *, "writing fort.21 file"
      OPEN(21+CK_LUN,FILE='fort.21')
      write(21+CK_LUN,'(A)') agrid
      do i=1, np
         write(21+CK_LUN,'(I8,E17.9)') i, fric(i)
      enddo
      close(21+CK_LUN)
C
      stop
      end


