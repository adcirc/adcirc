      program build12
c-------------------------------------------------------------------------
c  This utility constructs fort.12 and fort.21 from fort.13
c  vjp 8/1/2006
c-------------------------------------------------------------------------
      implicit none
      integer ne, np, ne2, np2, idumy, i, j, ntot, nattr, ncount, inode
      real(8) xdumy, ydumy, value
      real(8) sea_surface_hag_def, surface_substate_def
      real(8) surface_canopy_def, prim_weight_def
      real(8) mannings_def, surface_derl_def(12)
      character(24) agrid, agrid2, agrid3
      character(80) attr, inline
      character(42) :: prim_weight = "primitive_weighting_in_continuity_equation"
      character(23) :: mannings    = "mannings_n_at_sea_floor"
      character(25) :: surface_substate = "surface_submergence_state"
      character(26) :: surface_canopy = "surface_canopy_coefficient" 
      character(46) :: surface_derl = "surface_directional_effective_roughness_length"
      character(30) :: sea_surface_hag = "sea_surface_height_above_geoid"
      real(8),allocatable :: x(:), y(:), dp(:)
      real(8),allocatable :: startdry(:), vcanopy(:), z0land(:,:)
      real(8),allocatable :: fric(:)

      open(AMPI_LUN(14),file='fort.14')
      call AMPI_LUN_REGISTER(AMPI_LUN(14),file='fort.14')
      read(AMPI_LUN(14),'(A)') agrid
      read(AMPI_LUN(14),*) ne,np

      allocate( x(np), y(np), dp(np) )

      do i=1,np
         read(AMPI_LUN(14),*) idumy,x(i),y(i),dp(i)
      enddo
      close(AMPI_LUN(14))

      allocate( startdry(np), vcanopy(np), z0land(np,12), fric(np) )

c-----------------------------------------------------------------------
c   read variable descriptors
c-----------------------------------------------------------------------

      open(AMPI_LUN(13),file='fort.13')
      call AMPI_LUN_REGISTER(AMPI_LUN(13),file='fort.13')
      read(AMPI_LUN(13),'(A)') agrid2
      print *, agrid2
      read(AMPI_LUN(13),*)  np
      print *, "np = ",np
      read(AMPI_LUN(13),*)  nattr
      print *, "nattr = ", nattr
      do i=1, nattr
         read(AMPI_LUN(13),'(A)') attr
         print *, "found attribute = ",trim(attr)

         read(AMPI_LUN(13),'(A)') inline
         print *, "units = ", trim(inline)

         read(AMPI_LUN(13),*) idumy   
         print *, "values per node = ", idumy      

         if (trim(attr).eq. prim_weight) then
           read(AMPI_LUN(13),*) prim_weight_def 
           print *, "prim_weight default = ", prim_weight_def

         elseif (trim(attr).eq. mannings) then
           read(AMPI_LUN(13),*) mannings_def 
           print *, "mannings default = ", mannings_def

         elseif (trim(attr).eq. surface_substate) then
           read(AMPI_LUN(13),*) surface_substate_def 
           print *, "surface_substate default = ", surface_substate_def

         elseif (trim(attr).eq. surface_canopy) then
           read(AMPI_LUN(13),*) surface_canopy_def 
           print *, "surface_canopy default = ", surface_canopy_def

         elseif (trim(attr).eq. surface_derl) then
           read(AMPI_LUN(13),*) surface_derl_def(1:12) 
           print *, "surface_derl default = "
           print *, surface_derl(1:12)

         elseif (trim(attr).eq. sea_surface_hag) then
           read(AMPI_LUN(13),*) sea_surface_hag_def 
           print *, "surface_surface_hag default = ", sea_surface_hag_def
         endif
      enddo

c-----------------------------------------------------------------------
c  read variable values
c-----------------------------------------------------------------------
        
      do i=1, nattr
         read(AMPI_LUN(13),'(A)') attr
         print *, "processing: ",trim(attr)
         if (trim(attr).eq. prim_weight) then
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(AMPI_LUN(13),'(A)') inline
           enddo
         elseif (trim(attr).eq. mannings) then
           do j=1, np
              fric(j) = mannings_def
           enddo 
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(AMPI_LUN(13),*) inode, value
              fric(inode) = value
           enddo
         elseif (trim(attr).eq. surface_substate) then
           do j=1, np
              startdry(j) = surface_substate_def 
           enddo 
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(AMPI_LUN(13),*) inode, value
              startdry(inode) = 1.0  
           enddo
         elseif (trim(attr).eq. surface_canopy) then
           do j=1, np
              vcanopy(j) = surface_canopy_def   
           enddo 
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(AMPI_LUN(13),*) inode, value
              vcanopy(inode) = value
           enddo
         elseif (trim(attr).eq. surface_derl) then
           do j=1, np
              z0land(j,1:12) = surface_derl_def(1:12)
           enddo 
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           do j=1, ncount
              read(AMPI_LUN(13),*) inode, z0land(inode,1:12)
           enddo
         elseif (trim(attr).eq. sea_surface_hag) then
           read(AMPI_LUN(13),*) ncount           
           print *, "ncount = ", ncount
           if (ncount .gt. 0) then
             do j=1, ncount
              read(AMPI_LUN(13),'(A)') inline
             enddo
           endif
         endif
      enddo
      close(AMPI_LUN(13))

c Write fort.12 file

      print *, "writing fort.12 file"
      open(AMPI_LUN(12),file='fort.12')
      call AMPI_LUN_REGISTER(AMPI_LUN(12),file='fort.12')
      write(AMPI_LUN(12),'(A)') agrid2
      write(AMPI_LUN(12),*) ne,np
      do i=1,np
         write(AMPI_LUN(12),'(I8,16E17.9)') i,x(i),y(i),startdry(i),
     &    (z0land(i,j),j=1, 12),vcanopy(i)
      enddo  
      close(AMPI_LUN(12))

c Write fort.21 file

      print *, "writing fort.21 file"
      OPEN(AMPI_LUN(21),FILE='fort.21')
      call AMPI_LUN_REGISTER(AMPI_LUN(21),FILE='fort.21')
      write(AMPI_LUN(21),'(A)') agrid
      do i=1, np
         write(AMPI_LUN(21),'(I8,E17.9)') i, fric(i)
      enddo
      close(AMPI_LUN(21))
C
      stop
      end


