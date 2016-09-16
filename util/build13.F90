    program build13
    implicit none
    logical :: has_dry, has_fric
    character(80) agrid, agrid2, agrid3
    integer :: ne, np, ne2, np2, idumy, i, j, ntot
    real(8) xdumy, ydumy
    integer,allocatable :: iz0(:)
    real(8),allocatable :: dp(:), startdry(:), vcanopy(:), z0land(:,:)
    real(8),allocatable :: tau0var(:), fric(:)

!-------------------------------------------------------------------------
!  This utility constructs a fort.13 file for adcirc v46
!  vjp 9/1/2006
!-------------------------------------------------------------------------

    open(14,file='fort.14')
    read(14,'(A)') agrid
    read(14,*) ne,np
    allocate( dp(np) )
    do i=1,np
        read(14,*) idumy,xdumy,ydumy,dp(i)
    enddo
    close(14)

    inquire(file='fort.12', exist = has_dry)
    if (has_dry)  then
        open(12,file='fort.12')
        read(12,'(A)') agrid2
        read(12,*) ne2,np2
        allocate( startdry(np), vcanopy(np), z0land(np,12), iz0(np) )
        do i=1,np
            read(12,*) idumy,xdumy,ydumy,startdry(i), &
            (z0land(i,j),j=1, 12),vcanopy(i)
        enddo
        close(12)

        allocate( tau0var(np) )
        do i=1,np
            if(dp(i) <= 10.0d0) tau0var(i)=0.02d0
            if(dp(i) > 10.0d0) tau0var(i)=0.005d0
            if(startdry(i) < -60000.0d0) tau0var(i)=0.03d0
        enddo
    endif

! Read local fort.21 file

    inquire(file='fort.21', exist = has_fric)
    if (has_fric) then
        OPEN(21,FILE='fort.21')
        allocate( fric(np) )
        read(21,'(A)') agrid3
        do i=1, np
            read(21,*) idumy, fric(i)
        enddo
        close(21)
    endif

!-----------------------------------------------------------------------
!   write variable descriptors
!-----------------------------------------------------------------------

    open(13,file='fort.13')
    write(13,'(A)') trim(agrid2)
    write(13,'(I7)')  np
    if (has_dry .AND. has_fric) then
        write(13,'(I2)')  6
    elseif (has_dry .AND. .NOT. has_fric) then
        write(13,'(I2)')  5
    elseif (has_fric .AND. .NOT. has_dry) then
        write(13,'(I2)')  1
    else
        print *, "cannot find either a fort.12 or a for.21 file"
        stop
    endif


    if (has_dry) then
        write(13,'(A)')  "sea_surface_height_above_geoid"
        write(13,'(A2)')  " m"
        write(13,'(I2)')  1
        write(13,'(A7)')  "0.36576"

        write(13,'(A)')  "surface_submergence_state"
        write(13,'(I2)')  1
        write(13,'(I2)')  1
        write(13,'(A3)')  "0.0"

        write(13,'(A)') "surface_canopy_coefficient"
        write(13,'(I2)')  1
        write(13,'(I2)')  1
        write(13,'(A4)')  " 1.0"

        write(13,'(A)')  "primitive_weighting_in_continuity_equation"
        write(13,'(I2)')  1
        write(13,'(I2)')  1
        write(13,'(A4)')  "0.0"

        write(13,'(A)') "surface_directional_effective_roughness_length"
        write(13,'(A2)')  " m"
        write(13,'(I3)')  12
        write(13,'(A)') &
        "0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0  0.0"
    endif

    if (has_fric) then
        write(13,'(A)')  "mannings_n_at_sea_floor"
        write(13,'(A2)')  " m"
        write(13,'(I2)')  1
        write(13,'(A4)')  "0.0"
    endif

!-----------------------------------------------------------------------
!  write variable values
!-----------------------------------------------------------------------
            
    if (has_dry) then

    ! VARABLE 1
        write(13,'(A)')  "sea_surface_height_above_geoid"
        write(13,'(I1)')  0

    ! VARABLE 2
        write(13,'(A)')  "surface_submergence_state"
        ntot = 0
        do i=1, np
            if(startdry(i) == -88888.0d0) ntot = ntot+1
        enddo
        write(13,'(I7)') ntot
        do i=1, np
            if(startdry(i) == -88888.0d0) write(13,'(I7,2X,F3.1)')  i,1.0
        enddo

    ! VARABLE 3
        write(13,'(A)') "surface_canopy_coefficient"
        ntot = 0
        do i=1, np
            if(vcanopy(i) /= 1.0d0) ntot = ntot+1
        enddo
        write(13,'(I8)') ntot
        do i=1, np
            if(vcanopy(i) /= 1.0d0) write(13,'(I7,F4.1)') i,0.0
        enddo

    ! VARABLE 4
        write(13,'(A)')  "primitive_weighting_in_continuity_equation"
        write(13,'(I8)') np
        do i=1,np
            write(13,'(I8,F8.3)') i,tau0var(i)
        enddo

    ! VARABLE 5
        write(13,'(A)') "surface_directional_effective_roughness_length"
        ntot = 0
        do i=1, np
            iz0(i) = 0
            do j=1, 12
                if(z0land(i,j) /= 0.0d0) then
                    iz0(i) = 1
                    ntot = ntot+1
                    exit
                endif
            enddo
        enddo
        write(13,'(I8)') ntot
        do i=1, np
            if (iz0(i) == 1) then
                write(13,'(I7,12E16.8)') i,(z0land(i,j),j=1, 12)
            endif
        enddo

    endif

    if (has_fric) then

    ! VARABLE 6
        write(13,'(A)') "mannings_n_at_sea_floor"
        write(13,'(I8)') np
        do i=1, np
            write(13,'(I7,E16.8)') i, fric(i)
        enddo

    endif

    close(13)
    stop
    END PROGRAM
