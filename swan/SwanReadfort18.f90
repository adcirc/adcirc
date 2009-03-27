subroutine SwanReadfort18
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.95: Marcel Zijlema
!
!   Updates
!
!   40.95, June 2008: New subroutine
!
!   Purpose
!
!   Reads fort.18 to obtain the following ADCIRC variables
!
!   MNE, MNP, NSTAE, NSTAV, NSTAM and NSTAC
!
!   These variables are needed by the message-passing routines in module MESSENGER
!
!   Next, we also need the number of elements and nodes in global domain for SWAN computation
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwanGriddata, only: nvertsg, ncellsg
    use SIZES
    use GLOBAL, only: NSTAE, NSTAV, NSTAM, NSTAC
!
    implicit none
!
!   Local variables
!
    integer, save           :: ient = 0 ! number of entries in this subroutine
    integer                 :: idum1    ! dummy integer 1
    integer                 :: idum2    ! dummy integer 2
    integer                 :: idum3    ! dummy integer 3
    integer                 :: idum4    ! dummy integer 4
    integer                 :: iostat   ! I/O status in call FOR
    integer                 :: j        ! loop counter
    character(80)           :: msgfil   ! name of message-passing file including path
    integer                 :: ndsd     ! unit reference number of file
    logical                 :: stpnow   ! indicate whether program must be terminated or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanReadfort18')
    !
    ! open file fort.18
    !
    ndsd   = 0
    iostat = 0
    msgfil = trim(INPUTDIR)//DIRCH2//'fort.18'
    call for (ndsd, msgfil, 'OF', iostat)
    if (stpnow()) goto 900
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, idum3
    !
    read(ndsd,100, end=950, err=910) ncellsg, idum2, MNE   ! number of elements
    do j = 1, MNE
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,100, end=950, err=910) nvertsg, idum2, MNP   ! number of nodes
    do j = 1, MNP
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,'(8x,i8)', end=950, err=910) idum1
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, idum3
    do j = 1, idum3
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAE   ! number of elevation stations
    do j = 1, NSTAE
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAV   ! number of velocity stations
    do j = 1, NSTAV
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAM   ! number of meteorlogical stations
    do j = 1, NSTAM
       read(ndsd,'(i8)', end=950, err=910) idum4
    enddo
    !
    read(ndsd,100, end=950, err=910) idum1, idum2, NSTAC   ! number of concentration stations
    !
    ! close file fort.18
    !
    close(ndsd)
    !
 900 return
    !
 910 call msgerr (4, 'error reading data from grid file fort.18' )
    goto 900
 950 call msgerr (4, 'unexpected end of file in grid file fort.18' )
    goto 900
    !
 100 format((8x,3i8))
    !
end subroutine SwanReadfort18
