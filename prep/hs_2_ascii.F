!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
C***************************************************************************
C  PROGRAM TO READ A BINARY HOTSTART FILE AND CREATE A SIMILAR ASCII FILE  *
C                                                                          *
C                      rl 10/11/01                                         *
C jgf updated for v45.06  09/07/2005                                       *
C***************************************************************************


      PROGRAM hs_2_ascii

      IMPLICIT NONE

      INTEGER I,J,IHOTSTP,MNP,MNE
      INTEGER IMHSF,ITHSF
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW
      CHARACTER FNAME*60
      INTEGER,ALLOCATABLE  :: NODECODE(:), NOFF(:)  !tcm 20110510 v50.05 added NOFF
      REAL(8),ALLOCATABLE :: ETA1(:),ETA2(:),UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF

C--Determine which hotstart files to process

      WRITE(*,*)'Enter binary hotstart file name'
      READ(*,*) FNAME
      OPEN (99,FILE=FNAME,ACCESS='DIRECT',RECL=8)

C--Read stuff from each local file  

      WRITE(*,*) 'Enter the number of nodes in the file'
      READ(*,*) MNP
      WRITE(*,*) 'Enter the number of elements in the file'
      READ(*,*) MNE

C--Allocate local work arrays

      ALLOCATE ( ETA1(MNP),ETA2(MNP),UU2(MNP),
     &     VV2(MNP),NODECODE(MNP),CH1(MNP),NOFF(MNE) )           

      IHOTSTP=1
      READ(99,REC=IHOTSTP) IMHSF
      IHOTSTP=2
      READ(99,REC=IHOTSTP) TIMEHSF
      IHOTSTP=3
      READ(99,REC=IHOTSTP) ITHSF

      DO I=1,MNP
        READ(99,REC=IHOTSTP+1) ETA1(I)
        READ(99,REC=IHOTSTP+2) ETA2(I)
        READ(99,REC=IHOTSTP+3) UU2(I)
        READ(99,REC=IHOTSTP+4) VV2(I)
        IHOTSTP=IHOTSTP+4
c       IF(IM.EQ.10) THEN
c         READ(LOCHSF,REC=IHOTSTP+1) CH1(I)
c         IHOTSTP=IHOTSTP+1
c         ENDIF
        READ(99,REC=IHOTSTP+1) NODECODE(I)
        IHOTSTP=IHOTSTP+1
        ENDDO
      DO I=1,MNE
         READ(99,REC=IHOTSTP+1) NOFF(I)                     
         IHOTSTP=IHOTSTP+1                                    
      END DO


C--Read in more common info from higest processor hotstart file

      READ(99,REC=IHOTSTP+1 ) IESTP
      READ(99,REC=IHOTSTP+2 ) NSCOUE
      READ(99,REC=IHOTSTP+3 ) IVSTP
      READ(99,REC=IHOTSTP+4 ) NSCOUV
      READ(99,REC=IHOTSTP+5 ) ICSTP
      READ(99,REC=IHOTSTP+6 ) NSCOUC
      READ(99,REC=IHOTSTP+7 ) IPSTP
      READ(99,REC=IHOTSTP+8 ) IWSTP
      READ(99,REC=IHOTSTP+9 ) NSCOUM
      READ(99,REC=IHOTSTP+10) IGEP
      READ(99,REC=IHOTSTP+11) NSCOUGE
      READ(99,REC=IHOTSTP+12) IGVP
      READ(99,REC=IHOTSTP+13) NSCOUGV
      READ(99,REC=IHOTSTP+14) IGCP
      READ(99,REC=IHOTSTP+15) NSCOUGC
      READ(99,REC=IHOTSTP+16) IGPP
      READ(99,REC=IHOTSTP+17) IGWP
      READ(99,REC=IHOTSTP+18) NSCOUGW
      CLOSE(99)

C--Open ASCII Hot Start File output 

      OPEN(99,FILE='hs_2_ascii.out')

C--Write out info to global file

      WRITE(99,*) IMHSF
      WRITE(99,*) TIMEHSF
      WRITE(99,*) ITHSF
      DO I=1,MNP
        WRITE(99,*) I,ETA1(I),ETA2(I),UU2(I),VV2(I),NODECODE(I)
        END DO
      DO I=1,MNE
        WRITE(99,*) NOFF(I)
      ENDDO
      WRITE(99,*) IESTP
      WRITE(99,*) NSCOUE
      WRITE(99,*) IVSTP
      WRITE(99,*) NSCOUV
      WRITE(99,*) ICSTP
      WRITE(99,*) NSCOUC
      WRITE(99,*) IPSTP
      WRITE(99,*) IWSTP
      WRITE(99,*) NSCOUM
      WRITE(99,*) IGEP
      WRITE(99,*) NSCOUGE
      WRITE(99,*) IGVP
      WRITE(99,*) NSCOUGV
      WRITE(99,*) IGCP
      WRITE(99,*) NSCOUGC
      WRITE(99,*) IGPP
      WRITE(99,*) IGWP
      WRITE(99,*) NSCOUGW
 
      CLOSE(99)           

      END




