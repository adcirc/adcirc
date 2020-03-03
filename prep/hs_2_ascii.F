C***************************************************************************
C  PROGRAM TO READ A BINARY HOTSTART FILE AND CREATE A SIMILAR ASCII FILE  *
C                                                                          *
C                      rl 10/11/01                                         *
C jgf updated for v45.06  09/07/2005                                       *
C***************************************************************************


      PROGRAM hs_2_ascii

      IMPLICIT NONE

      INTEGER I,J,IHOTSTP,MNP,MNE
!$omp threadprivate(I,J,IHOTSTP,MNP,MNE)
      INTEGER IMHSF,ITHSF
!$omp threadprivate(IMHSF,ITHSF)
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW
!$omp threadprivate(IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM)
!$omp threadprivate(IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW)
      CHARACTER FNAME*60
!$omp threadprivate(FNAME)
      INTEGER,ALLOCATABLE  :: NODECODE(:), NOFF(:)  !tcm 20110510 v50.05 added NOFF
!$omp threadprivate(NODECODE,NOFF)
      REAL(8),ALLOCATABLE :: ETA1(:),ETA2(:),UU2(:),VV2(:),CH1(:)
!$omp threadprivate(ETA1,ETA2,UU2,VV2,CH1)
      REAL(8) TIMEHSF
!$omp threadprivate(TIMEHSF)

C--Determine which hotstart files to process

      WRITE(*,*)'Enter binary hotstart file name'
      READ(*,*) FNAME
      OPEN (99+CK_LUN,FILE=FNAME,ACCESS='DIRECT',RECL=8)

C--Read stuff from each local file  

      WRITE(*,*) 'Enter the number of nodes in the file'
      READ(*,*) MNP
      WRITE(*,*) 'Enter the number of elements in the file'
      READ(*,*) MNE

C--Allocate local work arrays

      ALLOCATE ( ETA1(MNP),ETA2(MNP),UU2(MNP),
     &     VV2(MNP),NODECODE(MNP),CH1(MNP),NOFF(MNE) )           

      IHOTSTP=1
      READ(99+CK_LUN,REC=IHOTSTP) IMHSF
      IHOTSTP=2
      READ(99+CK_LUN,REC=IHOTSTP) TIMEHSF
      IHOTSTP=3
      READ(99+CK_LUN,REC=IHOTSTP) ITHSF

      DO I=1,MNP
        READ(99+CK_LUN,REC=IHOTSTP+1) ETA1(I)
        READ(99+CK_LUN,REC=IHOTSTP+2) ETA2(I)
        READ(99+CK_LUN,REC=IHOTSTP+3) UU2(I)
        READ(99+CK_LUN,REC=IHOTSTP+4) VV2(I)
        IHOTSTP=IHOTSTP+4
c       IF(IM.EQ.10) THEN
c         READ(LOCHSF,REC=IHOTSTP+1) CH1(I)
c         IHOTSTP=IHOTSTP+1
c         ENDIF
        READ(99+CK_LUN,REC=IHOTSTP+1) NODECODE(I)
        IHOTSTP=IHOTSTP+1
        ENDDO
      DO I=1,MNE
         READ(99+CK_LUN,REC=IHOTSTP+1) NOFF(I)                     
         IHOTSTP=IHOTSTP+1                                    
      END DO


C--Read in more common info from higest processor hotstart file

      READ(99+CK_LUN,REC=IHOTSTP+1 ) IESTP
      READ(99+CK_LUN,REC=IHOTSTP+2 ) NSCOUE
      READ(99+CK_LUN,REC=IHOTSTP+3 ) IVSTP
      READ(99+CK_LUN,REC=IHOTSTP+4 ) NSCOUV
      READ(99+CK_LUN,REC=IHOTSTP+5 ) ICSTP
      READ(99+CK_LUN,REC=IHOTSTP+6 ) NSCOUC
      READ(99+CK_LUN,REC=IHOTSTP+7 ) IPSTP
      READ(99+CK_LUN,REC=IHOTSTP+8 ) IWSTP
      READ(99+CK_LUN,REC=IHOTSTP+9 ) NSCOUM
      READ(99+CK_LUN,REC=IHOTSTP+10) IGEP
      READ(99+CK_LUN,REC=IHOTSTP+11) NSCOUGE
      READ(99+CK_LUN,REC=IHOTSTP+12) IGVP
      READ(99+CK_LUN,REC=IHOTSTP+13) NSCOUGV
      READ(99+CK_LUN,REC=IHOTSTP+14) IGCP
      READ(99+CK_LUN,REC=IHOTSTP+15) NSCOUGC
      READ(99+CK_LUN,REC=IHOTSTP+16) IGPP
      READ(99+CK_LUN,REC=IHOTSTP+17) IGWP
      READ(99+CK_LUN,REC=IHOTSTP+18) NSCOUGW
      CLOSE(99+CK_LUN)

C--Open ASCII Hot Start File output 

      OPEN(99+CK_LUN,FILE='hs_2_ascii.out')

C--Write out info to global file

      WRITE(99+CK_LUN,*) IMHSF
      WRITE(99+CK_LUN,*) TIMEHSF
      WRITE(99+CK_LUN,*) ITHSF
      DO I=1,MNP
        WRITE(99+CK_LUN,*) I,ETA1(I),ETA2(I),UU2(I),VV2(I),NODECODE(I)
        END DO
      DO I=1,MNE
        WRITE(99+CK_LUN,*) NOFF(I)
      ENDDO
      WRITE(99+CK_LUN,*) IESTP
      WRITE(99+CK_LUN,*) NSCOUE
      WRITE(99+CK_LUN,*) IVSTP
      WRITE(99+CK_LUN,*) NSCOUV
      WRITE(99+CK_LUN,*) ICSTP
      WRITE(99+CK_LUN,*) NSCOUC
      WRITE(99+CK_LUN,*) IPSTP
      WRITE(99+CK_LUN,*) IWSTP
      WRITE(99+CK_LUN,*) NSCOUM
      WRITE(99+CK_LUN,*) IGEP
      WRITE(99+CK_LUN,*) NSCOUGE
      WRITE(99+CK_LUN,*) IGVP
      WRITE(99+CK_LUN,*) NSCOUGV
      WRITE(99+CK_LUN,*) IGCP
      WRITE(99+CK_LUN,*) NSCOUGC
      WRITE(99+CK_LUN,*) IGPP
      WRITE(99+CK_LUN,*) IGWP
      WRITE(99+CK_LUN,*) NSCOUGW
 
      CLOSE(99+CK_LUN)           

      END




