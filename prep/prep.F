C----------------------------------------------------------------------------
C
C                           MODULE PREP
C
C----------------------------------------------------------------------------
C      current for ADCIRC v46.44   10/07/2006
C
C      Version 1.1 5/04/99 vjp
c      jjw fixes 053100
C      Revisions by rl 10/12/01, MEB 3/03, rl 3/03, rl 5/21/03, rl 5/18/04, 
C                   vp 11/27/03 (by rl)
C      Revisions by MEB 4/04
C      jgf Created PREP11, rewrote PREP20 for ADCIRC v45.12 02/24/2006 
C      jgf Created PREP13 for ADCRIC v46.00 
C      vjp PREP13 updates for ADCIRC v46.44
C      vjp PREP67_68 localization updates for ADCIRC v46.44
C      vjp added Relocalize for relocalizing fort.13 and fort.15
C
C----------------------------------------------------------------------------

   
C---------------------------------------------------------------------------
C            S U B R O U T I N E   P R E P  N E T  C D F
C---------------------------------------------------------------------------
C
C     jgf48.08 This subroutine will initialize the fulldomain netcdf
C     output files for the parallel run. These files are initialized 
C     in the adcprep phase because they contain fulldomain data which
C     are not available to the subdomains during a parallel run.
C
C     Revision history:
C
C     Date      Programmer                       Description of change
C     ----      ----------                       ---------------------
C     03/1/08  Cristina Forbes, PSGS @ UNC-IMS  Added memory allocation, 
C                                               added boundary vector arrays,
C                                               corrected coordinate systems,
C                                               corrected stations' name variables,
C                                               removed variable initialization,
C                                               added logic for switching I/O mode,
C                                               added switch to subroutine calls,
C                                               corrected subroutine arguments,
C                                               code clean-up.
C      =================================================================
C
C---------------------------------------------------------------------------
      SUBROUTINE PREPNETCDF()
C---------------------------------------------------------------------------
#ifdef NETCDF
      USE PRESIZES, ONLY: NM, MNTIF, MNE, MNP
      USE VERSION
C      USE SIZES
      USE PRE_GLOBAL, ONLY: NSTAE, NOUTE, RUNDES, RUNID, AGRID,
     &    DTDP, TAU0, STATIM, REFTIM, RNDAY, DRAMP, A00, B00, C00,
     &    H0, SLAM0, SFEA0, CF, ESLM, CORI, DP, NBDV, NBVV,
     &    SLEL, SFEL, ICS, NT, NTRSPE, IHOT, NOLIBF, NOLIFA,
     &    XEV,YEV,XEM,YEM,XEL,YEL,X,Y,NELG,DT, NTIF,
     &    NOLICA, NOLICAT, NWP, NCOR, NTIP, NWS, NRAMP, NBFR,
     &    NHY, NOPE, NETA, NBOU, NVEL, RAD2DEG,
     &    NVDLL, NVELL, IBTYPE,
     &    NSTAV, NOUTV, SLEV, SFEV, NTRSPV, 
     &    NSTAM, NOUTM, SLEM, SFEM, NTRSPM, 
     &    NOUTGE, SLAM, SFEA, NDSETSE,
     &    NOUTGV, NDSETSV,
     &    NOUTGW, NDSETSW, 
     &    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &    nhstar,nhsinc,
     &    title, institution,
     &    source, history, references, comments, host, convention,
     &    contact, base_date
      USE GLOBAL, ONLY: STATNAME, STATNAMEV, STATNAMEM,
     &                  imhs, iestp,nscoue,ivstp,nscouv, icstp,
     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                  nscougw,
     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
     &                    nscouc, ncipstp, nciwstp, nccoum, 
     &                    ncigep, nccouge, ncigvp, nccougv, 
     &                    igcp, nscougc, ncigpp, ncigwp,
     &                    nccougw   
      use memory_usage
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype,
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22,
     &     eta11,etadd, nodecodenc, noffnc, uu2, vv2,
     &     writenc61, writenc62, writenc63, writenc64,
     &     writenc71, writenc72, writenc73, writenc74,
     &     writenc67,
     &     nbvvnc,
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc,
     &     ncmode
     
      IMPLICIT NONE      
      INTEGER :: nbytes = 0
      INTEGER lun
      REAL(8) time
      INTEGER num_records,num_recordse, it
      INTEGER i, k,j
      INTEGER nface
      INTEGER, PARAMETER :: SNAMELEN = 50
!      INTEGER FileFmtMajor, FileFmtMinor, FileFmtRev

      LOGICAL FOUND

!      time = 0.d0
      timenc_len=1
      nface=3
      ncmode=0
!         FileFmtMajor = 1
!         FileFmtMinor = 1
!         FileFmtRev   = 0


      IF(.NOT. ALLOCATED(timenc)) ALLOCATE(timenc(timenc_len))
!      timenc(timenc_len)=time

      IF (ABS(NOUTE).EQ.3) THEN  !fort.61
         lun = 61
         num_records = NSTAE
         IF(.NOT. ALLOCATED(et00nc))
     &            ALLOCATE(et00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xenc))
     &            ALLOCATE(xenc(num_records))
         IF(.NOT. ALLOCATED(yenc))
     &            ALLOCATE(yenc(num_records))
         IF (.NOT.ALLOCATED(STATNAME))
     &            ALLOCATE(STATNAME(num_records,SNAMELEN))
         
        IF(ICS.EQ.2) THEN
            DO I=1, num_records
               xenc(i)=SLEL(I)
               yenc(i)=SFEL(I)
            END DO
            ELSE
               DO I=1,num_records
             	  xenc(i)=XEL(I)
                  yenc(i)=YEL(I)
               END DO
         ENDIF

         CALL writenc61 (NOUTE,
     &          lun, timenc_len, num_records, timenc,
     &          et00nc, nt, NTRSPE, rundes, runid, agrid, irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xenc, yenc, dp,
     &          nope, neta, nbdv,
     &          nbou, nvel, nbvv,
     &          ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.61


      IF (ABS(NOUTV).EQ.3) THEN  !fort.62
         lun = 62
         num_records = NSTAV
         IF(.NOT. ALLOCATED(UU00nc))
     &      ALLOCATE(UU00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(VV00nc))
     &      ALLOCATE(VV00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xvnc))
     &      ALLOCATE(xvnc(num_records))
         IF(.NOT. ALLOCATED(yvnc))
     &      ALLOCATE(yvnc(num_records))
         IF (.NOT.ALLOCATED(STATNAMEV))
     &       ALLOCATE(STATNAMEV(num_records,50))

         IF(ICS.EQ.2) THEN
            DO I=1, num_records
                xvnc(i)=SLEV(I)
                yvnc(i)=SFEV(I)
            END DO
            ELSE
               DO I=1,num_records
              	  xvnc(i)=XEV(I)
                  yvnc(i)=YEV(I)
               END DO
         ENDIF

         CALL writenc62 (NOUTV,
     &          lun, timenc_len, num_records, timenc,
     &          nt,NTRSPV,UU00nc,VV00nc,rundes,runid,agrid,irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xvnc, yvnc, dp,
     &          nope, neta, nbdv,
     &          nbou, nvel, nbvv,
     &                      ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.62


      IF (ABS(NOUTGE).EQ.3) THEN  !fort.63
         lun = 63
         num_records = MNP
         IF(.NOT. ALLOCATED(eta22))
     &      ALLOCATE(eta22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xnc))
     &      ALLOCATE(xnc(num_records))
         IF(.NOT. ALLOCATED(ync))
     &      ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtype)) ALLOCATE(ibtype(nbou))

               DO I=1,num_records
             	     xnc(i)=X(i)
                     ync(i)=Y(i)
       	       END DO
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO

            WRITE(633,*)"**noutge, fileunit,timenc_len, np, timenc",
     &             noutge, fileunit,timenc_len, num_records, timenc,
     &             "**eta22",
     &             eta22, "** nt = ", nt, "** ndsetse = ", ndsetse,
     &             "** rundes = ", rundes, runid, agrid,
     &             "** irtype = ", irtype,
     &             "** dtdp = ", dtdp,  "** ihot = ", ihot,
     &             "** ics, nolibf, nolifa = ", ics, nolibf, nolifa,
     &             "** nolica, nolicat, nwp = ",nolica, nolicat, nwp,
     &             "** ncor, ntip, nws, nramp, tau0 = ", ncor, ntip,
     &             nws, nramp, tau0,
     &             "** statim,reftim, rnday, dramp = ",statim,
     &             reftim, rnday, dramp,
     &             "** a00, b00, c00, h0 = ",a00, b00, c00, h0,
     &             "** slam0, sfea0, cf,eslm, cori = ", slam0, sfea0,
     &             cf,eslm, cori,
     &             "**  ntif, nbfr, ne, nhy, nm = ", ntif, nbfr,
     &             nhy, nm,
     &             "**  x = ", xnc, "**  y = ", ync, "**  dp = ", dp,
     &             "nope, neta, nbdv = ", nope, neta, nbdv,
     &             "nbou, nvel, nbvv = ", nbou, nvel, nbvvnc, 
     &              "nbounc", nbounc,"nvellnc", nvellnc,
     &             "ibtypenc",ibtype,"netanc",netanc,"nopenc",nopenc,
     &             "nvdllnc",nvdllnc,
     &             title, institution, source, history,
     &             references, comments, host, convention, contact,
     &             base_date


         CALL writenc63 (NOUTGE,
     &          lun, timenc_len, num_records, timenc,
     &          eta22, nt, NDSETSE, rundes, runid, agrid, irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xnc, ync, dp,
     &          nope, neta, nbdv,nopenc,netanc,
     &          nvdllnc, nbou, nvel, nbvvnc,nbounc,nvelnc,
     &          nvellnc,ibtype, ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.63

      IF (ABS(NOUTGV).EQ.3) THEN  !fort.64
         lun = 64
         num_records = MNP
         IF(.NOT. ALLOCATED(uu22))
     &      ALLOCATE(uu22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(vv22))
     &      ALLOCATE(vv22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xnc))
     &      ALLOCATE(xnc(num_records))
         IF(.NOT. ALLOCATED(ync))
     &      ALLOCATE(ync(num_records))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtype)) ALLOCATE(ibtype(nbou))
              
               DO I=1,num_records
             	  xnc(i)=X(i)
                  ync(i)=Y(i)
               END DO
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO

         CALL writenc64(NOUTGV,
     &          lun, timenc_len, num_records, timenc,
     &          nt,NDSETSV,uu22,vv22,rundes,runid,agrid,irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xnc, ync, dp,
     &          nope, neta, nbdv,nopenc,netanc,
     &          nvdllnc, nbou, nvel, nbvvnc,nbounc,nvelnc,
     &          nvellnc,ibtype, ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.64

      IF (ABS(NOUTM).EQ.3) THEN  !fort.71
         lun = 71
         num_records = NSTAM
         IF(.NOT. ALLOCATED(RMP00nc))
     &      ALLOCATE(RMP00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xpnc))
     &      ALLOCATE(xpnc(num_records))
         IF(.NOT. ALLOCATED(ypnc))
     &      ALLOCATE(ypnc(num_records))
        IF (.NOT.ALLOCATED(STATNAMEM)) 
     &       ALLOCATE(STATNAMEM(num_records,SNAMELEN))
         
         IF(ICS.EQ.2) THEN
            DO I=1, num_records
                xpnc(i)=SLEM(I)
                ypnc(i)=SFEM(I)
                
            END DO
            ELSE
                DO I=1,num_records
          	   xpnc(i)=XEM(i)                       !maybe it should be xpnc(i)
                   ypnc(i)=YEM(i)
               END DO
         ENDIF

               DO I=1,nope
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO

         CALL writenc71 (NOUTM,
     &          lun, timenc_len, num_records, timenc,
     &          RMP00nc, nt, NTRSPM, rundes, runid, agrid, irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xpnc, ypnc, dp,
     &          nope, neta, nbdv,
     &          nbou, nvel, nbvvnc,
     &                      ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.71


      IF (ABS(NOUTM).EQ.3) THEN  !fort.72
         lun = 72
         num_records = NSTAM
         IF(.NOT. ALLOCATED(RMU00nc))
     &      ALLOCATE(RMU00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(RMV00nc))
     &      ALLOCATE(RMV00nc(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xwnc))
     &      ALLOCATE(xwnc(num_records))
         IF(.NOT. ALLOCATED(ywnc))
     &      ALLOCATE(ywnc(num_records))
         IF (.NOT.ALLOCATED(STATNAMEM)) 
     &    ALLOCATE(STATNAMEM(num_records,50))
         
         IF(ICS.EQ.2) THEN
            DO I=1, num_records
                xwnc(i)=SLEM(I)
                ywnc(i)=SFEM(I)
            END DO
            ELSE
               DO I=1,num_records
             	  xwnc(i)=XEM(i)
                  ywnc(i)=YEM(i)
               END DO
         ENDIF

               DO I=1,nbou
                  DO J=1,nvel
             	  nbvvnc(i,j)=nbvv(i,j)
                  END DO
               END DO

         CALL writenc72(NOUTM,
     &          lun, timenc_len, num_records, timenc,
     &          nt,NTRSPM,RMU00nc,RMV00nc,rundes,runid,agrid,irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xwnc, ywnc, dp,
     &          nope, neta, nbdv,
     &          nbou, nvel, nbvvnc,
     &                      ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.72

      
      IF (ABS(NOUTGW).EQ.3) THEN  !fort.73
         lun = 73
         num_records = MNP
         IF(.NOT. ALLOCATED(pr22))
     &      ALLOCATE(pr22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xnc))
     &      ALLOCATE(xnc(num_records))
         IF(.NOT. ALLOCATED(ync))
     &      ALLOCATE(ync(num_records))

            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtype)) ALLOCATE(ibtype(nbou))

               DO I=1,num_records
             	     xnc(i)=X(i)
                     ync(i)=Y(i)
       	       END DO
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO

         CALL writenc73 (NOUTGW,
     &          lun, timenc_len, num_records, timenc,
     &          pr22, nt, NDSETSW, rundes, runid, agrid, irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xnc, ync, dp,
     &          nope, neta, nbdv,nopenc,netanc,
     &             nvdllnc,
     &             nbou, nvel, nbvvnc,nbounc,nvelnc,
     &             nvellnc,ibtype,
     &                      ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.73

      IF (ABS(NOUTGW).EQ.3) THEN  !fort.74
         lun = 74
         num_records = MNP
         IF(.NOT. ALLOCATED(WNCX22))
     &      ALLOCATE(WNCX22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(WNCY22))
     &      ALLOCATE(WNCY22(num_records,timenc_len))
         IF(.NOT. ALLOCATED(xnc))
     &      ALLOCATE(xnc(num_records))
         IF(.NOT. ALLOCATED(ync))
     &      ALLOCATE(ync(num_records))
               DO I=1,num_records
             	        xnc(i)=X(i)
                	ync(i)=Y(i)
               END DO
         
               netanc=0
               DO i=1,nope
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
              END DO

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
               
               k=0
                  DO i=1,nope
                        DO j=1,nvdllnc(i)
                        k=k+1
             	        nopenc(k)=nbdv(i,j)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO


         CALL writenc74(NOUTGW,
     &          lun, timenc_len, num_records, timenc,
     &          nt,NDSETSW,WNCX22,WNCY22,rundes,runid,agrid,irtype,
     &          dt, ihot, ics, nolibf, nolifa,
     &          nolica, nolicat, nwp,
     &          ncor, ntip, nws, nramp, tau0,
     &          statim,reftim, rnday, dramp,
     &          a00, b00, c00, h0,
     &          slam0, sfea0, cf,eslm, cori,
     &          NTIF, nbfr, NELG, nface,nm,
     &          xnc, ync, dp,
     &          nope, neta, nbdv,nopenc,netanc,
     &          nvdllnc, nbou, nvel, nbvvnc,nbounc,nvelnc,
     &          nvellnc,ibtype, ncmode,
     &          title, institution, source, history,
     &          references, comments, host, convention, contact,
     &          base_date)

      ENDIF !fort.74
      
      PRINT *, "NHSTAR BEFORE ************** = ", NHSTAR
      
      IF (ABS(NHSTAR).EQ.3) THEN  !fort.67 or fort.68
      
      PRINT *, "======================================================="
      PRINT *, "******** PREPING NETCDF HOTSTART FILES ************* = "
!      PRINT *, "NHSTAR AFTER ************** = ", NHSTAR
      
!      PRINT *, "================================ DTMSG = ",DTMSG
!      PRINT *, "================================ DT = ",DT
!      PRINT *, "================================ DP = ",DP
!      PRINT *, "================================ DTDP = ",DTDP
!      PRINT *, "================================ statim = ",statim
!      PRINT *, "================================ reftim = ",reftim
!      PRINT *, "================================ statim = ",statim
!      PRINT *, "================================ rnday = ",rnday
!      PRINT *, "================================ slam0 = ",slam0
!      PRINT *, "================================ sfea0 = ",sfea0
!      PRINT *, "================================ rnday = ",rnday
!      PRINT *, "================================ eslm = ",eslm
!      PRINT *, "================================ cori = ",cori
!      PRINT *, "================================ ntif = ",ntif
      
!         if(IHOT.eq.0.or.IHOT.eq.-68)lun=67
!         if(IHOT.eq.-67)lun=68
!        num_records=np_g
!        num_recordse = NE_G
         num_records = MNP
         num_recordse = MNE
!        InputFileFmtVn=001001000

        IF(.NOT.ALLOCATED(eta11))ALLOCATE(eta11(num_records,timenc_len))
        IF(.NOT.ALLOCATED(eta22))ALLOCATE(eta22(num_records,timenc_len))
        IF(.NOT.ALLOCATED(etadd))ALLOCATE(etadd(num_records,timenc_len))
        IF(.NOT.ALLOCATED(uu22))ALLOCATE(uu22(num_records,timenc_len))
        IF(.NOT.ALLOCATED(vv22))ALLOCATE(vv22(num_records,timenc_len))
        IF(.NOT.ALLOCATED(nodecodenc))ALLOCATE(nodecodenc(num_records))
        IF(.NOT.ALLOCATED(noffnc))ALLOCATE (noffnc(num_recordse))
      
!             IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
!           IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
        IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(num_records))
        IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(num_records))
        
        IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
        IF(.NOT. ALLOCATED(NBDV)) ALLOCATE(NBDV(NOPE,NETA))
        IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
        IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
        IF(.NOT. ALLOCATED(ibtype)) ALLOCATE(ibtype(nbou))

!        PRINT *,"nope = ", nope
!        PRINT *,"neta = ", neta
        netanc=neta
!        PRINT *,"netanc = ", netanc
       
               DO I=1,num_records
             	     xnc(i)=X(i)
                     ync(i)=Y(i)
       	       END DO
!
!            elevation specified boundary forcing segments
!
               netanc=0
               DO i=1,nope
!        PRINT *,"nvdllnc b ",i, nvdllnc(i)
                  nvdllnc(i)=nvdll(i)
                  netanc=netanc+nvdllnc(i)
!        PRINT *,"nvdllnc a",i, nvdllnc(i)
              END DO
!        PRINT *,"netanc = ", netanc

            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
!        PRINT *,"allocate nopenc = "
!        PRINT *,"allocate nope = ", nope
               
               k=0
                  DO i=1,nope
!        PRINT *,"nope ",i, nope
!        PRINT *,"nvdllnc(i) ",i, nvdllnc(i)
                        DO j=1,nvdllnc(i)
                        k=k+1
!        PRINT *,"nopenc b",k, nopenc(k)
             	        nopenc(k)=nbdv(i,j)
!        PRINT *,"nopenc a",k, nopenc(k)
                     END DO
                  END DO
!
!              normal flow (discharge) specified boundary segments
!
               nvelnc=0
               DO i=1,nbou
                  nvellnc(i)=nvell(i)
                  nvelnc=nvelnc+nvellnc(i)
                  DO j=1,nvel
             	     nbvvnc(i,j)=nbvv(i,j)
                 END DO
               END DO

            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
               
               k=0
                  DO i=1,nbou
                        DO j=1,nvellnc(i)
                        k=k+1
             	        nbounc(k)=nbvv(i,j)
                     END DO
                  END DO
!        PRINT *,"nvelnc = ", nvelnc
!        PRINT *,"netanc = ", netanc
                  
!          PRINT *, "==== write_output: IM = ", IM, IHOTSTP
!          PRINT *, "==== write_output: TIME = ", TIME, IHOTSTP
!          PRINT *, "==== write_output: IM = ", IT, IHOTSTP
!          PRINT *, "==== write_output: NPX = ", MNP, IHOTSTP
!          PRINT *, "==== write_output: NEX = ", MNE, IHOTSTP


          lun=67
          
C         Determine if netcdf file exists
          INQUIRE(FILE="fort.67.nc",EXIST=FOUND)
          
          IF ( .NOT. FOUND ) THEN
!          PRINT *, "LUN = ", LUN
!     &                    iestp,nscoue,ivstp,nscouv, icstp,
!     &                    nscouc, ipstp, iwstp, nscoum, igep, nscouge,
!     &                    igvp, nscougv, igcp, nscougc, igpp, igwp,
!     &                    nscougw,
     
          call writenc67 (nhstar,abs(lun),timenc_len,num_records,timenc,
     &                    nt, nhsinc, uu22,vv22,
     &                    eta11, eta22, etadd, nodecodenc,noffnc,
     &                    imhs, it,
     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
     &                    nscouc, ncipstp, nciwstp, nccoum, 
     &                    ncigep, nccouge, ncigvp, nccougv, 
     &                    igcp, nscougc, ncigpp, ncigwp,
     &                    nccougw,     
     &                    rundes, runid, agrid, irtype,
     &                    dt, ihot, ics, nolibf, nolifa,
     &                    nolica, nolicat, nwp,
     &                    ncor, ntip, nws, nramp, tau0,
     &                    statim, reftim, rnday, dramp,
     &                    a00, b00, c00, h0,
     &                    slam0, sfea0, cf ,eslm, cori,
     &                    ntif, nbfr, NELG, nface, nm,
     &                    xnc, ync, dp,
     &                    nope, neta, nbdv,nopenc,netanc, nvdllnc,
     &                    nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                    nvellnc,ibtype, ncmode,
     &                    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &                    title, institution, source, history,
     &                    references, comments, host, convention, 
     &                    contact, base_date)
                  
          ENDIF
          
          
          INQUIRE(FILE="fort.68.nc",EXIST=FOUND)
          
          IF ( .NOT. FOUND ) THEN
          
          lun=68
!          PRINT *, "LUN = ", LUN
!     iestp,nscoue,ivstp,nscouv, icstp,
!     &                    nscouc, ipstp, iwstp, nscoum, igep, nscouge,
!     &                    igvp, nscougv, igcp, nscougc, igpp, igwp,
!     &                    nscougw,




          call writenc67 (nhstar,abs(lun),timenc_len,num_records,timenc,
     &                    nt, nhsinc, uu22,vv22,
     &                    eta11, eta22, etadd, nodecodenc,noffnc,
     &                    imhs, it,
     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
     &                    nscouc, ncipstp, nciwstp, nccoum, 
     &                    ncigep, nccouge, ncigvp, nccougv, 
     &                    igcp, nscougc, ncigpp, ncigwp,
     &                    nccougw,
     &                    rundes, runid, agrid, irtype,
     &                    dt, ihot, ics, nolibf, nolifa,
     &                    nolica, nolicat, nwp,
     &                    ncor, ntip, nws, nramp, tau0,
     &                    statim, reftim, rnday, dramp,
     &                    a00, b00, c00, h0,
     &                    slam0, sfea0, cf ,eslm, cori,
     &                    ntif, nbfr, NELG, nface, nm,
     &                    xnc, ync, dp,
     &                    nope, neta, nbdv,nopenc,netanc, nvdllnc,
     &                    nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                    nvellnc,ibtype, ncmode,
     &                    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &                    title, institution, source, history,
     &                    references, comments, host, convention, 
     &                    contact, base_date)
     
                  ENDIF
                  
                  ENDIF

      call memory_status()
      RETURN
#endif
C---------------------------------------------------------------------------
      END SUBROUTINE PREPNETCDF
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 0
C---------------------------------------------------------------------------
C
C     jgf46.28 from jgf45.16 This subroutine will break up the full
C     domain initial concentration file into subdomains. The fort.10
C     file may contain initial concentration for either 2D or 3D ADCIRC
C     runs.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP10()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL 
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J
      INTEGER iproc          ! subdomain index
      CHARACTER*80 header1, header2  ! header comments in unit 10 files
      CHARACTER*80 nvn_nvp           ! string representing nfen, np
      REAL(SZ) nvn           ! number of vertical nodes from unit 10 file
      REAL(SZ) nvp           ! number of horizontal nodes from unit 10 file
      INTEGER nhnn           ! horizontal nodes counter
      INTEGER nvnn           ! vertical nodes counter
      INTEGER sdu(nproc)     ! subdomain unit numbers for unit 10 files
      REAL(SZ), ALLOCATABLE :: fdData2D(:)   !(MNP)      full domain data
      REAL(SZ), ALLOCATABLE :: fdData3D(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: sdData2D(:)   !(MNP)      subdomain data
      REAL(SZ), ALLOCATABLE :: sdData3D(:,:) !(MNP,NFEN) subdomain data
      LOGICAL success        ! .true. if all files opened successfully 

      CALL OpenPrepFiles(10, 'initial concentration         ', 
     &     1, nproc, sdu, success)

      IF (.not.success) THEN
         WRITE(*,*) 'WARNING: Unit 10 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C     Read header information from full domain unit 10 file
      READ(10,80) header1
      READ(10,80) header2
C
C     Transcribe header information.
      DO iproc = 1, nproc
         WRITE(sdu(iproc),80)  header1
         WRITE(sdu(iproc),80)  header2
      ENDDO
C
C     Check node number data for consistency (paranoia).
      READ(10,80) nvn_nvp
      READ(nvn_nvp,*) nvn, nvp
      IF ( nvn .ne. nfen .or. nvp .ne. nnodg ) then
         WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
         WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
         WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
      ENDIF
C
C     Decompose concentration data
C
      IF (C2D_PTrans) THEN
c     read in the full domain data
         ALLOCATE ( fdData2D(MNP) )
         nbytes = 8*mnp
         call memory_alloc(nbytes)
         DO i=1, NNODG
            READ(10,*) nhnn, fdData2D(nhnn)
         ENDDO
c     write out subdomain data
         ALLOCATE ( sdData2D(MNP) )
         nbytes = 8*mnp
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nnodp(iproc)
            DO i=1, nnodp(iproc)
               sdData2D(i) = fdData2D(IMAP_NOD_LG(i,iproc))
               WRITE(sdu(iproc),*) i, sdData2D(i)
            ENDDO
         ENDDO
         DEALLOCATE ( fdData2D, sdData2D )
         nbytes = 16*mnp
         call memory_dealloc(nbytes)
      ENDIF
C
      IF (C3D_PTrans) THEN
c     read in the full domain data
         ALLOCATE ( fdData3D(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(10,*) nhnn, nvnn, fdData3D(nhnn,nvnn)
            ENDDO
         ENDDO
c     write out subdomain data
         ALLOCATE ( sdData3D(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData3D(i,j) = fdData3D(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData3D(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData3D, sdData3D )
         nbytes = 16*mnp*nfen
         call memory_dealloc(nbytes)
      ENDIF
C
C     Close full domain file and all subdomain files
      CLOSE(10)
      DO iproc=1, nproc
         CLOSE(sdu(iproc))
      ENDDO
C
  80  FORMAT(A80)
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP10
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 1
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will break up the full domain initial
C     density forcing file into subdomains. The fort.11 file may
C     contain initial density, temperature, and/or salinity depending on
C     the value of IDEN in the fort.15 file.
C
C     jgf45.12 This subroutine is designed to work for baroclinic 3D
C     runs only, not 2D runs.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP11()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL 
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J
      INTEGER iproc          ! subdomain index
      CHARACTER*80 header1, header2  ! header comments in unit 11 files
      CHARACTER*80 nvn_nvp           ! string representing nfen, np
      REAL(SZ) nvn           ! number of vertical nodes from unit 11 file
      REAL(SZ) nvp           ! number of horizontal nodes from unit 11 file
      INTEGER nhnn           ! horizontal nodes counter
      INTEGER nvnn           ! vertical nodes counter
      INTEGER sdu(nproc)     ! subdomain unit numbers for unit 11 files
      REAL(SZ), ALLOCATABLE :: fdData1(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: fdData2(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: sdData1(:,:) !(MNP,NFEN) subdomain data
      REAL(SZ), ALLOCATABLE :: sdData2(:,:) !(MNP,NFEN) subdomain data
      LOGICAL success        ! .true. if all files opened successfully 

      CALL OpenPrepFiles(11, 'initial density forcing       ', 
     &     1, nproc, sdu, success)

      IF (.not.success) THEN
         WRITE(*,*) 'WARNING: Unit 11 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C     Read header information from full domain unit 11 file
      READ(11,80) header1
      READ(11,80) header2
C
C     Transcribe header information.
      DO iproc = 1, nproc
         WRITE(sdu(iproc),80)  header1
         WRITE(sdu(iproc),80)  header2
      ENDDO
C
C     Check node number data for consistency (paranoia).
      READ(11,80) nvn_nvp
      READ(nvn_nvp,*) nvn, nvp
      IF ( nvn .ne. nfen .or. nvp .ne. nnodg ) then
         WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
         WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
         WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
      ENDIF
C
C     Decompose density forcing data; format based on value of IDEN.
C     jgf45.12 This is designed to work for baroclinic 3D runs only, not
C     2D runs.

c     read in the full domain data
      SELECT CASE (ABS(IDEN))
      
      CASE(1,2,3)
         ALLOCATE ( fdData1(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(11,*) nhnn, nvnn, fdData1(nhnn,nvnn)
            ENDDO
         ENDDO

      CASE(4)
         ALLOCATE ( fdData1(MNP,NFEN) )
         ALLOCATE ( fdData2(MNP,NFEN) )
         nbytes = 16*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(11,*) nhnn, nvnn, 
     &              fdData1(nhnn,nvnn),fdData2(nhnn,nvnn)
            ENDDO
         ENDDO

      END SELECT
         
c     write out subdomain data
      SELECT CASE (ABS(IDEN))
         
      CASE(1,2,3) 
         ALLOCATE ( sdData1(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData1(i,j) = fdData1(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData1(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData1, sdData1 )
         nbytes = 16*mnp*nfen
         call memory_dealloc(nbytes)
         
      CASE(4) 
         ALLOCATE ( sdData1(MNP,NFEN) )
         ALLOCATE ( sdData2(MNP,NFEN) )
         nbytes = 16*mnp*nfen
         call memory_alloc(nbytes)
          DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData1(i,j) = fdData1(IMAP_NOD_LG(i,iproc),j)
                  sdData2(i,j) = fdData2(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData1(i,j), sdData2(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData1, fdData2 )
         DEALLOCATE ( sdData1, sdData2 )
         nbytes = 32*mnp*nfen
         call memory_dealloc(nbytes)
         
      END SELECT
C
C     Close full domain file and all subdomain files
      CLOSE(11)
      DO iproc=1, nproc
         CLOSE(sdu(iproc))
      ENDDO
C
  80  FORMAT(A80)
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP11
C---------------------------------------------------------------------------




C---------------------------------------------------------------------------
C                S U B R O U T I N E   R E L O C A L I Z E
C---------------------------------------------------------------------------
C
C     This routine allows re-localizing selected files after a prepall 
C     operation.  vjp 10/2006
C
C---------------------------------------------------------------------------
      SUBROUTINE RELOCALIZE()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I, J, IPROC, INDEX, ITEMP, idumy
      CHARACTER(14) LOCFN
      CHARACTER(80) skipped
      INTEGER,ALLOCATABLE :: LUNP(:)
C
      print *, "entering relocalize"
!      print *, "nproc = ", nproc

      allocate( lunp(nproc) )  ! logical unit number for each subdomain
      do iproc=1, nproc
         lunp(iproc) = 105 + (iproc-1)
      enddo

      if (.not.allocated(nnodp)) then
        ALLOCATE(NNODP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(nelp)) then
        ALLOCATE(NELP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(netap)) then
        ALLOCATE(NETAP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(nfluxfp)) then
        ALLOCATE(NFLUXFP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
         LOCFN(1:14) = 'PE0000/fort.18'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (LUNP(IPROC),FILE=LOCFN)
      ENDDO


      print *, "from relocalize: reading local-to-global element maps"
      DO IPROC = 1,NPROC
         READ(LUNP(IPROC),'(A)') skipped    !  read past fileFmt header
         READ(LUNP(IPROC),'(8X,3I8)') NELG, MNEP, NELP(IPROC)
      ENDDO
!      print *, "nelg = ", nelg

      if (.not.allocated(imap_el_lg)) then
        ALLOCATE ( IMAP_EL_LG(MNEP, NPROC) )
        nbytes = 4*nproc*mnep
        call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
      DO I=1, NELP(IPROC)
         READ(LUNP(IPROC),'(I8)') idumy                 
         IMAP_EL_LG(I,IPROC) = abs(idumy)
      ENDDO
      ENDDO

      print *, "from relocalize: reading local-to-global node maps"
      DO IPROC = 1,NPROC
         READ(LUNP(IPROC),'(8X,3I8)') NNODG, MNPP, NNODP(IPROC)
      ENDDO
!      print *, "nnodg = ", nnodg

      if (.not.allocated(imap_nod_lg)) then
         ALLOCATE ( IMAP_NOD_LG(MNPP, NPROC) )
         nbytes = 4*nproc*mnpp
         call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
      DO I=1, NNODP(IPROC)
         READ(LUNP(IPROC),'(I8)') idumy                 
         IMAP_NOD_LG(I,IPROC) = abs(idumy)
      ENDDO
      ENDDO

C  This section for prep15
      IF (PREP_15) THEN
        print *, "from relocalize: reading nfluxf for each subdomain"
        DO IPROC = 1,NPROC
           READ(LUNP(IPROC),'(8X,I8)') NFLUXFP(IPROC)
        ENDDO

        print *, "from relocalize: reading neta for each subdomain"
        DO IPROC = 1,NPROC
           READ(LUNP(IPROC),'(8X,3I8)') idumy, NETA_MAX, NETAP(IPROC)
        ENDDO

        if (.not.allocated(obnode_lg)) then
           ALLOCATE ( OBNODE_LG(NETA_MAX, NPROC) )
           nbytes = 4*nproc*neta_max
           call memory_alloc(nbytes)
        endif

        print *, "from relocalize: reading open boundary table"
        DO IPROC = 1,NPROC
        DO I=1, NETAP(IPROC)
           READ(LUNP(IPROC),'(I8)') OBNODE_LG(I,IPROC)    
        ENDDO
        ENDDO
      ENDIF
      
! Build Global-to-Local Node Map
      if (.not.allocated(itotproc)) then
        ALLOCATE ( ITOTPROC(NNODG) )
        nbytes = 4*nnodg               
        call memory_alloc(nbytes)
      endif
      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO IPROC = 1,NPROC
         DO I = 1,NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            ITOTPROC(INDEX) = ITOTPROC(INDEX) + 1
         ENDDO
      ENDDO
      MNEI = 0    
      DO I = 1,NNODG
         IF (ITOTPROC(I) .gt. MNEI) MNEI = ITOTPROC(I)
      ENDDO
      print *, "MNEI = ", MNEI

      if (.not.allocated(imap_nod_gl2)) then
        ALLOCATE( IMAP_NOD_GL2( 2*MNEI, NNODG) )
        nbytes = 8*mnei*nnodg               
        call memory_alloc(nbytes)
      endif
      print *, "allocated imap_nod_GL2"

      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO IPROC = 1,NPROC
         DO J = 1,NNODP(IPROC)
            INDEX = IMAP_NOD_LG(J,IPROC)
            ITOTPROC(INDEX) = ITOTPROC(INDEX) + 1
            ITEMP = (ITOTPROC(INDEX)-1)*2 + 1
            IMAP_NOD_GL2(ITEMP,INDEX) = IPROC
            IMAP_NOD_GL2(ITEMP+1,INDEX) = J
         ENDDO
      ENDDO

      call memory_status()
      print *, "leaving relocalize"
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE RELOCALIZE
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 3
C---------------------------------------------------------------------------
C
C     jgf46.00 This subroutine will break up the full domain nodal
C     attributes file into subdomains.
C
C     jgf48.47 Rewritten to only open a smaller number of subdomains
C     at a time (256 by default)
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP13()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
C
      integer :: nbytes = 0
      INTEGER i              ! node loop counter
      INTEGER k              ! attribute loop counter
      INTEGER l              ! line loop counter
      INTEGER m              ! attribute default value counter
      INTEGER iproc          ! subdomain loop counter
      INTEGER sdu(nproc)     ! subdomain unit number for unit 13 files
      INTEGER NumOfNodes     ! must match the no. of nodes in grid (unit 14)
      INTEGER NAttr          ! number of nodal attributes in the file
      INTEGER NumNotDefault  ! number of nodes specified in the file
      CHARACTER(len=80) header       ! header comments in unit 13 files
      CHARACTER(len=80), ALLOCATABLE :: AttrName ! label for attribute
      CHARACTER(len=80), ALLOCATABLE :: Units ! label for physical units
      CHARACTER(len=80) Skipped      ! data we want to skip over
      REAL(SZ) DefaultVal(12)   ! default value of attribute
      INTEGER, ALLOCATABLE:: NoOfVals ! at each node for an attribute
      INTEGER Mode           !=0 to count, =1 to write
      LOGICAL success        ! .true. if all files opened successfully 
      INTEGER, ALLOCATABLE :: SDNumND(:,:)  ! subdomain # of nodes not default
C     jgf48.47 Do the decomposition for a max of 256 subdomains at a
C     time ... some platforms/compilers limit the number of files that 
C     can be open at any one time.
      INTEGER, PARAMETER :: maxOpenFiles = 256 
      INTEGER startProc  
      INTEGER endProc 
      INTEGER deltaProc
C
C     Perform decomposition over range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc
C
C        Open full domain and subdomain fort.13 files.
         CALL OpenPrepFiles(13, 'nodal attributes              ', 
     &       startProc, endProc, sdu, success)
C
         IF (.not.success) THEN
            WRITE(*,*) 'WARNING: Unit 13 files not preprocessed.'
            RETURN ! note early return
         ENDIF
C
C        Read header information from full domain unit 13 file
         READ(13,'(A80)') header
         READ(13,*) NumOfNodes     ! number of nodes according to unit 13
         READ(13,*) NAttr          ! number of attributes in the unit 13 file
C
C        Check to make sure that the number of nodes in the nodal
C        attributes file is the same as in the grid file (unit 14).
         IF (NumOfNodes.NE.NNODG) THEN
            WRITE(6,9900)
 9900       FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'The number of nodes in the grid file (unit 14) and'
     &        /,1X,'the nodal attributes file (unit 13) must match.',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP                   ! We're toast.
         ENDIF
C
C        Transcribe header information into subdomain unit 13 files
         DO iproc = startProc, endProc
            WRITE(sdu(iproc),'(A80)') header
            WRITE(sdu(iproc),*) NNODP(iproc) 
            WRITE(sdu(iproc),*) NAttr
         ENDDO
C
C        Transcribe attribute names from full domain file to subdomains.
         DO k=1, NAttr
            READ(13,'(A80)') AttrName
            READ(13,'(A80)') Units
            READ(13,*) NoOfVals
            READ(13,*) (DefaultVal(m),m=1,NoOfVals)
            DO iproc=startProc, endProc
               WRITE(sdu(IPROC),'(A80)') AttrName
               WRITE(sdu(IPROC),'(A80)') Units
               WRITE(sdu(IPROC),*) NoOfVals
               WRITE(sdu(IPROC),'(12(1x,e16.9))') 
     &          (DefaultVal(m),m=1,NoOfVals)
            END DO
         END DO
C
C        Allocate and initialize the matrix for the number of Non Default
C        nodes in each SubDomain for each nodal attribute
         ALLOCATE(SDNumND(nproc,NAttr))
         nbytes = 8*nproc*nattr
         call memory_alloc(nbytes)
         DO iproc=startProc,endProc
            DO k=1, NAttr
               SDNumND(iproc,k)=0
            END DO
         END DO
C        We need to figure out how many nodes go into each subdomain
C        for each attribute.
         CALL processNodalAttr(NAttr, 0, sdu, SDNumND,
     &      startProc, endProc)
C
C     Now rewind and advance to the beginning of the data again
         REWIND(13)
         DO l=1, 3
            READ(13,*) skipped    ! skip header, NumOfNodes, NAttr
         END DO
         DO k=1, NAttr
            DO l=1, 4
               READ(13,*) skipped ! skip AttrName,Units,NoOfVals,Default 
            END DO
         END DO
C
C        Now read each of the nodal attributes and transcribe them to the 
C        appropriate subdomain.
         CALL processNodalAttr(NAttr, 1, sdu, SDNumND,
     &      startProc, endProc)
         DEALLOCATE(SDNumND)
         nbytes = 8*nproc*nattr
         call memory_dealloc(nbytes)
C
C        Close full domain and subdomain files
         CLOSE (13)
         DO iproc=startProc, endProc
            CLOSE(sdu(iproc))
         ENDDO 
         startProc = endProc + 1
      END DO
C
      if (allocated(NodalAttributes)) then  
         DEALLOCATE(NodalAttributes)
         nbytes = 4*nwp
         call memory_dealloc(nbytes)
      endif
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP13
C---------------------------------------------------------------------------

C     ----------------------------------------------------------------
C        S U B R O U T I N E   P R O C E S S   N O D A L   A T T R 
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to support PREP13. This subroutine is called
C     twice; once to determine the number of nodes with non-default
C     values going into each subdomain, then a second time to actually
C     place the data in the subdomain files. 
C
C     This is necessary because the attributes in the subdomain files
C     must each have the number of non-default values at the top, and
C     this information cannot be known until we have processed the
C     entire fulldomain file. 
C     
C     ----------------------------------------------------------------
      SUBROUTINE ProcessNodalAttr(NAttr, mode, sdu, SDNumND,
     &   startProc, endProc)
      USE PRE_GLOBAL
      IMPLICIT NONE
C
      INTEGER,intent(in) :: NAttr  ! number of attributes in the file
      INTEGER,intent(in) :: Mode   !=0 to count and return, =1 to write
      INTEGER,intent(in),dimension(nproc) :: sdu !i/o unit number array 
      INTEGER,intent(inout),dimension(nproc,NAttr) :: SDNumND
      INTEGER,intent(in) :: startProc 
      INTEGER,intent(in) :: endProc
      INTEGER NumNotDefault      ! number of nodes specified in the file
      INTEGER NumCol             ! number of values per node for an attr
      INTEGER NodeNum            ! full domain node number
      INTEGER SDNode             ! subdomain node number
      INTEGER i                  ! node loop counter
      INTEGER j                  ! column loop counter
      INTEGER k                  ! attribute loop counter
      INTEGER m                  ! mapping loop counter
      INTEGER iproc              ! subdomain loop counter
      INTEGER iproc2             ! mapped subdomain
      CHARACTER(len=80) AttrName ! label for attribute
      REAL(SZ), ALLOCATABLE :: AttrData(:)  ! attribute data 
      CHARACTER(len=80) Skipped  ! data we want to skip over 
C
      DO k=1, NAttr
         READ(13,*) AttrName
         IF (Mode.eq.1) THEN
            DO iproc=startProc,endProc
               WRITE(sdu(iproc),'(A80)') AttrName
            END DO
         ENDIF
         SELECT CASE (AttrName)
         CASE("primitive_weighting_in_continuity_equation")
            NumCol=1
         CASE("surface_submergence_state")
            NumCol=1
         CASE("quadratic_friction_coefficient_at_sea_floor")
            NumCol=1
         CASE("mannings_n_at_sea_floor")
            NumCol=1
         CASE("chezy_friction_coefficient_at_sea_floor")
            NumCol=1
         CASE("sea_surface_height_above_geoid")
            NumCol=1
         CASE("surface_directional_effective_roughness_length") 
            NumCol=12
         CASE("surface_canopy_coefficient") 
            NumCol=1
         CASE("bridge_pilings_friction_parameters") 
            NumCol=4
         CASE
     &   ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth") 
            NumCol=1
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth") 
            NumCol=1
         CASE 
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            NumCol=2
         CASE
     &   ("subtract_from_depth_to_equal_mean_sea_level") 
            NumCol=1
         CASE DEFAULT
            NumCol=1
            IF (Mode.eq.0) THEN
               WRITE(6,1001)    ! Nodal Attributes file
               WRITE(6,1021) AttrName ! contains invalid name
            ELSE
               WRITE(6,1031) AttrName ! Process 1st column only
            ENDIF
         END SELECT
         ALLOCATE(AttrData(NumCol))
         READ(13,*) NumNotDefault
         IF (Mode.eq.1) THEN
            DO iproc=startProc,endProc
               WRITE(sdu(iproc),*) SDNumND(iproc,k)
            END DO
         ENDIF
         DO i=1, NumNotDefault
            READ(13,*) NodeNum, (AttrData(j),j=1,NumCol)
            IF (ITOTPROC(NodeNum).eq.1) THEN
               iproc = IMAP_NOD_GL(1,NodeNum)
               IF ( (iproc.lt.startProc) .or. (iproc.gt.endProc) ) THEN
                  CYCLE ! skip it if it does not map to our range of procs 
               ENDIF
               IF (Mode.eq.0) SDNumND(iproc,k) = SDNumND(iproc,k)+1
               IF (Mode.eq.1) THEN
                  SDNode = IMAP_NOD_GL(2,NodeNum)
                  WRITE(sdu(iproc),1100) SDNode,(AttrData(j),j=1,NumCol)
               ENDIF
            ELSE
               DO m=1, ITOTPROC(NodeNum)
                  iproc2 = IMAP_NOD_GL2(2*(m-1)+1,NodeNum)
                  DO iproc=startProc, endProc
                     IF (iproc.EQ.iproc2) THEN !f.d. node maps to this s.d. 
                        IF (Mode.eq.0) THEN
                           SDNumND(iproc,k)=SDNumND(iproc,k)+1
                        ENDIF
                        IF (Mode.eq.1) THEN
                           SDNode = IMAP_NOD_GL2(2*(m-1)+2,NodeNum)
                           WRITE(sdu(iproc),1100) SDNode,
     &                          (AttrData(j),j=1,NumCol)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            END IF
         END DO
         DEALLOCATE(AttrData)
         IF (Mode.eq.1) THEN
            WRITE(6,'(A25,A80)') '     Finished processing ', AttrName
            WRITE(6,*) 'for processor range ',startProc,' to ',endProc
         ENDIF
      END DO
C
C
 1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
 1021 FORMAT('contains invalid name: ',A80) 
 1031 FORMAT('WARNING: Processed only one column of unrecognized ',A80) 
 1100 FORMAT(I10,32000(2X,E16.8))
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ProcessNodalAttr
C     ----------------------------------------------------------------


      SUBROUTINE PREP14()
      USE PRE_GLOBAL 
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C 
C  This routine writes a Local Grid file "fort.14" file for each subdomain  C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER I,I1,J,K,L,ETYPE,ITEMP,ITEMP2,ILNODE,ILNODE2,ILNODE3
      INTEGER JD,JG,JP,IPROC,IPROC2,IPROC3,DISC,BBN,IBP
      INTEGER INDEX,INDEX2,ITOT,ITYPE,NUMS(10)
      CHARACTER LOCFN*14,PE*6
      CHARACTER*80 OUTMSG
C
      ETYPE = 3   ! The only Element-Type supported by ADCIRC is 3.
C
C--------------------------------------------------------------------------
C--MAIN LOOP:   Write a Local Grid File ( fort.14 ) for each PE
C--------------------------------------------------------------------------
C
      NETA_MAX = 0   ! max number of open boundary nodes on any subdomain

      DO 1000 IPROC = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.14'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (14,FILE=LOCFN)
C
C--------------------------------------------------------------------------
C--OPEN BOUNDARY NODES PROCESSING BEGINS HERE
C--------------------------------------------------------------------------
C
C--Partition the open boundary nodes between various processors
C
         NETAP(IPROC) = 0
         DO K=1, NOPE
            NVDLLP(K) = 0
            DO J=1, NETA
               OBNODE_LG(J,IPROC) = 0
               NBDVP(K,J) = 0
            ENDDO
         ENDDO
C
         ITOT = 0
         DO K = 1,NOPE
            DO I = 1,NVDLL(K)
               ITOT = ITOT + 1
               INDEX = NBDV(K,I)
               DO J = 1,ITOTPROC(INDEX)
                  ITEMP = (J-1)*2+1
                  IPROC2 = IMAP_NOD_GL2(ITEMP,INDEX)
                  ILNODE = IMAP_NOD_GL2(ITEMP+1,INDEX)
                  IF (IPROC.EQ.IPROC2) THEN
                    NETAP(IPROC) = NETAP(IPROC)+1
                    NVDLLP(K) = NVDLLP(K) + 1
                    NBDVP(K,NVDLLP(K)) = ILNODE
                    OBNODE_LG(NETAP(IPROC),IPROC)=ITOT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         IF (NETAP(IPROC) > NETA_MAX) NETA_MAX = NETAP(IPROC)
C
         NOPEP(IPROC) = 0
         DO K = 1,NOPE
            IF (NVDLLP(K).NE.0) THEN
              NOPEP(IPROC) = NOPEP(IPROC) + 1
            ENDIF
         ENDDO
C
C
C--------------------------------------------------------------------------
C--LAND BOUNDARY NODES PROCESSING BEGINS HERE
C--------------------------------------------------------------------------
C
C--Partition Land Boundary Segments between PEs
C
         NVELP(IPROC) = 0
         DO K = 1,NBOU
            NVELLP(K) = 0
            IBTYPEP(K,IPROC) = IBTYPE(K)
         ENDDO
C
         DO K = 1,NBOU
C
C--Weir Land Boundary Node-Pair Case
Cmod vjp 3/8/99
C  mod to allow that each of Weir-node pair might be ghosts nodes
C
         IF ((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24)) THEN
           DO I = 1,NVELL(K)
              INDEX = NBVV(K,I)
              INDEX2 = IBCONNR(K,I)
              DO J = 1,ITOTPROC(INDEX)
                 ITEMP = (J-1)*2 + 1
                 IPROC2  =  IMAP_NOD_GL2(ITEMP,INDEX)
                 ILNODE2 =  IMAP_NOD_GL2(ITEMP+1,INDEX)
                 IF (IPROC.EQ.IPROC2) THEN
                   DO JD = 1, ITOTPROC(INDEX2)
                      ITEMP2 = (JD-1)*2 + 1
                      IPROC3  = IMAP_NOD_GL2(ITEMP2,INDEX2)
                      ILNODE3 = IMAP_NOD_GL2(ITEMP2+1,INDEX2)
                      IF (IPROC.EQ.IPROC3) THEN 
                        NVELP(IPROC) = NVELP(IPROC) + 1
                        NVELLP(K) = NVELLP(K) + 1
                        LBINDEX_LG(K,NVELLP(K)) = I
                        NBVVP(K,NVELLP(K))   = ILNODE2
                        IBCONNRP(K,NVELLP(K)) = ILNODE3
                      ENDIF
                   ENDDO
                 ENDIF
              ENDDO
           ENDDO
C
C--All Other Land Boundary Node types
C
         ELSE
C
           DO I = 1,NVELL(K)
              INDEX = NBVV(K,I)
              DO J = 1,ITOTPROC(INDEX)
                 ITEMP = (J-1)*2 + 1
                 IPROC2 =  IMAP_NOD_GL2(ITEMP,INDEX)
                 ILNODE =  IMAP_NOD_GL2(ITEMP+1,INDEX)
                 IF (IPROC.EQ.IPROC2) THEN 
                   NVELP(IPROC) = NVELP(IPROC) + 1
                   NVELLP(K) = NVELLP(K) + 1
                   LBINDEX_LG(K,NVELLP(K)) = I
                   NBVVP(K,NVELLP(K)) = ILNODE
                 ENDIF
              ENDDO
           ENDDO
C
         ENDIF
C
         ENDDO
C
Cmod 05/18/2004 rl -- I don't think this next part is the correct
c  way to handle islands.  Rather, if an island is split by a domain, it 
c  should remain an island.  This will ensure that the boundary is 
c  closed.  The only error would occur in ghost node space, which is
c  not a problem since the answers are not used there anyway.

Cmod 12/18/98 vjp --this section re-written
C--If a PE has only part of a closed internal land boundary
C  modify its local IBTYPE to be an external land boundary segment
C  of the same type by decrementing its IBTYPE.
C  and remove a closing loop node if present

C
c        DO K=1, NBOU
c          IF (NVELLP(K).LT.NVELL(K)) THEN
c            IF (  (IBTYPEP(K,IPROC).EQ.1)
c    &         .OR.(IBTYPEP(K,IPROC).EQ.11)
c    &         .OR.(IBTYPEP(K,IPROC).EQ.21)) THEN
c decrement ibtype
c              IBTYPEP(K,IPROC) = IBTYPEP(K,IPROC)-1
C remove loop closing node
c              IF (NVELLP(K).GT.1.AND.
c    &           NBVVP(K,NVELLP(K)).EQ.NBVVP(K,1)) THEN
c                NVELLP(K) = NVELLP(K)-1
c              ENDIF
c            ENDIF
c          ENDIF
c        ENDDO

C If a segment contains only one node, remove the segment from the list     
C (NOTE: rl 5/18/04 I don't see how this could possibly happen, including 
C  ghost nodes)

         DO K=1, NBOU
           IF (NVELLP(K).EQ.1) NVELLP(K) = 0
         ENDDO

C    
C--Count the number of land boundary segments on PE IPROC.
C         
         NBOUP(IPROC) = 0
         DO K = 1,NBOU
            IF (NVELLP(K).NE.0) THEN
              NBOUP(IPROC) = NBOUP(IPROC) + 1
            ENDIF
         ENDDO
C
C--Count to check correctness of NVELP
C
         DISC=0  ! LB Nodes with non-zero normal discharge
         BBN=0   ! Mainland Barrier Boundary Nodes
         IBP=0   ! Internal Barrier Boundary Pairs
         ITEMP = 0
C
C     jgf46.21 Added support for IBTYPE=52.
         DO 400 K=1,NBOU
            IF (NVELLP(K).EQ.0) GOTO 400
            ITYPE = IBTYPEP(K,IPROC)
            IF ((ITYPE.EQ.2).OR.(ITYPE.EQ.12).OR.(ITYPE.EQ.22).OR.
     &           (ITYPE.EQ.52)) THEN
              DISC = DISC + NVELLP(K)
            ENDIF
            IF ((ITYPE.EQ.3).OR.(ITYPE.EQ.13).OR.(ITYPE.EQ.23)) THEN
              BBN = BBN + NVELLP(K)
            ENDIF
            IF ((ITYPE.EQ.4).OR.(ITYPE.EQ.24)) THEN
              IBP = IBP + NVELLP(K)
            ENDIF
            IF ((ITYPE.NE.2).AND.(ITYPE.NE.12).AND.(ITYPE.NE.22).AND.
     &           (ITYPE.NE.52).AND.
     &           (ITYPE.NE.3).AND.(ITYPE.NE.13).AND.(ITYPE.NE.23).AND.
     &           (ITYPE.NE.4).AND.(ITYPE.NE.24)) THEN
               ITEMP = ITEMP + NVELLP(K)
            ENDIF
            I1 = 0
            DO I=1,NVELLP(K)
               IF ((ITYPE.EQ.1).OR.(ITYPE.EQ.11).OR.
     &             (ITYPE.EQ.21)) THEN
                 IF ((I.EQ.NVELLP(K)).AND.(NBVVP(K,I).NE.I1)) THEN
                   ITEMP = ITEMP + 1
                 ENDIF
               ENDIF
               IF (I.EQ.1) I1 = NBVVP(K,I)
            ENDDO
 400     CONTINUE
C
c        print *, IPROC-1,ITEMP,DISC,BBN,2*IBP
         ITEMP  = ITEMP + DISC + BBN + 2*IBP
         IF (ITEMP.NE.NVELP(IPROC)) THEN
c          print *, "changed value from ",NVELP(IPROC)," to ",ITEMP
           NVELP(IPROC) = ITEMP
         ENDIF
         IF (NVELP(IPROC)+1.GT.MNVEL) THEN
           print *, "NVEL exceeds parameter value MNVEL on PE",IPROC
           print *, "local NVEL value = ",ITEMP
           stop
         ENDIF
C
C--Construct a LBCODE for each Land Boundary Node of this PE
C
         JP=0
         DO K = 1,NBOU
            DO I=1, NVELLP(K)
               JP = JP+1
               LBCODEP(JP,IPROC) = IBTYPEP(K,IPROC)
            ENDDO
         ENDDO
C
C--Determine whether there are any normal flow boundaries local to PE
C
         NFLUXFP(IPROC) = 0
         DO K=1, NBOU
            IF (NVELLP(K).GT.0) THEN
              IF ((IBTYPE(K).EQ.2).OR.(IBTYPE(K).EQ.12)
     &            .OR.(IBTYPE(K).EQ.22).OR.(IBTYPE(K).EQ.52)) THEN
                NFLUXFP(IPROC) = 1
              ENDIF
            ENDIF
         ENDDO

#ifdef ADCSWAN
C
C--------------------------------------------------------------------------
C--jgf48.17 UnSWAN BOUNDARY NODES PROCESSING BEGINS HERE
C--------------------------------------------------------------------------
C
C--Partition the UnSWAN boundary nodes between various processors
C
         numUnSWANBdryNodesTotalP(IPROC) = 0
         DO K=1, numUnSWANBdrySegs
            numUnSWANBdryNodesEachP(K) = 0
            DO J=1, numUnSWANBdryNodesTotal
               unswanBdryNodes_lG(J,IPROC) = 0
               unswanBdryNodesP(K,J) = 0
            ENDDO
         ENDDO
C
C
         ITOT = 0
         DO K = 1,numUnSWANBdrySegs
            DO I = 1,numUnSWANBdryNodesEach(K)
               ITOT = ITOT + 1
               INDEX = unswanBdryNodes(K,I)
               DO J = 1,ITOTPROC(INDEX)
                  ITEMP = (J-1)*2+1
                  IPROC2 = IMAP_NOD_GL2(ITEMP,INDEX)
                  ILNODE = IMAP_NOD_GL2(ITEMP+1,INDEX)
                  IF (IPROC.EQ.IPROC2) THEN
                    numUnSWANBdryNodesTotalP(IPROC)
     &                 = numUnSWANBdryNodesTotalP(IPROC)+1
                    numUnSWANBdryNodesEachP(K)
     &                 = numUnSWANBdryNodesEachP(K) + 1
                    unswanBdryNodesP(K,numUnSWANBdryNodesEachP(K)) 
     &                 = ILNODE
                    unswanBdryNodes_lG
     &                 (numUnSWANBdryNodesTotalP(IPROC),IPROC) = ITOT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
C
         numUnSWANBdrySegsP(IPROC) = 0
         DO K = 1,numUnSWANBdrySegs
            IF (numUnSWANBdryNodesEachP(K).NE.0) THEN
              numUnSWANBdrySegsP(IPROC) = numUnSWANBdrySegsP(IPROC) + 1
            ENDIF
         ENDDO
#endif

C
C--------------------------------------------------------------------------
C--BEGIN WRITING LOCAL GRID ( fort.14 ) FILE HERE        
C--------------------------------------------------------------------------
C
C--Write Mesh Data
C
          WRITE(14,80) AGRID
C
          NUMS(1) = NELP(IPROC)
          NUMS(2) = NNODP(IPROC)
c  
cjgf45.06    CALL INSERT(SIZEMSG,OUTMSG,NUMS,2)
cjgf45.06    WRITE(14,80) OUTMSG
          WRITE(14,43) NELP(IPROC),NNODP(IPROC) !jgf45.06
C
          DO J = 1,NNODP(IPROC)
             INDEX = IMAP_NOD_LG(J,IPROC)
             WRITE(14,44) J,X(INDEX),Y(INDEX),DP(INDEX)
          ENDDO
c
          DO J = 1,NELP(IPROC)
             WRITE(14,45) J,ETYPE,NNEP(1,J,IPROC),NNEP(2,J,IPROC),
     &                           NNEP(3,J,IPROC)
          ENDDO
 43       FORMAT(2I8)
 44       FORMAT(I8,3(E24.12))
 45       FORMAT(5I8)
C
C--Write Open Boundary Data
C
          CALL NEWINDEX(NOPEMSG,OUTMSG,NOPEP(IPROC))
          WRITE(14,80) OUTMSG
C
          CALL NEWINDEX(NETAMSG,OUTMSG,NETAP(IPROC))
          WRITE(14,80) OUTMSG
C
          ITOT = 0
          DO K = 1,NOPE
             IF (NVDLLP(K).GT.0)THEN
                ITOT = ITOT + 1
                CALL NEWINDEX(NVDLLMSG(K),OUTMSG,NVDLLP(K))
                WRITE(14,80) OUTMSG
                DO I = 1,NVDLLP(K)
                   WRITE(14,*) NBDVP(K,I)
                ENDDO
             ENDIF
          ENDDO
C
C--Write Land Boundary Data
C
          CALL NEWINDEX(NBOUMSG,OUTMSG,NBOUP(IPROC))
          WRITE(14,80) OUTMSG
C
          CALL NEWINDEX(NVELMSG,OUTMSG,NVELP(IPROC))
          WRITE(14,80) OUTMSG
C
          DO K = 1,NBOU
             IF(NVELLP(K).GT.0)THEN
                ITYPE = IBTYPEP(K,IPROC)
                NUMS(1) = NVELLP(K)
                NUMS(2) = ITYPE
                CALL INSERT(NVELLMSG(K),OUTMSG,NUMS,2)
                WRITE(14,80) OUTMSG
C
                IF ((ITYPE.NE.3).AND.(ITYPE.NE.13).AND.
     &             (ITYPE.NE.23).AND.(ITYPE.NE.4).AND.
     &             (ITYPE.NE.24)) THEN
                  DO I = 1,NVELLP(K)
                     WRITE(14,'(I8)') NBVVP(K,I)
                  ENDDO
                ELSEIF ((ITYPE.EQ.3).OR.(ITYPE.EQ.13).OR.
     &                  (ITYPE.EQ.23)) THEN 
                  DO I = 1,NVELLP(K)
                     INDEX = LBINDEX_LG(K,I)
                     WRITE(14,81) NBVVP(K,I),BAR1(K,INDEX),BAR2(K,INDEX)
                  ENDDO
C
                ELSEIF ((ITYPE.EQ.4).OR.(ITYPE.EQ.24)) THEN
                  DO I = 1,NVELLP(K)
                     INDEX = LBINDEX_LG(K,I)
                     WRITE(14,82) NBVVP(K,I),IBCONNRP(K,I),
     &                        BAR1(K,INDEX),BAR2(K,INDEX),BAR3(K,INDEX)
                  ENDDO
                ENDIF
             ENDIF
          ENDDO
          !jgf48.17 Added support for UnSWAN boundaries in fort.14
#ifdef ADCSWAN
C
C--Write UnSWAN Boundary Data
C
          CALL NEWINDEX(numUnSWANBdrySegsMSG,OUTMSG,
     &                  numUnSWANBdrySegsP(IPROC))
          WRITE(14,80) OUTMSG
C
          CALL NEWINDEX(numUnSWANBdryNodesTotalMSG,OUTMSG,
     &                  numUnSWANBdryNodesTotalP(IPROC))
          WRITE(14,80) OUTMSG
C
          ITOT = 0
          DO K = 1,numUnSWANBdrySegs
             IF (numUnSWANBdryNodesEachP(K).GT.0)THEN
                ITOT = ITOT + 1
                 WRITE(14,*) numUnSWANBdryNodesEachP(K),
     &                       unswanBdryMarkers(K)
C                CALL NEWINDEX(numSWANBdryNodesEachMSG(K),OUTMSG,
C     &                        numSWANBdryNodesEachP(K))
C                WRITE(14,80) OUTMSG
                DO I = 1,numUnSWANBdryNodesEachP(K)
                   WRITE(14,*) unswanBdryNodesP(K,I)
                ENDDO
             ENDIF
          ENDDO
#endif

C
          CLOSE(14)
C
1000  CONTINUE
C
C--Print Summary of Boundary Node Decomposition
C
      print *, " "
      print *, "Boundary Node Decomposition Data"
      print *, "DOMAIN      NOPE    NETA    NBOU  NVEL    NWEIR"
      WRITE(*,90)  "GLOBAL",NOPE, NETA, NBOU, NVEL, NWEIR
      DO IPROC=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,IPROC-1)
         WRITE(*,90)  PE,NOPEP(IPROC),NETAP(IPROC),
     &                NBOUP(IPROC),NVELP(IPROC),NWEIRP(IPROC)
      ENDDO
C
  80  FORMAT(A80)
  81  FORMAT(I8,E13.6,2X,E13.6)
  82  FORMAT(I8,2X,I8,2X,E13.6,2X,E13.6,2X,E13.6)
  90  FORMAT(1X,A6,5I8)
C
      RETURN
      END


      SUBROUTINE PREP15()
      USE PRE_GLOBAL 
      use memory_usage
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C 
C  This routine writes a Local Input file "fort.15" file for each subdomain C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C           Modifications by RL on 10/9/01 to accomodate NWS = -2           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,K,L,JG,JP,KK
      INTEGER INDEX,ITOT,ILNODE,IPROC,IPROC2,ITYPE,NUMS(10)
      CHARACTER LOCFN*14,PE*6
      CHARACTER*80 OUTMSG
C
C--Write a Local Input file ( fort.15 ) for each PE
C
      ! max number of stations in any subdomain
      NSTAE_MAX = 0; NSTAV_MAX = 0; NSTAM_MAX = 0; NSTAC_MAX = 0
C
      DO 1000 IPROC = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.15'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (15,FILE=LOCFN)
C
         WRITE(15,80) RUNDES
         WRITE(15,80) RUNID
         WRITE(15,80) OVERMSG
         WRITE(15,80) ABOUTMSG
         WRITE(15,80) SCREENMSG
         WRITE(15,80) HOTMSG
         WRITE(15,80) ICSMSG
         WRITE(15,80) IMMSG
         IF (CBaroclinic) THEN  !jgf46.28
            WRITE(15,80) IDENMSG  
         ENDIF
         WRITE(15,80) IBFMSG
         WRITE(15,80) IFAMSG
         WRITE(15,80) ICAMSG
         WRITE(15,80) ICATMSG
         WRITE(15,80) NWPMSG
         IF (NWP.gt.0) THEN !jgf46.00 write nodal attributes 
            DO I=1, NWP
               WRITE(15,80) NodalAttributes(I)
            ENDDO
         ENDIF
         WRITE(15,80) NCORMSG
         WRITE(15,80) NTIPMSG
C     jgfdebug46.02 Added check for NWS=45 to write NWS=5 
         IF (NWS.EQ.45) THEN
            WRITE(15,'(A1)') "5"
         ELSE
            WRITE(15,80) NWSMSG
         ENDIF
         WRITE(15,80) RAMPMSG
         WRITE(15,80) GMSG
         WRITE(15,80) TAU0MSG
C        jgf47.11 Added writing of limits for time varying tau0
         IF ( (TAU0.le.-5.d0).AND.(TAU0.gt.-6.d0) ) THEN
             WRITE(15,80) TAU0LIMMSG
         ENDIF
         WRITE(15,80) DTMSG
         WRITE(15,80) STATMSG
         WRITE(15,80) REFTMSG

         IF((NWS.EQ.0).AND.(NRS.GE.1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
         IF((NWS.EQ.1).AND.(NRS.GE.1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
C     jgfdebug46.02 Added check for NWS=45.
C     jgf46.02 Added NWS=8.
C     jgf46.16 Merged:
C     cf & cm added NWS=9: asymmetric hurricane wind model
C     sb46.28sb01 added NWS=12: OWI format
         IF ((ABS(NWS).EQ.2).OR.(ABS(NWS).EQ.4).OR.(ABS(NWS).EQ.45).OR.
     &        (ABS(NWS).EQ.5).OR.(ABS(NWS).EQ.6).OR.(ABS(NWS).EQ.8)
     &        .OR.(ABS(NWS).EQ.9).OR.(ABS(NWS).EQ.12))THEN
            WRITE(15,80) WSMSG1
         ENDIF
         IF (NWS.EQ.3) THEN
            WRITE(15,80) WSMSG1
            WRITE(15,80) WSMSG2
         ENDIF

         WRITE(15,80) RNDAYMSG
         WRITE(15,80) DRAMPMSG
         WRITE(15,80) COEFMSG
         WRITE(15,80) H0MSG
         WRITE(15,80) SLMSG        
         WRITE(15,80) TAUMSG        
         WRITE(15,80) ESLMSG        
         WRITE(15,80) CORIMSG
         WRITE(15,80) NTIFMSG
         DO I=1,NTIF
            WRITE(15,80)  TIPOTAG(I)
            WRITE(15,80)  TPKMSG(I)
         ENDDO
         
         WRITE(15,80) NBFRMSG
         DO I=1,NBFR
            WRITE(15,80) BOUNTAG(I)
            WRITE(15,80) AMIGMSG(I)
         ENDDO
         DO I=1,NBFR
            WRITE(15,80) ALPHA1(I)
            DO J=1,NETAP(IPROC)
               WRITE(15,80) EMOMSG(I,OBNODE_LG(J,IPROC))
            ENDDO 
         ENDDO   
         
         WRITE(15,80) ANGMSG
C
C--If there were any normal flow boundaries local to PE, process them
C
         IF (NFLUXFP(IPROC).EQ.1) THEN
C
            NFLBNP = 0
            DO I=1, NFLBN
               INDEX = FLBN(I)
               DO J=1, ITOTPROC(INDEX)
                  IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDEX)
                  IF (IPROC.EQ.IPROC2) THEN
                     NFLBNP = NFLBNP + 1
                     FLBNXP(NFLBNP) = FLBNX(I)
                  ENDIF
               ENDDO
            ENDDO
C     
            WRITE(15,80) NFFRMSG
            IF (NFFR.NE.0) THEN
               DO I=1,NFFR
                  WRITE(15,80) FBOUNTAG(I)
                  WRITE(15,80) FREQMSG(I)
               ENDDO
               DO I=1,NFFR
                  WRITE(15,80) ALPHA2(I)
                  DO J=1,NFLBNP
                     WRITE(15,80) QNMSG(I,FLBNXP(J))
cdebug               print *, "PE=",IPROC," FLUXNODE=",FLBNXP(J)
                  ENDDO
               ENDDO
            ENDIF
C     
         ENDIF
C
C--Write Local Elevation Station Info:
C--Create Local-to-Global element "ownership" of an elevation station
C
C     WRITE(15,80) STAEMSG !jgf45.07 we may have changed NOUTE in adcprep
         WRITE(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
C
         NSTAEP(IPROC) = 0
         DO K = 1,NSTAE
            DO J=1,NELP(IPROC)
               INDEX = abs(IMAP_EL_LG(J,IPROC))
               IF (INDEX.EQ.NNSEG(K)) THEN
                  NSTAEP(IPROC) = NSTAEP(IPROC) + 1
                  KK = K
                  if (STAE_SHARE(K) > -1) KK = -K
                  IMAP_STAE_LG(NSTAEP(IPROC),IPROC) = KK
                  STAE_SHARE(K) = IPROC
               ENDIF
            ENDDO  
         ENDDO        
         NSTAE_MAX = MAX(NSTAEP(IPROC),NSTAE_MAX) 
C     
         CALL INSERT(NSTAEMSG,OUTMSG,NSTAEP(IPROC),1)
         WRITE(15,80) OUTMSG       
C     
         DO K=1,NSTAEP(IPROC)
            INDEX = abs(IMAP_STAE_LG(K,IPROC))
            WRITE(15,80) STAELOC(INDEX)
         ENDDO             
C
C--Write Local Velocity Station Info:
C--Create Local-to-Global element "ownership" of an velocity station
C
         WRITE(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
C
         NSTAVP(IPROC) = 0
         DO K = 1,NSTAV
            DO J=1,NELP(IPROC)
               INDEX = abs(IMAP_EL_LG(J,IPROC))
               IF (INDEX.EQ.NNSVG(K)) THEN
                  NSTAVP(IPROC) = NSTAVP(IPROC) + 1
                  KK = K
                  if (STAV_SHARE(K) > -1) KK = -K
                  IMAP_STAV_LG(NSTAVP(IPROC),IPROC) = KK
                  STAV_SHARE(K) = IPROC
               ENDIF
            ENDDO 
         ENDDO   
         NSTAV_MAX = MAX(NSTAVP(IPROC),NSTAV_MAX) 
C
         CALL INSERT(NSTAVMSG,OUTMSG,NSTAVP(IPROC),1)
         WRITE(15,80) OUTMSG       
C
         DO K=1,NSTAVP(IPROC)
            INDEX = abs(IMAP_STAV_LG(K,IPROC))
            WRITE(15,80) STAVLOC(INDEX)
         ENDDO
C
C--If IM=10 Write Concentration Station Info:
C--Create Local-to-Global element "ownership" of an concentration station
C
         NSTACP(IPROC) = 0
         IF (C2D_PTrans.or.C3D_PTrans) THEN !jgf46.28
C
C     WRITE(15,80) STACMSG   !jgf45.07 we may have changed NOUTC in adcprep
            WRITE(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
C
            DO K = 1,NSTAC
               DO J=1,NELP(IPROC)
                  INDEX = abs(IMAP_EL_LG(J,IPROC))
                  IF (INDEX.EQ.NNSCG(K)) THEN
                     NSTACP(IPROC) = NSTACP(IPROC) + 1
                     KK = K
                     if (STAC_SHARE(K) > -1) KK = -K
                     IMAP_STAC_LG(NSTACP(IPROC),IPROC) = KK
                     STAC_SHARE(K) = IPROC
                  ENDIF
               ENDDO 
            ENDDO   
            NSTAC_MAX = MAX(NSTACP(IPROC),NSTAC_MAX) 
C
            CALL INSERT(NSTACMSG,OUTMSG,NSTACP(IPROC),1)
            WRITE(15,80) OUTMSG       
C     
            DO K=1,NSTACP(IPROC)
               INDEX = abs(IMAP_STAC_LG(K,IPROC))
               WRITE(15,80) STACLOC(INDEX)
            ENDDO
C     
         ENDIF
C
C--Write Local Meterological Station Info:
C--Create Local-to-Global element "ownership" of an elevation station
C
         NSTAMP(IPROC) = 0
         IF (NWS.NE.0) THEN
            WRITE(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
            DO K = 1,NSTAM
               DO J=1,NELP(IPROC)
                  INDEX = abs(IMAP_EL_LG(J,IPROC))
                  IF (INDEX.EQ.NNSMG(K)) THEN
                     NSTAMP(IPROC) = NSTAMP(IPROC) + 1
                     KK = K
                     if (STAM_SHARE(K) > -1) KK = -K
                     IMAP_STAM_LG(NSTAMP(IPROC),IPROC) = KK
                     STAM_SHARE(K) = IPROC
                  ENDIF
               ENDDO  
            ENDDO        
            NSTAM_MAX = MAX(NSTAMP(IPROC),NSTAM_MAX) 
C
            CALL INSERT(NSTAMMSG,OUTMSG,NSTAMP(IPROC),1)
            WRITE(15,80) OUTMSG       
C
            DO K=1,NSTAMP(IPROC)
               INDEX = abs(IMAP_STAM_LG(K,IPROC))
               WRITE(15,80) STAMLOC(INDEX)
            ENDDO             
         ENDIF

C
C--Write Local Elevation Data Output Info
C
C     WRITE(15,80) OUTGEMSG !jgf45.07 we may have changed NOUTGE in adcprep
      WRITE(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
C
C--Write Local Velocity Data Output Info
C
C     WRITE(15,80) OUTGVMSG !jgf45.07 we may have changed NOUTGV in adcprep
      WRITE(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
C
C     jgf45.07 write subdomain concentration data output info if necessary
C
      IF (IM.EQ.10) WRITE(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
C
C--Write Local Wind Velocity Data Output Info ( added 4/16/98 vjp )
C
C     jgf45.07 we may have changed NOUTGW in adcprep
C     IF (NWS.NE.0) WRITE(15,80) OUTGWMSG
      IF (NWS.NE.0) WRITE(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
C
C--Write Harmonic Analysis Data
C
         WRITE(15,80) HARFRMSG
         DO I=1,NHARFR
            WRITE(15,80) HAFNAM(I)
            WRITE(15,80) HAFREMSG(I)
c           WRITE(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         ENDDO
C
         WRITE(15,80) HARPARMSG
         WRITE(15,80) OUTHARMSG
C
C--Write Hot Start Info
C
         WRITE(15,80) HSTARMSG
C
C--Write Solver Info
C
         WRITE(15,80) SOLVMSG
C
C        jgf48.03 Write netCDF metadata, if necessary
C        TODO: Extend this to other file types that are netCDF-enabled
         IF ( (ABS(NOUTE).EQ.3) .OR. (ABS(NOUTV).EQ.3) .OR.
     &        (ABS(NOUTC).EQ.3) .OR. (ABS(NOUTM).EQ.3) .OR.
     &        (ABS(NOUTGE).EQ.3) .OR. (ABS(NOUTGV).EQ.3) .OR.
     &        (ABS(NOUTGC).EQ.3) .OR. (ABS(NOUTGW).EQ.3). OR.
     &        (ABS(NHSTAR).EQ.3))THEN                                   ! for netCDF hotstart option
            WRITE(15,80) title
            WRITE(15,80) institution
            WRITE(15,80) source
            WRITE(15,80) history
            WRITE(15,80) references
            WRITE(15,80) comments
            WRITE(15,80) host
            WRITE(15,80) convention
            WRITE(15,80) contact
            WRITE(15,80) base_date
         ENDIF
C
C--Write 3DVS Info
C
         IF(C3DVS) THEN
            CALL PREP15_3DVS(IPROC)
c        ELSEIF(C3DDSS) THEN
c           CALL PREP15_3DDSS(IPROC)
         ENDIF
C
         CLOSE(15)
C
 1000 CONTINUE
C     jgf45.12 Release memory that may have been allocated for vertical
C     grid specification.
      IF(C3DVS.and.(IGC.EQ.0)) THEN
         DEALLOCATE ( Sigma )
         nbytes = 8*nfen
         call memory_dealloc(nbytes)
      ENDIF
      IF(C3DVS.and.(IEVC.EQ.0)) THEN
         DEALLOCATE ( EVTot )
         nbytes = 8*nfen
         call memory_dealloc(nbytes)
      ENDIF

C
C--Print Summary of Stations
C
      print *, " "
      print *, "Station Data"
      print *, "DOMAIN      NSTAE   NSTAV    NSTAC    NSTAM"
      WRITE(*,92)  "GLOBAL",NSTAE,NSTAV,NSTAC,NSTAM
      DO IPROC=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,IPROC-1)
         WRITE(*,92)  PE,NSTAEP(IPROC),NSTAVP(IPROC),
     .        NSTACP(IPROC),NSTAMP(IPROC)
      ENDDO
C     
      RETURN
 80   FORMAT(A80)
 92   FORMAT(1X,A6,4I8)
      END

C
C---------------------------------------------------------------------------C
C                     (  Serial Version  6/24/02  )                         C 
C  This routine writes the 3DVS info in the Local Input file "fort.15" file C
C  for each subdomain using the domain decomposition of the ADCIRC grid     C
C  created by the routine DECOMP.                                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 41.11a                    C
C
C     jgf45.11 Updated to handle new format of 3D input files with stations
C     defined by coordinates rather than node numbers. 
C---------------------------------------------------------------------------C
      SUBROUTINE PREP15_3DVS(IPROC)
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER N          ! vertical grid layer counter
      INTEGER IPROC      ! subdomain counter
      INTEGER SDStation  ! subdomain station 
      INTEGER FDStation  ! full domain station 
      INTEGER SDEle      ! subdomain station element 
      INTEGER FDEle      ! full domain element 

c     jgf45.10 removed IDIAG
      WRITE(15,80) IDENMSG
      WRITE(15,80) SLIPMSG
      WRITE(15,80) Z0MSG
      WRITE(15,80) ALPMSG
      WRITE(15,80) FEMSG
C     jgf45.12 Added code to record thicknesses of vertical grid layers,
C     if necessary.
      IF(IGC.EQ.0) THEN
         DO N=1,NFEN
            WRITE(15,*) Sigma(N)
         ENDDO
      ENDIF
      WRITE(15,80) EVCMSG
C     jgf45.12 Add code to record vertical eddy viscosity profile.
      IF(IEVC.EQ.0) THEN
         DO N=1,NFEN
            WRITE(15,*) EVTot(N)
         ENDDO
      ENDIF
      IF(IEVC.EQ.50) WRITE(15,80) THETAMSG
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D density station
C     elements to corresponding elements in subdomains. Write out
C     subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DSD,TO3DSDS,TO3DSDF,NSPO3DSD
      IF(I3DSD.GT.0) THEN
         NNSTA3DDP(IPROC) = 0
         DO FDStation = 1, NSTA3DD
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DDG(FDStation) ) THEN
                  NNSTA3DDP(IPROC) = NNSTA3DDP(IPROC) + 1
                  IMAP_STA3DD_LG(NNSTA3DDP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO   
         WRITE(15,*) NNSTA3DDP(IPROC)
         DO SDStation = 1, NNSTA3DDP(IPROC)
            FDStation = IMAP_STA3DD_LG(SDStation,IPROC)
            WRITE(15,80) STA3DDLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DDMSG
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D velocity station
C     elements to corresponding elements in subdomains. Write out
C     velocity subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DSV,TO3DSVS,TO3DSVF,NSPO3DSV
      IF(I3DSV.GT.0) THEN
         NNSTA3DVP(IPROC) = 0
         DO FDStation = 1, NSTA3DV
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DVG(FDStation) ) THEN
                  NNSTA3DVP(IPROC) = NNSTA3DVP(IPROC) + 1
                  IMAP_STA3DV_LG(NNSTA3DVP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO   
         WRITE(15,*) NNSTA3DVP(IPROC)
         DO SDStation = 1, NNSTA3DVP(IPROC)
            FDStation = IMAP_STA3DV_LG(SDStation,IPROC)
            WRITE(15,80) STA3DVLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DVMSG
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D turbulence station
C     elements to corresponding elements in subdomains. Write out
C     turbulence subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DST,TO3DSTS,TO3DSTF,NSPO3DST
      IF(I3DST.GT.0) THEN
         NNSTA3DTP(IPROC) = 0
         DO FDStation = 1, NSTA3DT
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DTG(FDStation) ) THEN
                  NNSTA3DTP(IPROC) = NNSTA3DTP(IPROC) + 1
                  IMAP_STA3DT_LG(NNSTA3DTP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO   
         WRITE(15,*) NNSTA3DTP(IPROC)
         DO SDStation = 1, NNSTA3DTP(IPROC)
            FDStation = IMAP_STA3DT_LG(SDStation,IPROC)
            WRITE(15,80) STA3DTLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DTMSG
      ENDIF

      WRITE(15,80) DGDMSG
      WRITE(15,80) DGVMSG
      WRITE(15,80) DGTMSG
C
C     jgf45.12: Write out the parameters for the transport equation, if
C     necessary.
      IF (C3D_BTrans) THEN
C     Lateral and vertical diffusion coefficients.
         WRITE(15,*) NLSD, NVSD
         WRITE(15,*) NLTD, NVTD
C     Time stepping coefficient for the transport equation terms.
         WRITE(15,*) ALP4
C     Temperature boundary condition file type, if necessary
         IF ( IDEN .eq. 3 .or. IDEN .eq. 4 ) THEN
            WRITE(15,*) NTF
         ENDIF
      ENDIF

      RETURN
  80  FORMAT(A80)
  81  FORMAT(I8,2E15.8,2I8,A32)
  82  FORMAT(500I8)
      END
C-----------------------------------------------------------------------
C     End of subroutine PREP15_3DVS
C-----------------------------------------------------------------------


      SUBROUTINE PREP18()   
      USE PRE_GLOBAL 
      use memory_usage
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C 
C  This Routine writes a message-passing file "fort.18" for each subdomain  C
C  of the domain decomposition created by DECOMP.                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER N1, N2, N3, KMIN
      INTEGER I,J,K,L,ITEMP,IPR,IPR1
      INTEGER INDEX,ITOT,IEL,IELG,ILNODE,IPROC,ITYPE
      INTEGER,ALLOCATABLE :: RES_NODE(:)
      CHARACTER LOCFN*14,PE*6
C
C  Allocate local arrays
C
      ALLOCATE ( RES_NODE(MNPP) )
      nbytes = 4*mnpp
      call memory_alloc(nbytes)
C
C--Write Message-Passing File for each PE
C
      DO 1000 I = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.18'
         CALL IWRITE(LOCFN,3,6,I-1)
         OPEN (18,FILE=LOCFN)
 
         write(18,3050) FileFmtVersion, 0, 0
 
!vjp 9/17/06
C--Write the Global indexes of all local elements in local element order

         WRITE(18,3000) NELG, MNEP, NELP(I)  ! number of Global elements
         DO J = 1,NELP(I)
            INDEX = IMAP_EL_LG(J,I)
            WRITE(18,'(I8)') INDEX        ! Global index of local element
         ENDDO

C--Write the Global indexes of all local nodes in local node order
C  write global index as positive if a resident node and negative 
C  if a ghost node 

         WRITE(18,3001) NNODG, MNPP, NNODP(I)   ! number of Global nodes
         ITOT = 0
         DO J = 1,NNODP(I)
            INDEX = IMAP_NOD_LG(J,I)
            IPR = IMAP_NOD_GL(1,INDEX)
            IF (IPR.EQ.I)THEN
               ITOT = ITOT + 1
               RES_NODE(ITOT) = J
               WRITE(18,'(I8)') INDEX        ! Global index of resident node
            ELSE
               WRITE(18,'(I8)') -1*INDEX     ! Global index of ghost node
            ENDIF
         ENDDO
         IF (ITOT.NE.NOD_RES_TOT(I)) STOP 'ERROR IN # OF RES. NODES'

C--Write local normal flow boundary flag               
C--vjp This info is used only for relocalizing fort.15
         WRITE(18,3002) NFLUXFP(I) ! normal flow b.c. flag for subdomain

C--Write global and local total number of elevation boundary nodes
C--vjp This info is used only for relocalizing fort.15
         WRITE(18,3003) NETA, NETA_MAX, NETAP(I) ! number of global elevation b.c. nodes
         DO J = 1,NETAP(I)
            INDEX = OBNODE_LG(J,I)
            WRITE(18,'(I8)') INDEX           ! Global open boundary node index
         ENDDO

C--Write the Global indexes of all Elevation Stations in local node order
C  write global index as positive if a resident node and negative 
C  if a ghost node 

         WRITE(18,3004) NSTAE, NSTAE_MAX, NSTAEP(I) ! number of Global Elevation Stations
         DO J = 1,NSTAEP(I)
            INDEX = IMAP_STAE_LG(J,I)
            WRITE(18,'(I8)') INDEX           ! Global station number
         ENDDO

C--Write the Global indexes of all Velocity Stations in local node order
C  write global index as positive if a resident node and negative 
C  if a ghost node 

         WRITE(18,3005) NSTAV, NSTAV_MAX, NSTAVP(I) ! number of Global Velocity Stations
         DO J = 1,NSTAVP(I)
            INDEX = IMAP_STAV_LG(J,I)
            WRITE(18,'(I8)') INDEX           ! Global station number
         ENDDO

C--Write the Global indexes of all Elevation Stations in local node order
C  write global index as positive if a resident node and negative 
C  if a ghost node 

         WRITE(18,3006) NSTAM, NSTAM_MAX, NSTAMP(I) ! number of Global Meteorlogical Stations
         DO J = 1,NSTAMP(I)
            INDEX = IMAP_STAM_LG(J,I)
            WRITE(18,'(I8)') INDEX           ! Global station number
         ENDDO

C--Write the Global indexes of all Concentration Stations in local node order
C  write global index as positive if a resident node and negative 
C  if a ghost node 

         WRITE(18,3007) NSTAC, NSTAC_MAX, NSTACP(I) ! number of Global Concentration Stations
         DO J = 1,NSTACP(I)
            INDEX = IMAP_STAC_LG(J,I)
            WRITE(18,'(I8)') INDEX           ! Global station number
         ENDDO

C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C
C--Write the Resident Node List
C
         WRITE(18,3010) (I-1),NOD_RES_TOT(I)
         WRITE(18,1130) (RES_NODE(J),J=1,ITOT)
C
C--Write the Number of Communicating PEs
C
         WRITE(18,3020) NUM_COMM_PE(I)
C
C--Write the Receive List
C
         DO J = 1,NUM_COMM_PE(I)
            IPR = COMM_PE_NUM(J,I)
            IRECV_TOT(J,I) = 0
            DO K = 1,NNODP(I)
               INDEX = IMAP_NOD_LG(K,I)
               IF (IMAP_NOD_GL(1,INDEX).EQ.IPR) THEN
                  IRECV_TOT(J,I) = IRECV_TOT(J,I) + 1
                  IRECV(IRECV_TOT(J,I)) = K
c uncomment next line and comment preceding line for debugging
c                 IRECV(IRECV_TOT(J,I)) = INDEX
               ENDIF
            ENDDO
            WRITE(18,3030) (IPR-1), IRECV_TOT(J,I)
            WRITE(18,1130) (IRECV(K),K=1,IRECV_TOT(J,I))
         ENDDO
C              
C--write the send list
C
         DO J = 1,NUM_COMM_PE(I)
            IPR = COMM_PE_NUM(J,I)
            ISEND_TOT(J,I) = 0
            DO K = 1,NNODP(IPR)
               INDEX = IMAP_NOD_LG(K,IPR)
               IF (IMAP_NOD_GL(1,INDEX).EQ.I) THEN
                  ISEND_TOT(J,I) = ISEND_TOT(J,I) + 1
                  ISEND(ISEND_TOT(J,I)) = IMAP_NOD_GL(2,INDEX)
c uncomment next line and comment preceding line for debugging
c                 ISEND(ISEND_TOT(J,I)) = INDEX
               ENDIF
            ENDDO
            WRITE(18,3040)  IPR-1, ISEND_TOT(J,I)
            WRITE(18,1130) (ISEND(K),K=1,ISEND_TOT(J,I))
         ENDDO
         
         CLOSE(18)
C     
1000  CONTINUE
C
C--Compute the surface to volume ratio (in %)
C
      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NUM_COMM_PE(I)
            ITOT = ITOT + IRECV_TOT(J,I)
         ENDDO
         PROC_SV(I) = (ITOT/REAL(NOD_RES_TOT(I)))*100.0
c        WRITE(6,*) I-1,PROC_SV(I)
      ENDDO
C
      print *, " "
      print *, "Communication Data"
      print *, "DOMAIN  COMM_PE  %(SURF/VOL)"
      print *, "------  -------  -----------"
      DO I=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,I-1)
         WRITE(6,92) PE, NUM_COMM_PE(I),PROC_SV(I)
      ENDDO

      deallocate( res_node )
      nbytes = 4*mnpp
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
C
  92  FORMAT(1X,A6,2X,I7,2X,F8.2)
1130  FORMAT(8X,9I8)
3000  FORMAT('NELG    ',3I8)
3001  FORMAT('NNODG   ',3I8)
3002  FORMAT('NFLUXF  ',I8)
3003  FORMAT('NETA    ',3I8)
3004  FORMAT('NSTAE   ',3I8)
3005  FORMAT('NSTAV   ',3I8)
3006  FORMAT('NSTAM   ',3I8)
3007  FORMAT('NSTAC   ',3I8)
3010  FORMAT('RES NODE',2I8)    
3020  FORMAT('COMM PE ',2I8)    
3030  FORMAT('RECV PE ',2I8)    
3040  FORMAT('SEND PE ',2I8) 
3050  FORMAT('FileFmt ',3I8)
      END


      SUBROUTINE PREP19()
      USE PRE_GLOBAL 
      use memory_usage
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C 
C  This routine writes a Local "Aperiodic Elevation Boundary Condtions"     C
C  (fort.19) file for each subdomain using the domain decomposition of      C
C  the ADCIRC grid created by the routine DECOMP.                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C
C     jgf45.12 Added subroutine call to open files.
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC,ESBINP
      CHARACTER*40,ALLOCATABLE :: ESBIN(:)

C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 19 
C  Global Aperiodic Elevation Boundary Conditions file
C
C     Open full domain and subdomain fort.19 files
      CALL OpenPrepFiles(19, 'aperiodic elevation boundary  ', 
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 19 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C--Allocate local arrays
C
      ALLOCATE ( ESBIN(MNETA) )
      nbytes = 8*mneta
      call memory_alloc(nbytes)
C
      READ(19,40) ETIMINC
      DO IPROC = 1,NPROC
         WRITE(SDU(IPROC),40)  ETIMINC
      ENDDO
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
1000  CONTINUE
      DO I=1, NETA
         READ(19,40,END=9999)  ESBIN(I)
      ENDDO
C
      DO IPROC= 1,NPROC
         DO I=1, NETAP(IPROC)
            ESBINP = ESBIN(OBNODE_LG(I,IPROC))
            WRITE(SDU(IPROC),40) ESBINP
         ENDDO
      ENDDO
C
      GO TO 1000
C
C--Close Global file and all the Local Files     
C
 9999 CLOSE (19)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      deallocate( esbin )
      nbytes = 8*mneta
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
  40  FORMAT(A40)
      END


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 2 0
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will break up the full domain aperiodic
C     flux boundaries into subdomains using the domain decomposition of
C     the ADCIRC grid created by the routine DECOMP.
C
C     -Written by MEB 04/01/04
C     -Added by jgf to 45.06 10/07/2005
C     -jgf45.12 Rewritten to correct bugs in subdomain fort.20
C     formatting as well as the erroneous use of the GL mapping instead
C     of GL2. Also added subroutine call to open files.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP20()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER IPROC, FLUX_INC
      INTEGER INDEX14, I
      REAL(SZ)  FLUX_VAL
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER INDEX ! full domain node number for a flow boundary node
      INTEGER J     ! counter for subdomains that corrsp. to a single f.d. node
      INTEGER IPROC2! PE of a subdomain that matches a single full domain node

C     Open full domain and subdomain fort.20 files
      CALL OpenPrepFiles(20, 'aperiodic flux boundary       ', 
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 20 files not preprocessed.'
         RETURN ! note early return
      ENDIF
c
c     Write Increment into all flux files
c 
      READ(20,*) FLUX_INC
      DO IPROC=1,NPROC
         WRITE(SDU(IPROC),*) FLUX_INC
      ENDDO
C
C     jgf45.12 Write each full domain nodal flux value into each of the
C     subdomains that that full domain node maps to. The full domain
C     node may map to more than one subdomain node if it falls on a
C     boundary between subdomains (ghost nodes).
C
 33   DO I=1, EXIST_FLUX            ! loop through full domain flow nodes
         READ(20,*,END=40) FLUX_VAL ! get a flo val for this f.d. flow node
         INDEX = FLUX14_ARY(I)      ! get full domain flow node number 
         DO J=1, ITOTPROC(INDEX)    ! loop over subdomains for 1 f.d. node
            IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDEX) ! find next subdomain
            DO IPROC=1, NPROC
               IF (IPROC.EQ.IPROC2) THEN ! full domain node maps to this s.d. 
                  WRITE(SDU(IPROC),50) FLUX_VAL
               ENDIF
            END DO
         END DO
      END DO
      GOTO 33
 40   CLOSE (20)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      IF (allocated(FLUX14_ARY)) then
        DEALLOCATE (FLUX14_ARY)
        nbytes = 4*exist_flux
        call memory_dealloc(nbytes)
      ENDIF
      call memory_status()
      return
 50   FORMAT (F16.8,1x,I6,1x,I6,1x,I6)
c----------------------------------------------------------------------------
      END SUBROUTINE PREP20
c----------------------------------------------------------------------------




C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 2 2
C---------------------------------------------------------------------------
C
C                     (  Serial Version  2/28/98  )                         C 
C  This routine reads a global external meteorology file when NWS=1,+-2,3,  C
C  +-4,+-5,+-6.  In each case it wites a local meteorology file of the same C
C  format for each subdomain using the domain decomposition of the ADCIRC   C
C  grid created by the routine DECOMP.                                      C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C 
C  This version is compatible with ADCIRC version 43.03                     C
C
C     jgf46.02 Added subroutine call to open prep files; this provides
C     the user with the ability to skip the prepping of wind data files.
C
C     jgfdebug46.02 Added NWS=45 to imitate the behavior of the v42 (IPET)
C     code.
C
C     jgf46.02 Added NWS=8 to copy the wind files for the Holland model 
C     into the subdomains.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP22()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND,DONE
      INTEGER I,J,IPROC,IPROC2,ILNODE,INDEX,NHG,LINDEX
      CHARACTER*80 PBLJAGF
C      CHARACTER FNAME*60,LOCFN*14,CMD1*63,CMD2*7,CMD*70,INLINE*80
      CHARACTER FNAME*60,CMD1*63,CMD2*7,CMD*70
      CHARACTER*170 Line ! line of data from NWS=8 (Holland) file
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER,ALLOCATABLE  :: NG(:)
      REAL(SZ),ALLOCATABLE :: WVNXG(:),WVNYG(:),PRG(:)
      REAL(SZ),ALLOCATABLE :: WVNXL(:),WVNYL(:),PRL(:)
      REAL(SZ) U,V,PR
      REAL(SZ) RHOWAT     !jgfdebug46.02  
C     jgf48.47 Do the decomposition for a max of 256 subdomains at a
C     time ... some platforms/compilers limit the number of files that 
C     can be open at any one time.
      INTEGER, PARAMETER :: maxOpenFiles = 256 
      INTEGER startProc  
      INTEGER endProc 
      INTEGER deltaProc
C
C     Allocate local work arrays
C
      ALLOCATE ( NG(MNWP) )
      nbytes = 4*mnwp
      call memory_alloc(nbytes)
      ALLOCATE ( WVNXG(MNWP),WVNYG(MNWP),PRG(MNWP) )
      nbytes = 24*mnwp
      call memory_alloc(nbytes)
      ALLOCATE ( WVNXL(MNWP),WVNYL(MNWP),PRL(MNWP) )
      nbytes = 24*mnwp
      call memory_alloc(nbytes)
C
C     Perform decomposition over a range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc

C        Open full domain and all subdomain fort.22 files
         CALL OpenPrepFiles(22, 'wind information              ', 
     &      startProc, endProc, sdu, success)

         IF (.not.success) THEN
           WRITE(*,*) 'WARNING: Unit 22 files not preprocessed.'
           RETURN ! note early return
         ENDIF
C
C--Branch to Appropriate Code
C
         SELECT CASE(ABS(NWS))
C        -------------
         CASE(1,2,5,7)
C        -------------
C
C     MAIN LOOP FOR NWS = 1, +-2,+-5,+-7
C     (1)  Read a record from Global Wind Stress File
C     (2)  Use Decomp arrarys to Localize record to a subdomain
C     (3)  Write Local Wind Stress record in same format 
           DO                     ! loop forever (or until file ends)
              READ(22,*,END=9999)  
     &           (NG(I),WVNXG(I),WVNYG(I),PRG(I),I=1,NNODG)
              DO IPROC = STARTPROC, ENDPROC
                 DO I=1, NNODP(IPROC)
                    INDEX = IMAP_NOD_LG(I,IPROC)
                    WVNXL(I) = WVNXG(INDEX)
                    WVNYL(I) = WVNYG(INDEX)
                    PRL(I) = PRG(INDEX)
                 ENDDO
                 DO I=1, NNODP(IPROC)
                    WRITE(SDU(IPROC),1100)  I,WVNXL(I),WVNYL(I),PRL(I)
                 ENDDO
              ENDDO
           ENDDO
C
C        -------     
         CASE(4)
C        -------
C        MAIN LOOP FOR NWS = +- 4  ( PBL Format )
C        (1)  Read a record from Global Wind Stress File
C        (2)  Use Decomp arrarys to Localize record to a subdomain
C        (3)  Write out in PBL Format on subdomain
C
C--Read a wind field record from the global input file
C
         DO
            READ(22,'(A80)',END=9999) PBLJAGF
            IF(PBLJAGF(2:2).EQ.'#') THEN
               DO IPROC =  STARTPROC,ENDPROC
                  WRITE(SDU(IPROC),1101)
                  WRITE(SDU(IPROC),1100) 1,0.0,0.0,0.0 !victor didn't like this line 27/11/03
               ENDDO
            ELSE
!     vjp 27/11/03 
!     rewrote this section to handle ghost-nodes 
C              READ(PBLJAGF,'(I8,3E13.5)',END=9999) NHG,U,V,PR
               READ(PBLJAGF,*,END=9999) NHG,U,V,PR
               DO J=1, ITOTPROC(NHG)
                  IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
                  LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
                  WRITE(SDU(IPROC),1100) LINDEX,U,V,PR
               ENDDO
            ENDIF
         END DO
C
C        --------     
         CASE(45)
C        --------
C        jgf46.02 Convert NWS=4 winds to NWS=5 winds to imitate the Katrina
C        (IPET) version of the code.
C
C-- Read a wind field record from the global input file
C
         DO
            RHOWAT=1000.0d0
            CALL NWS4GET(WVNXG,WVNYG,PRG,G,RHOWAT,NNODG,DONE)
C     
            DO IPROC = STARTPROC,ENDPROC
               DO I=1, NNODP(IPROC)
                  INDEX = IMAP_NOD_LG(I,IPROC)
                  WVNXL(I) = WVNXG(INDEX)
                  WVNYL(I) = WVNYG(INDEX)
                  PRL(I) = PRG(INDEX)
               ENDDO
               DO I=1, NNODP(IPROC)
                  WRITE(SDU(IPROC),1100)  I,WVNXL(I),WVNYL(I),PRL(I)
               ENDDO
            ENDDO
C--   If reached EOF in NWS4GET last time go close files and return
C     
            IF (DONE) GOTO 9999
         ENDDO
C     
C        ---------     
         CASE(3,6,8,9)
C        ---------
C     jgf46.16 Added NWS=9, asymmetric hurricane wind model
C
C     jgf48.35 Moved NWS=3 and NWS=6 into this case, since it seems to 
C     work better than calling the shell to copy the files. Although
C     this approach seems to leave a lot of trailing spaces on each
C     fort.22 line, at least on the IBM platform.
C
C     Holland Wind Model (8) or asymetric hurricane wind model (9)
C
C     Hurricane track data are in a small file that describes the
C     hurricane characteristics and provides parameters for the wind
C     model. The actual wind velocity at each node will be calculated at
C     run time from these characteristics and parameters by ADCIRC
C     internally, the wind model. Therefore, each subdomain can use the
C     exact same information, and we can just rewrite the file into each
C     subdirectory. We will read it in and write it out to each
C     subdirectory rather than cp or ln it, because I (jgf) happen to
C     prefer to avoid calls to the platform-dependent shell from
C     Fortran.
         DO 
            READ(22,170,END=9999) Line
            DO IPROC = STARTPROC,ENDPROC
               WRITE(SDU(IPROC),170) Line
            ENDDO
         ENDDO
C
C        ------------
         CASE DEFAULT
C        ------------
         print *, "NWS=",NWS," has incorrect value in PREP22"
         RETURN

         END SELECT

C
C--Close Global file and all the Local Files     
C
 9999    CLOSE (22)
         DO IPROC=STARTPROC, ENDPROC
            CLOSE (SDU(IPROC))
         ENDDO
         startProc=endProc+1
      ENDDO
      

      DEALLOCATE ( NG,  WVNXG, WVNYG, PRG )
      DEALLOCATE ( WVNXL, WVNYL, PRL )
      nbytes = 52*mnwp
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
  60  FORMAT(A60)
 170  FORMAT(A170)
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1100 FORMAT(I8,3E13.5)
 1101 FORMAT(' #')
c----------------------------------------------------------------------------
      END SUBROUTINE PREP22
c----------------------------------------------------------------------------
      
      

      SUBROUTINE PREP23()
      USE PRE_GLOBAL
C
C---------------------------------------------------------------------------C
C                           (  add MEB 03/04/03  )                          C
C  This routine writes a Local Input file "fort.23" file for each subdomain C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER IPROC, NHG, J, LINDEX
      CHARACTER*80 PBLJAGF
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      REAL(SZ)                U,V
C
C--Open Global Wave Stress File ( UNIT 23 )
C 
C     Open full domain and subdomain fort.23 files
      CALL OpenPrepFiles(23, 'wave stress                   ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 23 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C--------------------------------------------------------------------------
C--MAIN LOOP
C   (1)  Read a record from Global Wave Stress File
C   (2)  Use Decomp arrays to Localize record to a subdomain
C   (3)  Write Local Wave Stress record in standard PBL format
C--------------------------------------------------------------------------
C 
C--Read a wave field record from the global input file                         
C--and write out to respective local fort.23 file.                                    
C
 170  READ(23,'(A80)',END=9999) PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') THEN
         DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),1101)
            WRITE(SDU(IPROC),1100) 1,0.0,0.0 !victor didn't like this line 27/11/03
         ENDDO
      ELSE
! vjp 27/11/03 
! rewrote this section to handle ghost-nodes 
! and changed if test from "and" to "or"
         READ(PBLJAGF,'(I8,2E13.5)',END=9999) NHG,U,V
         IF ((U.NE.0.).OR.(V.NE.0.)) THEN
            DO J=1, ITOTPROC(NHG)
               IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
               LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
               WRITE(SDU(IPROC),1100) LINDEX,U,V
            ENDDO
         ENDIF
      ENDIF
      
      GOTO 170

 9999 CLOSE(23)
      DO IPROC=1,NPROC
         CLOSE(SDU(IPROC))
      ENDDO
      
 1100 FORMAT(I8,2E13.5)
 1101 FORMAT (' #')
      
 99   RETURN
      END

      SUBROUTINE HOTLOCALIZE()
      USE PRE_GLOBAL 
      use presizes; use memory_usage
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C
C             started mods for harmonic analysis and 3D RL 5/22/03          C
C         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
C                                                                           C
C  This routine reads the global hot start file (either fort.67 or fort.68) C
C  and writes local hot start files of the same format.                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
#ifdef NETCDF
      USE GLOBAL, ONLY:  nciestp,nccoue,ncivstp,nccouv,
     &                    ncipstp, nciwstp, nccoum, 
     &                    ncigep, nccouge, ncigvp, nccougv, 
     &                    ncigpp, ncigwp,
     &                    nccougw
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype,
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22,
     &     eta11,etadd, nodecodenc, noffnc,
     &     writenc61, writenc62, writenc63, writenc64,
     &     writenc71, writenc72, writenc73, writenc74,
     &     writenc67, readnc67, nodecodenc,noffnc,
     &     nbvvnc,ibtypenc,
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc,
     &     ncmode


#endif
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND
      INTEGER I,J,IPROC,INDEX,IHOTSTP, not_active
      INTEGER IMHSF,ITHSF, NH, N
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW
      INTEGER N3DSD, I3DSDRec, N3DSV, I3DSVRec, N3DST, I3DSTRec
      INTEGER N3DGD, I3DGDRec, N3DGV, I3DGVRec, N3DGT, I3DGTRec
      CHARACTER FNAME*60,LOCFN*14
      CHARACTER*16 FNAME1
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME1)

      INTEGER,ALLOCATABLE  :: LOC(:),NODECODE(:),NOFF(:), domA(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:), 
     &    UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF
      integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

ckmd Added information for 3D hotstart
      REAL(SZ),ALLOCATABLE :: DUU(:),DUV(:),DVV(:)
      REAL(SZ),ALLOCATABLE :: UU(:),VV(:)
      REAL(SZ),ALLOCATABLE :: BSX(:),BSY(:)
      REAL(SZ),ALLOCATABLE :: WZ(:,:), q20(:,:), q20l(:,:)
      REAL(SZ),ALLOCATABLE :: SigT(:,:), Sal(:,:), Temp(:,:)
      REAL(SZ),ALLOCATABLE :: RealQ(:,:), ImagQ(:,:)
ckmd end of additions

      INTEGER INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV,
     &                                       IICALL,INFREQ,ITUD,NTSTEPS
      INTEGER IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
      CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
      REAL(8)  TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:),
     &                                         YVELAV(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)    
      REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
      REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
      REAL(SZ) TIME

C
C     ALLOCATE ARRAYS FOR NETCDF - MCF 
C
      REAL(SZ),ALLOCATABLE :: zeta1 (:,:)  
      REAL(SZ),ALLOCATABLE :: zeta2 (:,:)  
      REAL(SZ),ALLOCATABLE :: zetadisc (:,:) 
      REAL(SZ) DUMMY
      INTEGER IDUMMY
      INTEGER LUN
      INTEGER NHS
     
      integer, save       :: FileFmtMajor
      integer, save       :: FileFmtMinor
      integer, save       :: FileFmtRev     
!      REAL(SZ),ALLOCATABLE :: uu22 (:,:) 
!      REAL(SZ),ALLOCATABLE :: vv22 (:,:) 
!      INTEGER, ALLOCATABLE :: nodecodenc (:) 
!      INTEGER, ALLOCATABLE :: noffnc (:) 

C
C--   Open Appropriate Hot Start File based on the value of IHOT from
C--   the fort.15 file
C
      write(*,*) "enter IHOT: "
      read(*,*) IHOT
      IF(IHOT.EQ.67) FNAME='fort.67'
      IF(IHOT.EQ.68) FNAME='fort.68'

C     START ADDITION FOR NETCDF - MCF - 6/27/08  

      IF(IHOT.EQ.-67) FNAME='fort.67.nc'
      IF(IHOT.EQ.-68) FNAME='fort.68.nc'

C     END ADDITION FOR NETCDF - MCF - 6/27/08  

      INQUIRE(FILE=FNAME,EXIST=FOUND)
      IF (FOUND) THEN
         WRITE(*,1011) FNAME
         IF(IHOT.EQ.67.OR.IHOT.EQ.68)
     &      OPEN(IHOT,FILE=FNAME,ACCESS='DIRECT',RECL=8)
      ELSE
         WRITE(*,1010) FNAME
         STOP
      ENDIF


C     START ADDITION FOR NETCDF - MCF - 6/27/08  

C     If hotstart files are binary files then:

      IF(IHOT.EQ.67.or.IHOT.EQ.68) THEN

C--   Read in info from global hot start files

      IHOTSTP=1
      READ(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1

      if (.not. CMP_VERSION_NUMBERS(InputFileFmtVn, FileFmtVersion)) then
        write(*, 1001)
        write(*, 1012)
        write(*, 9973)
        stop
      endif

      READ(IHOT,REC=IHOTSTP) IMHSF        ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) TIMEHSF      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) ITHSF        ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NP_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NP_A_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_A_IN      ; IHOTSTP = IHOTSTP + 1
      if (nnodg == np_g_in) then
        MNP = nnodg
      else
        print *, "number global nodes does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nnodg
        write(*,'(A,I8)') "hotstart value = ", np_g_in
        stop 
      endif
      if (nelg ==  ne_g_in) then
        MNE = nelg
      else
        print *, "number global elements does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nelg  
        write(*,'(A,I8)') "hotstart value = ", ne_g_in
        stop 
      endif
C
C Allocate local work arrays
C
      nbytes = 4*nproc 
      call memory_alloc(nbytes)
      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      nbytes = 7*mnp*8   
      call memory_alloc(nbytes)
      ALLOCATE ( NOFF(MNE) )
      nbytes = 4*mne   
      call memory_alloc(nbytes)

      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      nbytes = 32*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      nbytes = 16*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      nbytes = 32*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      nbytes = 16*mnharf*mnstae
      call memory_alloc(nbytes)
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      nbytes = 16*mnharf*mnstav
      call memory_alloc(nbytes)
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
      nbytes = 16*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( XVELAV(MNP),XVELVA(MNP),YVELAV(MNP),YVELVA(MNP) )
      nbytes = 32*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      nbytes = 12*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( INAMEFR(MNHARF) )
      nbytes = 4*mnharf
      call memory_alloc(nbytes)

ckmd Added information for 3D hotstart
      ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
      nbytes = 3*mnp*8
      call memory_alloc(nbytes)
      ALLOCATE ( UU(MNP),VV(MNP))
      nbytes = 2*mnp*8
      call memory_alloc(nbytes)
      ALLOCATE ( BSX(MNP),BSY(MNP))
      nbytes = 2*mnp*8
      call memory_alloc(nbytes)
      ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_alloc(nbytes)
      ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_alloc(nbytes)
      ALLOCATE ( q20l(MNP,NFEN), SigT(MNP,NFEN))
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_alloc(nbytes)
      ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_alloc(nbytes)
ckmd end of additions

!  Continue reading global hot start file
      print *, "continuing to read global hotstart file"

      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) ETA1(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) ETA2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) EtaDisc(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) UU2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) VV2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      IF(IM.EQ.10) THEN
        DO I=1,MNP
           READ(IHOT,REC=IHOTSTP) CH1(I) ; IHOTSTP = IHOTSTP + 1
        END DO
      ENDIF
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) NODECODE(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNE
         READ(IHOT,REC=IHOTSTP) NOFF(I)  ; IHOTSTP = IHOTSTP + 1                    
      END DO

      READ(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
C
Ckmd Added information for 3D hotstart     
C     jgf46.02 Read in 3D hotstart data if appropriate
      IF (C3D) THEN
C     Start reading in the data
        READ(IHOT,REC=IHOTSTP) IDEN
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DSD
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DSDRec
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DSV
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DSVRec
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DST
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DSTRec
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DGD
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DGDRec
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DGV
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DGVRec
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) N3DGT
        IHOTSTP=IHOTSTP+1
        READ(IHOT,REC=IHOTSTP) I3DGTRec
        IHOTSTP=IHOTSTP+1

        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) DUU(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) DUV(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) DVV(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) UU(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) VV(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) BSX(NH)
          IHOTSTP=IHOTSTP+1
          READ(IHOT,REC=IHOTSTP) BSY(NH)
          IHOTSTP=IHOTSTP+1
        ENDDO

        DO NH=1,MNP
          DO N=1,NFEN
            READ(IHOT,REC=IHOTSTP) RealQ(NH,N)
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) ImagQ(NH,N)
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) WZ(NH,N)
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) q20(NH,N)
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) q20l(NH,N)
            IHOTSTP=IHOTSTP+1
            IF(ABS(IDen).EQ.1) THEN
              READ(IHOT,REC=IHOTSTP) SigT(NH,N)
              IHOTSTP=IHOTSTP+1
            ENDIF
            IF(ABS(IDen).EQ.2) THEN
              READ(IHOT,REC=IHOTSTP) Sal(NH,N)
              IHOTSTP=IHOTSTP+1
            ENDIF
            IF(ABS(IDen).EQ.3) THEN
              READ(IHOT,REC=IHOTSTP) Temp(NH,N)
              IHOTSTP=IHOTSTP+1
            ENDIF
            IF(ABS(IDen).EQ.4) THEN
              READ(IHOT,REC=IHOTSTP) Sal(NH,N)
              IHOTSTP=IHOTSTP+1
              READ(IHOT,REC=IHOTSTP) Temp(NH,N)
              IHOTSTP=IHOTSTP+1
            ENDIF
          ENDDO
        ENDDO
ckmd end of additions


      ENDIF
      
      
      ELSEIF(IHOT.EQ.-67.or.IHOT.EQ.-68) THEN

C   Start Addition for netCDF by MCF 6/27/08
C
C   READ NETCDF HOTSTART FILE
#ifdef NETCDF
!      PRINT *, "hstart: IHOT",IHOT
      IF (IHOT.EQ.-67.OR.IHOT.EQ.-68)THEN
      READ(14,*)AGRID
      READ(14,*)MNE,MNP
            DO I = 1,MNP
              READ(14,*)
      ENDDO
       DO I = 1,MNE
          READ(14,*)
       ENDDO
!      PRINT *, "prep: MNE,MNP",MNE,MNP
      READ(14,*) NOPE
!      WRITE(6,*) NOPE
      READ(14,*) NETA
!      WRITE(6,*) NETA
!      PRINT *, "prep: NOPE,NETA",NOPE,NETA
      
      DO J=1,NOPE
         READ(14,*) IDUMMY
         DO I=1,IDUMMY
            READ(14,*)
         ENDDO
      ENDDO
      READ(14,*) NBOU
      READ(14,*) NVEL
!      PRINT *, "prep: NBOU,NVEL",NBOU,NVEL

      
         timenc_len=1
!         timenc(timenc_len)=time
!         timenc(timenc_len)=TIMEHSF
         ncmode=0                                                      !serial version
         
!         FileFmtMajor = 1
!         FileFmtMinor = 2
!         FileFmtRev   = 0
          InputFileFmtVn = VERSION_NUMBER(FileFmtMajor, 
     &                                    FileFmtMinor,
     &                                    FileFmtRev)
         
         IF(.NOT. ALLOCATED(timenc))  ALLOCATE(timenc(timenc_len))
       IF(.NOT. ALLOCATED(zeta1 ))   ALLOCATE( zeta1(mnp, timenc_len))
       IF(.NOT. ALLOCATED(zeta2 ))   ALLOCATE( zeta2(mnp, timenc_len))
       IF(.NOT. ALLOCATED(zEtaDisc ))ALLOCATE(zEtaDisc(mnp, timenc_len))
       IF(.NOT. ALLOCATED(uu22 ))    ALLOCATE(uu22(mnp, timenc_len))
       IF(.NOT. ALLOCATED(vv22 ))    ALLOCATE(vv22(mnp, timenc_len))
       IF(.NOT. ALLOCATED(nodecodenc))ALLOCATE(nodecodenc(mnp))
       IF(.NOT. ALLOCATED(noffnc ))  ALLOCATE(noffnc (mne))
       IF(.NOT. ALLOCATED(nm ))  ALLOCATE(nm (mne,3))
       IF(.NOT. ALLOCATED(xnc))  ALLOCATE(xnc (mnp))
       IF(.NOT. ALLOCATED(ync))  ALLOCATE(ync (mnp))
       IF(.NOT. ALLOCATED(dp))  ALLOCATE(dp (mnp))
        IF(.NOT. ALLOCATED(NBDV)) ALLOCATE(NBDV(NOPE,NETA))
            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))
            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
            IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
            IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
            IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
            IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))

!     &   iestp,nscoue,ivstp,nscouv, icstp, 
!     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
!     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
!     &                  nscougw,


!     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
!     &                    nscouc, ncipstp, nciwstp, nccoum, 
!     &                    ncigep, nccouge, ncigvp, nccougv, 
!     &                    igcp, nscougc, ncigpp, ncigwp,
!     &                    nccougw,


         call readnc67 (nhsinc,ABS(IHOT),timenc_len,mnp,timenc,
     &                  nt, nhsinc, uu22,vv22,
     &                  zeta1, zeta2, zEtaDisc, nodecodenc,noffnc,
     &                  IMHSF, ITHSF,
     &   iestp,nscoue,ivstp,nscouv, icstp, 
     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                  nscougw,
     &                  rundes, runid, agrid, irtype,
     &                  dtdp, ihot, ics, nolibf, nolifa,
     &                  nolica, nolicat, nwp, 
     &                  ncor, ntip, nws, nramp, tau0,
     &                  statim, reftim, rnday, dramp,
     &                  a00, b00, c00, h0,
     &                  slam0, sfea0, cf ,eslm, cori,
     &                  ntif, nbfr, mne, nhy, nm,
     &                  xnc, ync, dp,
     &                  nope, neta, nbdv,nopenc,netanc,
     &                  nvdllnc,
     &                  nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                  nvellnc,ibtypenc,
     &                  ncmode,
     &                  title, institution, source, history,
     &                  references, comments, host, convention, contact,
     &                  base_date)
!       time= timenc(timenc_len)

         IF(IHOT.EQ.-67)THEN
         
          INQUIRE(FILE="fort.68.nc",EXIST=FOUND)
          
              IF ( .NOT. FOUND ) THEN
          
                 lun=68
                 nhs=3
          call writenc67 (nhs,abs(lun),timenc_len,mnp,timenc,
     &                    nt, nhsinc, uu22,vv22,
     &                  zeta1, zeta2, zEtaDisc, nodecodenc,noffnc,
     &                  IMHSF, ITHSF,
     &   iestp,nscoue,ivstp,nscouv, icstp, 
     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                  nscougw,
     &                  rundes, runid, agrid, irtype,
     &                  dtdp, ihot, ics, nolibf, nolifa,
     &                  nolica, nolicat, nwp, 
     &                  ncor, ntip, nws, nramp, tau0,
     &                  statim, reftim, rnday, dramp,
     &                  a00, b00, c00, h0,
     &                  slam0, sfea0, cf ,eslm, cori,
     &                  ntif, nbfr, mne, nhy, nm,
     &                  xnc, ync, dp,
     &                  nope, neta, nbdv,nopenc,netanc,
     &                  nvdllnc,
     &                  nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                  nvellnc,ibtypenc,
     &                  ncmode,
     &                    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &                    title, institution, source, history,
     &                    references, comments, host, convention, 
     &                    contact, base_date)
          
              ENDIF
          
         elseif(IHOT.EQ.-68)THEN
         
               INQUIRE(FILE="fort.67.nc",EXIST=FOUND)
          
               IF ( .NOT. FOUND ) THEN
          
                   lun=67
                   nhs=3
          call writenc67 (nhs,abs(lun),timenc_len,mnp,timenc,
     &                    nt, nhsinc, uu22,vv22,
     &                  zeta1, zeta2, zEtaDisc, nodecodenc,noffnc,
     &                  IMHSF, ITHSF,
     &   iestp,nscoue,ivstp,nscouv, icstp, 
     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                  nscougw,
     &                  rundes, runid, agrid, irtype,
     &                  dtdp, ihot, ics, nolibf, nolifa,
     &                  nolica, nolicat, nwp, 
     &                  ncor, ntip, nws, nramp, tau0,
     &                  statim, reftim, rnday, dramp,
     &                  a00, b00, c00, h0,
     &                  slam0, sfea0, cf ,eslm, cori,
     &                  ntif, nbfr, mne, nhy, nm,
     &                  xnc, ync, dp,
     &                  nope, neta, nbdv,nopenc,netanc,
     &                  nvdllnc,
     &                  nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                  nvellnc,ibtypenc,
     &                  ncmode,
     &                    FileFmtMajor, FileFmtMinor, FileFmtRev,
     &                    title, institution, source, history,
     &                    references, comments, host, convention, 
     &                    contact, base_date)
               ENDIF

         endif      
         
!          PRINT *, "LUN = ", LUN
!     iestp,nscoue,ivstp,nscouv, icstp,
!     &                    nscouc, ipstp, iwstp, nscoum, igep, nscouge,
!     &                    igvp, nscougv, igcp, nscougc, igpp, igwp,
!     &                    nscougw,


     
 !                 ENDIF


      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      ALLOCATE ( NOFF(MNE) )

!       PRINT *,"======== Print NC variables ================"
             DO I=1,MNP
                eta1 (I)=zeta1(I, 1)
!                PRINT *, "ETA1(I) = ",I, ETA1(I)
             END DO
             DO I=1,MNP
                eta2 (I)=zeta2(I, 1)
!                PRINT *, "ETA2(I) = ",I, ETA2(I)
             END DO
             DO I=1,MNP
                EtaDisc (I)=zEtaDisc(I, 1)
!                PRINT *, "EtaDisc(I) = ",I, EtaDisc(I)
             END DO
             DO I=1,MNP
                UU2 (I)=UU22(I, 1)
!                PRINT *, "UU2(I) = ",I, UU2(I)
             END DO
             DO I=1,MNP
                VV2 (I)=VV22(I, 1)
!                 PRINT *, "VV22(I) = ",I, VV2(I)
            END DO
             DO I=1,MNP
                nodecode (I)=nodecodenc(I)
!                  PRINT *, "nodecodenc(I) = ",I, nodecode(I)
            END DO
             DO I=1,MNE
                noff (I)=noffnc(I)
!                  PRINT *, "noff(I) = ",I, noff(I)
             END DO
             TIMEHSF=timenc(timenc_len)
!             NNODP=MNP
!             NELP=MNE
             
!                   PRINT *,"  ----------START--NETCDF-----------------"
!                   PRINT *,"  IMHSF = ", IMHSF
!                   PRINT *,"  TIMEHSF = ", TIMEHSF
!                   PRINT *,"  TIME = ", timenc(timenc_len)
!                   PRINT *,"  ITHSF = ", ITHSF
!                   PRINT *,"  ETA11 = ", ETA1(11)
!                   PRINT *,"  ETA22 = ", ETA2(11)
!                   PRINT *,"  ETADD = ", EtaDisc(11)
!                   PRINT *,"  uu22 = ", uu2(11)
!                   PRINT *,"  vv22 = ", vv2(11)
!                   PRINT *,"  nnodecode = ", nodecode(11)
!                   PRINT *,"  noff = ", noff(11)
!                   PRINT *,"  IESTP = ",IESTP
!                   PRINT *,"  NSCOUE = ",NSCOUE
!                   PRINT *,"  IVSTP = ",IVSTP
!                   PRINT *,"  NSCOUV = ",NSCOUV
!                   PRINT *,"  ICSTP = ",ICSTP
!                   PRINT *,"  NSCOUC = ",NSCOUC
!                   PRINT *,"  IPSTP = ",IPSTP
!                   PRINT *,"  IWSTP = ",IWSTP
!                   PRINT *,"  NSCOUM = ",NSCOUM
!                   PRINT *,"  IGEP = ",IGEP
!                   PRINT *,"  NSCOUGE = ",NSCOUGE
!                   PRINT *,"  IGVP = ",IGVP
!                   PRINT *,"  NSCOUGV = ",NSCOUGV
!                   PRINT *,"  IGCP = ",IGCP
!                   PRINT *,"  NSCOUGC = ",NSCOUGC
!                   PRINT *,"  IGPP = ",IGPP
!                   PRINT *,"  IGWP = ",IGWP
!                   PRINT *,"  NSCOUGW = ",NSCOUGW
!                   PRINT *,"  ---------END-------NETCDF---------------"
      ENDIF     


#endif

      
      ELSE
 
        PRINT *,"The hotstart files are neither binary nor netCDF files"

      ENDIF             
C   End Addition for netCDF by MCF 6/27/08
      
      
      
C     jgf48.03 harmonic analysis not supported yet
       IF(IHOT.GT.0) THEN
#if 0
C
C....DETERMINE HARMONIC ANALYSIS PARAMETERS

      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START 
C     INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS      

      IF(IHARIND.EQ.1) THEN
         ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
         ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
         ITMV = ITHAF - (ITHAF-ITHAS)*FMV
         IHABEG=ITHAS+NHAINC

C.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

         IF(ITHSF.GT.ITHAS) THEN
            READ(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
         ENDIF

         IF(ITHSF.GE.IHABEG) THEN
            READ(IHOT,REC=IHOTSTP) INZ ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INF ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IMM ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INP ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INSTAE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INSTAV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IISTAE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IISTAV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IIGLOE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IIGLOV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IICALL ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INFREQ ; IHOTSTP = IHOTSTP + 1
            
            DO I=1,INFREQ+INF
               READ(IHOT,REC=IHOTSTP) FNAM8(1) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) FNAM8(2) ; IHOTSTP = IHOTSTP + 1

               INAMEFR(I) = FNAME1
               READ(IHOT,REC=IHOTSTP) IFREQ(I) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) IFF(I) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) IFACE(I) ; IHOTSTP = IHOTSTP + 1
            ENDDO
            
            READ(IHOT,REC=IHOTSTP) TIMEUD ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) ITUD ; IHOTSTP = IHOTSTP + 1
            
            DO I=1,IMM
               DO J=1,IMM
                  READ(IHOT,REC=IHOTSTP) HA(I,J) ; IHOTSTP = IHOTSTP + 1
               ENDDO
            ENDDO
            
            IF(NHASE.EQ.1) THEN
               DO J=1,INSTAE
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) STAELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHASV.EQ.1) THEN
               DO J=1,INSTAV
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) STAULV(I,J) 
                     IHOTSTP = IHOTSTP + 1
                     READ(IHOT,REC=IHOTSTP) STAVLV(I,J) 
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHAGE.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) GLOELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHAGV.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) GLOULV(I,J)
                     IHOTSTP = IHOTSTP + 1
                     READ(IHOT,REC=IHOTSTP) GLOVLV(I,J)
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF
            
         ENDIF
         
         IF((FMV.GT.0.).AND.(INFREQ.GT.0).AND.(IM.EQ.0)) THEN !include means and variances
            IF(ITHSF.GT.ITMV) THEN
               READ(IHOT,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP) ELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) ELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
               IF(NHAGV.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP) XVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) YVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) XVELVA(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) YVELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
            ENDIF
         ENDIF    ! charmv
      ENDIF     ! HARIND
#endif
               ENDIF
C
C--Open All Local Hot Start files 
C

!                   PRINT *,"  ----------START-LOC-NETCDF---------------"
!                   PRINT *,"  IMHSF = ", IMHSF
!                   PRINT *,"  TIMEHSF = ", TIMEHSF
!                   PRINT *,"  TIME = ", time
!                   PRINT *,"  ITHSF = ", ITHSF
!                   PRINT *,"  ETA11 = ", ETA1(11)
!                   PRINT *,"  ETA22 = ", ETA2(11)
!                   PRINT *,"  ETADD = ", EtaDisc(11)
!                   PRINT *,"  uu22 = ", uu2(11)
!                   PRINT *,"  vv22 = ", vv2(11)
!                   PRINT *,"  nodecode = ", nodecode(11)
!                   PRINT *,"  noff = ", noff(11)
!                   PRINT *,"  IESTP = ",IESTP
!                   PRINT *,"  NSCOUE = ",NSCOUE
!                   PRINT *,"  IVSTP = ",IVSTP
!                   PRINT *,"  NSCOUV = ",NSCOUV
!                   PRINT *,"  ICSTP = ",ICSTP
!                   PRINT *,"  NSCOUC = ",NSCOUC
!                   PRINT *,"  IPSTP = ",IPSTP
!                   PRINT *,"  IWSTP = ",IWSTP
!                   PRINT *,"  NSCOUM = ",NSCOUM
!                   PRINT *,"  IGEP = ",IGEP
!                   PRINT *,"  NSCOUGE = ",NSCOUGE
!                   PRINT *,"  IGVP = ",IGVP
!                   PRINT *,"  NSCOUGV = ",NSCOUGV
!                   PRINT *,"  IGCP = ",IGCP
!                   PRINT *,"  NSCOUGC = ",NSCOUGC
!                   PRINT *,"  IGPP = ",IGPP
!                   PRINT *,"  IGWP = ",IGWP
!                   PRINT *,"  NSCOUGW = ",NSCOUGW
!                   DO I=1,MNP
!                     PRINT *, "ETA1(I) = ",I, ETA1(I)
!                   END DO
!
!                   PRINT *,"  ---------END----LOC---NETCDF------------"

      ALLOCATE ( LOC(NPROC) )
      DO IPROC = 1,NPROC
         LOC(IPROC) = 105 + (IPROC-1)
         LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (LOC(IPROC),FILE=LOCFN,ACCESS='DIRECT',RECL=8)
      ENDDO
C
C--Write out info to local hot start files
C
      DO IPROC = 1,NPROC
         IHOTSTP=1
         WRITE(LOC(IPROC),REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) IMHSF          ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) TIMEHSF        ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) ITHSF          ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NNODP(IPROC)   ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NELP(IPROC)    ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NNODP(IPROC)   ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NELP(IPROC)    ; IHOTSTP = IHOTSTP + 1

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) ETA1(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) ETA2(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) EtaDisc(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) UU2(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) VV2(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         IF(IM.EQ.10) THEN
            DO I=1, NNODP(IPROC)
               INDEX = ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOC(IPROC),REC=IHOTSTP) CH1(INDEX)
               IHOTSTP=IHOTSTP+1
            END DO
         ENDIF

         DO I=1, NNODP(IPROC)
            INDEX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) NODECODE(INDEX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1,NELP(IPROC)
            INDEX=ABS(IMAP_EL_LG(I,IPROC))
            WRITE(LOC(IPROC),REC=IHOTSTP) NOFF(INDEX)
            IHOTSTP=IHOTSTP+1                                    
         END DO
         
         WRITE(LOC(IPROC),REC=IHOTSTP) IESTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUE  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IVSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUV  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) ICSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUC  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IPSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) IWSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUM  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IGEP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUGE  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IGVP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUGV  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IGCP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUGC  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC(IPROC),REC=IHOTSTP) IGPP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) IGWP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC(IPROC),REC=IHOTSTP) NSCOUGW  ; IHOTSTP = IHOTSTP + 1

Ckmd Added information for 3D hotstart
C     jgf46.02 Write out 3D hotstart data if appropriate
         IF (C3D) THEN
C     Start writing out the 3D hotstart information
            WRITE(LOC(IPROC),REC=IHOTSTP) IDEN
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DSD
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DSDRec
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DSV
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DSVRec
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DST
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DSTRec
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DGD
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DGDRec
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DGV
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DGVRec
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) N3DGT
            IHOTSTP=IHOTSTP+1
            WRITE(LOC(IPROC),REC=IHOTSTP) I3DGTRec
            IHOTSTP=IHOTSTP+1

            DO I=1,NNODP(IPROC)
               INDEX=IMAP_NOD_LG(I,IPROC)
               WRITE(LOC(IPROC),REC=IHOTSTP) DUU(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) DUV(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) DVV(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) UU(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) VV(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) BSX(INDEX)
               IHOTSTP=IHOTSTP+1
               WRITE(LOC(IPROC),REC=IHOTSTP) BSY(INDEX)
               IHOTSTP=IHOTSTP+1
             ENDDO

             DO I=1,NNODP(IPROC)
               INDEX=IMAP_NOD_LG(I,IPROC)
               DO N=1,NFEN
                 WRITE(LOC(IPROC),REC=IHOTSTP) RealQ(INDEX,N)
                 IHOTSTP=IHOTSTP+1
                 WRITE(LOC(IPROC),REC=IHOTSTP) ImagQ(INDEX,N)
                 IHOTSTP=IHOTSTP+1
                 WRITE(LOC(IPROC),REC=IHOTSTP) WZ(INDEX,N)
                 IHOTSTP=IHOTSTP+1
                 WRITE(LOC(IPROC),REC=IHOTSTP) q20(INDEX,N)
                 IHOTSTP=IHOTSTP+1
                 WRITE(LOC(IPROC),REC=IHOTSTP) q20l(INDEX,N)
                 IF(ABS(IDen).EQ.1) THEN
                   IHOTSTP=IHOTSTP+1
                   WRITE(LOC(IPROC),REC=IHOTSTP) SigT(INDEX,N)
                 ENDIF
                 IF(ABS(IDen).EQ.2) THEN
                   IHOTSTP=IHOTSTP+1
                   WRITE(LOC(IPROC),REC=IHOTSTP) Sal(INDEX,N)
                 ENDIF
                 IF(ABS(IDen).EQ.3) THEN
                   IHOTSTP=IHOTSTP+1
                   WRITE(LOC(IPROC),REC=IHOTSTP) Temp(INDEX,N)
                 ENDIF
                 IF(ABS(IDen).EQ.4) THEN
                   IHOTSTP=IHOTSTP+1
                   WRITE(LOC(IPROC),REC=IHOTSTP) Sal(INDEX,N)
                   IHOTSTP=IHOTSTP+1
                   WRITE(LOC(IPROC),REC=IHOTSTP) Temp(INDEX,N)
                   IHOTSTP=IHOTSTP+1
                 ENDIF
               ENDDO
             ENDDO
ckmd end of additions
         ENDIF
#if 0
C 
C....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

c       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
c         WRITE(LOC(IPROC),REC=IHOTSTP+1) ICHA
c         IHOTSTP = IHOTSTP + 1
c         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
c    &                LOC(IPROC),IHOTSTP)
c
c         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC(IPROC),IHOTSTP)
c         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC(IPROC),IHOTSTP)
c         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC(IPROC),IHOTSTP)
c         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC(IPROC),IHOTSTP)
c         ENDIF
c
c       if(CHARMV) then
c         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
c           IHOTSTP=IHOTSTP+1
c           WRITE(LOC(IPROC),REC=IHOTSTP) NTSTEPS
c           IF(NHAGE.EQ.1) THEN
c             DO I=1, NNODP(IPROC)
c               INDEX = IMAP_NOD_LG(I,IPROC)
c               DO I=1,MNP
c                 WRITE(LOC(IPROC),REC=IHOTSTP+1) ELAV(INDEX)
c                 WRITE(LOC(IPROC),REC=IHOTSTP+2) ELVA(INDEX)
c                 IHOTSTP=IHOTSTP+2
c                 END DO
c             ENDIF
c           IF(NHAGV.EQ.1) THEN
c             DO I=1,NNODP(IPROC)
c               WRITE(LOC(IPROC),REC=IHOTSTP+1) XVELAV(INDEX)
c               WRITE(LOC(IPROC),REC=IHOTSTP+2) YVELAV(INDEX)
c               WRITE(LOC(IPROC),REC=IHOTSTP+3) XVELVA(INDEX)
c               WRITE(LOC(IPROC),REC=IHOTSTP+4) YVELVA(INDEX)
c               IHOTSTP=IHOTSTP+4
c               END DO
c             ENDIF
c           ENDIF
c         ENDIF               
#endif

      ENDDO
C
C--Close Global file and all the Local Files     
C
!            PRINT *, "BEFORE CLOSE GLOBAL = ", IHOT

      IF(IHOT.EQ.67.OR.IHOT.EQ.68)CLOSE (IHOT)
!            PRINT *, "BEFORE CLOSE LOCAL = ", IHOT
      DO IPROC=1, NPROC
         CLOSE (LOC(IPROC))
      ENDDO
!            PRINT *, "AFTER CLOSE LOCAL = ", IHOT
C
      IF(ALLOCATED(LOC)) DEALLOCATE ( LOC )
      nbytes = 4*nproc
      call memory_dealloc(nbytes)
      IF(ALLOCATED( ETA1 ))DEALLOCATE ( ETA1  )
      IF(ALLOCATED( ETA2 ))DEALLOCATE ( ETA2  )
      IF(ALLOCATED( EtaDisc ))DEALLOCATE ( EtaDisc )
      IF(ALLOCATED( UU2 ))DEALLOCATE ( UU2 )
      IF(ALLOCATED( VV2 ))DEALLOCATE ( VV2 )
      IF(ALLOCATED( NODECODE ))DEALLOCATE ( NODECODE )
      IF(ALLOCATED( CH1 ))DEALLOCATE ( CH1 )
      nbytes = 7*mnp*8  
      call memory_dealloc(nbytes)
      IF(ALLOCATED(NOFF))DEALLOCATE ( NOFF )
      nbytes = 6*mne  
      call memory_dealloc(nbytes)
      IF(ALLOCATED( DUU ))DEALLOCATE ( DUU  )
      IF(ALLOCATED( DUV ))DEALLOCATE ( DUV  )
      IF(ALLOCATED( DVV ))DEALLOCATE ( DVV )
      nbytes = 3*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED( UU )) DEALLOCATE ( UU )
      IF(ALLOCATED( VV )) DEALLOCATE ( VV )
      nbytes = 2*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( BSX )) DEALLOCATE ( BSX )
      IF(ALLOCATED  ( BSY )) DEALLOCATE ( BSY )
      nbytes = 2*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( WZ )) DEALLOCATE ( WZ )
      IF(ALLOCATED  ( q20 )) DEALLOCATE (q20 )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( RealQ ))  DEALLOCATE ( RealQ)
      IF(ALLOCATED  ( ImagQ ))  DEALLOCATE ( ImagQ)
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( q20l )) DEALLOCATE ( q20l )
      IF(ALLOCATED  ( SigT )) DEALLOCATE ( SigT )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED ( Sal)) DEALLOCATE ( Sal  )
      IF(ALLOCATED ( Temp )) DEALLOCATE (  Temp )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      call memory_status()
C
      RETURN
 1001 FORMAT('ERROR: The hot start file')
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1012 FORMAT('was a nonmatching version') 
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE HOTLOCALIZE


      SUBROUTINE HOTGLOBALIZE()
      USE PRE_GLOBAL 
      use presizes; use memory_usage
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C
C             started mods for harmonic analysis and 3D RL 5/22/03          C
C         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
C                                                                           C
C  This routine reads the global hot start file (either fort.67 or fort.68) C
C  and writes local hot start files of the same format.                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND
      INTEGER I,J,IPROC,INDEX,IHOTSTP, not_active
      INTEGER IMHSF,ITHSF,IVALUE,IDUMY, NH, N
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW
      CHARACTER FNAME*60,LOCFN*14
      CHARACTER*16 FNAME1
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME1)

      INTEGER,ALLOCATABLE  :: LOC(:),NODECODE(:),NOFF(:), domA(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:), 
     &    UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF, RVALUE
      integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

#if 0
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet 
      INTEGER INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV,
     &                                       IICALL,INFREQ,ITUD,NTSTEPS
      INTEGER IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
      CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
      REAL(8)  TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:),
     &                                         YVELAV(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)    
      REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
      REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
#endif
C
C--   Open Appropriate Hot Start File based on the value of IHOT from
C--   the fort.15 file
C
      write(*,*) "enter IHOT: "
      read(*,*) IHOT
      IF(IHOT.EQ.67) FNAME='fort.67'
      IF(IHOT.EQ.68) FNAME='fort.68'

C
C--Open All Local Hot Start files 
C
      ALLOCATE ( LOC(NPROC) )
      DO IPROC = 1,NPROC
         LOC(IPROC) = 105 + (IPROC-1)
         LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         INQUIRE(FILE=LOCFN,EXIST=FOUND)
         IF (FOUND) THEN
           WRITE(*,1011) LOCFN
           OPEN (LOC(IPROC),FILE=LOCFN,ACCESS='DIRECT',RECL=8)
         ELSE
          WRITE(*,1010) FNAME
          STOP
         ENDIF
      ENDDO
C
C Allocate local work arrays
C
      MNP  =  nnodg    !  global number of nodes    ( read from fort.18 )
      print *, "MNP =", MNP

      nbytes = 4*nproc 
      call memory_alloc(nbytes)
      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      nbytes = 7*mnp   

      MNE  =  nelg     !  global number of elements ( read from fort.18 )
      print *, "MNE =", MNE

      call memory_alloc(nbytes)
      ALLOCATE ( NOFF(MNE) )
      nbytes = 4*mne   
      call memory_alloc(nbytes)

#if HA
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet 
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      nbytes = 32*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      nbytes = 16*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      nbytes = 32*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      nbytes = 16*mnharf*mnstae
      call memory_alloc(nbytes)
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      nbytes = 16*mnharf*mnstav
      call memory_alloc(nbytes)
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
      nbytes = 16*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( XVELAV(MNP),XVELVA(MNP),YVELAV(MNP),YVELVA(MNP) )
      nbytes = 32*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      nbytes = 12*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( INAMEFR(MNHARF) )
      nbytes = 4*mnharf
      call memory_alloc(nbytes)
#endif

C
C--Read info from local hot start files
C
      DO IPROC = 1,NPROC
         IHOTSTP=1
         READ(LOC(IPROC),REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IMHSF          ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) TIMEHSF        ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) ITHSF          ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) ETA1(INDEX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) ETA2(INDEX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) EtaDisc(INDEX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) UU2(INDEX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) VV2(INDEX) = RVALUE
         END DO

         IF(IM.EQ.10) THEN
            DO I=1, NNODP(IPROC)
               INDEX = IMAP_NOD_LG(I,IPROC)
               READ(LOC(IPROC),REC=IHOTSTP) RVALUE           
               IHOTSTP=IHOTSTP+1
               IF (INDEX > 0) CH1(INDEX) = RVALUE
            END DO
         ENDIF

         DO I=1, NNODP(IPROC)
            INDEX = IMAP_NOD_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) IVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) NODECODE(INDEX) = IVALUE
         END DO

         DO I=1,NELP(IPROC)
            INDEX = IMAP_EL_LG(I,IPROC)
            READ(LOC(IPROC),REC=IHOTSTP) IVALUE           
            IHOTSTP=IHOTSTP+1
            IF (INDEX > 0) NOFF(INDEX) = IVALUE
         END DO
         
         READ(LOC(IPROC),REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

         READ(LOC(IPROC),REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
         READ(LOC(IPROC),REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
#if HA
C
C     jgf46.02 Write out 3D hotstart data if appropriate
         IF (C3D) THEN
            CALL WriteHotStart3D(LOC(IPROC),IHOTSTP)
         ENDIF
C 
C....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

c       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
c         READ(LOC(IPROC),REC=IHOTSTP+1) ICHA
c         IHOTSTP = IHOTSTP + 1
c         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
c    &                LOC(IPROC),IHOTSTP)
c
c         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC(IPROC),IHOTSTP)
c         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC(IPROC),IHOTSTP)
c         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC(IPROC),IHOTSTP)
c         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC(IPROC),IHOTSTP)
c         ENDIF
c
c       if(CHARMV) then
c         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
c           IHOTSTP=IHOTSTP+1
c           READ(LOC(IPROC),REC=IHOTSTP) NTSTEPS
c           IF(NHAGE.EQ.1) THEN
c             DO I=1, NNODP(IPROC)
c               INDEX = IMAP_NOD_LG(I,IPROC)
c               DO I=1,MNP
c                 READ(LOC(IPROC),REC=IHOTSTP+1) ELAV(INDEX)
c                 READ(LOC(IPROC),REC=IHOTSTP+2) ELVA(INDEX)
c                 IHOTSTP=IHOTSTP+2
c                 END DO
c             ENDIF
c           IF(NHAGV.EQ.1) THEN
c             DO I=1,NNODP(IPROC)
c               READ(LOC(IPROC),REC=IHOTSTP+1) XVELAV(INDEX)
c               READ(LOC(IPROC),REC=IHOTSTP+2) YVELAV(INDEX)
c               READ(LOC(IPROC),REC=IHOTSTP+3) XVELVA(INDEX)
c               READ(LOC(IPROC),REC=IHOTSTP+4) YVELVA(INDEX)
c               IHOTSTP=IHOTSTP+4
c               END DO
c             ENDIF
c           ENDIF
c         ENDIF               
#endif
          CLOSE (LOC(IPROC))

      ENDDO

C-----------------------------------------------------------------------
C--   Write info to global hot start files
C-----------------------------------------------------------------------

      OPEN(IHOT,FILE=trim(FNAME),ACCESS='DIRECT',RECL=8)
      print *, "opening global hotstart file"

      IHOTSTP=1
      WRITE(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IMHSF        ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) TIMEHSF      ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) ITHSF        ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNP          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNE          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNP          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNE          ; IHOTSTP = IHOTSTP + 1

      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) ETA1(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) ETA2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) EtaDisc(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) UU2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) VV2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      IF(IM.EQ.10) THEN
        DO I=1,MNP
           WRITE(IHOT,REC=IHOTSTP) CH1(I)
           IHOTSTP=IHOTSTP+1
        END DO
      ENDIF
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) NODECODE(I)
         IHOTSTP=IHOTSTP+1
      END DO
      DO I=1,MNE
         WRITE(IHOT,REC=IHOTSTP) NOFF(I)                     
         IHOTSTP=IHOTSTP+1                                    
      END DO

      WRITE(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
C
     
       ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet 
#if 0
C     jgf46.02 Read in 3D hotstart data if appropriate
      IF (C3D) THEN
         CALL ReadHotStart3D(IHOT,IHOTSTP)
      ENDIF
C
C.....DETERMINE HARMONIC ANALYSIS PARAMETERS

      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START INFORMATION FOR
C.....IN PROGRESS HARMONIC ANALYSIS      

      IF(IHARIND.EQ.1) THEN
         ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
         ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
         ITMV = ITHAF - (ITHAF-ITHAS)*FMV
         IHABEG=ITHAS+NHAINC

C.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

         IF(ITHSF.GT.ITHAS) THEN
            WRITE(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
         ENDIF

         IF(ITHSF.GE.IHABEG) THEN
            WRITE(IHOT,REC=IHOTSTP) INZ ; IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INF ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IMM ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INP ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INSTAE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INSTAV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IISTAE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IISTAV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IIGLOE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IIGLOV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IICALL ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INFREQ ;  IHOTSTP = IHOTSTP + 1
            
            DO I=1,INFREQ+INF
               WRITE(IHOT,REC=IHOTSTP) FNAM8(1) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) FNAM8(2) ;  IHOTSTP = IHOTSTP + 1
               INAMEFR(I) = FNAME1
               WRITE(IHOT,REC=IHOTSTP) IFREQ(I) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) IFF(I) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) IFACE(I) ;  IHOTSTP = IHOTSTP + 1
            ENDDO
            
            WRITE(IHOT,REC=IHOTSTP) TIMEUD ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) ITUD ;  IHOTSTP = IHOTSTP + 1
            
            DO I=1,IMM
               DO J=1,IMM
                  WRITE(IHOT,REC=IHOTSTP) HA(I,J)
                  IHOTSTP = IHOTSTP + 1
               ENDDO
            ENDDO
            
            IF(NHASE.EQ.1) THEN
               DO J=1,INSTAE
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) STAELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHASV.EQ.1) THEN
               DO J=1,INSTAV
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) STAULV(I,J)
                     IHOTSTP=IHOTSTP+1                    
                     WRITE(IHOT,REC=IHOTSTP) STAVLV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHAGE.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) GLOELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF
            
            IF(NHAGV.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) GLOULV(I,J)
                     IHOTSTP = IHOTSTP + 1
                     WRITE(IHOT,REC=IHOTSTP) GLOVLV(I,J)
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF
            
         ENDIF
         
         IF((FMV.GT.0.).AND.(INFREQ.GT.0).AND.(IM.EQ.0)) THEN !include means and variances
            IF(ITHSF.GT.ITMV) THEN
               WRITE(IHOT,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.EQ.1) THEN
                  DO I=1,INP
                     WRITE(IHOT,REC=IHOTSTP) ELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) ELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
               IF(NHAGV.EQ.1) THEN
                  DO I=1,INP
                     WRITE(IHOT,REC=IHOTSTP) XVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) YVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) XVELVA(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) YVELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
            ENDIF
         ENDIF    ! charmv
      ENDIF     ! HARIND
#endif

C
C--Close Global file and all the Local Files     
C
      CLOSE (IHOT)
C
      DEALLOCATE ( LOC )
      nbytes = 4*nproc
      call memory_dealloc(nbytes)
      DEALLOCATE ( ETA1, ETA2, EtaDisc, UU2, VV2, NODECODE, CH1 )
      nbytes = 7*mnp  
      call memory_dealloc(nbytes)
      DEALLOCATE ( NOFF )
      nbytes = 6*mne  
      call memory_dealloc(nbytes)
      call memory_status()
C
      RETURN
 1001 FORMAT('ERROR: The hot start file')
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1012 FORMAT('was a nonmatching version') 
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE HOTGLOBALIZE






C     ----------------------------------------------------------------------
C          S U B R O U T I N E     R E A D  H O T  S T A R T  3 D
C     ----------------------------------------------------------------------
C
C     jgf46.02 This subroutine supports PREP67_68. It reads in the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE ReadHotStart3D(UnitNumber,FilePosition)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL 
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
C
      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE ReadHotStart3D
C     ----------------------------------------------------------------------



C     ----------------------------------------------------------------------
C          S U B R O U T I N E     W R I T E   H O T  S T A R T  3 D
C     ----------------------------------------------------------------------
C
C     jgf46.02 This subroutine supports PREP67_68. It writes out the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE WriteHotStart3D(UnitNumber,FilePosition)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL 
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
C
      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE WriteHotStart3D
C     ----------------------------------------------------------------------



C---------------------------------------------------------------------------
C     S U B R O U T I N E    O P E N  F U L L  D O M A I N  F I L E 
C---------------------------------------------------------------------------
C
C     jgf47.02 This subroutine will open the full domain file
C
C---------------------------------------------------------------------------
      SUBROUTINE OpenFullDomainFile(UnitNumber, Description, Success) 
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open 
      CHARACTER(len=30), intent(in) :: Description ! description of file
      LOGICAL, intent(out):: Success     ! .true. if file opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=7) DefaultName! default name of full domain file
      INTEGER ErrorIO             ! zero if file opened successfully
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'

      DefaultName(1:5) = 'fort.'
      WRITE(DefaultName(6:7),2) UnitNumber
C
C     Determine the name of the file; if found, open it
 31   IF (USE_DEFAULT) THEN     
         FileName = DefaultName
      ELSE
         WRITE(*,850) ! type skip to bypass
         WRITE(*,900) Description
         WRITE(*,910) UnitNumber
         READ(*,'(A)') FileName
         FileName = trim(FILENAME)
      ENDIF
C
C     Determine if full domain file exists
      INQUIRE(FILE=FileName,EXIST=FOUND)
C
C     If it does exist, open it
      IF ( FOUND ) THEN
         WRITE(*,1011) FileName !found
         OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .false.
         ELSE 
            WRITE(*,*) "Successfully opened full domain file."
         ENDIF
      ELSE
C     Give the user a chance to opt out of prepping this file.
         IF (FileName .eq. skipstring) RETURN ! note the early RETURN
         WRITE(*,1010) FileName !not found
         GOTO 31
      ENDIF
 
 2    FORMAT(I2)
 30   FORMAT(A30)
 850  FORMAT(/,'Type ''skip'' to bypass preprocessing or')
 900  FORMAT('Enter the name of the ',A30)
 910  FORMAT('file (unit ',I3,'): ')
 1010 FORMAT('File ',A7,/,' WAS NOT FOUND! Try again or type "skip"',/)
 1011 FORMAT('File ',A7,/,' WAS FOUND!  Opening & Processing file.',/)
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE OpenFullDomainFile
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C     S U B R O U T I N E   O P E N  S U B D O M A I N  F I L E 
C---------------------------------------------------------------------------
C
C     jgf47.02 This subroutine will open a single subdomain file
C
C---------------------------------------------------------------------------
      SUBROUTINE OpenSubDomainFile(UnitNumber, IProc, sdu, Success) 
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number of full dom file 
      INTEGER, intent(in) :: iproc          ! subdomain number
      INTEGER, intent(out) :: sdu  ! i/o unit nunber that was opened
      LOGICAL, intent(out):: Success     ! .true. if files opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=7) DefaultName! default name of full domain file
      INTEGER ErrorIO             ! zero if file opened successfully
      CHARACTER*14 sdFileName     ! subdomain file name
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'

      DefaultName(1:5) = 'fort.'
      WRITE(DefaultName(6:7),2) UnitNumber

C     Open subdomain file
      sdu = 105 + (iproc-1)
      sdFileName(1:7) = 'PE0000/'
      sdFileName(8:14) = DefaultName
      CALL IWRITE(sdFileName, 3, 6, iproc-1)
      OPEN (UNIT=sdu, FILE=sdFileName, IOSTAT=ErrorIO)
      Success = .true.
      IF ( ErrorIO .GT. 0 ) THEN
         WRITE(*,*) "ERROR: Subdomain file cannot be opened."
         print *,sdu
         print *,sdFileName
         Success = .false.
      ENDIF
 2    FORMAT(I2)
      RETURN

C---------------------------------------------------------------------------
      END SUBROUTINE OpenSubDomainFile
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C           S U B R O U T I N E   O P E N  P R E P  F I L E S
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will open the full domain file and each
C     of the subdomain files in the range of subdomains provided in the
C     arguments. It assumes that all unit numbers are between 10 and 99.
C
C---------------------------------------------------------------------------
      SUBROUTINE OpenPrepFiles(UnitNumber, Description,
     &     startProc, endProc, SDU, Success) 
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open 
      CHARACTER(len=30), intent(in) :: Description ! description of file
      INTEGER, intent(in) :: startProc        ! subdomains to start with
      INTEGER, intent(in) :: endProc          ! subdomain to end on
      INTEGER, intent(out), dimension(nproc) :: SDU ! Subdomain unit numbers
      LOGICAL, intent(out):: Success     ! .true. if files opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=7) DefaultName! default name of full domain file
      INTEGER ErrorIO             ! zero if file opened successfully
      INTEGER iproc               ! subdomain index
      CHARACTER(len=14) sdFileName     ! subdomain file name
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'

      DefaultName(1:5) = 'fort.'
      WRITE(DefaultName(6:7),2) UnitNumber
C
C     Determine the name of the file; if found, open it
 31   IF (USE_DEFAULT) THEN     
         FileName = DefaultName
      ELSE
         WRITE(*,850) ! type skip to bypass
         WRITE(*,900) Description
         WRITE(*,910) UnitNumber
         READ(*,'(A)') FileName
         FileName = trim(FILENAME)
      ENDIF
C
C     Determine if full domain file exists
      INQUIRE(FILE=FileName,EXIST=FOUND)
C
C     If it does exist, open it
      IF ( FOUND ) THEN
         WRITE(*,1011) FileName !found
         OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .false.
         ENDIF
      ELSE
C     Give the user a chance to opt out of prepping this file.
         IF (FileName .eq. skipstring) RETURN ! note the early RETURN
         WRITE(*,1010) FileName !not found
         GOTO 31
      ENDIF
C
      If (.not.Success) RETURN ! failed to open full domain file
C
C     Open each of the subdomain files
      DO iproc = startProc, endProc
         sdu(iproc) = 105 + (iproc-1)
         sdFileName(1:7) = 'PE0000/'
         sdFileName(8:14) = DefaultName
#ifdef ADCSWAN
         sdFileName = 'PE0000/'//FileName
#endif
         CALL IWRITE(sdFileName, 3, 6, iproc-1)
         OPEN (UNIT=SDU(iproc), FILE=sdFileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .false.
            RETURN ! failed to open at least one subdomain file
         ENDIF
      ENDDO

 2    FORMAT(I2)
 30   FORMAT(A30)
 850  FORMAT(/,'Type ''skip'' to bypass preprocessing or')
 900  FORMAT('Enter the name of the ',A30)
 910  FORMAT('file (unit ',I3,'): ')
 1010 FORMAT('File ',A7,/,' WAS NOT FOUND! Try again or type "skip"',/)
 1011 FORMAT('File ',A7,/,' WAS FOUND!  Opening & Processing file.',/)
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE OpenPrepFiles
C---------------------------------------------------------------------------



c***********************************************************************
c   Subroutine to write out to the hotstart file (UNITS 67 and 68)     *
c   header information and the LHS matrix for the harmonic analysis    *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUT(NP,NSTAE,NSTAV,ISTAE,ISTAV,IGLOE,IGLOV,
c    &  IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NP,NSTAE,NSTAV,ISTAE,AE,ISTAV
c     INTEGER IGLOE,IGLOV,IOUNIT,IHOTSTP,I,J
c     CHARACTER*16 FNAME
c     CHARACTER*8 FNAM8(2)
c     EQUIVALENCE (FNAM8(1),FNAME)c
c
c
c***** Write Out various parameter values
c
c     WRITE(IOUNIT,REC=IHOTSTP+1) NZ
c     WRITE(IOUNIT,REC=IHOTSTP+2) NF
c     WRITE(IOUNIT,REC=IHOTSTP+3) MM
c     WRITE(IOUNIT,REC=IHOTSTP+4) NP
c     WRITE(IOUNIT,REC=IHOTSTP+5) NSTAE
c     WRITE(IOUNIT,REC=IHOTSTP+6) NSTAV
c     WRITE(IOUNIT,REC=IHOTSTP+7) ISTAE
c     WRITE(IOUNIT,REC=IHOTSTP+8) ISTAV
c     WRITE(IOUNIT,REC=IHOTSTP+9) IGLOE
c     WRITE(IOUNIT,REC=IHOTSTP+10) IGLOV
c     WRITE(IOUNIT,REC=IHOTSTP+11) ICALL
c     WRITE(IOUNIT,REC=IHOTSTP+12) NFREQ
c     IHOTSTP = IHOTSTP+12
c
c     do i=1,nfreq+nf
c        FNAME=NAMEFR(I)
c        WRITE(IOUNIT,REC=IHOTSTP+1) FNAM8(1)
c        WRITE(IOUNIT,REC=IHOTSTP+2) FNAM8(2)
c        IHOTSTP=IHOTSTP+2
c        WRITE(IOUNIT,REC=IHOTSTP+1) hafreq(i)
c        WRITE(IOUNIT,REC=IHOTSTP+2) HAFF(i)
c        WRITE(IOUNIT,REC=IHOTSTP+3) HAFACE(i)
c        IHOTSTP=IHOTSTP+3
c     end do
c
c
c***** Write Out time of most recent H.A. update
c
c     WRITE(IOUNIT,REC=IHOTSTP+1) TIMEUD
c     WRITE(IOUNIT,REC=IHOTSTP+2) ITUD
c     IHOTSTP=IHOTSTP+2
c
c***** Write Out LHS Matrix
c
c     do i=1,mm
c        do j=1,mm
c           IHOTSTP = IHOTSTP + 1
c           WRITE(IOUNIT,REC=IHOTSTP) HA(I,J)
c        END DO
c     END DO
c
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write global elevation harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTEG(NP,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER IOUNIT
c     INTEGER NP,IHOTSTP,N,I 
c
c***** Write Out Global Elevation RHS load vector
c
c     do n=1,np
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOELV(I,N)
c        end do
c     end do
c     
c     return
c     end subroutine

c***********************************************************************
c   Subroutine to write elevation station harmonic analysis RHS load   *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTES(NSTAE,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NSTAE,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Elevation RHS load vector
c
c     do n=1,NSTAE
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAELV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write global velocity harmonic analysis RHS load     *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTVG(NP,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NP,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Global Velocity RHS load vector
c
c     do n=1,np
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOULV(I,N)
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOVLV(I,N)
c        end do
c     end do
c     
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write velocity station harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTVS(NSTAV,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NSTAV,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Velocity LHS load vector
c
c     do N=1,NSTAV
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAULV(I,N)
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAVLV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine

CC    Addition by CF   8/2007
      SUBROUTINE PREP80()
      USE PRE_GLOBAL 
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C 
C  This routine writes the domain decomposition information into a file,    C
C  "fort.80".  This file is used by the ADCIRC post-processor ADCPOST.      C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C     jgf45.07 Added subdomain->fulldomain element mapping to handle the    C
C     processing of the NOFF array.
C     jgf45.11 Added IDEN for processing fort.44 file, added additional 
C     arrays to handle processing of 3D recording stations.
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER I,K
C
      OPEN(UNIT=80,FILE='fort.80')              ! output for ADCPOST
C
C--Write out the domain decomposition information into a file
C  which will later be used in post-processing the results
C
      WRITE(80,80) RUNDES
      WRITE(80,80) RUNID
      WRITE(80,80) AGRID
      WRITE(80,'(2I8,16x,A)') NELG,NNODG,'! Total # elements & nodes'
      WRITE(80,'(I8,24x,A)') NPROC,' ! Number of processors'
      WRITE(80,'(I8,24x,A)') MNPP,'! Max nodes on any processor'
      WRITE(80,'(I8,24x,A)') MNEP,'! Max elements on any processor'!jgf45.07
      WRITE(80,'(I8,24x,A)') IM,'! IM, run type'         !jgf46.02
      WRITE(80,'(I8,24x,A)') NWS,'! NWS, wind data type' !jgf46.02
      WRITE(80,'(I8,24x,A)') NSTAE,'! NSTAE'
      WRITE(80,'(I8,24x,A)') NSTAV,'! NSTAV'
      IF (IM.EQ.10) THEN
         WRITE(80,'(I8,24x,A)') NSTAC,' ! NSTAC' !jgf46.02
      ENDIF
      IF (NWS.NE.0) THEN
         WRITE(80,'(I8,24x,A)') NSTAM,'! NSTAM' !jgf46.02
      ENDIF
      WRITE(80,'(I8,24x,A)') MNHARF,'! MNHARF'
      WRITE(80,'(2I8,16x,A)') MNWLAT,MNWLON,'! NWLON, NWLAT'
C
      DO I = 1,NPROC
         WRITE(80,'(3I8,A33)') I-1, NNODP(I), NOD_RES_TOT(I), 
     &        '  ! PE, NNODP(PE), NOD_RES_TOT(PE)'
         WRITE(80,1130) (IMAP_NOD_LG(K,I),K=1,NNODP(I))
      ENDDO
C
      WRITE(80,*) "GLOBAL   PE     LOCAL   ( Global-to-Local Nodes )"
      DO I = 1,NNODG
         WRITE(80,1140) I, IMAP_NOD_GL(1,I)-1, IMAP_NOD_GL(2,I)
      ENDDO
C
C     jgf45.07 Add subdomain->fulldomain mapping to handle NOFF processing
C     IMAP_EL_LG(I,PE) = Global Element Number of Local Element I on PE
      DO I = 1,NPROC
         WRITE(80,'(2I8,A33)') I-1, NELP(I), '  ! PE, NELP(PE)'
         WRITE(80,1130) (IMAP_EL_LG(K,I),K=1,NELP(I))
      ENDDO
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTE,TOUTSE,TOUTFE,NSPOOLE, 
     &    '   ! NOUTE,TOUTSE,TOUTFE,NSPOOLE'
C
      DO I = 1,NPROC
         WRITE(80,*) I,NSTAEP(I)
         DO K = 1,NSTAEP(I)
            WRITE(80,*) IMAP_STAE_LG(K,I)
         ENDDO
      ENDDO
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTV,TOUTSV,TOUTFV,NSPOOLV,
     &    '   ! NOUTV,TOUTSV,TOUTFV,NSPOOLV'
C
      DO I = 1,NPROC
         WRITE(80,*) I,NSTAVP(I)
         DO K = 1,NSTAVP(I)
            WRITE(80,*) IMAP_STAV_LG(K,I)
         ENDDO
      ENDDO
C
      IF (IM.EQ.10) THEN ! jgf46.02
         WRITE(80,'(I8,2E15.8,I8,A32)') NOUTC,TOUTSC,TOUTFC,NSPOOLC,
     &        '   ! NOUTC,TOUTSC,TOUTFC,NSPOOLC'
         DO I = 1,NPROC
            WRITE(80,*) I,NSTACP(I)
            DO K = 1,NSTACP(I)
               WRITE(80,*) IMAP_STAC_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C
      IF (NWS.NE.0) THEN ! jgf46.02
         WRITE(80,'(I8,2E15.8,I8,A32)') NOUTM,TOUTSM,TOUTFM,NSPOOLM,
     &        '   ! NOUTM,TOUTSM,TOUTFM,NSPOOLM'
         DO I = 1,NPROC
            WRITE(80,*) I,NSTAMP(I)
            DO K = 1,NSTAMP(I)
               WRITE(80,*) IMAP_STAM_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGE, TOUTSGE,TOUTFGE,NSPOOLGE,
     &    '   ! NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGV, TOUTSGV,TOUTFGV,NSPOOLGV,
     &    '   ! NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGC, TOUTSGC,TOUTFGC,NSPOOLGC,
     &    '   ! NOUTGC, TOUTSGC, TOUTFGC, NSPOOLGC'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGW, TOUTSGW,TOUTFGW,NSPOOLGW,
     &    '   ! NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW'
C
      WRITE(80,'(4I4,A32)') NHASE,NHASV,NHAGE,NHAGV,
     &    '   ! NHASE, NHASV, NHAGE, NHAGV'

C     -------------------------------------------------------------
C     
C     S T A R T   3 D   D A T A 
C     
C     -------------------------------------------------------------

      WRITE(80,*) IDEN !jgf45.11 needed to post process the fort.44 file

C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D density stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD,
     &       '   ! I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD'
      IF(I3DSD.GT.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DDP(I)
            DO K = 1, NNSTA3DDP(I)
               WRITE(80,*) IMAP_STA3DD_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D velocity stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV,
     &       '   ! I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV'
      IF(I3DSV.GT.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DVP(I)
            DO K = 1, NNSTA3DVP(I)
               WRITE(80,*) IMAP_STA3DV_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D turbulence stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DST, TO3DSTS, TO3DSTF, NSPO3DST, NSTA3DT,
     &       '   ! I3DST, TO3DSST, TO3DFST, NSPO3DST, NSTA3DT'
      IF(I3DST.GT.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DTP(I)
            DO K = 1, NNSTA3DTP(I)
               WRITE(80,*) IMAP_STA3DT_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
      WRITE(80,82) I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD,
     &       '   ! I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD'
      WRITE(80,82) I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV,
     &       '   ! I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV'
      WRITE(80,82) I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT,
     &       '   ! I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT'
C
C     End 3D data
C
      WRITE(80,*) NBYTE
C
      CLOSE(80)
C
  80  FORMAT(A80)
  81  FORMAT(I8,2E15.8,2I8,A32)
  82  FORMAT(I8,2E15.8,I8,A32)
1130  FORMAT(8X,9I8)
1140  FORMAT(8X,3I8)
C
      RETURN
      END


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P   U N S W A N  
C---------------------------------------------------------------------------
C
C     jgf48.17 This subroutine just copies the UnSWAN input file to each
C     subdomain. Currently, it is only called if the compiler directive
C     ADCSWAN is defined. 
C
C---------------------------------------------------------------------------
      SUBROUTINE prepUnSWAN()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      INTEGER I,J,IPROC,IPROC2,ILNODE,INDEX,NHG,LINDEX
      CHARACTER*170 Line ! line of data from UnSWAN file
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      LOGICAL origUseDefault ! there is no default UnSWAN input file name
C
C     Open full domain and all subdomain fort.26 files (UnSWAN)
      origUseDefault = USE_DEFAULT
      USE_DEFAULT = .false. ! abusing the global variable to act as an argument
      CALL OpenPrepFiles(26, 'UnSWAN                        ', 
     &     1, nproc, sdu, success)
      USE_DEFAULT = origUseDefault

      IF (.not.success) THEN
         WRITE(*,*) 'WARNING: UnSWAN files not preprocessed.'
         RETURN ! note early return
      ELSE
         DO 
            READ(26,170,END=9999) Line
            DO IPROC = 1,NPROC
               WRITE(SDU(IPROC),60) trim(Line)
            ENDDO
         ENDDO
      ENDIF
C
C--Close fulldomain file and all the subdomain files     
C
 9999 CLOSE (26)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      RETURN
  60  FORMAT(A60)
 170  FORMAT(A170)
 1100 FORMAT(I8,3E13.5)
 1101 FORMAT(' #')
c----------------------------------------------------------------------------
      END SUBROUTINE prepUnSWAN
c----------------------------------------------------------------------------
      
