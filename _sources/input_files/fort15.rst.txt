.. _fort15:

Fort.15: Model Parameter and Periodic Boundary Condition File
=============================================================

The fort.15 file contains the majority of the parameters required to run both the 2DDI and 3D versions of ADCIRC and the information to drive the model with harmonic boundary conditions (either elevation or flux). This file is required to run the ADCIRC model.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input. Conditional input is indicated by a clause following the variable name(s).

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`
   :ref:`RUNID <RUNID>`
   :ref:`NFOVER <NFOVER>`
   :ref:`NABOUT <NABOUT>`
   :ref:`NSCREEN <NSCREEN>`
   :ref:`IHOT <IHOT>`
   :ref:`ICS <ICS>`
   :ref:`IM <IM>`
   :ref:`IDEN <IDEN>`      if IM = 21
   :ref:`NOLIBF <NOLIBF>`
   :ref:`NOLIFA <NOLIFA>`
   :ref:`NOLICA <NOLICA>`
   :ref:`NOLICAT <NOLICAT>`
   :ref:`NWP <NWP>`
   for j=1 to :ref:`NWP <NWP>`
      :ref:`AttrName(j) <AttrName>`
   end j loop
   :ref:`NCOR <NCOR>`
   :ref:`NTIP <NTIP>`
   :ref:`NWS <NWS>`
   :ref:`NRAMP <NRAMP>`
   :ref:`G <G>`
   :ref:`TAU0 <TAU0>`
   :ref:`Tau0FullDomainMin <Tau0FullDomainMin>`, :ref:`Tau0FullDomainMax <Tau0FullDomainMax>`      if TAU0 = -5.0
   :ref:`DTDP <DTDP>`
   :ref:`STATIM <STATIM>`
   :ref:`REFTIM <REFTIM>`
   :ref:`WTIMINC <WTIMINC>`
   :ref:`RNDAY <RNDAY>`
   :ref:`DRAMP <DRAMP>`      if NRAMP = 0 or 1
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`      if NRAMP = 2
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`      if NRAMP = 3
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`, :ref:`DRAMPElev <DRAMPElev>`      if NRAMP = 4
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`, :ref:`DRAMPElev <DRAMPElev>`, :ref:`DRAMPTip <DRAMPTip>`      if NRAMP = 5
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`, :ref:`DRAMPElev <DRAMPElev>`, :ref:`DRAMPTip <DRAMPTip>`, :ref:`DRAMPMete <DRAMPMete>`      if NRAMP = 6
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`, :ref:`DRAMPElev <DRAMPElev>`, :ref:`DRAMPTip <DRAMPTip>`, :ref:`DRAMPMete <DRAMPMete>`, :ref:`DRAMPWRad <DRAMPWRad>`      if NRAMP = 7
   :ref:`DRAMP <DRAMP>`, :ref:`DRAMPExtFlux <DRAMPExtFlux>`, :ref:`FluxSettlingTime <FluxSettlingTime>`, :ref:`DRAMPIntFlux <DRAMPIntFlux>`, :ref:`DRAMPElev <DRAMPElev>`, :ref:`DRAMPTip <DRAMPTip>`, :ref:`DRAMPMete <DRAMPMete>`, :ref:`DRAMPWRad <DRAMPWRad>`, :ref:`DUnRampMete <DUnRampMete>`      if NRAMP = 8
   :ref:`A00 <A00>`, :ref:`B00 <B00>`, :ref:`C00 <C00>`
   :ref:`H0 <H0>`      if NOLIFA = 0, 1
   :ref:`H0 <H0>` :ref:`INTEGER <INTEGER>` :ref:`INTEGER <INTEGER>` :ref:`VELMIN <VELMIN>`      if NOLIFA = 2 or 3
   :ref:`SLAM0 <SLAM0>`, :ref:`SFEA0 <SFEA0>`
   :ref:`TAU <TAU>`      if NOLIBF = 0
   :ref:`CF <CF>`      if NOLIBF = 1
   :ref:`CF <CF>`, :ref:`HBREAK <HBREAK>`, :ref:`FTHETA <FTHETA>`, :ref:`FGAMMA <FGAMMA>`      if NOLIBF = 2
   :ref:`ESLM <ESLM>`      if IM = 0, 1, 2
   :ref:`ESLM <ESLM>`, :ref:`ESLC <ESLC>`      if IM = 10
   :ref:`CORI <CORI>`      if NCOR = 1
   :ref:`NTIF <NTIF>`
   for k=1 to :ref:`NTIF <NTIF>`
      :ref:`TIPOTAG(k) <TIPOTAG>`
      :ref:`TPK(k) <TPK>`, :ref:`AMIGT(k) <AMIGT>`, :ref:`ETRF(k) <ETRF>`, :ref:`FFT(k) <FFT>`, :ref:`FACET(k) <FACET>`
   end k loop
   :ref:`NBFR <NBFR>`
   for k=1 to :ref:`NBFR <NBFR>`
      :ref:`BOUNTAG(k) <BOUNTAG>`
      :ref:`AMIG(k) <AMIG>`, :ref:`FF(k) <FF>`, :ref:`FACE(k) <FACE>`
   end k loop
   for k=1 to :ref:`NBFR <NBFR>`
      :ref:`ALPHAE(k) <ALPHAE>`
      for j=1 to :ref:`NETA <NETA>`
         :ref:`EMO(k,j) <EMO>`, :ref:`EFA(k,j) <EFA>`
      end j loop
   end k loop
   :ref:`ANGINN <ANGINN>`
   :ref:`NFFR <NFFR>`      include this line only if IBTYPE = 2, 12, 22, 32 or 52 in the Grid and Boundary Information File
   for k=1 to :ref:`NFFR <NFFR>`
      :ref:`FBOUNTAG(k) <FBOUNTAG>`
      :ref:`FAMIGT(k) <FAMIGT>`, :ref:`FFF(k) <FFF>`, :ref:`FFACE(k) <FFACE>`
   end k loop
   for k=1 to :ref:`NFFR <NFFR>`
      :ref:`ALPHAQ(k) <ALPHAQ>`
      for j=1 to :ref:`NVEL <NVEL>`
         :ref:`QNAM(k,j) <QNAM>`, :ref:`QNPH(k,j) <QNPH>`      use this line if IBTYPE = 2, 12, 22 in the Grid and Boundary Information File
         :ref:`QNAM(k,j) <QNAM>`, :ref:`QNPH(k,j) <QNPH>`, :ref:`ENAM(k,j) <ENAM>`, :ref:`ENPH(k,j) <ENPH>`      use this line if IBTYPE = 32 in the Grid and Boundary Information File
      end j loop
   end k loop
   :ref:`NOUTE <NOUTE>`, :ref:`TOUTSE <TOUTSE>`, :ref:`TOUTFE <TOUTFE>`, :ref:`NSPOOLE <NSPOOLE>`
   :ref:`NSTAE <NSTAE>`
   for k=1 to :ref:`NSTAE <NSTAE>`
      :ref:`XEL(k) <XEL>`, :ref:`YEL(k) <YEL>`
   end k loop
   :ref:`NOUTV <NOUTV>`, :ref:`TOUTSV <TOUTSV>`, :ref:`TOUTFV <TOUTFV>`, :ref:`NSPOOLV <NSPOOLV>`
   :ref:`NSTAV <NSTAV>`
   for k=1 to :ref:`NSTAV <NSTAV>`
      :ref:`XEV(k) <XEV>`, :ref:`YEV(k) <YEV>`
   end k loop
   :ref:`NOUTC <NOUTC>`, :ref:`TOUTSC <TOUTSC>`, :ref:`TOUTFC <TOUTFC>`, :ref:`NSPOOLC <NSPOOLC>`      include this line only if IM =10
   :ref:`NSTAC <NSTAC>`      include this line only if IM =10
   for k=1 to :ref:`NSTAC <NSTAC>`
      :ref:`XEC(k) <XEC>`, :ref:`YEC(k) <YEC>`
   end k loop
   :ref:`NOUTM <NOUTM>`, :ref:`TOUTSM <TOUTSM>`, :ref:`TOUTFM <TOUTFM>`, :ref:`NSPOOLM <NSPOOLM>`      include this line only if NWS is not equal to zero
   :ref:`NSTAM <NSTAM>`      include this line only if NWS is not equal to zero
   for k=1 to :ref:`NSTAM <NSTAM>`
      :ref:`XEM(k) <XEM>`, :ref:`YEM(k) <YEM>`
   end k loop
   :ref:`NOUTGE <NOUTGE>`, :ref:`TOUTSGE <TOUTSGE>`, :ref:`TOUTFGE <TOUTFGE>`, :ref:`NSPOOLGE <NSPOOLGE>`
   :ref:`NOUTGV <NOUTGV>`, :ref:`TOUTSGV <TOUTSGV>`, :ref:`TOUTFGV <TOUTFGV>`, :ref:`NSPOOLGV <NSPOOLGV>`
   :ref:`NOUTGC <NOUTGC>`, :ref:`TOUTSGC <TOUTSGC>`, :ref:`TOUTFGC <TOUTFGC>`, :ref:`NSPOOLGC <NSPOOLGC>`      include this line only if IM =10
   :ref:`NOUTGW <NOUTGW>`, :ref:`TOUTSGW <TOUTSGW>`, :ref:`TOUTFGW <TOUTFGW>`, :ref:`NSPOOLGW <NSPOOLGW>`      include this line only if NWS is not equal to zero
   :ref:`NFREQ <NFREQ>`
   for k=1 to :ref:`NFREQ <NFREQ>`
   :ref:`NAMEFR(k) <NAMEFR>`
   :ref:`HAFREQ(k) <HAFREQ>`, :ref:`HAFF(k) <HAFF>`, :ref:`HAFACE(k) <HAFACE>`
   end k loop
   :ref:`THAS <THAS>`, :ref:`THAF <THAF>`, :ref:`NHAINC <NHAINC>`, :ref:`FMV <FMV>`
   :ref:`NHASE <NHASE>`, :ref:`NHASV <NHASV>`, :ref:`NHAGE <NHAGE>`, :ref:`NHAGV <NHAGV>`
   :ref:`NHSTAR <NHSTAR>`, :ref:`NHSINC <NHSINC>`
   :ref:`ITITER <ITITER>`, :ref:`ISLDIA <ISLDIA>`, :ref:`CONVCR <CONVCR>`, :ref:`ITMAX <ITMAX>`
   For a 2DDI ADCIRC run that does not use NetCDF, the file ends here. For any ADCIRC run that uses NetCDF, the lines NCPROJ through NCDATE (described at the end of this file format) are required metadata and must be added at the end of the fort.15 file.
   The following information is only included for a 3D run:
   :ref:`IDEN <IDEN>`
   :ref:`ISLIP <ISLIP>`, :ref:`KP <KP>`
   :ref:`Z0S <Z0S>`, :ref:`Z0B <Z0B>`
   :ref:`ALP1 <ALP1>`, :ref:`ALP2 <ALP2>`, :ref:`ALP3 <ALP3>`
   :ref:`IGC <IGC>`, :ref:`NFEN <NFEN>`
   for k=1 to :ref:`NFEN <NFEN>`      include this loop only if IGC = 0, k=1 at bottom, k= NFEN at surface
      :ref:`SIGMA(k) <SIGMA>`
   end k loop
   :ref:`IEVC <IEVC>`, :ref:`EVMIN <EVMIN>`, :ref:`EVCON <EVCON>`
   for k=1 to :ref:`NFEN <NFEN>`      include this loop only if IEVC = 0, k=1 at bottom, k= NFEN at surface
      :ref:`EVTOT(k) <EVTOT>`
   end k loop
   :ref:`THETA1 <THETA1>`, :ref:`THETA2 <THETA2>`      include this line only if IEVC = 50 or 51
   :ref:`I3DSD <I3DSD>`, :ref:`TO3DSDS <TO3DSDS>`, :ref:`TO3DSDF <TO3DSDF>`, :ref:`NSPO3DSD <NSPO3DSD>`
   :ref:`NSTA3DD <NSTA3DD>`
   for k=1 to :ref:`NSTA3DD <NSTA3DD>`
      :ref:`X3DS(k) <X3DS>`, :ref:`Y3DS(k) <Y3DS>`
   end k loop
   :ref:`I3DSV <I3DSV>`, :ref:`TO3DSVS <TO3DSVS>`, :ref:`TO3DFVF <TO3DFVF>`, :ref:`NSPO3DSV <NSPO3DSV>`
   :ref:`NSTA3DV <NSTA3DV>`
   for k=1 to :ref:`NSTA3DV <NSTA3DV>`
      :ref:`X3DS(k) <X3DS>`, :ref:`Y3DS(k) <Y3DS>`
   end k loop
   :ref:`I3DST <I3DST>`, :ref:`TO3DSTS <TO3DSTS>`, :ref:`TO3DSTF <TO3DSTF>`, :ref:`NSPO3DST <NSPO3DST>`
   :ref:`NSTA3DT <NSTA3DT>`
   for k=1 to :ref:`NSTA3DT <NSTA3DT>`
      :ref:`X3DS(k) <X3DS>`, :ref:`Y3DS(k) <Y3DS>`
   end k loop
   :ref:`I3DGD <I3DGD>`, :ref:`TO3DGDS <TO3DGDS>`, :ref:`TO3DGDF <TO3DGDF>`, :ref:`NSPO3DGD <NSPO3DGD>`
   :ref:`I3DGV <I3DGV>`, :ref:`TO3DGVS <TO3DGVS>`, :ref:`TO3DGVF <TO3DGVF>`, :ref:`NSPO3DGV <NSPO3DGV>`
   :ref:`I3DGT <I3DGT>`, :ref:`TO3DGTS <TO3DGTS>`, :ref:`TO3DGTF <TO3DGTF>`, :ref:`NSPO3DGT <NSPO3DGT>`
   The following line will be read in if IM is 21 or 31.
   :ref:`RES_BC_FLAG <RES_BC_FLAG>`, :ref:`BCFLAG_LNM <BCFLAG_LNM>`, :ref:`BCFLAG_TEMP <BCFLAG_TEMP>`
   The following two lines will be read in if RES_BC_FLAG is negative.
   :ref:`RBCTIMEINC <RBCTIMEINC>`
   :ref:`BCSTATIM <BCSTATIM>`
   The following two lines will be read in if RES_BC_FLAG = 2.
   :ref:`RBCTIMEINC <RBCTIMEINC>`, :ref:`SBCTIMEINC <SBCTIMEINC>`
   :ref:`BCSTATIM <BCSTATIM>`, :ref:`SBCSTATIM <SBCSTATIM>`
   The following two lines will be read in if RES_BC_FLAG = 3.
   :ref:`RBCTIMEINC <RBCTIMEINC>`, :ref:`TBCTIMEINC <TBCTIMEINC>`
   :ref:`BCSTATIM <BCSTATIM>`, :ref:`TBCSTATIM <TBCSTATIM>`
   The following two lines will be read in if RES_BC_FLAG = 4.
   :ref:`RBCTIMEINC <RBCTIMEINC>`, :ref:`SBCTIMEINC <SBCTIMEINC>`, :ref:`TBCTIMEINC <TBCTIMEINC>`
   :ref:`BCSTATIM <BCSTATIM>`, :ref:`SBCSTATIM <SBCSTATIM>`, :ref:`TBCSTATIM <TBCSTATIM>`
   The following two lines will be read in if RES_BC_FLAG = 3 or 4 and BCFLAG_TEMP is not equal to 0.
   :ref:`TTBCTIMEINC <TTBCTIMEINC>`, :ref:`TTBCSTATIM <TTBCSTATIM>`
   :ref:`TTBCTIMEINC <TTBCTIMEINC>`
   The following two lines will be read in only if IM is 21 or 31.
   :ref:`SPONGEDIST <SPONGEDIST>`
   :ref:`EQNSTATE <EQNSTATE>`
   The following lines will be read in only if IDEN is > 0.
   :ref:`NLSD <NLSD>`, :ref:`NVSD <NVSD>`
   :ref:`NLTD <NLTD>`, :ref:`NVTD <NVTD>`
   :ref:`ALP4 <ALP4>`
   The following line will be read in only if IDEN = 3 or 4.
   :ref:`NTF <NTF>`
   The following lines will be read in only if the NetCDF output or hotstart format is chosen
   :ref:`NCPROJ <NCPROJ>`
   :ref:`NCINST <NCINST>`
   :ref:`NCSOUR <NCSOUR>`
   :ref:`NCHIST <NCHIST>`
   :ref:`NCREF <NCREF>`
   :ref:`NCCOM <NCCOM>`
   :ref:`NCHOST <NCHOST>`
   :ref:`NCCONV <NCCONV>`
   :ref:`NCCONT <NCCONT>`
   :ref:`NCDATE <NCDATE>`
   The following Fortran namelist lines are optional, but if they appear, they must appear at the very end of the fort.15 file.
   :ref:`metControl <metControl>` :ref:`WindDragLimit`=floatValue, :ref:`DragLawString`='stringValue', :ref:`rhoAir`=floatValue 
   :ref:`timeBathyControl <timeBathyControl>` NDDT=integerValue, BTIMINC=floatValue, BCHGTIMINC=floatValue 
   :ref:`waveCoupling <waveCoupling>` WindWaveMultiplier=floatValue 
   :ref:`SWANOutputControl <SWANOutputControl>` SWAN_OutputHS=logicalValue, SWAN_OutputDIR=logicalValue, SWAN_OutputTM01=logicalValue, SWAN_OutputTPS=logicalValue, SWAN_OutputWIND=logicalValue, SWAN_OutputTM02=logicalValue, SWAN_OutputTMM10=logicalValue 
   :ref:`subdomainModeling <subdomainModeling>` subdomainOn=logicalValue
   :ref:`wetDryControl <wetDryControl>` outputNodeCode=logicalValue, outputNOFF=logicalValue, noffActive=logicalValue 
   :ref:`inundationOutputControl <inundationOutputControl>` inundationOutput=logicalValue0, inunThresh =floatValue 
   :ref:`TVWControl <TVWControl>` use_TVW=logicalValue, TVW_file='stringValue', nout_TVW =integerValue, touts_TVW =floatValue, toutf_TVW=floatValue, nspool_TVW =integerValue

Example
-------

The following is a simple example of the beginning of a fort.15 file:

.. code-block:: none

   Sample ADCIRC Run
   Sample Run ID
   1                ! NFOVER
   0                ! NABOUT
   1                ! NSCREEN
   0                ! IHOT
   2                ! ICS
   1                ! IM
   1                ! NOLIBF
   1                ! NOLIFA
   1                ! NOLICA
   1                ! NOLICAT
   3                ! NWP
   mannings_n_at_sea_floor
   surface_directional_effective_roughness_length
   average_horizontal_eddy_viscosity_in_sea_water_wrt_depth
   1                ! NCOR
   1                ! NTIP
   0                ! NWS
   1                ! NRAMP
   9.81             ! G
   0.0001           ! TAU0
   600.0            ! DTDP
   0.0              ! STATIM
   0.0              ! REFTIM
   60.0 30.0 3      ! WTIMINC
   5.0              ! RNDAY
   43200.0          ! DRAMP
   0.0 0.0 0.0      ! A00, B00, C00
   0.1              ! H0 