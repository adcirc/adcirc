VPATH = $(top_srcdir)/autotools/lib/libswan_serial_source

AM_FFLAGS = -I../libswan_serial1 -I../libadcswan
AM_FCFLAGS = -I../libswan_serial1 -I../libadcswan
ARFLAGS = cr

if BUILD_SWAN

noinst_LIBRARIES = libswan_serial2.a

libswan_serial2_a_SOURCES = \
      swanmain.f swanpre1.f swanpre2.f swancom1.f swancom2.f swancom3.f \
      swancom4.f swancom5.f swanout1.f swanout2.f swanser.f swanparll.f \
      SwanReadGrid.f90 SwanReadADCGrid.f90 SwanReadTriangleGrid.f90 \
      SwanReadEasymeshGrid.f90 SwanInitCompGrid.f90 SwanCheckGrid.f90 \
      SwanCreateEdges.f90 SwanGridTopology.f90 SwanGridVert.f90 \
      SwanGridCell.f90 SwanGridFace.f90 SwanPrintGridInfo.f90 SwanFindPoint.f90 \
      SwanPointinMesh.f90 SwanBpntlist.f90 SwanPrepComp.f90 SwanVertlist.f90 \
      SwanCompUnstruc.f90 SwanDispParm.f90 SwanPropvelX.f90 SwanSweepSel.f90 \
      SwanPropvelS.f90 SwanTranspAc.f90 SwanTranspX.f90 SwanDiffPar.f90 \
      SwanGSECorr.f90 SwanGradDepthorK.f90 SwanGradVel.f90 SwanInterpolatePoint.f90 \
      SwanInterpolateAc.f90 SwanInterpolateOutput.f90 SwanConvAccur.f90 SwanConvStopc.f90 \
      SwanFindObstacles.f90 SwanCrossObstacle.f90 SwanComputeForce.f90 SwanIntgratSpc.f90 \
      SwanBndStruc.f90 SwanReadfort18.f90 SwanPunCollect.f90 SwanSumOverNodes.f90 \
	  SwanMinOverNodes.f90 SwanMaxOverNodes.f90 SwanVTKWriteHeader.f90 SwanVTKWriteData.f90 \
      SwanVTKPDataSets.f90 ocpids.f ocpcre.f ocpmix.f fftpack51.f90

endif
