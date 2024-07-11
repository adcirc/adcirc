###########################################################################
#                 CMake Build File for GRIB2 (within ADCIRC)
#
#    Written By: Zach Cobell
#
###########################################################################

###########################################################################
#  GENERAL OPTIONS
###########################################################################
cmake_minimum_required(VERSION 3.5)
PROJECT(wgrib2)
###########################################################################

###########################################################################
#  COMPILERS/LANGUAGE
###########################################################################

###########################################################################
#...Language Specifications
ENABLE_LANGUAGE(Fortran)
ENABLE_LANGUAGE(C)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
###########################################################################

###########################################################################
# COMPILER DEPENDENT FLAGS
###########################################################################
IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    SET(IPOLATES_R8_FLAG "-r8")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    SET(IPOLATES_R8_FLAG "-fdefault-real-8 -fdefault-double-8")
ENDIF()
###########################################################################
    
###########################################################################
# JASPER
###########################################################################
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/jasper EXCLUDE_FROM_ALL )
###########################################################################

###########################################################################
# config.h file generation
###########################################################################
SET( USE_REGEX TRUE )
SET( USE_TIGGE FALSE )
SET( DISABLE_ALARM FALSE )
SET( DISABLE_TIMEZONE TRUE )
SET( USE_UDF FALSE )
SET( USE_IPOLATES TRUE )
SET( USE_IPOLATES_TYPE 3 )
SET( IPOLATES_LIB TRUE )
SET( IPOLATES_LIB_STRING "ip2lib_d" )
SET( USE_G2CLIB FALSE )
SET( USE_PROJ4 FALSE )
SET( USE_JASPER TRUE )
SET( USE_AEC FALSE )
SET( USE_MYSQL FALSE )
SET( USE_PNG FALSE )
SET( WMO_VALIDATION FALSE )
IF(NETCDF_WORKING)
    IF(NETCDF4_WORKING)
        SET(USE_NETCDF3 FALSE)
        SET(USE_NETCDF4 TRUE)
    ELSE()
        SET(USE_NETCDF3 TRUE)
        SET(USE_NETCDF4 FALSE)
    ENDIF()
ELSE()
    SET(USE_NETCDF3 FALSE)
    SET(USE_NETCDF4 FALSE)
ENDIF()
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)
###########################################################################

###########################################################################
# JPEG
###########################################################################
FIND_LIBRARY( JPEG TRUE )
###########################################################################

###########################################################################
# libgeo
###########################################################################
ADD_LIBRARY( geo STATIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/gctp.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/alberfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/alberinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/alconfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/alconinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/azimfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/aziminv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/eqconfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/eqconinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/equifor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/equiinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/for_init.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/gnomfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/gnominv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/goodfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/goodinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/gvnspfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/gvnspinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/hamfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/haminv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/imolwfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/imolwinv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/inv_init.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/lamazfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/lamazinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/lamccfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/lamccinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/merfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/merinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/millfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/millinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/molwfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/molwinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/obleqfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/obleqinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/omerfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/omerinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/orthfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/orthinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/paksz.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/polyinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/polyfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/psinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/psfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/robinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/robfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sininv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sinfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sominv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/somfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sphdz.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sterinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/sterfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/stplninv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/stplnfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/tminv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/tmfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/utminv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/utmfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/untfz.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/vandginv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/vandgfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/wivinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/wivfor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/wviiinv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/wviifor.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/cproj.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source/report.c )
###########################################################################

###########################################################################
# IPOLATES 
###########################################################################
ADD_LIBRARY( ipolates STATIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/earth_radius.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_c.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_equid_cylind_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_gaussian_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_lambert_conf_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_mercator_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_polar_stereo_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_rot_equid_cylind_egrid_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/gdswzd_rot_equid_cylind_mod.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ijkgds0.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ijkgds1.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipolates.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipolatev.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipxetas.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipxwafs2.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipxwafs3.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/ipxwafs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/movect.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates0.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates1.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates2.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates3.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polates6.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev0.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev1.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev2.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev3.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev4.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polatev6.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polfixs.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/polfixv.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ip2lib_d/splat.f90 )
TARGET_INCLUDE_DIRECTORIES( ipolates PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} )
SET_TARGET_PROPERTIES(ipolates PROPERTIES COMPILE_FLAGS ${IPOLATES_R8_FLAG})
SET_TARGET_PROPERTIES(ipolates PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/ipolates )
###########################################################################

###########################################################################
# libgrib2 
###########################################################################
ADD_LIBRARY( grib2 STATIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/AAIG.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/AAIGlong.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/addtime.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/aec_pk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Aerosol.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Alarm.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Alias.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ave_test.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ave_var.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Bbox.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/bitstream.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Box_ave.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Checksum.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/cname.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/codetable_4_230.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/CodeTable.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Code_Values.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/complex_pk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Config.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/copy.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/crc32.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Cress_lola.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Csv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Csv_long.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Cyclic.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Data.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/dec_png_clone.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Dump.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Earth.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/enc_jpeg2000_clone.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/End.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ensemble.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ens_processing.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/EOF.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ExtName.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/fatal_error.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fcst_ave.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ffopen.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fi.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/File.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_CFSv2_fcst.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_ncep_2.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_ncep_3.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_ncep_4.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_ncep.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/FlagTable.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/flt2ieee.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/flt2ieee_nan.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Flush.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/fnlist.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/For.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ftime1.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ftime2.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ftime.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ftn_api.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Gctpc.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc_ll2xy.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/GDT.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/geo.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grads.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib_out_irr2.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib_out_irr.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gribtab.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grid_def.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grid_id.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/grid_ident.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Gridout.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Help.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Hybrid.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee2flt.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee2flt_nan.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee_pk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/If_rec.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_netcdf.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/init.c wgrib2/int8.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/intpower.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Inv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Inv_no.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Irr_grids.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/itoshort_a.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/JMA.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/jpeg_pk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Last.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/lat2ij.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Latlon.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Level.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Limit.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Lola.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Macro.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/manage_inv_out.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Match.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Match_fs.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mem_buffer.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Merge.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Misc.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/missing.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/mk_gdt.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/mk_kgds.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Model_version_date.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mod_grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mysql.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mysql_dump.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mysql_speed.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ncep_grids.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/NCEP_norm.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/NCEP_uv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ncpu.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ndate.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ndates.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Netcdf.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Netcdf_sup.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/New_grid.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/new_grid_lambertc.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/openmp_util.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/parse_loop.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/parse_msg.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Precision.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Prob.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/prod_def_temp_size.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Proj4.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/proj4_initialize.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/proj4_transform.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Range.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rdieee_file.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rd_inventory.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rd_seq_grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Reduced_gaussian_grid.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rotll.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/RPN.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/save_string.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Scaling_0001.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Scan.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec0.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec1.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec3.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/sec3_grids.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec4.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec5.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec6.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Sec_len.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_byte.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_date.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ensm_derived_fcst.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ens_num.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ftime1.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ftime2.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_metadata.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_pdt.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_percentile.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_prob.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_sec_size.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ts_dates.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/setup_user_gribtable.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_val.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/simple_pk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Small_grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/space_view2ij.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/space_view.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/SpatialProcessing.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Spectral_bands.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/stagger.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Status.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Submsg_uv.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/subtime.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Summary.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/swap_buffer.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/test_sec.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Tigge.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Time_processing.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Tosubmsg.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/UDF.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Undefine.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/units.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Unix_time.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_0.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_complex.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_run_length.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/update_sec4.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/VerfTime.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Warn_old_g2lib.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Waves.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wgrib2_api.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wind_dir.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wind_speed.c
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wrt_grib.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wrtieee.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wxtext.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/ftn_api/fort_wgrib2.c 
    ${CMAKE_CURRENT_SOURCE_DIR}/ftn_api/wgrib2api.f90  
    ${CMAKE_CURRENT_SOURCE_DIR}/ftn_api/wgrib2lowapi.f90 
    ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wgrib2.c ) 

ADD_DEPENDENCIES( grib2 ipolates geo libjasper ) 
TARGET_INCLUDE_DIRECTORIES( grib2 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/ipolates 
                                          ${CMAKE_CURRENT_SOURCE_DIR}/gctpc/source
                                          ${CMAKE_CURRENT_SOURCE_DIR} )
SET_TARGET_PROPERTIES( grib2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/grib2 )
TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE CALLABLE_WGRIB2 USE_IPOLATES=3 )
IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE "IFORT" )
ELSEIF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE "GFORTRAN" )
ENDIF()
TARGET_LINK_LIBRARIES( grib2 PRIVATE ipolates geo libjasper ${JPEG_LIBRARIES} )
TARGET_INCLUDE_DIRECTORIES( grib2 PRIVATE ${JPEG_INCLUDE_DIR} )
TARGET_INCLUDE_DIRECTORIES( grib2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/jasper/src/libjasper/include )
TARGET_INCLUDE_DIRECTORIES( grib2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jasper/src/libjasper/include )
IF(NETCDF_WORKING)
    TARGET_INCLUDE_DIRECTORIES(grib2 PRIVATE ${NETCDF_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(grib2 PRIVATE ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})
ENDIF()
###########################################################################
