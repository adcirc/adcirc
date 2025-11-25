###########################################################################
#                 CMake Build File for GRIB2 (within ADCIRC)
#
#    Written By: Zach Cobell
#
###########################################################################

###########################################################################
#  GENERAL OPTIONS
###########################################################################
cmake_minimum_required(VERSION 3.14)
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
# JASPER
###########################################################################
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/jasper EXCLUDE_FROM_ALL )
###########################################################################

###########################################################################
# config.h file generation
###########################################################################
SET( USE_NETCDF FALSE)
SET( USE_REGEX TRUE)
SET( USE_TIGGE FALSE)
SET( USE_MYSQL FALSE)
SET( USE_IPOLATES FALSE)
SET( USE_UDF FALSE)
SET( USE_OPENMP FALSE)
SET( USE_PROJ4 FALSE)
SET( USE_WMO_VALIDATION FALSE)
SET( DISABLE_TIMEZONE TRUE)
SET( DISABLE_ALARM FALSE)
SET( USE_NAMES NCEP)
SET( USE_G2CLIB_HIGH FALSE)
SET( USE_G2CLIB_LOW TRUE)
SET( USE_AEC FALSE)
SET( BUILD_EXTRA FALSE)
SET( FTP_TEST_FILES FALSE)
SET( FTP_EXTRA_TEST_FILES FALSE)
SET( MAKE_FTN_API TRUE)
SET( DISABLE_STAT TRUE)
SET( BUILD_LIB TRUE)
SET( BUILD_SHARED_LIB FALSE)
SET( BUILD_WGRIB FALSE)
SET( TEST_FILE_DIR ".")
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)
###########################################################################

###########################################################################
# JPEG
###########################################################################
FIND_LIBRARY( JPEG TRUE )
###########################################################################

###########################################################################
# libgeo - Geographic projection library
###########################################################################
add_library(geo STATIC
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/alberfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/alberinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/alconfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/alconinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/azimfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/aziminv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/br_gctp.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/cproj.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/eqconfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/eqconinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/equifor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/equiinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/for_init.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/gctp.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/gnomfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/gnominv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/goodfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/goodinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/gvnspfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/gvnspinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/hamfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/haminv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/imolwfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/imolwinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/inv_init.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/lamazfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/lamazinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/lamccfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/lamccinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/merfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/merinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/millfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/millinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/molwfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/molwinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/obleqfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/obleqinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/omerfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/omerinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/orthfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/orthinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/paksz.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/polyfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/polyinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/psfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/psinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/report.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/robfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/robinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sinfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sininv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/somfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sominv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sphdz.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sterfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/sterinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/stplnfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/stplninv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/tmfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/tminv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/untfz.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/utmfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/utminv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/vandgfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/vandginv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/wivfor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/wivinv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/wviifor.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/wviiinv.c)
###########################################################################

###########################################################################
# G2C 
###########################################################################
ADD_LIBRARY(g2c STATIC 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/grib2_int.h 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/cmplxpack.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/compack.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/comunpack.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/drstemplates.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_addfield.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_addgrid.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_addlocal.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_create.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_free.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_getfld.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_gribend.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_info.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_miss.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack1.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack2.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack3.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack4.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack5.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack6.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/g2_unpack7.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/gbits.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/getdim.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/getpoly.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/gridtemplates.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/int_power.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/misspack.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/mkieee.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/pack_gp.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/pdstemplates.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/rdieee.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/reduce.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/seekgb.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/simpack.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/simunpack.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/specpack.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/specunpack.c 
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/util.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/decenc_jpeg2000.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/jpcpack.c
            ${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/jpcunpack.c
          )
target_compile_definitions(g2c PRIVATE _LARGEFILE64_SOURCE JPEG USE_JPEG2000)
target_link_libraries(g2c PUBLIC libjasper)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/g2c/src/grib2.h.in" "${CMAKE_CURRENT_BINARY_DIR}/grib2.h" @ONLY)
target_include_directories(g2c PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                                       ${CMAKE_CURRENT_SOURCE_DIR}/jasper/src/libjasper/include
                                       ${CMAKE_CURRENT_BINARY_DIR}/jasper/src/libjasper/include)
add_dependencies(g2c libjasper)
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Check_pdt_size.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Checksum.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Cluster.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/cname.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/codetable_4_230.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/CodeTable.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Code_Values.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Code_Values_JMA.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/complex_pk.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Config.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/copy.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/crc32.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Cress_lola.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Csv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Csv_long.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/cubed_sphere2ll.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/CubeFace2global.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Cyclic.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Data.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Dump.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Earth.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Else.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/End.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Endif.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ensemble.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ens_processing.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Ens_qc.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/EOF.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Export_lonlat.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Fix_undef.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gaussian2ij.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Gctpc.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc_ll2xy.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/GDT.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/geo.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Geolocation.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/get_unixtime.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grads.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib_out_irr2.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grib_out_irr.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gribtab.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Gribtable_used.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grid_def.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Grid_id.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/grid_ident.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Gridout.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Help.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Hybrid.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee2flt.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee2flt_nan.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ieee_pk.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/If_delayed_error.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/If_rec.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_grib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_grib_fs.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_lonlat.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Import_netcdf.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/init.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/int8.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Lvl.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Macro.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/manage_inv_out.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Match.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Match_fs.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Match_inv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mem_buffer.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Merge.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Misc.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/missing.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/mk_gdt.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Model_version_date.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Mod_grib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Names.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/New_grid_order.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/openmp_util.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/parse_loop.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/parse_msg.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/pdt_len.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Precision.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Prob.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Proj4.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/proj4_initialize.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/proj4_transform.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/PyInv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Range.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rdieee_file.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rd_inventory.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/rd_seq_grib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/read_latlon.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Reduced_gaussian_grid.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Reset_delayed_error.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_ftime_mode.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Set_gds.c
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
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/string2time_unit.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Submsg_uv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/subtime.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Summary.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/swap_buffer.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/test_sec.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Tigge.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Time_processing.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Timer.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Tosubmsg.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Type_reftime.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/UDF.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Undefine.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/units.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Unix_time.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Unmerge_fcst.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_0.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_complex.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/unpk_run_length.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/update_sec3.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/update_sec4.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/v1_v2_if.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/VerfTime.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Warn_old_g2lib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Waves.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wgrib2_api.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wgrib2.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wgrib2_main.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wind_dir.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wind_speed.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wind_uv.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Write_sec.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/Wrt_grib.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wrtieee.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/wxtext.c
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ftn_api/fort_wgrib2.c 
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ftn_api/wgrib2api.f90  
             ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/ftn_api/wgrib2lowapi.f90)

ADD_DEPENDENCIES( grib2 geo libjasper g2c )
TARGET_INCLUDE_DIRECTORIES( geo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/include)
TARGET_INCLUDE_DIRECTORIES( grib2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                                          ${CMAKE_CURRENT_SOURCE_DIR}/wgrib2/gctpc/source/include
                                          ${CMAKE_CURRENT_SOURCE_DIR}
                                          ${CMAKE_BINARY_DIR}/CMakeFiles/mod/grib2
                                          ${JPEG_INCLUDE_DIR}
                                          ${CMAKE_CURRENT_BINARY_DIR}/jasper/src/libjasper/include
                                          ${CMAKE_CURRENT_SOURCE_DIR}/jasper/src/libjasper/include )
SET_TARGET_PROPERTIES( grib2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/grib2 )
TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE CALLABLE_WGRIB2 G2_JPEG2000_ENABLED)
IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
    TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE "IFORT" )
ELSEIF(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    TARGET_COMPILE_DEFINITIONS( grib2 PRIVATE "GFORTRAN" )
ENDIF()
TARGET_LINK_LIBRARIES( grib2 PUBLIC geo g2c libjasper ${JPEG_LIBRARIES} )
IF(NETCDF_WORKING)
    TARGET_INCLUDE_DIRECTORIES(grib2 PRIVATE ${NETCDF_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(grib2 PRIVATE ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})
ENDIF()
###########################################################################
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
    set(ADDITIONAL_WGRIB2_COMPILER_FLAGS "-Wno-incompatible-pointer-types;-Wno-stringop-overflow")
    set(ADDITIONAL_GCTPC_COMPILER_FLAGS "-Wno-old-style-definition")
    message(STATUS "Adding additional compiler flags to gctpc: ${ADDITIONAL_GCTPC_COMPILER_FLAGS}")
    message(STATUS "Adding additional compiler flags to wgrib2: ${ADDITIONAL_WGRIB2_COMPILER_FLAGS}")
    target_compile_options(geo PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_GCTPC_COMPILER_FLAGS}>)
    target_compile_options(grib2 PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_WGRIB2_COMPILER_FLAGS}>)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "IntelLLVM")
    set(ADDITIONAL_G2C_COMPILER_FLAGS "-Wno-absolute-value")
    set(ADDITIONAL_WGRIB2_COMPILER_FLAGS "-Wno-deprecated-non-prototype;-Wno-fortify-source;-Wno-enum-conversion;-Wno-incompatible-pointer-types-discards-qualifiers")
    set(ADDITIONAL_GCTPC_COMPILER_FLAGS "-Wno-deprecated-non-prototype;-Wno-absolute-value;-Wno-return-type")
    message(STATUS "Adding additional compiler flags to g2c: ${ADDITIONAL_G2C_COMPILER_FLAGS}")
    message(STATUS "Adding additional compiler flags to gctpc: ${ADDITIONAL_GCTPC_COMPILER_FLAGS}")
    message(STATUS "Adding additional compiler flags to wgrib2: ${ADDITIONAL_WGRIB2_COMPILER_FLAGS}")
    target_compile_options(g2c PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_G2C_COMPILER_FLAGS}>)
    target_compile_options(geo PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_GCTPC_COMPILER_FLAGS}>)
    target_compile_options(grib2 PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_WGRIB2_COMPILER_FLAGS}>)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
    message(STATUS "No specific wgrib2 compiler flags set for Intel classic compiler")
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "PGI" OR ${CMAKE_C_COMPILER_ID} STREQUAL "NVHPC")
    set(ADDITIONAL_WGRIB2_COMPILER_FLAGS "--diag_suppress code_is_unreachable,integer_sign_change,declared_but_not_referenced,set_but_not_used")
    set(ADDITIONAL_GCTPC_COMPILER_FLAGS "--diag_suppress set_but_not_used")
    message(STATUS "Adding additional compiler flags to gctpc: ${ADDITIONAL_GCTPC_COMPILER_FLAGS}")
    message(STATUS "Adding additional compiler flags to wgrib2: ${ADDITIONAL_WGRIB2_COMPILER_FLAGS}")
    target_compile_options(geo PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_GCTPC_COMPILER_FLAGS}>)
    target_compile_options(grib2 PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_WGRIB2_COMPILER_FLAGS}>)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Clang")
    set(ADDITIONAL_G2C_COMPILER_FLAGS "-Wno-absolute-value")
    set(ADDITIONAL_GCTPC_COMPILER_FLAGS "-Wno-absolute-value")
    message(STATUS "Adding additional compiler flags to g2c: ${ADDITIONAL_G2C_COMPILER_FLAGS}")
    message(STATUS "Adding additional compiler flags to gctpc: ${ADDITIONAL_GCTPC_COMPILER_FLAGS}")
    target_compile_options(g2c PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_G2C_COMPILER_FLAGS}>)
    target_compile_options(geo PRIVATE $<$<COMPILE_LANGUAGE:C>:${ADDITIONAL_GCTPC_COMPILER_FLAGS}>)
endif()

###########################################################################
# Modern ALIAS targets with namespacing
###########################################################################
# Create modern namespaced alias targets for external use
add_library(GRIB2::grib2 ALIAS grib2)
add_library(GRIB2::g2c ALIAS g2c)
add_library(GRIB2::geo ALIAS geo)
###########################################################################
