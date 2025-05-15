/** @file
 * @brief Pack and unpack an array of float/double using JPEG2000.
 * @author Stephen Gilbert @date 2003-08-17
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2003-08-17 | Gilbert | Initial.
 * 2004-11-92 | Gilbert | Fixed bug packing a near constant field.
 * 2004-07-19 | Gilbert | If jpeg2000 encoding fails, try again with different encoder options.
 * 2005-05-10 | Gilbert | Imposed minimum size on cpack.
 * 2022-08-12 | Hartnett | Now handle doubles too.
 */
#include "grib2_int.h"
#include <math.h>
#include <stdlib.h>

/**
 * Packs a float or double array into a JPEG2000 code stream.
 *
 * This function is used by jpcpack(), g2c_jpcpackf(), and
 * g2c_jpcpackd().
 *
 * @param fld Pointer to the float or double data values to pack.
 * @param fld_is_double If non-zero, then fld points to array of
 * doubles, otherwise an array of floats.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000. May be modified in this function.
 * @param cpack A pointer that will get the packed data field. Must be
 * allocated before this function is called. Pass the allocated size
 * in the lcpack parameter.
 * @param lcpack Pointer that gets the length of packed field in
 * cpack. This must also be set by the calling function to the size
 * available in cpack.
 * @param verbose If non-zero, error messages will be printed in case
 * of error. Otherwise, error codes will be return but no error
 * messages printed. Calls to the original g2c API may cause error
 * messages to be printed in case of error. For the new g2c_ API, no
 * error messages will be printed - instead an error code will be
 * returned. Call g2c_strerror() to get the error message for any
 * error code.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EJPEG Error encoding/decoding JPEG data.
 *
 * @author Stephen Gilbert, Ed Hartnett
 */
static int
jpcpack_int(void *fld, int fld_is_double, g2int width, g2int height, g2int *idrstmpl,
            unsigned char *cpack, g2int *lcpack, int verbose)
{
    g2int *ifld = NULL;
    static float alog2 = ALOG2; /*  ln(2.0) */
    g2int j, nbits, imin, imax, maxdif;
    g2int ndpts, nbytes, nsize, retry;
    float bscale, dscale, rmax, rmin, temp;
    double rmaxd, rmind;
    unsigned char *ctemp;
    float *ffld = fld;
    double *dfld = fld;
    int ret = G2C_NOERROR;

    LOG((2, "jpcpack_int() fld_is_double %d width %ld height %ld idrstmpl[1] %d *lcpack %ld",
         fld_is_double, width, height, idrstmpl[1], *lcpack));
    LOG((3, "idrstmpl: %ld %ld %ld %ld %ld %ld %ld", idrstmpl[0], idrstmpl[1], idrstmpl[2],
         idrstmpl[3], idrstmpl[4], idrstmpl[5], idrstmpl[6]));

    ndpts = width * height;
    bscale = int_power(2.0, -idrstmpl[1]);
    dscale = int_power(10.0, idrstmpl[2]);
    LOG((3, "ndpts %ld bscale %g dscale %g", ndpts, bscale, dscale));

    /* Find max and min values in the data. */
    rmaxd = dfld[0];
    rmind = dfld[0];
    rmax = ffld[0];
    rmin = ffld[0];
    if (fld_is_double)
    {
        for (j = 1; j < ndpts; j++)
        {
            if (dfld[j] > rmaxd)
                rmaxd = dfld[j];
            if (dfld[j] < rmind)
                rmind = dfld[j];
        }
        if (idrstmpl[1] == 0)
            maxdif = (g2int)(rint(rmaxd * dscale) - rint(rmind * dscale));
        else
            maxdif = (g2int)rint((rmaxd - rmind) * dscale * bscale);
    }
    else
    {
        for (j = 1; j < ndpts; j++)
        {
            if (ffld[j] > rmax)
                rmax = ffld[j];
            if (ffld[j] < rmin)
                rmin = ffld[j];
        }
        if (idrstmpl[1] == 0)
            maxdif = (g2int)(rint(rmax * dscale) - rint(rmin * dscale));
        else
            maxdif = (g2int)rint((rmax - rmin) * dscale * bscale);
    }
    LOG((3, "rmax %g rmaxd %g rmin %g rmind %g", rmax, rmaxd, rmin, rmind));

    /* If max and min values are not equal, pack up field. If they are
     * equal, we have a constant field, and the reference value (rmin)
     * is the value for each point in the field and set nbits to 0. */
    if (((fld_is_double && rmind != rmaxd) || (!fld_is_double && rmin != rmax)) && maxdif != 0)
    {
        ifld = malloc(ndpts * sizeof(g2int));

        /* Determine which algorithm to use based on user-supplied
         * binary scale factor and number of bits. */
        if (idrstmpl[1] == 0)
        {
            /*  No binary scaling and calculate minumum number of bits
             *  in which the data will fit. */
            imin = (g2int)rint((fld_is_double ? rmind : rmin) * dscale);
            imax = (g2int)rint((fld_is_double ? rmaxd : rmax) * dscale);
            maxdif = imax - imin;
            temp = log((double)(maxdif + 1)) / alog2;
            nbits = (g2int)ceil(temp);
            /*   scale data */
            if (fld_is_double)
            {
                rmind = (float)imin;
                for (j = 0; j < ndpts; j++)
                    ifld[j] = (g2int)rint(dfld[j] * dscale) - imin;
            }
            else
            {
                rmin = (float)imin;
                for (j = 0; j < ndpts; j++)
                    ifld[j] = (g2int)rint(ffld[j] * dscale) - imin;
            }
        }
        else
        {
            /* Use binary scaling factor and calculate minumum number
             * of bits in which the data will fit. */
            if (fld_is_double)
            {
                rmind = rmind * dscale;
                rmaxd = rmaxd * dscale;
                maxdif = (g2int)rint((rmaxd - rmind) * bscale);
            }
            else
            {
                rmin = rmin * dscale;
                rmax = rmax * dscale;
                maxdif = (g2int)rint((rmax - rmin) * bscale);
            }

            temp = log((double)(maxdif + 1)) / alog2;
            nbits = (g2int)ceil(temp);
            /*   scale data */
            if (fld_is_double)
            {
                for (j = 0; j < ndpts; j++)
                    ifld[j] = (g2int)rint(((dfld[j] * dscale) - rmind) * bscale);
            }
            else
            {
                for (j = 0; j < ndpts; j++)
                    ifld[j] = (g2int)rint(((ffld[j] * dscale) - rmin) * bscale);
            }
        }

        /* Pack data into full octets, then do JPEG 2000 encode and
         * calculate the length of the packed data in bytes. */
        retry = 0;
        nbytes = (nbits + 7) / 8;
        nsize = *lcpack; /* needed for input to enc_jpeg2000 */
        ctemp = calloc(ndpts, nbytes);
        sbits(ctemp, ifld, 0, nbytes * 8, 0, ndpts);
        if ((*lcpack = (g2int)enc_jpeg2000(ctemp, width, height, nbits, idrstmpl[5],
                                           idrstmpl[6], retry, (char *)cpack, nsize)) <= 0)
        {
            if (verbose)
                printf("jpcpack: ERROR Packing JPC = %d\n", (int)*lcpack);
            ret = G2C_EJPEG;
            if (*lcpack == -3)
            {
                retry = 1;
                if ((*lcpack = (g2int)enc_jpeg2000(ctemp, width, height, nbits, idrstmpl[5],
                                                   idrstmpl[6], retry, (char *)cpack, nsize)) <= 0)
                {
                    if (verbose)
                        printf("jpcpack: Retry Failed.\n");
                    ret = G2C_EJPEG;
                }
                else
                {
                    if (verbose)
                        printf("jpcpack: Retry Successful.\n");
                    ret = G2C_NOERROR;
                }
            }
        }
        free(ctemp);
    }
    else
    {
        nbits = 0;
        *lcpack = 0;
    }

    /* Fill in ref value and number of bits in Template 5.0. */
    if (fld_is_double)
        rmin = (float)rmind;
    mkieee(&rmin, idrstmpl, 1); /* ensure reference value is IEEE format. */
    idrstmpl[3] = nbits;
    idrstmpl[4] = 0; /* original data were reals */
    if (idrstmpl[5] == 0)
        idrstmpl[6] = 255; /* lossy not used */
    if (ifld)
        free(ifld);

    return ret;
}

/**
 * This function packs up a float array into a JPEG2000 code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are treated as a grayscale image and passed to a
 * JPEG2000 encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.40
 * or 5.40000 with the appropriate values.
 *
 * @param fld Pointer to the float data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * - 0 Reference value - ignored on input, set by jpcpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 if 0 use lossless compression, if 1 use lossy compression.
 * - 6 Desired compression ratio, if idrstmpl[5]=1. Set to 255, if
 idrstmpl[5]=0.
 * May be modified in this function.
 * @param cpack A pointer that will get the packed data field. Must be
 * allocated before this function is called. Pass the allocated size
 * in the lcpack parameter.
 * @param lcpack Pointer that gets the length of packed field in
 * cpack. This must be set by the calling function to the size
 * available in cpack.
 *
 * @author Stephen Gilbert, Ed Hartnett
 */
void
jpcpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
        unsigned char *cpack, g2int *lcpack)
{
    jpcpack_int(fld, 0, width, height, idrstmpl, cpack, lcpack, 1);
}

/**
 * This function packs up a float array into a JPEG2000 code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are treated as a grayscale image and passed to a
 * JPEG2000 encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.40
 * or 5.40000 with the appropriate values.
 *
 * This function is the V2 API version of jpcpack() for floats.
 *
 * @param fld Pointer to the float data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * - 0 Reference value - ignored on input, set by jpcpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 if 0 use lossless compression, if 1 use lossy compression.
 * - 6 Desired compression ratio, if idrstmpl[5]=1. Set to 255, if
 idrstmpl[5]=0.
 * May be modified in this function.
 * @param cpack A pointer that will get the packed data field. Must be
 * allocated before this function is called. Pass the allocated size
 * in the lcpack parameter.
 * @param lcpack Pointer that gets the length of packed field in
 * cpack. This must be set by the calling function to the size
 * available in cpack.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EJPEG Error encoding/decoding JPEG data.
 *
 * @author Ed Hartnett
 */
int
g2c_jpcpackf(float *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, size_t *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_JPEG_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = jpcpack_int(fld, 0, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}

/**
 * This function packs up a double array into a JPEG2000 code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are treated as a grayscale image and passed to a
 * JPEG2000 encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.40
 * or 5.40000 with the appropriate values.
 *
 * This function is the V2 API version of jpcpack() for doubles.
 *
 * @param fld Pointer to the double data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * - 0 Reference value - ignored on input, set by jpcpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by jpcpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 if 0 use lossless compression, if 1 use lossy compression.
 * - 6 Desired compression ratio, if idrstmpl[5]=1. Set to 255, if
 idrstmpl[5]=0.
 * May be modified in this function.
 * @param cpack A pointer that will get the packed data field. Must be
 * allocated before this function is called. Pass the allocated size
 * in the lcpack parameter.
 * @param lcpack Pointer that gets the length of packed field in
 * cpack. This must be set by the calling function to the size
 * available in cpack.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EJPEG Error encoding/decoding JPEG data.
 *
 * @author Ed Hartnett
 */
int
g2c_jpcpackd(double *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, size_t *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_JPEG_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = jpcpack_int(fld, 1, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}
