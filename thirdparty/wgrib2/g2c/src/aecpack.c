/** @file
 * @brief Pack array of float/double using AEC/CCSDS compression.
 * @author Eric Engle @date 2023-09-10
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2023-09-10 | Engle | Initial; Adapted from aecpack.c.
 * 2023-10-16 | Engle | Added include libaec to set default values
 *                      for CCSDS parameters.
 */

#include "grib2_int.h"
#include <libaec.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * Pack a float or double array into a AEC/CCSDS code stream.
 *
 * @param fld Pointer to the float or double data values to pack.
 * @param fld_is_double If non-zero, then fld points to array of
 * doubles, otherwise an array of floats.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * May be modified in this function.
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
 * - ::G2C_EAEC Error encoding/decoding AEC data.
 *
 * @author Eric Engle
 */
static int
aecpack_int(void *fld, int fld_is_double, g2int width, g2int height, g2int *idrstmpl,
            unsigned char *cpack, g2int *lcpack, int verbose)
{
    g2int ctemplen;
    g2int *ifld = NULL;
    g2int j;
    g2int imin, imax;
    g2int maxdif, nbits, ndpts, nbytes;
    g2int ccsds_flags, ccsds_block_size, ccsds_rsi;
    float bscale, dscale, rmax, rmin, temp;
    double rmaxd, rmind;
    static float alog2 = ALOG2; /* ln(2.0) */
    unsigned char *ctemp;
    float *ffld = fld;
    double *dfld = fld;
    int ret = G2C_NOERROR;

    LOG((2, "aecpack_int() fld_is_double %d width %ld height %ld idrstmpl[1] %d *lcpack %ld",
         fld_is_double, width, height, idrstmpl[1], *lcpack));
    LOG((3, "idrstmpl: %ld %ld %ld %ld %ld %ld %ld %ld", idrstmpl[0], idrstmpl[1], idrstmpl[2],
         idrstmpl[3], idrstmpl[4], idrstmpl[5], idrstmpl[6], idrstmpl[7]));

    ctemplen = 0;
    ccsds_flags = 0;
    ccsds_block_size = 0;
    ccsds_rsi = 0;
    nbits = 0;

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
        /* Allocate memory for scaled data. */
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

        /* Define AEC compression options */
        if (idrstmpl[3] <= 0)
        {
            nbits = pow(2, ceil(log(nbits) / log(2))); // Round to nearest base 2 int
        }
        else
        {
            nbits = idrstmpl[3];
            nbits = pow(2, ceil(log(nbits) / log(2))); // Round to nearest base 2 int
        }
        nbits = nbits < 8 ? 8 : nbits;

        if (idrstmpl[5] == 0)
        {
            ccsds_flags = AEC_DATA_SIGNED | AEC_DATA_PREPROCESS | AEC_DATA_MSB;
        }
        else
        {
            ccsds_flags = idrstmpl[5];
        }
        if (idrstmpl[6] == 0)
        {
            ccsds_block_size = 16;
        }
        else
        {
            ccsds_block_size = idrstmpl[6];
        }
        if (idrstmpl[7] == 0)
        {
            ccsds_rsi = 128;
        }
        else
        {
            ccsds_rsi = idrstmpl[7];
        }

        /* Pack data into full octets, then do AEC encode and
         * calculate the length of the packed data in bytes. */
        nbytes = (nbits + 7) / 8;
        ctemp = calloc(ndpts, nbytes);
        ctemplen = ndpts * nbytes;
        sbits(ctemp, ifld, 0, nbytes * 8, 0, ndpts);

        *lcpack = enc_aec(ctemp, ctemplen, nbits, ccsds_flags, ccsds_block_size, ccsds_rsi, cpack, lcpack);
        if (*lcpack < 0)
        {
            if (verbose)
                printf("aecpack: ERROR Packing AEC = %d\n", ret);
            nbits = 0;
            *lcpack = 0;
            ret = G2C_EAEC;
        }

        free(ctemp);
    }
    else
    {
        nbits = 0;
        *lcpack = 0;
    }

    /* Fill in values for template 5.42. */
    if (fld_is_double)
        rmin = (float)rmind;
    mkieee(&rmin, idrstmpl, 1); /* ensure reference value is IEEE format. */
    idrstmpl[3] = nbits;
    idrstmpl[5] = ccsds_flags;
    idrstmpl[6] = ccsds_block_size;
    idrstmpl[7] = ccsds_rsi;
    if (ifld)
        free(ifld);

    return ret;
}

/**
 * This function packs up a float array into a AEC code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are passed to the AEC encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.42
 * with the appropriate values.
 *
 * @param fld Pointer to the float data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * - 0 Reference value - ignored on input, set by aecpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 CCSDS compression options mask.
 * - 6 Block size.
 * - 7 Reference sample interval.
 * May be modified in this function.
 * @param cpack A pointer that will get the packed data field. Must be
 * allocated before this function is called. Pass the allocated size
 * in the lcpack parameter.
 * @param lcpack Pointer that gets the length of packed field in
 * cpack. This must be set by the calling function to the size
 * available in cpack.
 *
 * @author Eric Engle (adapted from jpcpack)
 */
void
aecpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
        unsigned char *cpack, g2int *lcpack)
{
    aecpack_int(fld, 0, width, height, idrstmpl, cpack, lcpack, 1);
}

/**
 * This function packs up a float array into a AEC code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are passed to the AEC encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.42
 * with the appropriate values.
 *
 * This function is the V2 API version of aecpack() for floats.
 *
 * @param fld Pointer to the float data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * - 0 Reference value - ignored on input, set by aecpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 CCSDS compression options mask.
 * - 6 Block size.
 * - 7 Reference sample interval.
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
 * - ::G2C_EAEC Error encoding/decoding AEC data.
 *
 * @author Eric Engle (adapted from jpcpack)
 */
int
g2c_aecpackf(float *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, size_t *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_AEC_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = aecpack_int(fld, 0, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}

/**
 * This function packs up a double array into a AEC code stream.
 *
 * After the data are scaled, and the reference value is subtracted
 * out, the data are passed to the AEC encoder.
 *
 * This function also fills in GRIB2 Data Representation Template 5.42
 * with the appropriate values.
 *
 * This function is the V2 API version of aecpack() for floats.
 *
 * @param fld Pointer to the float data values to pack.
 * @param width The number of points in the x direction.
 * @param height The number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template [Table
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * - 0 Reference value - ignored on input, set by aecpack routine.
 * - 1 Binary Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 2 Decimal Scale Factor - used on input, unchanged by aecpack
 routine.
 * - 3 number of bits for each data value - ignored on input
 * - 4 Original field type - currently ignored on input Data values
 assumed to be reals. Set to 0 on output.
 * - 5 CCSDS compression options mask.
 * - 6 Block size.
 * - 7 Reference sample interval.
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
 * - ::G2C_EAEC Error encoding/decoding AEC data.
 *
 * @author Eric Engle (adapted from jpcpack)
 */
int
g2c_aecpackd(double *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, size_t *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_AEC_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = aecpack_int(fld, 1, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}
