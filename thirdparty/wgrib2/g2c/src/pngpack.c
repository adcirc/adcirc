/** @file
 * @brief Pack data with PNG compression.
 * @author Stephen Gilbert @date 2003-08-27
 */
#include "grib2_int.h"
#include <math.h>
#include <stdlib.h>

/**
 * Packs float or double data into PNG image
 * format. This is called by pngpack() and pngpackd().
 *
 * After the data field is scaled, and the reference value is
 * subtracted out, it is treated as a grayscale image and passed to a
 * PNG encoder. It also fills in GRIB2 Data Representation Template
 * 5.41 or 5.40010 with the appropriate values.
 *
 * @param fld Pointer to array of float or double that contains the
 * data values to pack.
 * @param fld_is_double If non-zero, then fld is double, otherwise float.
 * @param width Number of points in the x direction. This is passed to
 * the PNG layer as a uint32.
 * @param height Number of points in the y direction. This is passed
 * to the PNG layer as a uint32.
 * @param idrstmpl Contains the array of values for Data
 * Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml)
 * or 5.40010.
 * - 0 Reference value - ignored on input, set by pngpack routine.
 * - 1 Binary Scale Factor - used on input.
 * - 2 Decimal Scale Factor - used on input.
 * - 3 number of bits for each grayscale pixel value - ignored on
 input.
 * - 4 Original field type - currently ignored on input, set = 0 on
 output. Data values assumed to be reals.
 * @param cpack The packed data field.
 * @param lcpack length of packed field cpack.
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
 * - ::G2C_EPNG Error encoding/decoding PNG data.
 *
 * @author Ed Hartnett @date Aug 8, 2022
 */
static int
pngpack_int(void *fld, int fld_is_double, g2int width, g2int height, g2int *idrstmpl,
            unsigned char *cpack, g2int *lcpack, int verbose)
{
    g2int *ifld = NULL;
    static float alog2 = ALOG2; /*  ln(2.0) */
    g2int j, nbits, imin, imax, maxdif;
    g2int ndpts, nbytes;
    float bscale, dscale, rmax, rmin, temp;
    double rmaxd, rmind;
    unsigned char *ctemp;
    float *ffld = fld;
    double *dfld = fld;
    int ret = G2C_NOERROR;

    LOG((2, "pngpack_int fld_is_double %d width %ld height %ld idrstmpl[1] %d",
         fld_is_double, width, height, idrstmpl[1]));

    ndpts = width * height;
    bscale = int_power(2.0, -idrstmpl[1]);
    dscale = int_power(10.0, idrstmpl[2]);
    LOG((3, "ndpts %d bscale %g dscale %g", ndpts, bscale, dscale));

    /* Find max and min values in the data. Either rmax and rmin will
     * be used (if fld_is_double is not true), or rmaxd and rmind will
     * be used (if fld_is_double is true). */
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
            /* No binary scaling and calculate minumum number of bits
             * in which the data will fit. */
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

        /* Pack data into full octets, then do PNG encode and
         * calculate the length of the packed data in bytes. */
        if (nbits <= 1)
            nbits = 1;
        else if (nbits <= 2)
            nbits = 2;
        else if (nbits <= 4)
            nbits = 4;
        else if (nbits <= 8)
            nbits = 8;
        else if (nbits <= 16)
            nbits = 16;
        else if (nbits <= 24)
            nbits = 24;
        else
            nbits = 32;

        int bytes_per_row = (nbits * width) / 8;
        if ((width * nbits) % 8 != 0)
        {
            bytes_per_row++;
        }
        nbytes = bytes_per_row * height;
        ctemp = calloc(nbytes, 1);
        for (j = 0; j < height; j++)
            sbits(ctemp + (j * bytes_per_row), ifld + (j * width), 0, nbits, 0, width);

        /* Encode data into PNG Format. */
        if ((*lcpack = (g2int)enc_png(ctemp, width, height, nbits, cpack)) <= 0)
        {
            if (verbose)
                printf("pngpack: ERROR Packing PNG = %d\n", (int)*lcpack);
            ret = G2C_EPNG;
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
    mkieee(&rmin, idrstmpl, 1); /* ensure reference value is IEEE format */
    idrstmpl[3] = nbits;
    idrstmpl[4] = 0; /* original data were reals */

    if (ifld)
        free(ifld);

    return ret;
}

/**
 * This subroutine packs up a float data field into PNG image format.
 *
 * After the data field is scaled, and the reference value is
 * subtracted out, it is treated as a grayscale image and passed to a
 * PNG encoder. It also fills in GRIB2 Data Representation Template
 * 5.41 or 5.40010 with the appropriate values.
 *
 * @param fld Pointer to array of float that contains the data values
 * to pack.
 * @param width Number of points in the x direction.
 * @param height Number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml)
 * or 5.40010.
 * - 0 Reference value - ignored on input, set by pngpack routine.
 * - 1 Binary Scale Factor - used on input.
 * - 2 Decimal Scale Factor - used on input.
 * - 3 number of bits for each grayscale pixel value - ignored on
 input.
 * - 4 Original field type - currently ignored on input, set = 0 on
 output. Data values assumed to be reals.
 * @param cpack The packed data field.
 * @param lcpack length of packed field cpack.
 *
 * @author Stephen Gilbert @date 2003-08-27
 * @author Ed Hartnett
 */
void
pngpack(float *fld, g2int width, g2int height, g2int *idrstmpl,
        unsigned char *cpack, g2int *lcpack)
{
    /* Ignore the return value. */
    pngpack_int(fld, 0, width, height, idrstmpl, cpack, lcpack, 1);
}

/**
 * This subroutine packs up a float data field into PNG image format.
 *
 * After the data field is scaled, and the reference value is
 * subtracted out, it is treated as a grayscale image and passed to a
 * PNG encoder. It also fills in GRIB2 Data Representation Template
 * 5.41 or 5.40010 with the appropriate values.
 *
 * @param fld Pointer to array of float that contains the data values
 * to pack.
 * @param width Number of points in the x direction.
 * @param height Number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml)
 * or 5.40010.
 * - 0 Reference value - ignored on input, set by pngpack routine.
 * - 1 Binary Scale Factor - used on input.
 * - 2 Decimal Scale Factor - used on input.
 * - 3 number of bits for each grayscale pixel value - ignored on
 input.
 * - 4 Original field type - currently ignored on input, set = 0 on
 output. Data values assumed to be reals.
 * @param cpack The packed data field.
 * @param lcpack length of packed field cpack.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EPNG Error encoding/decoding PNG data.
 *
 * @author Ed Hartnett
 */
int
g2c_pngpackf(float *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, int *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_PNG_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = pngpack_int(fld, 0, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}

/**
 * This subroutine packs up a double data field into PNG image format.
 *
 * After the data field is scaled, and the reference value is
 * subtracted out, it is treated as a grayscale image and passed to a
 * PNG encoder. It also fills in GRIB2 Data Representation Template
 * 5.41 or 5.40010 with the appropriate values.
 *
 * @param fld Pointer to array of double that contains the data values
 * to pack.
 * @param width Number of points in the x direction.
 * @param height Number of points in the y direction.
 * @param idrstmpl Contains the array of values for Data
 * Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml)
 * or 5.40010.
 * - 0 Reference value - ignored on input, set by pngpack routine.
 * - 1 Binary Scale Factor - used on input.
 * - 2 Decimal Scale Factor - used on input.
 * - 3 number of bits for each grayscale pixel value - ignored on
 input.
 * - 4 Original field type - currently ignored on input, set = 0 on
 output. Data values assumed to be reals.
 * @param cpack The packed data field.
 * @param lcpack length of packed field cpack.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EPNG Error encoding/decoding PNG data.
 *
 * @author Ed Hartnett @date Aug 8, 2022
 */
int
g2c_pngpackd(double *fld, size_t width, size_t height, int *idrstmpl,
             unsigned char *cpack, int *lcpack)
{
    g2int width8 = width, height8 = height, lcpack8 = *lcpack;
    g2int idrstmpl8[G2C_PNG_DRS_TEMPLATE_LEN];
    int i, ret;

    for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    ret = pngpack_int(fld, 1, width8, height8, idrstmpl8, cpack, &lcpack8, 0);

    if (!ret)
    {
        for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
            idrstmpl[i] = (int)idrstmpl8[i];
        *lcpack = (g2int)lcpack8;
    }
    return ret;
}
