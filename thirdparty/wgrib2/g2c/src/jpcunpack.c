/** @file
 * @brief Unpack a data field that was packed with JPEG2000.
 * stream
 * @author Stephem Gilbert @date 2003-08-27
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack JPEG2000 compressed data into an array of floats or doubles,
 * using info from the GRIB2 Data Representation [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 *
 * This function is used by jpcunpack(), g2c_jpcunpackf(), and
 * g2c_jpcunpackd().
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values.
 * @param fld_is_double Non-zero if the data are to be unpacked into a
 * double array, otherwise data will be unpacked into a float array.
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
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date 2022-09-06
 */
static int
jpcunpack_int(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
              void *fld, int fld_is_double, int verbose)
{
    g2int *ifld;
    g2int j, nbits;
    float ref, bscale, dscale;
    float *ffld = fld;
    double *dfld = fld;

    LOG((2, "jpcunpack_int len %ld ndpts %ld fld_is_double %d", len, ndpts, fld_is_double));

    rdieee(idrstmpl, &ref, 1);
    bscale = int_power(2.0, idrstmpl[1]);
    dscale = int_power(10.0, -idrstmpl[2]);
    nbits = idrstmpl[3];

    /* If nbits equals 0, we have a constant field where the reference
     * value is the data value at each gridpoint. */
    if (nbits != 0)
    {
        if (!(ifld = calloc(ndpts, sizeof(g2int))))
        {
            if (verbose)
                fprintf(stderr, "Could not allocate space in jpcunpack.\n  Data field NOT upacked.\n");
            return G2C_ENOMEM;
        }
        dec_jpeg2000((char *)cpack, len, ifld);
        if (fld_is_double)
        {
            for (j = 0; j < ndpts; j++)
                dfld[j] = (((float)ifld[j] * bscale) + ref) * dscale;
        }
        else
        {
            for (j = 0; j < ndpts; j++)
                ffld[j] = (((float)ifld[j] * bscale) + ref) * dscale;
        }
        free(ifld);
    }
    else
    {
        if (fld_is_double)
        {
            for (j = 0; j < ndpts; j++)
                dfld[j] = ref;
        }
        else
        {
            for (j = 0; j < ndpts; j++)
                ffld[j] = ref;
        }
    }

    return G2C_NOERROR;
}

/**
 * Unpack JPEG2000 compressed data into an array of floats, using info
 * from the GRIB2 Data Representation [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of float.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2_JPCUNPACK_MEM Out of memory.
 *
 * @author Stephem Gilbert @date 2003-08-27
 */
g2int
jpcunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
          float *fld)
{
    int ret;

    LOG((2, "g2c_jpcunpack len %lld ndpts %lld", len, ndpts));

    if ((ret = jpcunpack_int(cpack, len, idrstmpl, ndpts, fld, 0, 1)) == G2_JPCUNPACK_MEM)
        return G2_JPCUNPACK_MEM;

    return ret;
}

/**
 * Unpack JPEG2000 compressed data into an array of floats, using info
 * from the GRIB2 Data Representation [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of float.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date 2022-09-08
 */
int
g2c_jpcunpackf(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               float *fld)
{
    g2int idrstmpl8[G2C_JPEG_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    LOG((2, "g2c_jpcunpackf len %d ndpts %lld", len, ndpts));

    for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return jpcunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 0, 0);
}

/**
 * Unpack JPEG2000 compressed data into an array of doubles, using info
 * from the GRIB2 Data Representation [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 *
 * This function is the V2 API version of jpcunpack() for doubles.
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.40](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml)
 * or 5.40000.
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of double.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date 2022-08-12
 */
int
g2c_jpcunpackd(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               double *fld)
{
    g2int idrstmpl8[G2C_JPEG_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    LOG((2, "g2c_jpcunpackd len %lld ndpts %lld", len, ndpts));

    for (i = 0; i < G2C_JPEG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return jpcunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 1, 0);
}
