/** @file
 * @brief Unpack data packed with PNG compression.
 * @author Stephen Gilbert @date 2003-08-27
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack a data field that was packed into a PNG
 * image format using info from the GRIB2 Data Representation Template
 * 5.41 or 5.40010.
 *
 * @param cpack The packed data field (character*1 array).
 * @param len The length of packed field cpack().
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml)
 * or 5.40010.
 * @param ndpts The number of data values to unpack.
 * @param fld Pointer that will get the unpacked data values.
 * @param fld_is_double If non-zero, then fld will get data as double,
 * otherwise float.
 * @param verbose If non-zero, error messages will be printed in case
 * of error. Otherwise, error codes will be return but no error
 * messages printed. Calls to the original g2c API may cause error
 * messages to be printed in case of error. For the new g2c_ API, no
 * error messages will be printed - instead an error code will be
 * returned. Call g2c_strerror() to get the error message for any
 * error code.
 *
 * @return
 * - ::G2C_NOERROR No Error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Stephen Gilbert, Ed Hartnett @date Aug 8, 2022
 */
static int
pngunpack_int(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
              void *fld, int fld_is_double, int verbose)
{
    g2int *ifld;
    g2int j, nbits, width, height;
    float ref, bscale, dscale;
    unsigned char *ctemp;
    float *ffld = fld;
    double *dfld = fld;

    LOG((2, "pngunpack_int len %ld ndpts %ld fld_is_double %d", len, ndpts, fld_is_double));

    rdieee(idrstmpl, &ref, 1);
    bscale = int_power(2.0, idrstmpl[1]);
    dscale = int_power(10.0, -idrstmpl[2]);
    nbits = idrstmpl[3];
    LOG((2, "bscale %g dscale %g nbits %ld", bscale, dscale, nbits));

    /* If nbits equals 0, we have a constant field where the reference
     * value is the data value at each gridpoint. */
    if (nbits != 0)
    {
        ifld = calloc(ndpts, sizeof(g2int));
        ctemp = calloc(ndpts * 4, 1);
        if (!ifld || !ctemp)
        {
            if (verbose)
                fprintf(stderr, "Could not allocate space in jpcunpack.\n  Data field NOT upacked.\n");
            return G2C_ENOMEM;
        }
        dec_png(cpack, &width, &height, ctemp);

        int bytes_per_row = (nbits * width) / 8;
        if ((width * nbits) % 8 != 0)
        {
            bytes_per_row++;
        }
        for (j = 0; j < height; j++)
        {
            gbits(ctemp + (j * bytes_per_row), ifld + (j * width), 0, nbits, 0, width);
        }

        for (j = 0; j < ndpts; j++)
        {
            if (fld_is_double)
                dfld[j] = (((double)ifld[j] * bscale) + ref) * dscale;
            else
                ffld[j] = (((float)ifld[j] * bscale) + ref) * dscale;
        }
        free(ctemp);
        free(ifld);
    }
    else
    {
        for (j = 0; j < ndpts; j++)
        {
            if (fld_is_double)
                dfld[j] = ref;
            else
                ffld[j] = ref;
        }
    }

    return 0;
}

/**
 * This subroutine unpacks a data field that was packed into a PNG
 * image format using info from the GRIB2 Data Representation Template
 * 5.41 or 5.40010.
 *
 * @param cpack The packed data field (character*1 array).
 * @param len length of packed field cpack().
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml) or 5.40010.
 * @param ndpts The number of data values to unpack.
 * @param fld Contains the unpacked data values.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2_JPCUNPACK_MEM Out of memory.
 *
 * @author Stephen Gilbert @date 2003-08-27
 * @author Ed Hartnett
 */
g2int
pngunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
          float *fld)
{
    int ret;

    if ((ret = pngunpack_int(cpack, len, idrstmpl, ndpts, fld, 0, 1)) == G2C_ENOMEM)
        return G2_JPCUNPACK_MEM;

    return ret;
}

/**
 * This subroutine unpacks a data field that was packed into a PNG
 * image format using info from the GRIB2 Data Representation Template
 * 5.41 or 5.40010.
 *
 * @param cpack The packed data field (character*1 array).
 * @param len length of packed field cpack().
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml) or 5.40010.
 * @param ndpts The number of data values to unpack.
 * @param fld Contains the unpacked data values.
 *
 * @return
 * - ::G2C_NOERROR No Error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date Sep 8, 2022
*/
int
g2c_pngunpackf(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               float *fld)
{
    g2int idrstmpl8[G2C_PNG_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return pngunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 0, 0);
}

/**
 * This subroutine unpacks a data field that was packed into a PNG
 * image format using info from the GRIB2 Data Representation Template
 * 5.41 or 5.40010.
 *
 * @param cpack The packed data field (character*1 array).
 * @param len length of packed field cpack().
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template 5.41](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml) or 5.40010.
 * @param ndpts The number of data values to unpack.
 * @param fld Contains the unpacked data values.
 *
 * @return
 * - ::G2C_NOERROR No Error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date Aug 8, 2022
 */
int
g2c_pngunpackd(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               double *fld)
{
    g2int idrstmpl8[G2C_PNG_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    for (i = 0; i < G2C_PNG_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return pngunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 1, 0);
}
