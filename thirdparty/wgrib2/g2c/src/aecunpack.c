/** @file
 * @brief Unpack a data field that was packed with AEC compression.
 * stream
 * @author Eric Engle @date 2023-10-16
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack AEC compressed data into an
 * array of floats or doubles, using info from the GRIB2 Data
 * Representation [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 *
 * This function is used by aecunpack(), g2c_aecunpackf(), and
 * g2c_aecunpackd().
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
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
 * @author Eric Engle @date 2023-10-16
 */
static int
aecunpack_int(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
              void *fld, int fld_is_double, int verbose)
{
    g2int *ifld;
    g2int j, ctemplen = 0, nbits;
    g2int ccsds_flags, ccsds_block_size, ccsds_rsi;
    //g2int ifld1 = 0;
    int ret = 0;
    float ref, bscale, dscale;
    float *ffld = fld;
    double *dfld = fld;
    unsigned char *ctemp;
    size_t nbytes = 0;

    LOG((2, "aecunpack_int len %ld ndpts %ld fld_is_double %d", len, ndpts, fld_is_double));

    /* Get compression parameters from data representation template array. */
    rdieee(idrstmpl, &ref, 1);
    bscale = int_power(2.0, idrstmpl[1]);
    dscale = int_power(10.0, -idrstmpl[2]);
    nbits = idrstmpl[3];
    ccsds_flags = idrstmpl[5];
    ccsds_block_size = idrstmpl[6];
    ccsds_rsi = idrstmpl[7];

    /* If nbits equals 0, we have a constant field where the reference
     * value is the data value at each gridpoint. */
    if (nbits != 0)
    {
        if (!(ifld = calloc(ndpts, sizeof(g2int))))
        {
            if (verbose)
                fprintf(stderr, "Could not allocate space in aecunpack.\n  Data field NOT upacked.\n");
            return G2C_ENOMEM;
        }

        /* Determine the number of bytes needed for each value, then allocate the
         * buffer for the decoded AEC stream. */
        nbytes = (nbits + 7) / 8;
        if (nbytes == 3)
            nbytes = 4;
        ctemplen = nbytes * ndpts;
        if ((ctemp = (unsigned char *)malloc(ctemplen)) == NULL)
        {
            if (verbose)
                fprintf(stderr, "Allocation error.\n");
            return G2C_ENOMEM;
        }

        /* Decode the AEC stream. */
        ret = dec_aec(cpack, len, nbits, ccsds_flags, ccsds_block_size, ccsds_rsi, ctemp, ctemplen);
        if (ret < 0)
            return ret;

        /* IMPORTANT: The decoded AEC stream is byte-aligned (not bit), so when extracting the data
         * values from the buffer via gbits, we need to pass the byte size in bits, not nbits. */
        gbits(ctemp, ifld, 0, nbytes * 8, 0, ndpts);

        /* Zero out all higher-order bits, preserving nbits. NOTE: This might be unecessary. */
        //for (j = 0; j < ndpts; j++)
        //    ifld1 = ifld[j];
        //    ifld1 = ifld1 << (sizeof(ifld1)*8-nbits);
        //    ifld1 = ifld1 >> (sizeof(ifld1)*8-nbits);
        //    ifld[j] = ifld1;

        /* Unscale data. */
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

        /* Clean up. */
        free(ctemp);
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
 * Unpack AEC compressed data into an array of floats, using info
 * from the GRIB2 Data Representation [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of float.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Eric Engle @date 2023-10-16
 */
g2int
aecunpack(unsigned char *cpack, g2int len, g2int *idrstmpl, g2int ndpts,
          float *fld)
{
    int ret;

    LOG((2, "g2c_aecunpack len %lld ndpts %lld", len, ndpts));

    ret = aecunpack_int(cpack, len, idrstmpl, ndpts, fld, 0, 1);

    return ret;
}

/**
 * Unpack AEC compressed data into an array of floats, using info
 * from the GRIB2 Data Representation [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of float.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Eric Engle @date 2022-10-16
 */
int
g2c_aecunpackf(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               float *fld)
{
    g2int idrstmpl8[G2C_AEC_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    LOG((2, "g2c_aecunpackf len %d ndpts %lld", len, ndpts));

    for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return aecunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 0, 0);
}

/**
 * Unpack AEC compressed data into an array of doubles, using info
 * from the GRIB2 Data Representation [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 *
 * This function is the V2 API version of aecunpack() for doubles.
 *
 * @param cpack The packed data.
 * @param len The length of the packed data.
 * @param idrstmpl Pointer to array of values for Data Representation
 * [Template
 * 5.42](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml).
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values as an array
 * of double.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Eric Engle @date 2023-10-16
 */
int
g2c_aecunpackd(unsigned char *cpack, size_t len, int *idrstmpl, size_t ndpts,
               double *fld)
{
    g2int idrstmpl8[G2C_AEC_DRS_TEMPLATE_LEN];
    g2int len8 = len, ndpts8 = ndpts;
    int i;

    LOG((2, "g2c_aecunpackd len %lld ndpts %lld", len, ndpts));

    for (i = 0; i < G2C_AEC_DRS_TEMPLATE_LEN; i++)
        idrstmpl8[i] = idrstmpl[i];

    return aecunpack_int(cpack, len8, idrstmpl8, ndpts8, fld, 1, 0);
}
