/** @file
 * @brief Unpack Section 7 (Data Section) of a GRIB2 message.
 * @author Stephen Gilbert @date 2002-10-31
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-31 | Gilbert | Initial
 * 2002-12-20 | Gilbert | Added GDT info to arguments and added 5.51 processing.
 * 2003-08-29 | Gilbert | New templates using PNG and JPEG2000 algorithms/templates.
 * 2004-11-29 | Gilbert | JPEG2000 now allowed to use WMO Template 5.40 PNG allowed to use 5.41
 * 2004-12-16 | Taylor | Added check on comunpack return code.
 * 2008-12-23 | Wesley | Initialize Number of data points unpacked.
 * 2022-10-04 | Hartnett | Added g2c_unpack7().
 * 2023-10-16 | Engle | Added support for DRT 5.42, AEC compression.
 *
 */
#include "grib2_int.h"
#include <memory.h>
#include <string.h>

/**
 * Unpacks Section 7 (Data Section) of a GRIB2 message.
 *
 * This function is the internal function called by both g2_unpack7()
 * and g2c_unpack7().
 *
 * @param cgrib char array containing Section 7 of the GRIB2 message
 * @param iofst Pointer to a bit offset of the beginning of Section 7
 * in cgrib. This is updated by this function to reflect the data read
 * in this function. After this function is successfully called, the
 * value pointed to by iofst will be the number of bits to the end of
 * section 7 in cbuf.
 * @param igdsnum Grid Definition Template Number (see Code Table
 * 3.0). (Only used for DRS Template 5.51. May be 0 for other
 * templates.)
 * @param igdstmpl Pointer to an integer array containing the data
 * values for the specified Grid Definition Template (N=igdsnum). Each
 * element of this integer array contains an entry (in the order
 * specified) of Grid Definition Template 3.N. (Only used for DRS
 * Template 5.51, may be NULL for other templates).
 * @param idrsnum Data Representation Template Number (see Code Table
 * 5.0).
 * @param idrstmpl Pointer to an integer array containing the data
 * values for the specified Data Representation Template
 * (N=idrsnum). Each element of this integer array contains an entry
 * (in the order specified) of Data Representation Template 5.N
 * @param ndpts Number of data points to be unpacked and returned.
 * @param v1 If non-zero, then act like the V1 G2C API. This includes:
 * - printing error messages to stderr in the event of error.
 * - returning V1 error codes.
 * @param fld Pointer to a float pointer which gets a pointer to an
 * array allocated by this function to hold the unpacked data. This
 * memory must be freed by the caller.
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK7_BAD_DRT Unrecognized Data Representation Template.
 * - ::G2_UNPACK7_WRONG_GDT need one of GDT 3.50 through 3.53 to decode DRT 5.51
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 * - ::G2_UNPACK7_CORRUPT_SEC Corrupt section 7.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
static g2int
g2c_unpack7_int(unsigned char *cgrib, g2int *iofst, g2int igdsnum, g2int *igdstmpl,
                g2int idrsnum, g2int *idrstmpl, g2int ndpts, int v1, float **fld)
{
    g2int isecnum;
    g2int ipos, lensec;
    float *lfld;

    assert(cgrib && iofst && idrstmpl && fld);

    LOG((2, "g2c_unpack7_int *iofst %ld igdsnum %ld idrsnum %ld ndpts %ld v1 %d",
         *iofst, igdsnum, idrsnum, ndpts, v1));

    /* Get Length of Section */
    gbit(cgrib, &lensec, *iofst, 32);
    *iofst = *iofst + 32;
    LOG((3, "lensec %ld", lensec));

    /* Get Section Number */
    gbit(cgrib, &isecnum, *iofst, 8);
    *iofst = *iofst + 8;

    if (isecnum != 7)
        return G2_UNPACK_BAD_SEC;

    ipos = *iofst / 8;

    /* If we're using the V1 API, allocate memory for the unpacked
     * data. If we are using the v2 API, the caller is responsible for
     * allocating memory to hold the data. */
    if (v1)
    {
        if (!(lfld = calloc(ndpts ? ndpts : 1, sizeof(float))))
            return G2_UNPACK_NO_MEM;
        *fld = lfld;
    }
    else
        lfld = *fld;

    if (idrsnum == 0)
    {
        simunpack(cgrib + ipos, idrstmpl, ndpts, *fld);
    }
    else if (idrsnum == 2 || idrsnum == 3)
    {
        if (comunpack(cgrib + ipos, lensec, idrsnum, idrstmpl, ndpts, *fld))
            return G2_UNPACK7_CORRUPT_SEC;
    }
    else if (idrsnum == 50)
    {
        /* Spectral Simple */
        simunpack(cgrib + ipos, idrstmpl, ndpts - 1, lfld + 1);
        rdieee(idrstmpl + 4, lfld, 1);
    }
    else if (idrsnum == 51)
    {
        /* Spectral complex */
        if (igdsnum >= 50 && igdsnum <= 53)
            specunpack(cgrib + ipos, idrstmpl, ndpts, igdstmpl[0], igdstmpl[2],
                       igdstmpl[2], lfld);
        else
        {
            if (v1)
                fprintf(stderr, "g2_unpack7: Cannot use GDT 3.%d to unpack Data Section 5.51.\n",
                        (int)igdsnum);
            if (lfld)
                free(lfld);
            *fld = NULL;
            return G2_UNPACK7_WRONG_GDT;
        }
    }
#if defined USE_JPEG2000 || defined USE_OPENJPEG
    else if (idrsnum == 40 || idrsnum == 40000)
    {
        jpcunpack(cgrib + ipos, lensec - 5, idrstmpl, ndpts, *fld);
    }
#endif /* USE_JPEG2000 */
#ifdef USE_PNG
    else if (idrsnum == 41 || idrsnum == 40010)
    {
        pngunpack(cgrib + ipos, lensec - 5, idrstmpl, ndpts, *fld);
    }
#endif /* USE_PNG */
#ifdef USE_AEC
    else if (idrsnum == 42)
    {
        aecunpack(cgrib + ipos, lensec - 5, idrstmpl, ndpts, *fld);
    }
#endif /* USE_AEC */
    else
    {
        if (v1)
            fprintf(stderr, "g2_unpack7: Data Representation Template 5.%d not yet "
                            "implemented.\n",
                    (int)idrsnum);
        if (lfld)
            free(lfld);
        *fld = NULL;
        return G2_UNPACK7_BAD_DRT;
    }

    *iofst = *iofst + (8 * lensec);

    return G2_NO_ERROR;
}

/**
 * This subroutine unpacks Section 7 (Data Section) of a GRIB2
 * message.
 *
 * This function is maintained for backward compatibility. Users may
 * wish to use the newer g2c_unpack7() function instead.
 *
 * @param cgrib char array containing Section 7 of the GRIB2 message
 * @param iofst Pointer to a bit offset of the beginning of Section 7
 * in cgrib. This is updated by this function to reflect the data read
 * in this function. After this function is successfully called, the
 * value pointed to by iofst will be the number of bits to the end of
 * section 7 in cbuf.
 * @param igdsnum Grid Definition Template Number (see Code Table
 * 3.0). (Only used for DRS Template 5.51.)
 * @param igdstmpl Pointer to an integer array containing the data
 * values for the specified Grid Definition Template (N=igdsnum). Each
 * element of this integer array contains an entry (in the order
 * specified) of Grid Definition Template 3.N. (Only used for DRS
 * Template 5.51).
 * @param idrsnum Data Representation Template Number (see Code Table
 * 5.0).
 * @param idrstmpl Pointer to an integer array containing the data
 * values for the specified Data Representation Template
 * (N=idrsnum). Each element of this integer array contains an entry
 * (in the order specified) of Data Representation Template 5.N
 * @param ndpts Number of data points to be unpacked and returned.
 * @param fld Pointer to a float pointer which gets a pointer to an
 * array allocated by this function to hold the unpacked data. This
 * memory must be freed by the caller.
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK7_BAD_DRT Unrecognized Data Representation Template.
 * - ::G2_UNPACK7_WRONG_GDT need one of GDT 3.50 through 3.53 to decode DRT 5.51
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 * - ::G2_UNPACK7_CORRUPT_SEC Corrupt section 7.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
g2int
g2_unpack7(unsigned char *cgrib, g2int *iofst, g2int igdsnum, g2int *igdstmpl,
           g2int idrsnum, g2int *idrstmpl, g2int ndpts, float **fld)
{
    return g2c_unpack7_int(cgrib, iofst, igdsnum, igdstmpl, idrsnum, idrstmpl,
                           ndpts, 1, fld);
}

/**
 * This subroutine unpacks Section 7 (Data Section) of a GRIB2
 * message.
 *
 * This function is the newer version of g2_unpack7().
 *
 * @param cgrib char array containing Section 7 of the GRIB2 message
 * @param igdsnum Grid Definition Template Number (see Code Table
 * 3.0). (Only used for DRS Template 5.51.) May be zero for other
 * templates.
 * @param gds_tmpl_len Number of elements in the GDS template.
 * @param gdstmpl Pointer to an integer array containing the data
 * values for the specified Grid Definition Template (N=igdsnum). Each
 * element of this integer array contains an entry (in the order
 * specified) of Grid Definition Template 3.N. (Only used for DRS
 * Template 5.51). May be NULL.
 * @param idrsnum Data Representation Template Number (see Code Table
 * 5.0).
 * @param drs_tmpl_len Number of elements in the DRS template.
 * @param drstmpl Pointer to an integer array containing the data
 * values for the specified Data Representation Template
 * (N=idrsnum). Each element of this integer array contains an entry
 * (in the order specified) of Data Representation Template 5.N
 * @param ndpts Number of data points to be unpacked and returned.
 * @param fld Pointer which the data. Memory must be allocated in
 * advance by caller.
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK7_BAD_DRT Unrecognized Data Representation Template.
 * - ::G2_UNPACK7_WRONG_GDT need one of GDT 3.50 through 3.53 to decode DRT 5.51
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 * - ::G2_UNPACK7_CORRUPT_SEC Corrupt section 7.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
int
g2c_unpack7(unsigned char *cgrib, int igdsnum, int gds_tmpl_len, long long int *gdstmpl,
            int idrsnum, int drs_tmpl_len, long long int *drstmpl, int ndpts, float *fld)
{
    g2int iofst = 0;
    g2int *igdstmpl = NULL, *idrstmpl;
    int i;
    int ret;

    LOG((1, "g2c_unpack7 igdsnum %d gds_tmpl_len %d idrsnum %d drs_tmpl_len %d ndpts %d",
         igdsnum, gds_tmpl_len, idrsnum, drs_tmpl_len, ndpts));

    /* Check inputs. */
    assert(cgrib && drstmpl && fld);
    if (gds_tmpl_len && !gdstmpl)
        return G2C_EINVAL;

    /* Allocate memory to hold the g2int versions of the DRS and GDS
     * template arrays. */
    if (gds_tmpl_len)
        if (!(igdstmpl = malloc(gds_tmpl_len * sizeof(g2int))))
            return G2C_ENOMEM;
    if (!(idrstmpl = malloc(drs_tmpl_len * sizeof(g2int))))
        return G2C_ENOMEM;

    /* Copy the templates. */
    if (gds_tmpl_len)
        for (i = 0; i < gds_tmpl_len; i++)
            igdstmpl[i] = gdstmpl[i];
    for (i = 0; i < drs_tmpl_len; i++)
        idrstmpl[i] = drstmpl[i];

    /* Call the internal function that does the work. */
    ret = g2c_unpack7_int(cgrib, &iofst, igdsnum, igdstmpl, idrsnum, idrstmpl,
                          ndpts, 0, &fld);

    /* Free the g2int versions of the templates. */
    if (igdstmpl)
        free(igdstmpl);
    free(idrstmpl);

    return ret;
}
