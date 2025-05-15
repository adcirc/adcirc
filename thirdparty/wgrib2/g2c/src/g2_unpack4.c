/** @file
 * @brief Unpack Section 4 (Product Definition Section) of a GRIB2 message.
 * @author Stephen Gilbert @date 2002-10-31
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack [Section 4 (Product Definition
 * Section)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect4.shtml)
 * of a GRIB2 message.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-31 | Gilbert | Initial
 * 2009-01-14 | Vuong | Changed structure name template to gtemplate
 *
 * @param cgrib Array containing Section 4 of the GRIB2 message.
 * @param iofst Bit offset of the beginning of Section 4 in
 * cgrib. Returned with updated bit offset.
 * @param ipdsnum Product Definition Template Number (see [Table
 * 4.0](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-0.shtml)).
 * @param ipdstmpl Pointer that gets an integer array containing the data
 * values for the Product Definition Template specified by ipdsnum.
 * @param mappdslen Number of elements in ipdstmpl - i.e. number of
 * entries in Product Defintion Template specified by ipdsnum.
 * @param coordlist Pointer that gets an array containing floating
 * point values intended to document the vertical discretisation
 * associated to model data on hybrid coordinate vertical levels.
 * @param numcoord number of values in array coordlist.
 *
 * @returns
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK4_BAD_PDT message contains an undefined Product
 Definition Template.
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
g2int
g2_unpack4(unsigned char *cgrib, g2int *iofst, g2int *ipdsnum, g2int **ipdstmpl,
           g2int *mappdslen, float **coordlist, g2int *numcoord)
{

    g2int needext, i, j, nbits, isecnum;
    g2int lensec, isign, newlen;
    g2int *coordieee;
    g2int *lipdstmpl = 0;
    float *lcoordlist;
    gtemplate *mappds;

    *ipdstmpl = NULL;
    *coordlist = NULL;

    gbit(cgrib, &lensec, *iofst, 32); /* Get Length of Section */
    *iofst = *iofst + 32;
    gbit(cgrib, &isecnum, *iofst, 8); /* Get Section Number */
    *iofst = *iofst + 8;

    if (isecnum != 4)
    {
        *numcoord = 0;
        *mappdslen = 0;
        return G2_UNPACK_BAD_SEC;
    }

    gbit(cgrib, numcoord, *iofst, 16); /* Get num of coordinate values */
    *iofst = *iofst + 16;
    gbit(cgrib, ipdsnum, *iofst, 16); /* Get Prod. Def Template num. */
    *iofst = *iofst + 16;

    /* Get Product Definition Template */
    if (!(mappds = getpdstemplate(*ipdsnum)))
    {
        *mappdslen = 0;
        return G2_UNPACK4_BAD_PDT;
    }
    *mappdslen = mappds->maplen;
    needext = mappds->needext;

    /* Unpack each value into array ipdstmpl from the the
     * appropriate number of octets, which are specified in
     * corresponding entries in array mappds. */
    if (*mappdslen > 0)
        lipdstmpl = calloc(*mappdslen, sizeof(g2int));
    if (!lipdstmpl)
    {
        *mappdslen = 0;
        *ipdstmpl = NULL;
        if (mappds)
            free(mappds);
        return G2_UNPACK_NO_MEM;
    }
    *ipdstmpl = lipdstmpl;

    for (i = 0; i < mappds->maplen; i++)
    {
        nbits = abs(mappds->map[i]) * 8;
        if (mappds->map[i] >= 0)
        {
            gbit(cgrib, lipdstmpl + i, *iofst, nbits);
        }
        else
        {
            gbit(cgrib, &isign, *iofst, 1);
            gbit(cgrib, lipdstmpl + i, *iofst + 1, nbits - 1);
            if (isign == 1)
                lipdstmpl[i] = -1 * lipdstmpl[i];
        }
        *iofst = *iofst + nbits;
    }

    /* Check to see if the Product Definition Template needs to be
     * extended. The number of values in a specific template may
     * vary depending on data specified in the "static" part of the
     * gtemplate. */
    if (needext == 1)
    {
        free(mappds);
        mappds = extpdstemplate(*ipdsnum, lipdstmpl);
        newlen = mappds->maplen + mappds->extlen;
        lipdstmpl = realloc(lipdstmpl, newlen * sizeof(g2int));
        *ipdstmpl = lipdstmpl;
        /*   Unpack the rest of the Product Definition Template */
        j = 0;
        for (i = *mappdslen; i < newlen; i++)
        {
            nbits = abs(mappds->ext[j]) * 8;
            if (mappds->ext[j] >= 0)
            {
                gbit(cgrib, lipdstmpl + i, *iofst, nbits);
            }
            else
            {
                gbit(cgrib, &isign, *iofst, 1);
                gbit(cgrib, lipdstmpl + i, *iofst + 1, nbits - 1);
                if (isign == 1)
                    lipdstmpl[i] = -1 * lipdstmpl[i];
            }
            *iofst = *iofst + nbits;
            j++;
        }
        *mappdslen = newlen;
    }
    if (mappds->ext)
        free(mappds->ext);
    if (mappds)
        free(mappds);

    /* Get Optional list of vertical coordinate values after the
     * Product Definition Template, if necessary. */
    *coordlist = NULL;
    if (*numcoord != 0)
    {
        coordieee = calloc(*numcoord, sizeof(g2int));
        lcoordlist = calloc(*numcoord, sizeof(float));
        if (coordieee == 0 || lcoordlist == 0)
        {
            *numcoord = 0;
            *coordlist = NULL;
            if (coordieee)
                free(coordieee);
            if (lcoordlist)
                free(lcoordlist);
            return G2_UNPACK_NO_MEM;
        }
        else
        {
            *coordlist = lcoordlist;
        }
        gbits(cgrib, coordieee, *iofst, 32, 0, *numcoord);
        rdieee(coordieee, *coordlist, *numcoord);
        free(coordieee);
        *iofst = *iofst + (32 * (*numcoord));
    }

    return G2_NO_ERROR;
}
