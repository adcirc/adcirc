/** @file
 * @brief Unpack Section 5 (Data Representation Section) of a GRIB2 message.
 * @author Stephen Gilbert @date 2002-10-31
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack Section 5 (Data Representation Section) of a GRIB2 message.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-31 | Gilbert | Initial
 * 2009-01-14 | Vuong | Changed structure name template to gtemplate
 *
 * @param cgrib char array containing Section 5 of the GRIB2 message.
 * @param iofst Bit offset for the beginning of Section 5 in
 * cgrib. Returned with bit offset at the end of Section 5.
 * @param ndpts Number of data points unpacked and returned.
 * @param idrsnum Data Representation Template Number (see Code Table 5.0).
 * @param idrstmpl Pointer to an integer array containing the data
 * values for the specified Data Representation Template
 * (N=idrsnum). Each element of this integer array contains an entry
 * (in the order specified) of Data Representation Template 5.N.
 * @param mapdrslen- Number of elements in idrstmpl. i.e. number of
 * entries in Data Representation Template 5.N (N=idrsnum).
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 * - ::G2_UNPACK5_BAD_DRT "GRIB" message contains an undefined Data Representation Template.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
g2int
g2_unpack5(unsigned char *cgrib, g2int *iofst, g2int *ndpts, g2int *idrsnum,
           g2int **idrstmpl, g2int *mapdrslen)
{
    g2int needext, i, j, nbits, isecnum;
    g2int lensec, isign, newlen;
    g2int *lidrstmpl = 0;
    gtemplate *mapdrs;

    *idrstmpl = 0; /* NULL*/

    gbit(cgrib, &lensec, *iofst, 32); /* Get Length of Section */
    *iofst = *iofst + 32;
    gbit(cgrib, &isecnum, *iofst, 8); /* Get Section Number */
    *iofst = *iofst + 8;

    if (isecnum != 5)
    {
        *ndpts = 0;
        *mapdrslen = 0;
        return G2_UNPACK_BAD_SEC;
    }

    gbit(cgrib, ndpts, *iofst, 32); /* Get num of data points */
    *iofst = *iofst + 32;
    gbit(cgrib, idrsnum, *iofst, 16); /* Get Data Rep Template Num. */
    *iofst = *iofst + 16;

    /*   Gen Data Representation Template */
    if (!(mapdrs = getdrstemplate(*idrsnum)))
    {
        *mapdrslen = 0;
        return G2_UNPACK5_BAD_DRT;
    }
    *mapdrslen = mapdrs->maplen;
    needext = mapdrs->needext;

    /* Unpack each value into array ipdstmpl from the appropriate
     * number of octets, which are specified in corresponding
     * entries in array mapdrs. */
    if (*mapdrslen > 0)
        lidrstmpl = calloc(*mapdrslen, sizeof(g2int));
    if (!lidrstmpl)
    {
        *mapdrslen = 0;
        *idrstmpl = NULL;
        if (mapdrs)
            free(mapdrs);
        return G2_UNPACK_NO_MEM;
    }
    else
    {
        *idrstmpl = lidrstmpl;
    }
    for (i = 0; i < mapdrs->maplen; i++)
    {
        nbits = abs(mapdrs->map[i]) * 8;
        if (mapdrs->map[i] >= 0)
        {
            gbit(cgrib, lidrstmpl + i, *iofst, nbits);
        }
        else
        {
            gbit(cgrib, &isign, *iofst, 1);
            gbit(cgrib, lidrstmpl + i, *iofst + 1, nbits - 1);
            if (isign == 1)
                lidrstmpl[i] = -1 * lidrstmpl[i];
        }
        *iofst = *iofst + nbits;
    }

    /* Check to see if the Data Representation Template needs to be
     * extended. The number of values in a specific gtemplate may
     * vary depending on data specified in the "static" part of the
     * gtemplate. */
    if (needext == 1)
    {
        free(mapdrs);
        mapdrs = extdrstemplate(*idrsnum, lidrstmpl);
        newlen = mapdrs->maplen + mapdrs->extlen;
        lidrstmpl = realloc(lidrstmpl, newlen * sizeof(g2int));
        *idrstmpl = lidrstmpl;

        /*   Unpack the rest of the Data Representation Template */
        j = 0;
        for (i = *mapdrslen; i < newlen; i++)
        {
            nbits = abs(mapdrs->ext[j]) * 8;
            if (mapdrs->ext[j] >= 0)
            {
                gbit(cgrib, lidrstmpl + i, *iofst, nbits);
            }
            else
            {
                gbit(cgrib, &isign, *iofst, 1);
                gbit(cgrib, lidrstmpl + i, *iofst + 1, nbits - 1);
                if (isign == 1)
                    lidrstmpl[i] = -1 * lidrstmpl[i];
            }
            *iofst = *iofst + nbits;
            j++;
        }
        *mapdrslen = newlen;
    }
    if (mapdrs->ext)
        free(mapdrs->ext);
    if (mapdrs)
        free(mapdrs);

    return G2_NO_ERROR;
}
