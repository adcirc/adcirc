/** @file
 * @brief Unpack a data field that was packed using a
 * complex packing algorithm as defined in the GRIB2 documention.
 * @author Stephen Gilbert @date 2002-10-29
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack a data field that was packed using a complex packing
 * algorithm, using info from the GRIB2 Data Representation [Template
 * 5.2](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-2.shtml)
 * or [Template
 * 5.3](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-3.shtml).
 * Supports GRIB2 complex packing templates with or without spatial
 * differences (i.e. DRTs 5.2 and 5.3).
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-29 | Gilbert | Initial
 * 2004-12-16 | Gilbert | Added test (from Arthur Taylor/MDL) verifying group widths/lengths
 *
 * @param cpack pointer to the packed data field.
 * @param lensec length of section 7 (used for error checking).
 * @param idrsnum Data Representation Template number. Must equal 2
 * or 3.
 * @param idrstmpl pointer to the array of values for Data
 * Representation Template 5.2 or 5.3
 * @param ndpts The number of data values to unpack
 * @param fld Contains the unpacked data values. Must be allocated
 * with at least ndpts * sizeof(float) bytes before calling this
 * routine.
 *
 * @return 0 for success, error code otherwise.
 *
 * @author Stephen Gilbert @date 2002-10-29
 */
int
comunpack(unsigned char *cpack, g2int lensec, g2int idrsnum,
          g2int *idrstmpl, g2int ndpts, float *fld)
{
    g2int nbitsd = 0, isign;
    g2int j, iofst, ival1, ival2, minsd, itemp, l, k, n, non = 0;
    g2int *ifld, *ifldmiss = 0;
    g2int *gref, *gwidth, *glen;
    g2int itype, ngroups, nbitsgref, nbitsgwidth, nbitsglen;
    g2int msng1, msng2;
    float ref, bscale, dscale, rmiss1, rmiss2;
    g2int totBit, totLen;

    LOG((3, "comunpack lensec %ld idrsnum %ld ndpts %ld", lensec, idrsnum, ndpts));

    rdieee(idrstmpl, &ref, 1);
    bscale = (float)int_power(2.0, idrstmpl[1]);
    dscale = (float)int_power(10.0, -idrstmpl[2]);
    nbitsgref = idrstmpl[3];
    itype = idrstmpl[4];
    ngroups = idrstmpl[9];
    nbitsgwidth = idrstmpl[11];
    nbitsglen = idrstmpl[15];
    if (idrsnum == 3)
        nbitsd = idrstmpl[17] * 8;

    /*   Constant field */
    if (ngroups == 0)
    {
        for (j = 0; j < ndpts; j++)
            fld[j] = ref;
        return (0);
    }

    iofst = 0;
    ifld = (g2int *)calloc(ndpts, sizeof(g2int));
    gref = (g2int *)calloc(ngroups, sizeof(g2int));
    gwidth = (g2int *)calloc(ngroups, sizeof(g2int));

    /*  Get missing values, if supplied */
    if (idrstmpl[6] == 1)
    {
        if (itype == 0)
            rdieee(idrstmpl + 7, &rmiss1, 1);
        else
            rmiss1 = (float)idrstmpl[7];
    }
    if (idrstmpl[6] == 2)
    {
        if (itype == 0)
        {
            rdieee(idrstmpl + 7, &rmiss1, 1);
            rdieee(idrstmpl + 8, &rmiss2, 1);
        }
        else
        {
            rmiss1 = (float)idrstmpl[7];
            rmiss2 = (float)idrstmpl[8];
        }
    }

    /*  Extract Spatial differencing values, if using DRS Template 5.3 */
    if (idrsnum == 3)
    {
        if (nbitsd != 0)
        {
            /* wne mistake here shoujld be unsigned int */
            gbit(cpack, &ival1, iofst, nbitsd);
            iofst = iofst + nbitsd;
            if (idrstmpl[16] == 2)
            {
                /* wne mistake here shoujld be unsigned int */
                gbit(cpack, &ival2, iofst, nbitsd);
                iofst = iofst + nbitsd;
            }
            gbit(cpack, &isign, iofst, 1);
            iofst = iofst + 1;
            gbit(cpack, &minsd, iofst, nbitsd - 1);
            iofst = iofst + nbitsd - 1;
            if (isign == 1)
                minsd = -minsd;
        }
        else
        {
            ival1 = 0;
            ival2 = 0;
            minsd = 0;
        }
    }

    /*  Extract Each Group's reference value */
    if (nbitsgref != 0)
    {
        gbits(cpack, gref, iofst, nbitsgref, 0, ngroups);
        itemp = nbitsgref * ngroups;
        iofst = iofst + itemp;
        if (itemp % 8 != 0)
            iofst = iofst + (8 - (itemp % 8));
    }
    else
    {
        for (j = 0; j < ngroups; j++)
            gref[j] = 0;
    }

    /*  Extract Each Group's bit width */
    if (nbitsgwidth != 0)
    {
        gbits(cpack, gwidth, iofst, nbitsgwidth, 0, ngroups);
        itemp = nbitsgwidth * ngroups;
        iofst = iofst + itemp;
        if (itemp % 8 != 0)
            iofst = iofst + (8 - (itemp % 8));
    }
    else
    {
        for (j = 0; j < ngroups; j++)
            gwidth[j] = 0;
    }

    for (j = 0; j < ngroups; j++)
        gwidth[j] = gwidth[j] + idrstmpl[10];

    /*  Extract Each Group's length (number of values in each group) */
    glen = calloc(ngroups, sizeof(g2int));
    if (nbitsglen != 0)
    {
        gbits(cpack, glen, iofst, nbitsglen, 0, ngroups);
        itemp = nbitsglen * ngroups;
        iofst = iofst + itemp;
        if (itemp % 8 != 0)
            iofst = iofst + (8 - (itemp % 8));
    }
    else
    {
        for (j = 0; j < ngroups; j++)
            glen[j] = 0;
    }
    for (j = 0; j < ngroups; j++)
        glen[j] = (glen[j] * idrstmpl[13]) + idrstmpl[12];
    glen[ngroups - 1] = idrstmpl[14];

    /*  Test to see if the group widths and lengths are consistent
     *  with number of values, and length of section 7. */
    totBit = 0;
    totLen = 0;
    for (j = 0; j < ngroups; j++)
    {
        totBit += (gwidth[j] * glen[j]);
        totLen += glen[j];
    }
    if (totLen != ndpts)
        return 1;
    if (totBit / 8. > lensec)
        return 1;

    /*  For each group, unpack data values */
    if (idrstmpl[6] == 0)
    { /* no missing values */
        n = 0;
        for (j = 0; j < ngroups; j++)
        {
            if (gwidth[j] != 0)
            {
                gbits(cpack, ifld + n, iofst, gwidth[j], 0, glen[j]);
                for (k = 0; k < glen[j]; k++)
                {
                    ifld[n] = ifld[n] + gref[j];
                    n = n + 1;
                }
            }
            else
            {
                for (l = n; l < n + glen[j]; l++)
                    ifld[l] = gref[j];
                n = n + glen[j];
            }
            iofst = iofst + (gwidth[j] * glen[j]);
        }
    }
    else if (idrstmpl[6] == 1 || idrstmpl[6] == 2)
    {
        /* missing values included */
        ifldmiss = malloc(ndpts * sizeof(g2int));
        for (j = 0; j < ndpts; j++)
            ifldmiss[j] = 0;
        n = 0;
        non = 0;
        for (j = 0; j < ngroups; j++)
        {
            if (gwidth[j] != 0)
            {
                msng1 = (g2int)int_power(2.0, gwidth[j]) - 1;
                msng2 = msng1 - 1;
                gbits(cpack, ifld + n, iofst, gwidth[j], 0, glen[j]);
                iofst = iofst + (gwidth[j] * glen[j]);
                for (k = 0; k < glen[j]; k++)
                {
                    if (ifld[n] == msng1)
                        ifldmiss[n] = 1;
                    else if (idrstmpl[6] == 2 && ifld[n] == msng2)
                        ifldmiss[n] = 2;
                    else
                    {
                        ifldmiss[n] = 0;
                        ifld[non++] = ifld[n] + gref[j];
                    }
                    n++;
                }
            }
            else
            {
                msng1 = (g2int)int_power(2.0, nbitsgref) - 1;
                msng2 = msng1 - 1;
                if (gref[j] == msng1)
                {
                    for (l = n; l < n + glen[j]; l++)
                        ifldmiss[l] = 1;
                }
                else if (idrstmpl[6] == 2 && gref[j] == msng2)
                {
                    for (l = n; l < n + glen[j]; l++)
                        ifldmiss[l] = 2;
                }
                else
                {
                    for (l = n; l < n + glen[j]; l++)
                        ifldmiss[l] = 0;
                    for (l = non; l < non + glen[j]; l++)
                        ifld[l] = gref[j];
                    non += glen[j];
                }
                n = n + glen[j];
            }
        }
    }

    if (gref)
        free(gref);
    if (gwidth)
        free(gwidth);
    if (glen)
        free(glen);

    /*  If using spatial differences, add overall min value, and sum up recursively */
    if (idrsnum == 3)
    { /* spatial differencing */
        if (idrstmpl[16] == 1)
        { /* first order */
            ifld[0] = ival1;
            if (idrstmpl[6] == 0)
                itemp = ndpts; /* no missing values */
            else
                itemp = non;
            for (n = 1; n < itemp; n++)
            {
                ifld[n] = ifld[n] + minsd;
                ifld[n] = ifld[n] + ifld[n - 1];
            }
        }
        else if (idrstmpl[16] == 2)
        { /* second order */
            ifld[0] = ival1;
            ifld[1] = ival2;
            if (idrstmpl[6] == 0)
                itemp = ndpts; /* no missing values */
            else
                itemp = non;
            for (n = 2; n < itemp; n++)
            {
                ifld[n] = ifld[n] + minsd;
                ifld[n] = ifld[n] + (2 * ifld[n - 1]) - ifld[n - 2];
            }
        }
    }

    /*  Scale data back to original form */
    if (idrstmpl[6] == 0)
    { /* no missing values */
        for (n = 0; n < ndpts; n++)
        {
            fld[n] = (((float)ifld[n] * bscale) + ref) * dscale;
        }
    }
    else if (idrstmpl[6] == 1 || idrstmpl[6] == 2)
    {
        /* missing values included */
        non = 0;
        for (n = 0; n < ndpts; n++)
        {
            if (ifldmiss[n] == 0)
            {
                fld[n] = (((float)ifld[non++] * bscale) + ref) * dscale;
            }
            else if (ifldmiss[n] == 1)
                fld[n] = rmiss1;
            else if (ifldmiss[n] == 2)
                fld[n] = rmiss2;
        }
        if (ifldmiss)
            free(ifldmiss);
    }

    if (ifld)
        free(ifld);

    return 0;
}
