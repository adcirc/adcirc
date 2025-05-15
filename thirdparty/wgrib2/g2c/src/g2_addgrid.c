/**
 * @file
 * @brief Pack a Grid Definition Section (Section 3) and adds it to
 * a GRIB2 message.
 * @author Stephen Gilbeert @date 2002-11-01
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Packs a [Grid Definition Section (Section
 * 3)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml)
 * and adds it to a GRIB2 message.
 *
 * This function is used with routines
 * g2_create(), g2_addlocal(), g2_addfield(), and g2_gribend() to
 * create a complete GRIB2 message. g2_create() must be called first
 * to initialize a new GRIB2 message.
 *
 * @param cgrib Char array that contains the GRIB2 message to which
 * section should be added. Must be allocated large enough to store
 * the entire GRIB2 message.
 * @param igds Contains information needed for GRIB Grid Definition
 * Section 3. Must be dimensioned >= 5.
 * - igds[0] Source of grid definition (see [Code Table
 * 3.0](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-0.shtml)).
 * - igds[1] Number of grid points in the defined grid.
 * - igds[2] Number of octets needed for each additional grid points
 * definition. Used to define number of points in each row (or column)
 * for non-regular grids. = 0, if using regular grid.
 * - igds[3] Interpretation of list for optional points
 * definition. (See [Code Table
 * 3.11](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-11.shtml)).
 * - igds[4] Grid Definition Template Number (See [Code Table
 * 3.1](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-1.shtml)).
 * @param igdstmpl Contains the data values for the specified Grid
 * Definition Template (igds[4]). Each element of this integer
 * array contains an entry (in the order specified) of Grid Defintion
 * Template.
 * @param ideflist (Used if igds[2] != 0) This array contains the
 * number of grid points contained in each row (or column).
 * @param idefnum (Used if igds[2] != 0) The number of entries in
 * array ideflist.  i.e. number of rows (or columns) for which
 * optional grid points are defined.
 *
 * @returns
 * - > 0 Current size of updated GRIB2 message
 * - ::G2_ADD_MSG_INIT GRIB message was not initialized. Need to
 *   call routine gribcreate first.
 * - ::G2_ADD_MSG_COMPLETE GRIB message already complete. Cannot
 *   add new section.
 * - ::G2_BAD_SEC_COUNTS Sum of Section byte counts doesn't add to
 *   total byte count
 * - ::G2_BAD_SEC Previous Section was not 1, 2 or 7.
 * - ::G2_ADDGRID_BAD_GDT Could not find requested Grid Definition
 *   Template.
 *
 * @note The Grid Def Section (Section 3) can only follow Section 1,
 * 2 or Section 7 in a GRIB2 message.
 *
 * @author Stephen Gilbeert @date 2002-11-01
 */
g2int
g2_addgrid(unsigned char *cgrib, g2int *igds, g2int *igdstmpl, g2int *ideflist,
           g2int idefnum)
{
    static g2int one = 1, three = 3, miss = 65535;
    g2int lensec3, iofst, ibeg, lencurr, len;
    g2int i, j, temp, ilen, isecnum, nbits;
    gtemplate *mapgrid = 0;
    int ret;

    /* Check for GRIB header and terminator. Translate the error codes
     * to the legacy G2 error codes. */
    if ((ret = g2c_check_msg(cgrib, &lencurr, 1)))
    {
        if (ret == G2C_ENOTGRIB)
            return G2_ADD_MSG_INIT;
        if (ret == G2C_EMSGCOMPLETE)
            return G2_ADD_MSG_COMPLETE;
    }

    /* Loop through all current sections of the GRIB message to find
     * the last section number. */
    len = 16; /* length of Section 0 */
    for (;;)
    {
        /* Get section number and length of next section. */
        iofst = len * 8;
        gbit(cgrib, &ilen, iofst, 32);
        iofst = iofst + 32;
        gbit(cgrib, &isecnum, iofst, 8);
        len = len + ilen;

        /* Exit loop if last section reached. */
        if (len == lencurr)
            break;

        /* If byte count for each section doesn't match current total
         * length, then there is a problem. */
        if (len > lencurr)
        {
            printf("g2_addgrid: Section byte counts don''t add to total.\n");
            printf("g2_addgrid: Sum of section byte counts = %ld\n", len);
            printf("g2_addgrid: Total byte count in Section 0 = %ld\n", lencurr);
            return G2_BAD_SEC_COUNTS;
        }
    }

    /* Section 3 can only be added after sections 1, 2 and 7. */
    if (isecnum != 1 && isecnum != 2 && isecnum != 7)
    {
        printf("g2_addgrid: Section 3 can only be added after Section 1, 2 or 7.\n");
        printf("g2_addgrid: Section ',isecnum,' was the last found in given GRIB message.\n");
        return G2_BAD_SEC;
    }

    /* Add Section 3  - Grid Definition Section. */
    ibeg = lencurr * 8;            /*   Calculate offset for beginning of section 3 */
    iofst = ibeg + 32;             /*   leave space for length of section */
    sbit(cgrib, &three, iofst, 8); /* Store section number (3) */
    iofst = iofst + 8;
    sbit(cgrib, igds + 0, iofst, 8); /* Store source of Grid def. */
    iofst = iofst + 8;
    sbit(cgrib, igds + 1, iofst, 32); /* Store number of data pts. */
    iofst = iofst + 32;
    sbit(cgrib, igds + 2, iofst, 8); /* Store number of extra octets. */
    iofst = iofst + 8;
    sbit(cgrib, igds + 3, iofst, 8); /* Store interp. of extra octets. */
    iofst = iofst + 8;

    /* if Octet 6 is not equal to zero, Grid Definition Template may
     * not be supplied. */
    if (igds[0] == 0)
        sbit(cgrib, igds + 4, iofst, 16); /* Store Grid Def Template num. */
    else
        sbit(cgrib, &miss, iofst, 16); /* Store missing value as Grid Def Template num. */
    iofst = iofst + 16;

    /* Get Grid Definition Template. */
    if (igds[0] == 0)
    {
        if (!(mapgrid = getgridtemplate(igds[4])))
            return G2_ADDGRID_BAD_GDT;

        /* Extend the Grid Definition Template, if necessary. The
         *  number of values in a specific template may vary depending
         *  on data specified in the "static" part of the template. */
        if (mapgrid->needext)
        {
            free(mapgrid);
            mapgrid = extgridtemplate(igds[4], igdstmpl);
        }
    }

    /* Pack up each input value in array igdstmpl into the
     * appropriate number of octets, which are specified in
     * corresponding entries in array mapgrid. */
    for (i = 0; i < mapgrid->maplen; i++)
    {
        nbits = abs(mapgrid->map[i]) * 8;
        if ((mapgrid->map[i] >= 0) || (igdstmpl[i] >= 0))
            sbit(cgrib, igdstmpl + i, iofst, nbits);
        else
        {
            sbit(cgrib, &one, iofst, 1);
            temp = abs(igdstmpl[i]);
            sbit(cgrib, &temp, iofst + 1, nbits - 1);
        }
        iofst = iofst + nbits;
    }

    /* Pack template extension, if appropriate. */
    j = mapgrid->maplen;
    if (mapgrid->needext && mapgrid->extlen > 0)
    {
        for (i = 0; i < mapgrid->extlen; i++)
        {
            nbits = abs(mapgrid->ext[i]) * 8;
            if (mapgrid->ext[i] >= 0 || igdstmpl[j] >= 0)
                sbit(cgrib, igdstmpl + j, iofst, nbits);
            else
            {
                sbit(cgrib, &one, iofst, 1);
                temp = abs(igdstmpl[j]);
                sbit(cgrib, &temp, iofst + 1, nbits - 1);
            }
            iofst = iofst + nbits;
            j++;
        }
    }
    if (mapgrid->ext)
        free(mapgrid->ext);
    free(mapgrid);

    /* If requested, insert optional list of numbers defining number
     * of points in each row or column. This is used for non regular
     * grids. */
    if (igds[2] != 0)
    {
        nbits = igds[2] * 8;
        sbits(cgrib, ideflist, iofst, nbits, 0, idefnum);
        iofst = iofst + (nbits * idefnum);
    }

    /* Calculate length of section 3 and store it in octets 1-4 of section 3. */
    lensec3 = (iofst - ibeg) / 8;
    sbit(cgrib, &lensec3, ibeg, 32);

    /* Update current byte total of message in Section 0. */
    lencurr += lensec3;
    sbit(cgrib, &lencurr, 96, 32);

    return lencurr;
}
