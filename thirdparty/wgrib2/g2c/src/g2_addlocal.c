/** @file
 * @brief Add a Local Use Section (Section 2) to a GRIB2 message.
 * @author Stephen Gilbeert @date 2002-11-01
 */

#include "grib2_int.h"
#include <stdio.h>

/**
 * Adds a [Local Use Section (Section
 * 2)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect2.shtml)
 * to a GRIB2 message.
 *
 * This function is used with routines g2_create(), g2_addgrid(),
 * g2_addfield(), and g2_gribend() to create a complete GRIB2 message.
 *
 * @param cgrib Char array that contains the GRIB2 message to which
 * section 2 should be added. Must be allocated large enough to store
 * the entire GRIB2 message.
 * @param csec2 Character array containing information to be added in
 * Section 2.
 * @param lcsec2 Number of bytes of character array csec2 to be added
 * to Section 2.
 *
 * @returns > 0 = Current size of updated GRIB2 message.

 * - ::G2_ADD_MSG_INIT GRIB message was not initialized. Need to call
 * routine gribcreate first.
 * - ::G2_ADD_MSG_COMPLETE GRIB message already complete. Cannot add
 *   new section.
 * - ::G2_BAD_SEC_COUNTS Sum of Section byte counts doesn't add to
 *   total byte count.
 * - ::G2_BAD_SEC Previous Section was not 1 or 7.
 *
 * @note The Local Use Section (Section 2) can only follow Section 1
 * or Section 7 in a GRIB2 message.
 *
 * @author Stephen Gilbeert @date 2002-11-01
 */
g2int
g2_addlocal(unsigned char *cgrib, unsigned char *csec2, g2int lcsec2)
{
    static g2int two = 2;
    g2int j, k, lensec2, iofst, ibeg, lencurr, ilen, len, istart;
    g2int isecnum;
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

    /*  Loop through all current sections of the GRIB message to find
     *  the last section number. */
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
            printf("g2_addlocal: Section byte counts don't add to total.\n");
            printf("g2_addlocal: Sum of section byte counts = %ld\n", len);
            printf("g2_addlocal: Total byte count in Section 0 = %ld\n", lencurr);
            return G2_BAD_SEC_COUNTS;
        }
    }

    /* Section 2 can only be added after sections 1 and 7. */
    if (isecnum != 1 && isecnum != 7)
    {
        printf("g2_addlocal: Section 2 can only be added after Section 1 or Section 7.\n");
        printf("g2_addlocal: Section %ld was the last found in given GRIB message.\n", isecnum);
        return G2_BAD_SEC;
    }

    /* Add Section 2  - Local Use Section. */
    ibeg = lencurr * 8;          /*   Calculate offset for beginning of section 2 */
    iofst = ibeg + 32;           /*   leave space for length of section */
    sbit(cgrib, &two, iofst, 8); /* Store section number (2) */
    istart = lencurr + 5;
    k = 0;
    for (j = istart; j < istart + lcsec2; j++)
        cgrib[j] = csec2[k++];

    /* Calculate length of section 2 and store it in octets 1-4 of
     * section 2. */
    lensec2 = lcsec2 + 5; /* bytes */
    sbit(cgrib, &lensec2, ibeg, 32);

    /* Update current byte total of message in Section 0. */
    lencurr += lensec2;
    sbit(cgrib, &lencurr, 96, 32);

    return lencurr;
}
