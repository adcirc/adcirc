/** @file
 * @brief Finalize a GRIB2 message after all grids and fields have
 * been added.
 * @author Stephen Gilbert @date 2002-10-31
 */
#include "grib2_int.h"
#include <stdio.h>

/**
 * Finalize a GRIB2 message after all grids and fields
 * have been added.
 *
 * This function adds the End Section ("7777") to the end of
 * the GRIB message and calculates the length and stores it in the
 * appropriate place in Section 0. This routine is used with routines
 * g2_create(), g2_addlocal(), g2_addgrid(), and g2_addfield() to
 * create a complete GRIB2 message.
 *
 * @param cgrib Char array containing all the data sections added be
 * previous calls to g2_create(), g2_addlocal(), g2_addgrid(), and
 * g2_addfield(). After function is called, contains the finalized
 * GRIB2 message.
 *
 * @return
 * - > 0 Length of the final GRIB2 message in bytes.
 * - ::G2_GRIBEND_MSG_INIT GRIB message was not initialized - call
 g2_create() first.
 * - ::G2_BAD_SEC_COUNTS Sum of Section byte counts doesn't add to
 total byte count.

 * - ::G2_BAD_SEC Previous Section was not 7.
 *
 * @note This routine is intended for use with routines g2_create(),
 * g2_addlocal(), g2_addgrid(), and g2_addfield() to create a complete
 * GRIB2 message.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
g2int
g2_gribend(unsigned char *cgrib)
{
    g2int iofst, lencurr, len, ilen, isecnum;
    g2int lengrib;
    unsigned char seven = 0x37; /* '7' */
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
    len = 16; /* Length of Section 0. */
    for (;;)
    {
        /*    Get number and length of next section. */
        iofst = len * 8;
        gbit(cgrib, &ilen, iofst, 32);
        iofst = iofst + 32;
        gbit(cgrib, &isecnum, iofst, 8);
        len = len + ilen;

        /*    Exit loop if last section reached. */
        if (len == lencurr)
            break;

        /*    If byte count for each section doesn't match current
         *    total length, then there is a problem. */
        if (len > lencurr)
        {
            printf("g2_gribend: Section byte counts don''t add to total.\n");
            printf("g2_gribend: Sum of section byte counts  =  %d\n", (int)len);
            printf("g2_gribend: Total byte count in Section 0 = %d\n", (int)lencurr);
            return G2_BAD_SEC_COUNTS;
        }
    }

    /* Can only add End Section (Section 8) after Section 7. */
    if (isecnum != 7)
    {
        printf("g2_gribend: Section 8 can only be added after Section 7.\n");
        printf("g2_gribend: Section %ld was the last found in given GRIB message.\n", isecnum);
        return G2_BAD_SEC;
    }

    /* Add Section 8  - End Section */
    cgrib[lencurr] = seven;
    cgrib[lencurr + 1] = seven;
    cgrib[lencurr + 2] = seven;
    cgrib[lencurr + 3] = seven;

    /* Update current byte total of message in Section 0. */
    lengrib = lencurr + 4;
    sbit(cgrib, &lengrib, 96, 32);

    /* Return the length of the message. */
    return lengrib;
}
