/** @file
 *
 * @brief Search through a GRIB2 message and returns the number of
 * gridded fields found in the message and the number (and maximum
 * size) of Local Use Sections.
 *
 * @author Stephen Gilbeert @date 2002-10-28
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Search through a GRIB2 message and return the
 * number of gridded fields found in the message and the number (and
 * maximum size) of Local Use Sections.
 *
 * Other checks are performed to see if the message is a valid GRIB2
 * message.
 *
 * @param cgrib Pointer to a buffer containing the GRIB2 message.
 * @param listsec0 Pointer to an array that gets the information
 * decoded from GRIB Indicator Section 0. Must be allocated with >= 3
 * elements (see ::G2C_SECTION0_LEN).
 * - listsec0(0) Discipline-GRIB Master Table Number ([Code Table 0.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table0-0.shtml)).
 * - listsec0[1] GRIB Edition Number (currently 2).
 * - listsec0[2] Length of GRIB message.
 * @param listsec1 Pointer to an array that gets the information read
 * from GRIB Identification Section 1. Must be allocated with >= 13
 * elements (see ::G2C_SECTION1_LEN).
 * - listsec1[0] Id of orginating centre ([Table 0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/on388/table0.html)).
 * - listsec1[1] Id of orginating sub-centre ([Table C]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/on388/tablec.html)).
 * - listsec1[2] GRIB Master Tables Version Number ([Table 1.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-0.shtml)).
 * - listsec1[3] GRIB Local Tables Version Number ([Table 1.1]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-1.shtml)).
 * - listsec1[4] Significance of Reference Time ([Table 1.2]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-1.shtml))
 * - listsec1[5] Reference Time - Year (4 digits)
 * - listsec1[6] Reference Time - Month
 * - listsec1[7] Reference Time - Day
 * - listsec1[8] Reference Time - Hour
 * - listsec1[9] Reference Time - Minute
 * - listsec1[10] Reference Time - Second
 * - listsec1[11] Production status of data ([Table 1.3]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-3.shtml)).
 * - listsec1[12] Type of processed data ([Table 1.4]
 *  (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-4.shtml)).
 * @param numfields A pointer that gets the number of gridded fields
 * found in the GRIB message. That is, the number of occurences of
 * Sections 4 - 7.
 * @param numlocal A pointer that gets the number of Local Use
 * Sections (section 2) found in the GRIB message.
 *
 * @returns 0 for success, otherwise:
 * - ::G2_INFO_NO_GRIB Beginning characters "GRIB" not found.
 * - ::G2_INFO_GRIB_VERSION GRIB message is not Edition 2.
 * - ::G2_INFO_NO_SEC1 Could not find Section 1, where expected.
 * - ::G2_INFO_WRONG_END End string "7777" found, but not where expected.
 * - ::G2_INFO_BAD_END End string "7777" not found at end of message.
 * - ::G2_INFO_INVAL_SEC Invalid section number found.
 *
 * @author Stephen Gilbeert @date 2002-10-28
 */
g2int
g2_info(unsigned char *cgrib, g2int *listsec0, g2int *listsec1,
        g2int *numfields, g2int *numlocal)
{
    g2int mapsec1len = 13;
    g2int mapsec1[13] = {2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1};
    g2int i, j, istart, iofst, lengrib, lensec0, lensec1;
    g2int ipos, isecnum, nbits, lensec;

    *numlocal = 0;
    *numfields = 0;

    LOG((2, "g2_info"));

    /*  Check for beginning of GRIB message in the first 100 bytes. */
    istart = -1;
    for (j = 0; j < 100; j++)
    {
        if (cgrib[j] == 'G' && cgrib[j + 1] == 'R' && cgrib[j + 2] == 'I' &&
            cgrib[j + 3] == 'B')
        {
            istart = j;
            break;
        }
    }
    if (istart == -1)
    {
        printf("g2_info:  Beginning characters GRIB not found.");
        return G2_INFO_NO_GRIB;
    }

    LOG((3, "msg found at byte %ld", istart));

    /*  Unpack Section 0 - Indicator Section. */
    iofst = 8 * (istart + 6);
    gbit(cgrib, listsec0, iofst, 8); /* Discipline */
    iofst = iofst + 8;
    gbit(cgrib, &listsec0[1], iofst, 8); /* GRIB edition number */
    iofst = iofst + 8;
    /* iofst = iofst + 32; */
    /* gbit(cgrib, &lengrib, iofst, 32);        /\* Length of GRIB message *\/ */
    gbit(cgrib, &lengrib, iofst, 64); /* Length of GRIB message */
    /* iofst = iofst + 32; */
    iofst = iofst + 64;
    listsec0[2] = lengrib;
    lensec0 = 16;
    ipos = istart + lensec0;
    LOG((3, "unpacked section 0, lengrib %ld now at byte %ld", lengrib, ipos));

    /*  Currently handles only GRIB Edition 2. */
    if (listsec0[1] != 2)
    {
        printf("g2_info: can only decode GRIB edition 2.");
        return G2_INFO_GRIB_VERSION;
    }

    /*  Unpack Section 1 - Identification Section */
    gbit(cgrib, &lensec1, iofst, 32); /* Length of Section 1 */
    iofst = iofst + 32;
    gbit(cgrib, &isecnum, iofst, 8); /* Section number (1) */
    iofst = iofst + 8;
    if (isecnum != 1)
    {
        printf("g2_info: Could not find section 1.");
        return G2_INFO_NO_SEC1;
    }

    /* Unpack each input value in array listsec1 into the
       appropriate number of octets, which are specified in
       corresponding entries in array mapsec1. */
    for (i = 0; i < mapsec1len; i++)
    {
        nbits = mapsec1[i] * 8;
        gbit(cgrib, &listsec1[i], iofst, nbits);
        iofst = iofst + nbits;
    }
    ipos = ipos + lensec1;
    LOG((3, "unpacked section 1, now at byte %ld", ipos));

    /*  Loop through the remaining sections to see if they are
     *  valid. Also count the number of times Section 2 and Section
     *  4 appear. */
    for (;;)
    {
        if (cgrib[ipos] == '7' && cgrib[ipos + 1] == '7' && cgrib[ipos + 2] == '7' &&
            cgrib[ipos + 3] == '7')
        {
            LOG((3, "found 7777 at byte %ld", ipos));
            ipos = ipos + 4;
            if (ipos != (istart + lengrib))
            {
                printf("g2_info: '7777' found, but not where expected.\n");
                return G2_INFO_WRONG_END;
            }
            break;
        }

        iofst = ipos * 8;
        gbit(cgrib, &lensec, iofst, 32); /* Get Length of Section */
        iofst = iofst + 32;
        gbit(cgrib, &isecnum, iofst, 8); /* Get Section number */
        LOG((3, "found section number %ld of length %ld", isecnum, lensec));
        iofst = iofst + 8;
        ipos = ipos + lensec; /* Update beginning of section pointer */
        if (ipos > (istart + lengrib))
        {
            printf("g2_info: '7777'  not found at end of GRIB message.\n");
            return G2_INFO_BAD_END;
        }
        if (isecnum >= 2 && isecnum <= 7)
        {
            /* Increment counter for total number of local sections
             * or fields found. */
            if (isecnum == 2)
                (*numlocal)++;
            else if (isecnum == 4)
                (*numfields)++;
        }
        else
        {
            printf("g2_info: Invalid section number found in GRIB message: %ld\n", isecnum);
            return G2_INFO_INVAL_SEC;
        }
    }

    return 0;
}
