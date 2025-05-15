/** @file
 * @brief Initialize a new GRIB2 message and pack GRIB2 sections 0
 * (Indicator Section) and 1 (Identification Section).
 * @author Stephen Gilbert @date 2002-10-31
 */

#include "grib2_int.h"
#include <stdio.h>

#define MAPSEC1LEN 13 /**< Length of Map Section 1. */
#define LENSEC0 16    /**< Length of GRIB Section 0. */

/**
 * Initialize a new GRIB2 message and pack GRIB2
 * [Section 0 (Indicator
 * Section)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect0.shtml)
 * and [Section 1 (Identification
 * Section)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect1.shtml).
 *
 * This routine is used with routines g2_addlocal(), g2_addgrid(),
 * g2_addfield(), and g2_gribend() to create a complete GRIB2 message.
 * g2_create() must be called first to initialize a new GRIB2 message.
 * A call to g2_gribend() is required to complete GRIB2 message after
 * all fields have been added.
 *
 * @param[in] cgrib Character array to contain the GRIB2 message. Must
 * be allocated large enough to store the entire GRIB2 message.
 * @param[in] listsec0 Contains information needed for GRIB Indicator
 * Section 0. Must be dimensioned >= 2.
 * - listsec0[0] Discipline-GRIB Master Table Number ([Code Table 0.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table0-0.shtml)).
 * - listsec0[1] GRIB Edition Number (currently 2).
 * @param[in] listsec1 Contains information needed for GRIB
 * Identification Section 1. Must be dimensioned >= 13.
 * - listsec1[0] Id of orginating centre ([Table 0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/on388/table0.html)).
 * - listsec1[1] Id of orginating sub-centre ([Table C]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/on388/tablec.html)).
 * - listsec1[2] GRIB Master Tables Version Number ([Table 1.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-0.shtml)).
 * - listsec1[3] GRIB Local Tables Version Number ([Table 1.1]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-1.shtml)).
 * - listsec1[4] Significance of Reference Time ([Table 1.2]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table1-2.shtml))
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
 *
 * @return
 * - > 0 Current size of new GRIB2 message
 * - ::G2_CREATE_GRIB_VERSION Tried to use for version other than GRIB
 * Edition 2
 *
 * This routine is intended for use with routines g2_addlocal(),
 * g2_addgrid(), g2_addfield(), and g2_gribend() to create a complete
 * GRIB2 message.
 *
 * @author Stephen Gilbeert @date 2002-10-31
 */
g2int
g2_create(unsigned char *cgrib, g2int *listsec0, g2int *listsec1)
{
    g2int zero = 0, one = 1;

    /* The mapsec1 array tells us how many bytes are used in the GRIB
     * message by each element of the listsec1 array. For example the
     * first two elements of listsec1 are identifcation of originating
     * section and sub-section - these are each 2-byte entries in the
     * GRIB section 1 table, the IDENTIFICATION SECTION. (See
     * https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect1.shtml). */
    g2int mapsec1[MAPSEC1LEN] = {2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1};
    g2int i, lensec1, iofst, ibeg, nbits, len;

    /* Only GRIB Edition 2 is acceptable. */
    if (listsec0[1] != 2)
    {
        printf("g2_create: can only code GRIB edition 2.");
        return G2_CREATE_GRIB_VERSION;
    }

    /* Pack Section 0 - Indicator Section (except for total length of
     * GRIB message). */
    cgrib[0] = 0x47;                  /* 'G' */
    cgrib[1] = 0x52;                  /* 'R' */
    cgrib[2] = 0x49;                  /* 'I' */
    cgrib[3] = 0x42;                  /* 'B' */
    sbit(cgrib, &zero, 32, 16);       /* reserved for future use */
    sbit(cgrib, listsec0 + 0, 48, 8); /* Discipline */
    sbit(cgrib, listsec0 + 1, 56, 8); /* GRIB edition number */

    /* Pack Section 1 - Identification Section. */
    ibeg = LENSEC0 * 8;          /* Calculate offset for beginning of section 1. */
    iofst = ibeg + 32;           /* Leave space for length of section. */
    sbit(cgrib, &one, iofst, 8); /* Store section number (1). */
    iofst = iofst + 8;

    /* Pack up each input value in array listsec1 into the the
     * appropriate number of octets, which are specified in
     * corresponding entries in array mapsec1. */
    for (i = 0; i < MAPSEC1LEN; i++)
    {
        nbits = mapsec1[i] * 8;
        sbit(cgrib, listsec1 + i, iofst, nbits);
        iofst = iofst + nbits;
    }

    /* Calculate length of section 1 and store it in octets 1-4 of
     * section 1. */
    lensec1 = (iofst - ibeg) / 8;
    sbit(cgrib, &lensec1, ibeg, 32);

    /* Put current byte total of message into Section 0. */
    sbit(cgrib, &zero, 64, 32);
    len = LENSEC0 + lensec1;
    sbit(cgrib, &len, 96, 32);
    return (len);
}
