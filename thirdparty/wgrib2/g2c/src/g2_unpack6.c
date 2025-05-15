/** @file
 * @brief Unpack Section 6 (Bit-Map Section) of a GRIB2 message.
 * @author Stephen Gilbert @date 2002-10-31
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack Section 6 (Bit-Map Section) of a GRIB2 message.
 *
 * @param cgrib char array containing Section 6 of the GRIB2 message.
 * @param iofst Bit offset of the beginning of Section 6 in cgrib.
 * @param ngpts Number of grid points specified in the bit-map
 * @param ibmap Bitmap indicator (see Code Table 6.0)
 * - 0 bitmap applies and is included in Section 6.
 * - 1-253 Predefined bitmap applies
 * - 254 Previously defined bitmap applies to this field
 * - 255 Bit map does not apply to this product.
 * @param bmap Pointer to an integer array containing decoded
 * bitmap. (if ibmap=0)
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK6_BAD_BITMAP Unrecognized pre-defined bit-map.
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */
g2int
g2_unpack6(unsigned char *cgrib, g2int *iofst, g2int ngpts, g2int *ibmap,
           g2int **bmap)
{
    g2int j, isecnum;
    g2int *lbmap = 0;
    g2int *intbmap;

    *bmap = NULL;

    *iofst = *iofst + 32;             /* skip Length of Section */
    gbit(cgrib, &isecnum, *iofst, 8); /* Get Section Number */
    *iofst = *iofst + 8;

    if (isecnum != 6)
    {
        fprintf(stderr, "g2_unpack6: Not Section 6 data.\n");
        return G2_UNPACK_BAD_SEC;
    }

    gbit(cgrib, ibmap, *iofst, 8); /* Get bit-map indicator */
    *iofst = *iofst + 8;

    if (*ibmap == 0)
    { /* Unpack bitmap */
        if (ngpts > 0)
            lbmap = calloc(ngpts, sizeof(g2int));
        if (!lbmap)
            return G2_UNPACK_NO_MEM;

        *bmap = lbmap;
        intbmap = calloc(ngpts, sizeof(g2int));
        gbits(cgrib, intbmap, *iofst, 1, 0, ngpts);
        *iofst = *iofst + ngpts;
        for (j = 0; j < ngpts; j++)
            lbmap[j] = (g2int)intbmap[j];
        free(intbmap);
    }

    return G2_NO_ERROR;
}
