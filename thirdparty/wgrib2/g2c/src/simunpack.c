/** @file
 * @brief Unpack a data field that was packed using a simple packing
 * algorithm.
 * @author Stephen Gilbert @date 2002-10-29
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack a data field that was packed using a simple packing
 * algorithm, using info from the GRIB2 Data Representation Template
 * 5.0.
 *
 * @param cpack pointer to the packed data field.
 * @param idrstmpl pointer to the array of values for Data
 * Representation Template 5.0.
 * @param ndpts The number of data values to unpack.
 * @param fld A pointer that gets the unpacked data values. fld must be
`* allocated with at least ndpts * sizeof(float) bytes before calling
 * this routine.
 *
 * @return 0 for success, error code otherwise.
 *
 * @author Stephen Gilbert @date 2002-10-29
 */
g2int
simunpack(unsigned char *cpack, g2int *idrstmpl, g2int ndpts, float *fld)
{
    g2int *ifld;
    g2int j, nbits;
    float ref, bscale, dscale;

    assert(cpack && idrstmpl && fld);

    LOG((3, "simunpack ndpts %ld idrstmpl: %ld %ld %ld %ld %ld", ndpts, idrstmpl[0],
         idrstmpl[1], idrstmpl[2], idrstmpl[3], idrstmpl[4]));

    rdieee(idrstmpl, &ref, 1);
    bscale = int_power(2.0, idrstmpl[1]);
    dscale = int_power(10.0, -idrstmpl[2]);
    nbits = idrstmpl[3];

    if (!(ifld = calloc(ndpts, sizeof(g2int))))
    {
        fprintf(stderr, "Could not allocate space in simunpack.\n  "
                        "Data field NOT upacked.\n");
        return G2_JPCUNPACK_MEM;
    }

    /* If nbits equals 0, we have a constant field where the reference
     * value is the data value at each gridpoint. */
    if (nbits != 0)
    {
        gbits(cpack, ifld, 0, nbits, 0, ndpts);
        for (j = 0; j < ndpts; j++)
            fld[j] = (((float)ifld[j] * bscale) + ref) * dscale;
    }
    else
    {
        for (j = 0; j < ndpts; j++)
            fld[j] = ref;
    }

    free(ifld);
    return G2_NO_ERROR;
}
