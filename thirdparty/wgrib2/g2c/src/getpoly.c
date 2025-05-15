/** @file
 * @brief Return the J, K, and M pentagonal resolution parameters
 * specified in a GRIB Grid Definition Section used spherical harmonic
 * coefficients using GDT 5.50 through 5.53
 * @author Stephen Gilbert @date 2002-12-11
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * This subroutine returns the J, K, and M pentagonal resolution
 * parameters specified in a GRIB Grid Definition Section (GDS) used
 * spherical harmonic coefficients using GDT 5.50 through 5.53.
 *
 * If 51 - complex data spectral packing is used as the data
 * representation template number (see [Table
 * 5.0](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table5-0.shtml)),
 * then the Grid Definition Template Number in section 3 should be one
 * of:
 * - 50: Spherical Harmonic Coefficients (See Template 3.50)
 * - 51: Rotated Spherical Harmonic Coefficients (See Template 3.51)
 * - 52: Stretched Spherical Harmonic Coefficients (See Template 3.52)
 * - 53: Rotated and Stretched Spherical Harmonic Coefficients (See
 * Template 3.53)
 *
 * @param csec3 Character array that contains the packed GRIB2 GDS.
 * @param jj J pentagonal resolution parameter.
 * @param kk K pentagonal resolution parameter.
 * @param mm M pentagonal resolution parameter.
 *
 * @return always returns 0.
 *
 * @note Returns jj, kk, and mm set to zero, if grid template not
 * recognized.
 *
 * @author Stephen Gilbert @date 2002-12-11
 */
g2int
getpoly(unsigned char *csec3, g2int *jj, g2int *kk, g2int *mm)
{

    g2int *igdstmpl, *list_opt;
    g2int *igds;
    g2int iofst, igdtlen, num_opt;

    iofst = 0; /* set offset to beginning of section */
    if (!g2_unpack3(csec3, &iofst, &igds, &igdstmpl, &igdtlen,
                    &list_opt, &num_opt))
    {
        switch (igds[4]) /*  Template number */
        {
        case 50: /* Spherical harmonic coefficients */
        case 51:
        case 52:
        case 53: {
            *jj = igdstmpl[0];
            *kk = igdstmpl[1];
            *mm = igdstmpl[2];
            break;
        }
        default: {
            *jj = 0;
            *kk = 0;
            *mm = 0;
            break;
        }
        } /* end switch */
    }
    else
    {
        *jj = 0;
        *kk = 0;
        *mm = 0;
    }

    if (igds)
        free(igds);
    if (igdstmpl)
        free(igdstmpl);
    if (list_opt)
        free(list_opt);

    return 0;
}
