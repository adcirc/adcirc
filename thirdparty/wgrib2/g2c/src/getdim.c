/** @file
 * @brief Return the dimensions and scanning mode of a grid definition.
 * @author Stephen Gilbert @date 2002-12-11
 */

#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Return the dimensions and scanning mode of a grid definition.
 *
 * @param csec3 Character array that contains the packed GRIB2 GDS.
 * @param width x (or i) dimension of the grid. 0 if grid is not
 * recognized.
 * @param height y (or j) dimension of the grid. 0 if grid is not
 * recognized.
 * @param iscan Scanning mode (see [Table
 * 3.4](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-4.shtml)). 0
 * of grid not recognized.
 *
 * @return Always returns 0.
 *
 * @author Stephen Gilbert @date 2002-12-11
 */
g2int
getdim(unsigned char *csec3, g2int *width, g2int *height, g2int *iscan)
{
    g2int *igdstmpl = NULL, *list_opt = NULL;
    g2int *igds = NULL;
    g2int iofst, igdtlen, num_opt, jerr;

    iofst = 0; /* set offset to beginning of section */
    if (!(jerr = g2_unpack3(csec3, &iofst, &igds, &igdstmpl, &igdtlen,
                            &list_opt, &num_opt)))
    {
        switch (igds[4]) /*  Template number */
        {
        case 0: /* Lat/Lon */
        case 1:
        case 2:
        case 3: {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[18];
            break;
        }
        case 10: /* Mercator */
        {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[15];
            break;
        }
        case 20: /* Polar Stereographic */
        {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[17];
            break;
        }
        case 30: /* Lambert Conformal */
        {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[17];
            break;
        }
        case 40: /* Gaussian */
        case 41:
        case 42:
        case 43: {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[18];
            break;
        }
        case 90: /* Space View/Orthographic */
        {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[16];
            break;
        }
        case 110: /* Equatorial Azimuthal */
        {
            *width = igdstmpl[7];
            *height = igdstmpl[8];
            *iscan = igdstmpl[15];
            break;
        }
        default: {
            *width = 0;
            *height = 0;
            *iscan = 0;
            break;
        }
        } /* end switch */
    }
    else
    {
        *width = 0;
        *height = 0;
    }

    if (igds)
        free(igds);
    if (igdstmpl)
        free(igdstmpl);
    if (list_opt)
        free(list_opt);

    return 0;
}
