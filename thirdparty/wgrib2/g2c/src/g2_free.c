/**
 * @file
 * @brief Free up memory that was allocated for struct ::gribfield.
 * @author Stephen Gilbeert @date 2002-10-28
 */

#include "grib2_int.h"
#include <stdlib.h>

/**
 * Free memory that was allocated for struct ::gribfield.
 *
 * @param gfld pointer to ::gribfield structure (defined in include file
 * grib2.h) returned from routine g2_getfld().
 *
 * @note This routine must be called to free up memory used by the
 * decode routine, g2_getfld(), when user no longer needs to reference
 * this data.
 *
 * @author Stephen Gilbeert @date 2002-10-28
 */
void
g2_free(gribfield *gfld)
{
    if (gfld->idsect)
        free(gfld->idsect);
    if (gfld->local)
        free(gfld->local);
    if (gfld->list_opt)
        free(gfld->list_opt);
    if (gfld->igdtmpl)
        free(gfld->igdtmpl);
    if (gfld->ipdtmpl)
        free(gfld->ipdtmpl);
    if (gfld->coord_list)
        free(gfld->coord_list);
    if (gfld->idrtmpl)
        free(gfld->idrtmpl);
    if (gfld->bmap)
        free(gfld->bmap);
    if (gfld->fld)
        free(gfld->fld);
    free(gfld);

    return;
}
