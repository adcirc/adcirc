/** @file
 * @brief Check the Data Representation Template to see if missing
 * value management is used, and returns the missing value(s) in the
 * data field.
 * @author Stephen Gilbeert @date 2004-12-16
 */

#include "grib2_int.h"

/**
 * Check the Data Representation Template to see if
 * missing value management is used, and return the missing value(s)
 * in the data field.
 *
 * @param gfld pointer to ::gribfield structure.
 * @param rmiss List of the missing values used.
 * @param nmiss NUmber of the missing values included in the
 * field. rmiss must be allocated in the calling program with enough
 * space hold all the missing values.
 *
 * @author Stephen Gilbeert @date 2004-12-16
 */
void
g2_miss(gribfield *gfld, float *rmiss, int *nmiss)
{
    g2int itype;

    /* Missing value management currnetly only used in DRT's 5.2 and
     * 5.3. */
    if (gfld->idrtnum != 2 && gfld->idrtnum != 3)
    {
        *nmiss = 0;
        return;
    }

    itype = gfld->idrtmpl[4];
    if (gfld->idrtmpl[6] == 1)
    {
        *nmiss = 1;
        if (itype == 0)
            rdieee(gfld->idrtmpl + 7, rmiss + 0, 1);
        else
            rmiss[0] = (float)gfld->idrtmpl[7];
    }
    else if (gfld->idrtmpl[6] == 2)
    {
        *nmiss = 2;
        if (itype == 0)
        {
            rdieee(gfld->idrtmpl + 7, rmiss, 1);
            rdieee(gfld->idrtmpl + 8, rmiss + 1, 1);
        }
        else
        {
            rmiss[0] = (float)gfld->idrtmpl[7];
            rmiss[1] = (float)gfld->idrtmpl[8];
        }
    }
    else
    {
        *nmiss = 0;
    }
}
