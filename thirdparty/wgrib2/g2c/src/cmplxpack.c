/** @file
 * @brief Pack a data field using a complex packing algorithm.
 * @author Stephen Gilbert @date 2004-08-27
 */

#include "grib2_int.h"

/**
 * Pack up a data field using a complex packing
 * algorithm.
 *
 * This function supports GRIB2
 * complex packing templates with or without spatial differences
 * (i.e. DRTs 5.2 and 5.3). It also fills in GRIB2 Data Representation
 * Template 5.2 or 5.3 with the appropriate values.
 *
 * @param fld Contains the data values to pack.
 * @param ndpts The number of data values in array fld
 * @param idrsnum Data Representation Template number. Must equal 2 or
 * 3.
 * @param idrstmpl Contains the array of values for Data
 * Representation Template 5.2 or 5.3
 * - 0 Reference value - ignored on input, set by compack routine.
 * - 1 Binary Scale Factor
 * - 2 Decimal Scale Factor
 * - 6 Missing value management
 * - 7 Primary missing value
 * - 8 Secondary missing value
 * - 16 Order of Spatial Differencing  ( 1 or 2 )
 * @param cpack The packed data field.
 * @param lcpack length of packed field cpack. Will be set to -1 if
 * missing value management field is not 1 or 2.
 *
 * @author Stephen Gilbert @date 2004-08-27
 */
void
cmplxpack(float *fld, g2int ndpts, g2int idrsnum, g2int *idrstmpl,
          unsigned char *cpack, g2int *lcpack)
{
    if (idrstmpl[6] == 0) /* No internal missing values */
        compack(fld, ndpts, idrsnum, idrstmpl, cpack, lcpack);
    else if (idrstmpl[6] == 1 || idrstmpl[6] == 2)
        misspack(fld, ndpts, idrsnum, idrstmpl, cpack, lcpack);
    else
    {
        printf("cmplxpack: Don:t recognize Missing value option.");
        *lcpack = -1;
    }
}
