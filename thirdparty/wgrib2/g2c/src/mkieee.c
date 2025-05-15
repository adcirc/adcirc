/** @file
 * @brief Store a list of real values in 32-bit IEEE floating point format.
 * @author Stephen Gilbert @date 2002-10-29
 */

#include "grib2_int.h"
#include <math.h>
#include <stdlib.h>

/**
 * Store a list of real values in 32-bit IEEE floating point format.
 *
 * @param a Input array of floating point values.
 * @param num Number of floating point values to convert.
 * @param rieee Output array of data values in 32-bit IEEE format
 * stored in g2int integer array. rieee must be allocated with at
 * least 4*num bytes of memory before calling this function.
 *
 * @author Stephen Gilbert @date 2002-10-29
 */
void
mkieee(float *a, g2int *rieee, g2int num)
{
    g2int j, n, ieee, iexp, imant;
    double atemp;
    static double two23, two126;
    static g2int test = 0;

    if (test == 0)
    {
        two23 = (double)int_power(2.0, 23);
        two126 = (double)int_power(2.0, 126);
        test = 1;
    }

    for (j = 0; j < num; j++)
    {

        ieee = 0;

        if (a[j] == 0.0)
        {
            rieee[j] = ieee;
            continue;
        }

        //  Set Sign bit (bit 31 - leftmost bit).
        if (a[j] < 0.0)
        {
            ieee = 1 << 31;
            atemp = -1.0 * a[j];
        }
        else
        {
            ieee = 0 << 31;
            atemp = a[j];
        }

        //  Determine exponent n with base 2.
        if (atemp >= 1.0)
        {
            n = 0;
            while (int_power(2.0, n + 1) <= atemp)
            {
                n++;
            }
        }
        else
        {
            n = -1;
            while (int_power(2.0, n) > atemp)
                n--;
        }
        iexp = n + 127;
        if (n > 127)
            iexp = 255; // overflow
        if (n < -127)
            iexp = 0;

        //      set exponent bits ( bits 30-23 )
        ieee = ieee | (iexp << 23);

        //  Determine Mantissa
        if (iexp != 255)
        {
            if (iexp != 0)
                atemp = (atemp / int_power(2.0, n)) - 1.0;
            else
                atemp = atemp * two126;
            imant = (g2int)rint(atemp * two23);
        }
        else
        {
            imant = 0;
        }

        //      set mantissa bits ( bits 22-0 )
        ieee = ieee | imant;

        //  Transfer IEEE bit string to rieee array
        rieee[j] = ieee;
    }

    return;
}
