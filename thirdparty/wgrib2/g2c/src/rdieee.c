/** @file
 * @brief Read a list of real values in 32-bit IEEE floating
 * point format.
 * @author Stephen Gilbert @date 2002-10-25
 */
#include "grib2_int.h"

/**
 * Read a list of real values in 32-bit IEEE floating point format.
 *
 * @param rieee g2int array of floating point values in 32-bit IEEE
 * format.
 * @param num Number of floating point values to convert.
 * @param a float array of real values.  a must be allocated with at
 * least 4*num bytes of memory before calling this function.
 *
 * @author Stephen Gilbert @date 2002-10-25
 */
void
rdieee(g2int *rieee, float *a, g2int num)
{

    g2int j;
    g2int isign, iexp, imant;

    float sign, temp;
    static float two23, two126;
    static g2int test = 0;
    uint64_t msk1 = 0x80000000; /* 10000000000000000000000000000000 binary */
    g2int msk2 = 0x7F800000;    /* 01111111100000000000000000000000 binary */
    g2int msk3 = 0x007FFFFF;    /* 00000000011111111111111111111111 binary */

    if (test == 0)
    {
        two23 = (float)int_power(2.0, -23);
        two126 = (float)int_power(2.0, -126);
        test = 1;
    }

    for (j = 0; j < num; j++)
    {
        /*  Extract sign bit, exponent, and mantissa */
        isign = (rieee[j] & msk1) >> 31;
        iexp = (rieee[j] & msk2) >> 23;
        imant = (rieee[j] & msk3);
        /*printf("SAGieee= %ld %ld %ld\n",isign,iexp,imant); */

        sign = 1.0;
        if (isign == 1)
            sign = -1.0;

        if ((iexp > 0) && (iexp < 255))
        {
            temp = (float)int_power(2.0, (iexp - 127));
            a[j] = sign * temp * (1.0 + (two23 * (float)imant));
        }
        else if (iexp == 0)
        {
            if (imant != 0)
                a[j] = sign * two126 * two23 * (float)imant;
            else
                a[j] = sign * 0.0;
        }
        else if (iexp == 255)
            a[j] = sign * (1E+37);
    }
}
