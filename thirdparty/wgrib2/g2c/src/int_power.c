/** @file
 * @brief Provide function similar to C pow() power function.
 * @author Wesley Ebisuzaki
 */
#include "grib2_int.h"

/**
 * Function similar to C pow() power function.
 *
 * @param x The base value whose power is to be calculated.
 * @param y The power value.
 *
 * @return x**y
 *
 * @author Wesley Ebisuzaki
 */
double
int_power(double x, g2int y)
{
    double value;

    if (y < 0)
    {
        y = -y;
        x = 1.0 / x;
    }
    value = 1.0;

    while (y)
    {
        if (y & 1)
        {
            value *= x;
        }
        x = x * x;
        y >>= 1;
    }
    return value;
}
