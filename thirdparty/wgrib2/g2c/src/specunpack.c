/** @file
 * @brief Unpack a spectral data field that was packed using the
 * complex packing algorithm for spherical harmonic data
 * @author Stephen Gilbert @date 2000-06-21
 */
#include "grib2_int.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack a spectral data field that was packed using the complex
 * packing algorithm for spherical harmonic data as defined in the
 * GRIB2 documention, using info from the GRIB2 Data Representation
 * [Template
 * 5.51](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-51.shtml).
 *
 * @param cpack pointer to the packed data field.
 * @param idrstmpl pointer to the array of values for Data
 * Representation Template 5.51.
 * @param ndpts The number of data values to unpack (real and
 * imaginary parts).
 * @param JJ pentagonal resolution parameter.
 * @param KK pentagonal resolution parameter.
 * @param MM pentagonal resolution parameter.
 * @param fld Contains the unpacked data values. fld must be allocated
 * with at least ndpts * sizeof(float) bytes before calling this
 * routine.
 *
 * @return 0 for success, -3 for wrong type.
 *
 * @author Stephen Gilbert @date 2000-06-21
 */
g2int
specunpack(unsigned char *cpack, g2int *idrstmpl, g2int ndpts, g2int JJ,
           g2int KK, g2int MM, float *fld)
{
    g2int *ifld, j, iofst, nbits;
    float ref, bscale, dscale, *unpk;
    float *pscale, tscale;
    g2int Js, Ks, Ms, Ts, Ns, Nm, n, m;
    g2int inc, incu, incp;

    rdieee(idrstmpl + 0, &ref, 1);
    bscale = int_power(2.0, idrstmpl[1]);
    dscale = int_power(10.0, -idrstmpl[2]);
    nbits = idrstmpl[3];
    Js = idrstmpl[5];
    Ks = idrstmpl[6];
    Ms = idrstmpl[7];
    Ts = idrstmpl[8];

    if (idrstmpl[9] == 1)
    { /* unpacked floats are 32-bit IEEE */

        unpk = malloc(ndpts * sizeof(float));
        ifld = malloc(ndpts * sizeof(g2int));

        gbits(cpack, ifld, 0, 32, 0, Ts);
        iofst = 32 * Ts;
        rdieee(ifld, unpk, Ts);                          /* read IEEE unpacked floats */
        gbits(cpack, ifld, iofst, nbits, 0, ndpts - Ts); /* unpack scaled data */

        /* Calculate Laplacian scaling factors for each possible wave
         * number. */
        pscale = malloc((JJ + MM + 1) * sizeof(float));
        tscale = idrstmpl[4] * 1E-6;
        for (n = Js; n <= JJ + MM; n++)
            pscale[n] = pow((float)(n * (n + 1)), -tscale);

        /* Assemble spectral coeffs back to original order. */
        inc = 0;
        incu = 0;
        incp = 0;
        for (m = 0; m <= MM; m++)
        {
            Nm = JJ; /* triangular or trapezoidal */
            if (KK == JJ + MM)
                Nm = JJ + m; /* rhombodial */
            Ns = Js;         /* triangular or trapezoidal */
            if (Ks == Js + Ms)
                Ns = Js + m; /* rhombodial */
            for (n = m; n <= Nm; n++)
            {
                if (n <= Ns && m <= Ms)
                {                              /* grab unpacked value */
                    fld[inc++] = unpk[incu++]; /* real part */
                    fld[inc++] = unpk[incu++]; /* imaginary part */
                }
                else
                { /* Calc coeff from packed value */
                    fld[inc++] = (((float)ifld[incp++] * bscale) + ref) *
                                 dscale * pscale[n]; /* real part */
                    fld[inc++] = (((float)ifld[incp++] * bscale) + ref) *
                                 dscale * pscale[n]; /* imaginary part */
                }
            }
        }

        free(pscale);
        free(unpk);
        free(ifld);
    }
    else
    {
        printf("specunpack: Cannot handle 64 or 128-bit floats.\n");
        for (j = 0; j < ndpts; j++)
            fld[j] = 0.0;
        return G2_SPECUNPACK_TYPE;
    }

    return G2_NO_ERROR;
}
