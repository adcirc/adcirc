/** @file
 * @brief Search a file for the next GRIB message.
 * @author Stephen Gilbert @date 2002-10-28
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-28 | GILBERT | Modified from Iredell's skgb subroutine
 * 2009-01-16 | VUONG | Changed  lskip to 4 instead of sizof(g2int)
 * 2022-09-11 | Hartnett | Added g2c_seekgb() function.
 *
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Search a file for the next GRIB Message.
 *
 * The search is done starting at byte offset iseek of the file
 * referenced by lugb for mseek bytes at a time. If found, the
 * starting position and length of the message are returned in lskip
 * and lgrib, respectively. The search is terminated when an EOF or
 * I/O error is encountered.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-28 | GILBERT | Modified from Iredell's skgb subroutine
 * 2009-01-16 | VUONG | Changed  lskip to 4 instead of sizof(g2int)
 *
 * @param lugb FILE pointer for the file to search. File must be
 * opened before this routine is called.
 * @param iseek The number of bytes in the file to skip before search.
 * @param mseek The maximum number of bytes to search at a time (must
 * be at least 16, but larger numbers like 4092 will result in better
 * perfomance).
 * @param lskip Pointer that gets the number of bytes to skip from the
 * beggining of the file to where the GRIB message starts.
 * @param lgrib Pointer that gets the number of bytes in message (set
 * to 0, if no message found).
 *
 * @author Stephen Gilbert @date 2002-10-28
 */
void
seekgb(FILE *lugb, g2int iseek, g2int mseek, g2int *lskip, g2int *lgrib)
{
    g2int k, k4, ipos, nread, lim, start, vers, lengrib;
    int end;
    unsigned char *cbuf;

    LOG((3, "seekgb iseek %ld mseek %ld", iseek, mseek));

    *lgrib = 0;
    cbuf = (unsigned char *)malloc(mseek);
    nread = mseek;
    ipos = iseek;

    /* Loop until grib message is found. */
    while (*lgrib == 0 && nread == mseek)
    {
        /* Read partial section. */
        fseek(lugb, ipos, SEEK_SET);
        nread = fread(cbuf, sizeof(unsigned char), mseek, lugb);
        lim = nread - 8;

        /* Look for 'grib...' in partial section. */
        for (k = 0; k < lim; k++)
        {
            /* Look at the first 4 bytes - should be 'GRIB'. */
            gbit(cbuf, &start, k * BYTE, 4 * BYTE);

            /* Look at the 8th byte, it has the GRIB version. */
            gbit(cbuf, &vers, (k + 7) * BYTE, 1 * BYTE);

            /* If the message starts with 'GRIB', and is version 1 or
             * 2, then this is a GRIB message. */
            if (start == 1196575042 && (vers == 1 || vers == 2))
            {
                /* Find the length of the message. */
                if (vers == 1)
                    gbit(cbuf, &lengrib, (k + 4) * BYTE, 3 * BYTE);
                if (vers == 2)
                    gbit(cbuf, &lengrib, (k + 12) * BYTE, 4 * BYTE);
                LOG((4, "lengrib %ld", lengrib));

                /* Read the last 4 bytesof the message. */
                fseek(lugb, ipos + k + lengrib - 4, SEEK_SET);
                k4 = fread(&end, 4, 1, lugb);

                /* Look for '7777' at end of grib message. */
                if (k4 == 1 && end == 926365495)
                {
                    /* GRIB message found. */
                    *lskip = ipos + k;
                    *lgrib = lengrib;
                    LOG((4, "found end of message lengrib %ld", lengrib));
                    break;
                }
            }
        }
        ipos = ipos + lim;
    }

    free(cbuf);
}
