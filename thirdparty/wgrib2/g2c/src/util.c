/**
 * @file
 * @brief Internal utility functions for the g2c library.
 * @author Ed Hartnett @date Nov 11, 2021
 */

#include "grib2_int.h"
#include <stdarg.h>

/**
 * Check for 'GRIB' at the beginning of a GRIB message, and check to
 * see if the message is already terminated with '7777'.
 *
 * @param cgrib Buffer that contains the GRIB message.
 * @param lencurr Pointer that gets the length of the GRIB message.
 * @param verbose If non-zero, print any error messages to stdout.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOTGRIB GRIB header not found.
 * - ::G2C_EMSGCOMPLETE GRIB message already complete.
 *
 * @author Ed Hartnett @date Nov 11, 2021
 */
int
g2c_check_msg(unsigned char *cgrib, g2int *lencurr, int verbose)
{
    unsigned char G = 0x47;     /* 'G' */
    unsigned char R = 0x52;     /* 'R' */
    unsigned char I = 0x49;     /* 'I' */
    unsigned char B = 0x42;     /* 'B' */
    unsigned char seven = 0x37; /* '7' */

    assert(cgrib && lencurr);

    /* Check to see if beginning of GRIB message exists. */
    if (cgrib[0] != G || cgrib[1] != R || cgrib[2] != I || cgrib[3] != B)
    {
        if (verbose)
            printf("GRIB not found in given message. A call to routine g2_create() "
                   "is required to to initialize GRIB messge.\n");
        return G2C_ENOTGRIB;
    }

    /* Get current length of GRIB message. */
    gbit(cgrib, lencurr, 96, 32);

    /* Check to see if GRIB message is already complete. */
    if (cgrib[*lencurr - 4] == seven && cgrib[*lencurr - 3] == seven &&
        cgrib[*lencurr - 2] == seven && cgrib[*lencurr - 1] == seven)
    {
        if (verbose)
            printf("GRIB message already complete.  Cannot add new section.\n");
        return G2C_EMSGCOMPLETE;
    }

    return G2C_NOERROR;
}

#ifdef LOGGING
/* This is the severity level of messages which will be logged. Use
   severity 0 for errors, 1 for important log messages, 2 for less
   important, etc. */
int g2_log_level = -1;

/* This function prints out a message, if the severity of
 * the message is lower than the global g2_log_level. To use it, do
 * something like this:
 *
 * g2_log(0, "this computer will explode in %d seconds", i);
 *
 * After the first arg (the severity), use the rest like a normal
 * printf statement. Output will appear on stderr.
 *
 * This function is not included in the build unless NCEPLIBS-g2c was
 * built with -DLOGGING.
 *
 * Ed Hartnett
 */
void
g2_log(int severity, const char *fmt, ...)
{
    va_list argp;
    int t;
    FILE *f = stderr;

    /* If the severity is greater than the log level, we don't print
     * this message. */
    if (severity > g2_log_level)
        return;

    /* If the severity is zero, this is an error. Otherwise insert that
       many tabs before the message. */
    if (!severity)
        fprintf(f, "ERROR: ");
    for (t = 0; t < severity; t++)
        fprintf(f, "\t");

    /* Print out the variable list of args with vprintf. */
    va_start(argp, fmt);
    vfprintf(f, fmt, argp);
    va_end(argp);

    /* Put on a final linefeed. */
    fprintf(f, "\n");
    fflush(f);
}
#endif /* LOGGING */

/**
 * Use this to set the global log level.
 *
 * Settings:
 * - -1 turn off all logging.
 * - 0 show only errors.
 * - 1 output useful as verbose to utilities.
 * - 2 or 3 shows some/all calls to top-level functions.
 * - 4+ ever greater levels of detail.
 *
 * If logging is not enabled when building NCEPLIBS-g2c, this function
 * will do nothing.
 *
 * @param new_level The new logging level.
 *
 * @return ::G2C_NOERROR No error.
 * @author Ed Hartnett
 */
int
g2c_set_log_level(int new_level)
{
#ifdef LOGGING
    /* Remember the new level. */
    g2_log_level = new_level;

    LOG((1, "log_level changed to %d", g2_log_level));
#endif
    return G2C_NOERROR;
}
