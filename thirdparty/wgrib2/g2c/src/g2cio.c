/**
 * @file
 * @brief File I/O functions for the g2c library.
 * @author Ed Hartnett @date Nov 11, 2022
 */

#include "grib2_int.h"

#define BITSHIFT_7 7   /**< 7 bits. */
#define BITSHIFT_15 15 /**< 15 bits. */
#define BITSHIFT_31 31 /**< 31 bits. */
#define BITSHIFT_63 63 /**< 63 bits. */

/**
 * Read or write a big-endian integer type to an open file, with
 * conversion between native and big-endian format.
 *
 * GRIB2 handles negative numbers in a special way. Instead of storing
 * two-compliments, like every other programmer and computing
 * organization in the world, GRIB2 flips the first bit, then stores
 * the rest of the int as an unsigned number in the remaining 31
 * bits. How exciting!
 *
 * This function takes the excitement out of GRIB2 negative numbers.
 *
 * @param f Pointer to the open FILE.
 * @param write Non-zero if function should write, otherwise function
 * will read.
 * @param g2ctype The type to be read or written.
 * @param var Pointer to the int to be written, or pointer to the
 * storage that gets the int read.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett @date 11/7/22
 */
int
g2c_file_io(FILE *f, int write, int g2ctype, void *var)
{
    void *void_be;
    signed char *bvar = NULL;
    short *svar = NULL;
    int *ivar = NULL;
    long long int *i64var = NULL;
    unsigned long long int int64_be, int64_tmp;
    unsigned int int_be, int_tmp;
    unsigned short short_be, short_tmp;
    unsigned char byte_be, byte_tmp;
    int type_len;

    /* Check inputs. */
    if (!f || !var || g2ctype < G2C_BYTE || g2ctype > G2C_UINT64)
        return G2C_EINVAL;

    switch (g2ctype)
    {
    case G2C_BYTE:
    case G2C_UBYTE:
        type_len = ONE_BYTE;
        bvar = var;
        void_be = &byte_be;
        if (write)
        {
            /* Are we writing a negative number? */
            if (g2ctype == G2C_BYTE && *bvar < 0)
            {
                byte_tmp = -1 * *bvar;         /* Store as positive. */
                byte_tmp |= 1UL << BITSHIFT_7; /* Set sign bit. */
            }
            else
                byte_tmp = *bvar;

            /* Convert result to big-endian. */
            byte_be = byte_tmp;
        }
        break;
    case G2C_SHORT:
    case G2C_USHORT:
        type_len = TWO_BYTES;
        svar = var;
        void_be = &short_be;
        if (write)
        {
            /* Are we writing a negative number? */
            if (g2ctype == G2C_SHORT && *svar < 0)
            {
                short_tmp = -1 * *svar;          /* Store as positive. */
                short_tmp |= 1UL << BITSHIFT_15; /* Set sign bit. */
            }
            else
                short_tmp = *svar;

            /* Convert result to big-endian. */
            short_be = ntohs(short_tmp);
        }
        break;
    case G2C_INT:
    case G2C_UINT:
        type_len = FOUR_BYTES;
        ivar = var;
        void_be = &int_be;
        if (write)
        {
            /* Are we writing a negative number? */
            if (g2ctype == G2C_INT && *ivar < 0)
            {
                int_tmp = -1 * *ivar;          /* Store as positive. */
                int_tmp |= 1UL << BITSHIFT_31; /* Set sign bit. */
            }
            else
                int_tmp = *ivar;

            /* Convert result to big-endian. */
            int_be = ntohl(int_tmp);
        }
        break;
    case G2C_INT64:
    case G2C_UINT64:
        type_len = EIGHT_BYTES;
        i64var = var;
        void_be = &int64_be;
        if (write)
        {
            /* Are we writing a negative number? */
            if (g2ctype == G2C_INT64 && *i64var < 0)
            {
                int64_tmp = -1 * *i64var;         /* Store as positive. */
                int64_tmp |= 1ULL << BITSHIFT_63; /* Set sign bit. */
            }
            else
                int64_tmp = *i64var;

            /* Convert result to big-endian. */
            int64_be = ntoh64(int64_tmp);
        }
        break;
    default:
        return G2C_EBADTYPE;
    }

    if (write)
    {
        if ((fwrite(void_be, type_len, 1, f)) != 1)
            return G2C_EFILE;
    }
    else
    {
        /* Read from the file. */
        if ((fread(void_be, type_len, 1, f)) != 1)
            return G2C_EFILE;

        switch (g2ctype)
        {
        case G2C_BYTE:
        case G2C_UBYTE:
            /* No conversion needed for one-byte values. */
            *bvar = byte_be;

            /* Did we read a negative number? Check the sign bit... */
            if (g2ctype == G2C_BYTE && *bvar & 1 << BITSHIFT_7)
            {
                *bvar &= ~(1UL << BITSHIFT_7); /* Clear sign bit. */
                *bvar *= -1;                   /* Make it negative. */
            }
            break;
        case G2C_SHORT:
        case G2C_USHORT:
            /* Convert from big-endian. */
            *svar = htons(short_be);

            /* Did we read a negative number? Check the sign bit... */
            if (g2ctype == G2C_SHORT && *svar & 1 << BITSHIFT_15)
            {
                *svar &= ~(1UL << BITSHIFT_15); /* Clear sign bit. */
                *svar *= -1;                    /* Make it negative. */
            }
            break;
        case G2C_INT:
        case G2C_UINT:
            /* Convert from big-endian. */
            *ivar = htonl(int_be);

            /* Did we read a negative number? Check the sign bit... */
            if (g2ctype == G2C_INT && *ivar & 1 << BITSHIFT_31)
            {
                *ivar &= ~(1UL << BITSHIFT_31); /* Clear sign bit. */
                *ivar *= -1;                    /* Make it negative. */
            }
            break;
        case G2C_INT64:
        case G2C_UINT64:
            /* Convert from big-endian. */
            *i64var = hton64(int64_be);

            /* Did we read a negative number? Check the sign bit... */
            if (g2ctype == G2C_INT64 && *i64var & 1ULL << BITSHIFT_63)
            {
                *i64var &= ~(1ULL << BITSHIFT_63); /* Clear sign bit. */
                *i64var *= -1;                     /* Make it negative. */
            }
            break;
        default:
            return G2C_EBADTYPE;
        }
    }

    return G2C_NOERROR;
}

/**
 * Read or write a big-endian 4-byte signed int to an open GRIB2 file,
 * with conversion between native and big-endian format, and special
 * GRIB2 handling of negative numbers.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the int.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/7/22
 */
int
g2c_file_io_int(FILE *f, int write, int *var)
{
    return g2c_file_io(f, write, G2C_INT, var);
}

/**
 * Read or write a big-endian 4-byte unsigned int to an open GRIB2
 * file, with conversion between native and big-endian format.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the unsigned int.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/7/22
 */
int
g2c_file_io_uint(FILE *f, int write, unsigned int *var)
{
    return g2c_file_io(f, write, G2C_UINT, var);
}

/**
 * Read or write a big-endian signed short to an open GRIB2 file, with
 * conversion between native and big-endian format, and special GRIB2
 * handling of negative numbers.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the short.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_short(FILE *f, int write, short *var)
{
    return g2c_file_io(f, write, G2C_SHORT, var);
}

/**
 * Read or write a big-endian unsigned short to an open GRIB2 file,
 * with conversion between native and big-endian format.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the unsigned short.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_ushort(FILE *f, int write, unsigned short *var)
{
    return g2c_file_io(f, write, G2C_USHORT, var);
}

/**
 * Read or write a big-endian signed byte to an open GRIB2 file, with
 * conversion between native and big-endian format, and special GRIB2
 * handling of negative numbers.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the byte.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_byte(FILE *f, int write, char *var)
{
    return g2c_file_io(f, write, G2C_BYTE, var);
}

/**
 * Read or write a big-endian unsigned byte to an open GRIB2 file,
 * with conversion between native and big-endian format.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the unsigned byte.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_ubyte(FILE *f, int write, unsigned char *var)
{
    return g2c_file_io(f, write, G2C_UBYTE, var);
}

/**
 * Read or write a big-endian signed long long to an open GRIB2 file, with
 * conversion between native and big-endian format, and special GRIB2
 * handling of negative numbers.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the long long.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_longlong(FILE *f, int write, long long *var)
{
    return g2c_file_io(f, write, G2C_INT64, var);
}

/**
 * Read or write a big-endian unsigned long long to an open GRIB2 file,
 * with conversion between native and big-endian format.
 *
 * @param f Pointer to the open GRIB2 FILE.
 * @param write Non-zero to write, zero to read.
 * @param var Pointer to the unsigned long long.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/13/22
 */
int
g2c_file_io_ulonglong(FILE *f, int write, unsigned long long *var)
{
    return g2c_file_io(f, write, G2C_UINT64, var);
}

/**
 * Read or write a big-endian 4-byte int or unsigned int from or to an
 * open file, with conversion between native and big-endian format,
 * and handling of GRIB negative numbers. This is for template values.
 *
 * With template values, if the map value is negative, then the
 * template value may be negative.
 *
 * @param f Pointer to the open FILE.
 * @param rw_flag Non-zero if function should write, otherwise function
 * will read.
 * @param map The map value for this template item.
 * @param template_value Pointer to the template value to be written,
 * or pointer to the storage that gets the templage value read.
 *
 * @return
 * - :: G2C_NOERROR No error.
 * - :: G2C_EINVAL Invalid input.
 * - :: G2C_EFILE Error reading/writing file.
 *
 * @author Ed Hartnett 11/7/22
 */
int
g2c_file_io_template(FILE *f, int rw_flag, int map, long long int *template_value)
{
    int ret;

    /* Take the absolute value of map[t] because some of the
     * numbers are negative - used to indicate that the
     * cooresponding fields can contain negative data (needed for
     * unpacking). */
    switch (abs(map))
    {
    case ONE_BYTE:
        if (map < 0)
        {
            char my_byte;
            if (rw_flag)
                my_byte = *template_value;
            if ((ret = g2c_file_io_byte(f, rw_flag, &my_byte)))
                return ret;
            if (!rw_flag)
                *template_value = my_byte;
        }
        else
        {
            unsigned char my_ubyte;
            if (rw_flag)
                my_ubyte = *template_value;
            if ((ret = g2c_file_io_ubyte(f, rw_flag, &my_ubyte)))
                return ret;
            if (!rw_flag)
                *template_value = my_ubyte;
        }
        /* FILE_BE_INT1P(f, rw_flag, template_value); */
        break;
    case TWO_BYTES:
        if (map < 0)
        {
            short my_short;
            if (rw_flag)
                my_short = *template_value;
            if ((ret = g2c_file_io_short(f, rw_flag, &my_short)))
                return ret;
            if (!rw_flag)
                *template_value = my_short;
        }
        else
        {
            unsigned short my_ushort;
            if (rw_flag)
                my_ushort = *template_value;
            if ((ret = g2c_file_io_ushort(f, rw_flag, &my_ushort)))
                return ret;
            if (!rw_flag)
                *template_value = my_ushort;
        }
        /* FILE_BE_INT2P(f, rw_flag, template_value); */
        break;
    case FOUR_BYTES:
        if (map < 0)
        {
            int my_int;
            if (rw_flag)
                my_int = *template_value;
            if ((ret = g2c_file_io_int(f, rw_flag, &my_int)))
                return ret;
            if (!rw_flag)
                *template_value = my_int;
        }
        else
        {
            unsigned int my_uint;
            if (rw_flag)
                my_uint = *template_value;
            if ((ret = g2c_file_io_uint(f, rw_flag, &my_uint)))
                return ret;
            if (!rw_flag)
                *template_value = my_uint;
        }
        break;
    default:
        return G2C_EBADTEMPLATE;
    }

    return G2C_NOERROR;
}
