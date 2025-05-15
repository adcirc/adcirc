/**
 * @file
 * @brief Write summary output to a file, as is done by utility
 * grb2index.
 *
 * @author Ed Hartnett @date Sep 17, 2022
 */
#include "grib2_int.h"
#include <math.h>
#include <stdarg.h>

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** Pointer to the list of code tables. */
extern G2C_CODE_TABLE_T *g2c_table;

/** Maximum length of "Type of first fixed surface" string. */
#define G2C_MAX_TYPE_OF_FIXED_SURFACE_LEN 40

#define G2C_DATE_TIME_LEN 25 /**< Size of the date/time string. */

/**
 * Convert date and time from GRIB2 info to string output.
 *
 * @note This function is based on the prvtime() function in
 * NCEPLIBS-grib_utils.
 *
 * @param[in] ipdtn Product Definition Template Number ([Code Table
 * 4.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-0.shtml)).
 * @param ipdtmpl Array of data values for the Product Definition
 * Template specified by ipdtn.
 * @param year year
 * @param month month
 * @param day day
 * @param hour hour
 * @param minute minute
 * @param second second
 * @param tabbrev Character array that will get the date and time
 * string. Must be of length 100.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date Sep 28, 2022
 */
int
g2c_get_datetime(int ipdtn, long long int *ipdtmpl, short year, unsigned char month, unsigned char day,
                 unsigned char hour, unsigned char minute, unsigned char second, char *tabbrev)
{
    int iutpos, iutpos2, iunit, iunit2;
    char tunit[G2C_DATE_TIME_LEN], tunit2[G2C_DATE_TIME_LEN];
    char reftime[G2C_DATE_TIME_LEN], endtime[G2C_DATE_TIME_LEN], tmpval2[G2C_DATE_TIME_LEN];
    int itemp, itemp2, is;
    /* character(len = 16) :: reftime, endtime */
    /* character(len = 10) :: tmpval, tmpval2 */
    /* character(len = 10) :: tunit, tunit2 */
    /* integer, dimension(200) :: ipos, ipos2 */
    /* integer :: is, itemp, itemp2, iunit, iuni2t2, iunit2, iutpos, iutpos2, j */

    /* data ipos  /7*0, 16, 23, 17, 19, 18, 32, 31, 27*0, 17, 20, 0, 0, 22,  & */
    /*      25, 43*0, 23, 109*0/ */

    int ipos[200] = {
        0, 0, 0, 0, 0, 0, 0, 16, 23, 17, 19, 18, 32, 31, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 17, 20, 0, 0, 22, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    /* data ipos2 /7*0, 26, 33, 27, 29, 28, 42, 41, 27*0, 22, 30, 0, 0, 32,  & */
    /*      35, 43*0, 33, 109*0/ */

    int ipos2[200] = {
        0, 0, 0, 0, 0, 0, 0, 26, 33, 27,
        29, 28, 42, 41, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 22, 30, 0, 0, 32, 35, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        33, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    /* tabbrev(1:100) = " " */

    /* Check inputs. */
    assert(ipdtn >= 0 && ipdtmpl && tabbrev);

    /* The unit of time range is stored in the product template, but
     * stored at different places for different template numbers. Find
     * where it is stored for this product template number. */
    if (ipdtn <= 15 || ipdtn == 32 || ipdtn == 50 || ipdtn == 51 || ipdtn == 91)
        iutpos = 7;
    else if (ipdtn >= 40 && ipdtn <= 43)
        iutpos = 8;
    else if (ipdtn >= 44 && ipdtn <= 47)
        iutpos = 13;
    else if (ipdtn == 48)
        iutpos = 18;
    else if (ipdtn == 52)
        iutpos = 10;
    else
        iutpos = 7;

    /* Determine first unit of time range. */
    switch (ipdtmpl[iutpos])
    {
    case 0:
        strcpy(tunit, "minute");
        iunit = 1;
        break;
    case 1:
        strcpy(tunit, "hour");
        iunit = 1;
        break;
    case 2:
        strcpy(tunit, "day");
        iunit = 1;
        break;
    case 3:
        strcpy(tunit, "month");
        iunit = 1;
        break;
    case 4:
        strcpy(tunit, "year");
        iunit = 1;
        break;
    case 10:
        strcpy(tunit, "hour");
        iunit = 3;
        break;
    case 11:
        strcpy(tunit, "hour");
        iunit = 6;
        break;
    default:
        strcpy(tunit, "hour");
        iunit = 1;
    }

    /* Determine second unit of time range. */
    if (ipdtn > 0)
        iutpos2 = ipos2[ipdtn - 1];
    else
        iutpos2 = 0;
    switch (ipdtmpl[iutpos2])
    {
    case 0:
        strcpy(tunit2, "minute");
        iunit2 = 1;
        break;
    case 1:
        strcpy(tunit2, "hour");
        iunit2 = 1;
        break;
    case 2:
        strcpy(tunit2, "day");
        iunit2 = 1;
        break;
    case 3:
        strcpy(tunit2, "month");
        iunit2 = 1;
        break;
    case 4:
        strcpy(tunit2, "year");
        iunit2 = 1;
        break;
    case 10:
        strcpy(tunit2, "hour");
        iunit2 = 3;
        break;
    case 11:
        strcpy(tunit2, "hour");
        iunit2 = 6;
        break;
    default:
        strcpy(tunit2, "hour");
        iunit2 = 1;
    }

    /* Write a string with the date and time from section 1 of the message. */
    sprintf(reftime, "%4.4d%2.2d%2.2d%2.2d:%2.2d:%2.2d", year, month, day, hour, minute, second);

    itemp = abs((int)ipdtmpl[iutpos + 1]) * iunit;

    /* write(tmpval, '(I0)') itemp */

    sprintf(tabbrev, "valid at  %lld", ipdtmpl[iutpos + 1]);
    /* write(tabbrev, fmt = '("valid at  ", i4)') ipdtmpl(iutpos + 1) */

    /* Determine Reference Time: Year, Month, Day, Hour, Minute, Second. */

    if ((ipdtn >= 0 && ipdtn <= 7) || ipdtn == 15 || ipdtn == 20 || (ipdtn >= 30 && ipdtn <= 32) || ipdtn == 40 ||
        ipdtn == 41 || ipdtn == 44 || ipdtn == 45 || ipdtn == 48 || (ipdtn >= 50 && ipdtn <= 52)) /*  Point in time. */
    {
        sprintf(tabbrev, "valid  %d %s after %s", itemp, tunit, reftime);
    }
    else
    {
        is = ipos[ipdtn - 1] - 1; /* Continuous time interval. */
        sprintf(endtime, "%4.4d%2.2d%2.2d%2.2d:%2.2d:%2.2d", (int)ipdtmpl[is], (int)ipdtmpl[is + 1],
                (int)ipdtmpl[is + 2], (int)ipdtmpl[is + 3], (int)ipdtmpl[is + 4], (int)ipdtmpl[is + 5]);
        /*    write(endtime, fmt = '(i4,3i2.2,":",i2.2,":",i2.2)') (ipdtmpl(j), j = is, is + 5) */
        itemp2 = abs((int)ipdtmpl[iutpos2]) * iunit2;
        itemp2 = itemp + itemp2;
        sprintf(tmpval2, "%d", itemp2);
        if (ipdtn == 8 && ipdtmpl[8] < 0)
        {
            /*       tabbrev = "(" // trim(tmpval) // " -" // trim(tmpval2) // ") valid  " // trim(tmpval) // " " // trim(tunit) // " before " // reftime // " to " //endtime */
            sprintf(tabbrev, "(%d -%d) valid  %d %s before %s to %s", itemp, itemp2, itemp, tunit, reftime, endtime);
        }
        else if ((ipdtn >= 8 && ipdtn <= 14) || (ipdtn >= 42 && ipdtn <= 47) ||
                 ipdtn == 91) /* Continuous time interval */
        {
            /*       write(tmpval2, '(I0)') itemp2 */
            /*       tabbrev = "(" // trim(tmpval) // " -" // trim(tmpval2) // " hr) valid  " // trim(tmpval) // " " // trim(tunit) // " after " // reftime // " to " // endtime */
            sprintf(tabbrev, "(%d -%d hr) valid  %d %s after %s to %s", itemp, itemp2, itemp, tunit, reftime, endtime);
        }
    }

    return G2C_NOERROR;
}

/**
 * Format the level string for degrib2 output.
 *
 * @param cbuf Character buffer of length 11 which will get the
 * formatted level description.
 * @param ival Integer value which may need scaling.
 * @param scale Scaling factor information.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 11/16/22
 */
static int
format_level(char *cbuf, int ival, int scale)
{
    char tmpcbuf[37];

    assert(cbuf);

    if (!scale)
        sprintf(cbuf, "%d", ival);
    else
    {
        float rval;

        rval = (float)ival * pow(10, -1 * scale);
        if (rval == (int)rval)
        {
            /* This adds a space in front of the number, to match
             * behavor of the Fortran degrib2 utility. */
            sprintf(cbuf, " %d", (int)rval);
        }
        else
        {
            char fmt[37];
            sprintf(fmt, "%s%d%s", "%.", (int)abs(scale), "f");
            sprintf(tmpcbuf, fmt, rval);
            strcpy(cbuf, &tmpcbuf[1]);
        }
    }
    return G2C_NOERROR;
}

/**
 * Determine the string that describes the level information, given
 * the GRIB2 Product Definition Template information.
 *
 * This subroutine finds the "Type of first fixed surface" (see [Code
 * table
 * 4.5](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp4-0.shtml))
 * using the product definition template value and product template
 * array. The location in the template array varies depending on the
 * template number.
 *
 * @param ipdtn Product Definition Template Number ([Code Table
 * 4.0]
 * (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-0.shtml)).
 * @param ipdtmpl Array of data values for the Product Definition
 * Template specified by ipdtn.
 * @param level_desc Character array which gets the string
 * describing the level. Must be of length 40.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date Sep 28, 2022
 */
int
g2c_get_level_desc(int ipdtn, long long int *ipdtmpl, char *level_desc)
{
    int ipos;
    char tmpval1[37];
    char tmpval2[37];
    int ret;

    /* Check inputs. */
    assert(ipdtn >= 0 && ipdtmpl && level_desc);

    /* Use the template number to determine which element of the
     * product template array holds the "Type of first fixed surface"
     * value (which will be stored in ipos). */
    if (ipdtn <= 15)
        ipos = 9;
    else if (ipdtn >= 40 && ipdtn <= 43)
        ipos = 10;
    else if (ipdtn >= 44 && ipdtn <= 47)
        ipos = 15;
    else if (ipdtn == 48)
        ipos = 20;
    else if (ipdtn >= 50 && ipdtn <= 51)
        ipos = 9;
    else if (ipdtn == 52)
        ipos = 12;
    else if (ipdtn == 91)
        ipos = 9;
    else
        ipos = 9;

    /* Pressure Level. */
    if (ipdtmpl[ipos] == 100 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1] + 2)))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "mb");
    }
    /* Pressure Layer. */
    else if (ipdtmpl[ipos] == 100 && ipdtmpl[ipos + 3] == 100)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1] + 2)))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4] + 2)))
            return ret;
        sprintf(level_desc, "%s - %s mb", tmpval1, tmpval2);
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1) + 2) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4) + 2) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" mb" */
    }
    else if (ipdtmpl[ipos] == 101)
    {
        strcpy(level_desc, " Mean Sea Level");
    }
    /* Altitude above MSL. */
    else if (ipdtmpl[ipos] == 102 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "m above MSL");
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* level_desc = trim(tmpval1)//" m above MSL" */
    }
    /* Height above Ground. */
    else if (ipdtmpl[ipos] == 103 && ipdtmpl[ipos + 3] == 255)
    {
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* level_desc = trim(tmpval1)//" m above ground" */
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "m above ground");
    }
    /* Height above Ground. */
    else if (ipdtmpl[ipos] == 103 && ipdtmpl[ipos + 3] == 103)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4])))
            return ret;
        sprintf(level_desc, "%s - %s m AGL", tmpval1, tmpval2);
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4)) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" m AGL" */
    }
    /* Sigma Level. */
    else if (ipdtmpl[ipos] == 104 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "sigma");
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* level_desc = trim(tmpval1)//" sigma" */
    }
    /* Sigma Layer. */
    else if (ipdtmpl[ipos] == 104 && ipdtmpl[ipos + 3] == 104)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4])))
            return ret;
        sprintf(level_desc, "%s - %s sigma", tmpval1, tmpval2);
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4)) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" sigma" */
    }
    /* Hybrid Level. */
    else if (ipdtmpl[ipos] == 105 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "hybrid lvl");
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* level_desc = trim(tmpval1)//" hybrid lvl" */
    }
    /* Hybrid Level. */
    else if (ipdtmpl[ipos] == 105 && ipdtmpl[ipos + 3] == 105)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4])))
            return ret;
        sprintf(level_desc, "%s - %s hybrid lvl", tmpval1, tmpval2);
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4)) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" hybrid lvl" */
    }
    /* Depth Below Land Sfc. */
    else if (ipdtmpl[ipos] == 106 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "m below land");
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* level_desc = trim(tmpval1)//" m below land" */
    }
    /* Depth Below Land Sfc Layer. */
    else if (ipdtmpl[ipos] == 106 && ipdtmpl[ipos + 3] == 106)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4])))
            return ret;
        sprintf(level_desc, "%s - %s m DBLY", tmpval1, tmpval2);
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4)) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" m DBLY" */
    }
    else if (ipdtmpl[ipos] == 107)
    {
        strcpy(level_desc, " Isentropic level");
    }
    /* Press Diff from Ground Layer. */
    else if (ipdtmpl[ipos] == 108 && ipdtmpl[ipos + 3] == 108)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1] + 2)))
            return ret;
        if ((ret = format_level(tmpval2, ipdtmpl[ipos + 5], ipdtmpl[ipos + 4] + 2)))
            return ret;
        sprintf(level_desc, "%s - %s mb SPDY", tmpval1, tmpval2);
        /* write(tmpval1, *) ipdtmpl(ipos + 2)/100. */
        /* write(tmpval2, *) ipdtmpl(ipos + 5)/100. */
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1) + 2) */
        /* call frmt(tmpval2, ipdtmpl(ipos + 5), ipdtmpl(ipos + 4) + 2) */
        /* level_desc = trim(tmpval1)//" - "//trim(tmpval2)//" mb SPDY" */
    }
    else if (ipdtmpl[ipos] == 110)
    {
        strcpy(level_desc, " Layer bet 2-hyb lvl");
    }
    /* Potential Vorticity Sfc. */
    else if (ipdtmpl[ipos] == 109 && ipdtmpl[ipos + 3] == 255)
    {
        if ((ret = format_level(tmpval1, ipdtmpl[ipos + 2], ipdtmpl[ipos + 1])))
            return ret;
        sprintf(level_desc, "%s %s", tmpval1, "pv surface");
        /* write(tmpval1, *) ipdtmpl(ipos + 2). */
        /* call frmt(tmpval1, ipdtmpl(ipos + 2), ipdtmpl(ipos + 1)-6) */
        /*     level_desc = trim(tmpval1)//" pv surface" */
    }
    else if (ipdtmpl[ipos] == 111)
        strcpy(level_desc, " Eta level");
    else if (ipdtmpl[ipos] == 114)
        strcpy(level_desc, " Layer bet. 2-isent.");
    else if (ipdtmpl[ipos] == 117)
        strcpy(level_desc, " Mixed layer depth");
    else if (ipdtmpl[ipos] == 120)
        strcpy(level_desc, " Layer bet. 2-Eta lvl");
    else if (ipdtmpl[ipos] == 121)
        strcpy(level_desc, " Layer bet. 2-isob.");
    else if (ipdtmpl[ipos] == 125)
        strcpy(level_desc, " Specified height lvl");
    else if (ipdtmpl[ipos] == 126)
        strcpy(level_desc, " Isobaric level");
    else if (ipdtmpl[ipos] == 160)
        strcpy(level_desc, " Depth below sea lvl");
    else if (ipdtmpl[ipos] == 170)
        strcpy(level_desc, " Ionospheric D-region lvl");
    else if (ipdtmpl[ipos] == 1)
        strcpy(level_desc, " Surface");
    else if (ipdtmpl[ipos] == 2)
        strcpy(level_desc, " Cloud base lvl");
    else if (ipdtmpl[ipos] == 3)
        strcpy(level_desc, " Cloud top lvl");
    else if (ipdtmpl[ipos] == 4)
        strcpy(level_desc, " 0 Deg Isotherm");
    else if (ipdtmpl[ipos] == 5) /* from the surface. */
        strcpy(level_desc, " Level of adiabatic");
    else if (ipdtmpl[ipos] == 6)
        strcpy(level_desc, " Max wind lvl");
    else if (ipdtmpl[ipos] == 7)
        strcpy(level_desc, " Tropopause");
    else if (ipdtmpl[ipos] == 8)
        strcpy(level_desc, " Nom. top");
    else if (ipdtmpl[ipos] == 9)
        strcpy(level_desc, " Sea Bottom");
    else if (ipdtmpl[ipos] == 10)
        strcpy(level_desc, " Entire Atmosphere");
    else if (ipdtmpl[ipos] == 11)
        strcpy(level_desc, " Cumulonimbus Base");
    else if (ipdtmpl[ipos] == 12)
        strcpy(level_desc, " Cumulonimbus Top");
    else if (ipdtmpl[ipos] == 20)
        strcpy(level_desc, " Isothermal level");
    else if (ipdtmpl[ipos] == 200)
        strcpy(level_desc, " Entire Atmosphere");
    else if (ipdtmpl[ipos] == 201)
        strcpy(level_desc, " Entire ocean");
    else if (ipdtmpl[ipos] == 204)
        strcpy(level_desc, " Highest Frz. lvl");
    else if (ipdtmpl[ipos] == 206)
        strcpy(level_desc, " Grid scale cloud bl");
    else if (ipdtmpl[ipos] == 207)
        strcpy(level_desc, " Grid scale cloud tl");
    else if (ipdtmpl[ipos] == 209)
        strcpy(level_desc, " Boundary layer cbl");
    else if (ipdtmpl[ipos] == 210)
        strcpy(level_desc, " Boundary layer ctl");
    else if (ipdtmpl[ipos] == 211)
        strcpy(level_desc, " Boundary layer cl");
    else if (ipdtmpl[ipos] == 212)
        strcpy(level_desc, " Low cloud bot. lvl");
    else if (ipdtmpl[ipos] == 213)
        strcpy(level_desc, " Low cloud top lvl");
    else if (ipdtmpl[ipos] == 214)
        strcpy(level_desc, " Low cloud layer");
    else if (ipdtmpl[ipos] == 215)
        strcpy(level_desc, " Cloud ceiling");
    else if (ipdtmpl[ipos] == 220)
        strcpy(level_desc, " Planetary boundary");
    else if (ipdtmpl[ipos] == 221)
        strcpy(level_desc, " Layer 2 Hybrid lvl");
    else if (ipdtmpl[ipos] == 222)
        strcpy(level_desc, " Mid. cloud bot. lvl");
    else if (ipdtmpl[ipos] == 223)
        strcpy(level_desc, " Mid. cloud top lvl");
    else if (ipdtmpl[ipos] == 224)
        strcpy(level_desc, " Middle cloud layer");
    else if (ipdtmpl[ipos] == 232)
        strcpy(level_desc, " High cloud bot. lvl");
    else if (ipdtmpl[ipos] == 233)
        strcpy(level_desc, " High cloud top lvl");
    else if (ipdtmpl[ipos] == 234)
        strcpy(level_desc, " High cloud layer");
    else if (ipdtmpl[ipos] == 235)
        strcpy(level_desc, " Ocean Isotherm lvl");
    else if (ipdtmpl[ipos] == 236)
        strcpy(level_desc, " Layer 2-depth below");
    else if (ipdtmpl[ipos] == 237)
        strcpy(level_desc, " Bot. Ocean mix. lyr");
    else if (ipdtmpl[ipos] == 238)
        strcpy(level_desc, " Bot. Ocean iso. lyr");
    else if (ipdtmpl[ipos] == 239) /* isothermal level (S26CY). */
        strcpy(level_desc, " layer ocean sfc 26C");
    else if (ipdtmpl[ipos] == 240)
        strcpy(level_desc, " Ocean Mixed Layer");
    else if (ipdtmpl[ipos] == 241)
        strcpy(level_desc, " Order Seq. Of Data");
    else if (ipdtmpl[ipos] == 242)
        strcpy(level_desc, " Con. cloud bot. lvl");
    else if (ipdtmpl[ipos] == 243)
        strcpy(level_desc, " Con. cloud top lvl");
    else if (ipdtmpl[ipos] == 244)
        strcpy(level_desc, " Conv. cloud layer");
    else if (ipdtmpl[ipos] == 245)
        strcpy(level_desc, " Lowest lvl wet bulb");
    else if (ipdtmpl[ipos] == 246)
        strcpy(level_desc, " Max. equi. potential");
    else if (ipdtmpl[ipos] == 247)
        strcpy(level_desc, " Equilibrium level");
    else if (ipdtmpl[ipos] == 248)
        strcpy(level_desc, " Shallow con. cld bl");
    else if (ipdtmpl[ipos] == 249)
        strcpy(level_desc, " Shallow con. cld tl");
    else if (ipdtmpl[ipos] == 251)
        strcpy(level_desc, " Deep conv. cld bl");
    else if (ipdtmpl[ipos] == 252)
        strcpy(level_desc, " Deep conv. cld tl");
    else if (ipdtmpl[ipos] == 253) /* liquid water layer (LBLSW). */
        strcpy(level_desc, " Lowest bot. lvl sup");
    else if (ipdtmpl[ipos] == 254) /* liquid water layer (HBLSW). */
        strcpy(level_desc, " highest top lvl sup");
    else
    {
        sprintf(level_desc, " %4d (Unknown Lvl)", (int)ipdtmpl[ipos]);
        /* write(level_desc, fmt = '(1x,I4," (Unknown Lvl)")') ipdtmpl[ipos] */
    }

    return G2C_NOERROR;
}

/**
 * Write a summary file like the degrib2 utility.
 *
 * @param g2cid Indentifier for the file, returned by g2c_open() or
 * g2c_create().
 * @param fileout Path of output file. Any existing file of this name
 * will be overwritten.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid parameters.
 * - ::G2C_EFILE File I/O error.
 *
 * @author Ed Hartnett @date Sep 17, 2022
 */
int
g2c_degrib2(int g2cid, const char *fileout)
{
    FILE *f;
    G2C_MESSAGE_INFO_T *msg;
    int total_fields = 0;
    int i;
    int ret;

    /* Check inputs. */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES || g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;
    if (!fileout)
        return G2C_EINVAL;

    LOG((2, "g2c_degrib2 fileout %s", fileout));

    /* Open output file. */
    if (!(f = fopen(fileout, "w")))
        return G2C_EFILE;

    /* Print a summary for each message in the file. */
    for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
    {
        int fld;

        fprintf(f, "\n");
        fprintf(f, " GRIB MESSAGE  %ld  starts at %zu\n", msg->msg_num + 1, msg->bytes_to_msg + 1);
        fprintf(f, "\n");
        fprintf(f, "  SECTION 0:  %d 2 %zu\n", msg->discipline, msg->bytes_in_msg);
        fprintf(f, "  SECTION 1:  %d %d %d %d %d %d %d %d %d %d %d %d %d\n", msg->center, msg->subcenter,
                msg->master_version, msg->local_version, msg->sig_ref_time, msg->year, msg->month, msg->day,
                msg->hour, msg->minute, msg->second, msg->status, msg->type);
        fprintf(f, "  Contains  %d  Local Sections  and  %d  data fields.\n", msg->num_local, msg->num_fields);
        LOG((5, "  SECTION 1:  %d %d %d %d %d %d %d %d %d %d %d %d %d\n", msg->center, msg->subcenter,
             msg->master_version, msg->local_version, msg->sig_ref_time, msg->year, msg->month, msg->day,
             msg->hour, msg->minute, msg->second, msg->status, msg->type));

        /* For each field, print info. */
        for (fld = 0; fld < msg->num_fields; fld++)
        {
            G2C_SECTION_INFO_T *sec, *sec2, *sec3, *sec5, *sec6;
            G2C_SECTION3_INFO_T *sec3_info;
            G2C_SECTION4_INFO_T *sec4_info;
            G2C_SECTION5_INFO_T *sec5_info;
            G2C_SECTION6_INFO_T *sec6_info;
            char abbrev[G2C_MAX_NOAA_ABBREV_LEN + 1];
            char level_desc[G2C_MAX_TYPE_OF_FIXED_SURFACE_LEN + 1];
            char date_time[100 + 1];
            float *data = NULL;
            float total = 0.0, max = 0.0, min = 0.0, avg = 0;
            int t;

            fprintf(f, "\n");
            fprintf(f, "  FIELD  %d\n", fld + 1);
            /* Only print section 0 and 1 data for the first field. */
            if (fld == 0)
            {
                fprintf(f, "  SECTION 0:  %d 2\n", msg->discipline);
                fprintf(f, "  SECTION 1:  %d %d %d %d %d %d %d %d %d %d %d %d %d\n", msg->center, msg->subcenter,
                        msg->master_version, msg->local_version, msg->sig_ref_time, msg->year, msg->month, msg->day,
                        msg->hour, msg->minute, msg->second, msg->status, msg->type);
            }

            /* Find this field (a.k.a. product, a.k.a. section 4). */
            for (sec = msg->sec; sec; sec = sec->next)
                if (sec->sec_num == 4 && ((G2C_SECTION4_INFO_T *)(sec->sec_info))->field_num == fld)
                    break;
            if (!sec)
                return G2C_ENOSECTION;

            /* Is there a local section? */
            for (sec2 = sec; sec2; sec2 = sec2->prev)
                if (sec2->sec_num == 2)
                    break;
            if (sec2)
                fprintf(f, "  SECTION 2: %d bytes\n", sec2->sec_len - 5); /* Subtract 4 byte seclen and 1 byte secnum. */

            /* Find the sec3 that applies to this field. */
            for (sec3 = sec; sec3; sec3 = sec3->prev)
                if (sec3->sec_num == 3)
                    break;
            if (!sec3)
                return G2C_ENOSECTION;

            /* Print the section 3 (grid) info. */
            sec3_info = (G2C_SECTION3_INFO_T *)sec3->sec_info;
            fprintf(f, "  SECTION 3:  %d %d %d %d %d\n", sec3_info->source_grid_def, sec3_info->num_data_points, sec3_info->num_opt,
                    sec3_info->interp_list, sec3_info->grid_def);
            fprintf(f, "  GRID TEMPLATE 3. %d : ", sec3_info->grid_def);
            /* Cast to int to match behavior of  degrib2.F90. */
            for (t = 0; t < sec3->template_len; t++)
                fprintf(f, " %d", (int)sec3->template[t]);
            fprintf(f, "\n");
            if (!sec3_info->optional)
                fprintf(f, "  NO Optional List Defining Number of Data Points.\n");

            /* Print the section 4 (product) info. */
            sec4_info = (G2C_SECTION4_INFO_T *)sec->sec_info;
            fprintf(f, "  PRODUCT TEMPLATE 4. %d: ", sec4_info->prod_def);

            /* We can look up the parameter abbreviation with the
             * discipline and the first two elements of the product
             * template, which will be the category and product
             * number. */
            if ((ret = g2c_param_abbrev(msg->discipline, sec->template[0], sec->template[1], abbrev)))
                return ret;
            fprintf(f, "( PARAMETER = %-8s %d %lld %lld ) ", abbrev, sec->msg->discipline, sec->template[0],
                    sec->template[1]);
            /* Cast to int to match behavior of  degrib2.F90. */
            for (t = 0; t < sec->template_len; t++)
                fprintf(f, " %d", (int)sec->template[t]);
            fprintf(f, "\n");

            /* Using the product template number, and the template
             * values, we can figure out a description for the
             * horizontal level description. */
            if ((ret = g2c_get_level_desc(sec4_info->prod_def, sec->template, level_desc)))
                return ret;

            /* The Fortran code for level descriptions puts a space in
             * front of level descriptions that don't start with a
             * number. It's annoying, but do it here to match the
             * behavior of Fortran degrib2. */
            /* if (!isdigit(level_desc[0])) */
            /* { */
            /*     char tmp_level_desc[G2C_MAX_TYPE_OF_FIXED_SURFACE_LEN + 2]; */
            /*     sprintf(tmp_level_desc, " %s", level_desc); */
            /*     strcpy(level_desc, tmp_level_desc); */
            /* } */

            /* Put the date/time in a formatted string. */
            if ((ret = g2c_get_datetime(sec4_info->prod_def, sec->template, msg->year, msg->month, msg->day,
                                        msg->hour, msg->minute, msg->second, date_time)))
                return ret;
            fprintf(f, "  FIELD: %-8s%s %s\n", abbrev, level_desc, date_time);
            if (!sec4_info->num_coord)
                fprintf(f, "  NO Optional Vertical Coordinate List.\n");

            /* Find the sec5 that applies to this field. */
            for (sec5 = sec; sec5; sec5 = sec5->next)
                if (sec5->sec_num == 5)
                    break;
            if (!sec5)
                return G2C_ENOSECTION;

            /* Find the sec6 that applies to this field, if any. */
            for (sec6 = sec; sec6; sec6 = sec6->next)
                if (sec6->sec_num == 6)
                    break;

            /* Section 5 and 6 info. */
            sec5_info = (G2C_SECTION5_INFO_T *)sec5->sec_info;
            if (sec6)
                sec6_info = (G2C_SECTION6_INFO_T *)sec6->sec_info;
            if (sec6 && sec6_info->indicator != 255)
                fprintf(f, "  Num. of Data Points =  %d    with BIT-MAP  0\n", sec5_info->num_data_points);
            else
                fprintf(f, "  Num. of Data Points =  %d     NO BIT-MAP \n", sec5_info->num_data_points);
            fprintf(f, "  DRS TEMPLATE 5. %d : ", sec5_info->data_def);
            /* Cast to int to match behavior of  degrib2.F90. */
            for (t = 0; t < sec5->template_len; t++)
                fprintf(f, " %d", (int)sec5->template[t]);
            fprintf(f, "\n");
            fprintf(f, "  Data Values:\n");
            fprintf(f, "  Num. of Data Points =  %d   Num. of Data Undefined = 0\n", sec5_info->num_data_points);

            /* Now read the data and find the min, max, and average values. */

            LOG((5, "sec5_info->num_data_points %ld", sec5_info->num_data_points));
            if (sec5_info->num_data_points)
            {
                /* Allocate storage for the data. */
                if (!(data = malloc(sec5_info->num_data_points * sizeof(float))))
                    return G2C_ERROR;
            }

            /* Get the data from message 0, product 0. */
            if ((ret = g2c_get_prod(g2cid, msg->msg_num, fld, NULL, data)))
                return ret;

            /* Find min/max/avg. */
            if (sec5_info->num_data_points)
            {
                max = data[0];
                min = data[0];
                total = data[0];
                for (i = 1; i < sec5_info->num_data_points; i++)
                {
                    total += data[i];
                    if (data[i] > max)
                        max = data[i];
                    if (data[i] < min)
                        min = data[i];
                }
                avg = total / sec5_info->num_data_points;
            }
            fprintf(f, "( PARM= %s ) :  MIN=%25.8f AVE=%25.8f MAX=%25.8f\n",
                    abbrev, min, avg, max);

            /* Free the data. */
            if (sec5_info->num_data_points)
                free(data);

            total_fields++;
        }
    }

    fprintf(f, "  \n  Total Number of Fields Found =  %d\n", total_fields);

    /* Close output file. */
    if (fclose(f))
        return G2C_EFILE;

    return G2C_NOERROR;
}
