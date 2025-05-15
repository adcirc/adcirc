/**
 * @file
 * @brief Internal utility functions for the NCEPLIBS-g2c library for
 * the new g2c API.
 *
 * @author Ed Hartnett @date Oct 11, 2023
 */

#include "grib2_int.h"
#include <stdarg.h>

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** Pointer to the list of code tables. */
extern G2C_CODE_TABLE_T *g2c_table;

/**
 * Given an error code, return an error message.
 *
 * @param g2cerr An error number returned by one of the g2c_*
 * functions.
 *
 * @return The error message.
 *
 * @author Ed Hartnett
 */
const char *
g2c_strerror(int g2cerr)
{
    switch (g2cerr)
    {
    case G2C_NOERROR:
        return "No error";
    case G2C_ENOTGRIB:
        return "GRIB2 header not found";
    case G2C_EMSGCOMPLETE:
        return "GRIB message is already complete.";
    case G2C_ENAMETOOLONG:
        return "Name too long";
    case G2C_EINVAL:
        return "Invalid input";
    case G2C_EFILE:
        return "Error reading file";
    case G2C_EBADID:
        return "Bad ID";
    case G2C_ETOOMANYFILES:
        return "Too many files open";
    case G2C_ENOMEM:
        return "Out of memory";
    case G2C_EMSG:
        return "Error decoding message";
    case G2C_ENOMSG:
        return "No GRIB message found";
    case G2C_EXML:
        return "Error parsing XML";
    case G2C_ENOTFOUND:
        return "Table or entry not found";
    case G2C_ENOTGRIB2:
        return "Not GRIB 2";
    case G2C_ENOSECTION:
        return "Cannot find section";
    case G2C_ENOEND:
        return "Cannot find end of GRIB message";
    case G2C_EBADEND:
        return "End of message in wrong place";
    case G2C_EBADSECTION:
        return "Invalid section number";
    case G2C_EJPEG:
        return "Error encoding/decoding JPEG data";
    case G2C_EPNG:
        return "Error encoding/decoding PNG data";
    case G2C_ENOTEMPLATE:
        return "Template not found";
    case G2C_EBADTEMPLATE:
        return "Template problem";
    case G2C_ENOPARAM:
        return "Parameter not found";
    case G2C_ENOPRODUCT:
        return "Product not found";
    case G2C_EBADTYPE:
        return "Bad type";

    default:
        return "Unknown Error";
    }
}

/**
 * Log section 0 information.
 *
 * @param msg Pointer to ::G2C_MESSAGE_INFO_T which contains section 0
 * information.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 10/16/22
 */
int
g2c_log_section1(G2C_MESSAGE_INFO_T *msg)
{
#ifdef LOGGING
    char desc[G2C_MAX_GRIB_DESC_LEN + 1];
    int ret;

    /* Read in the CSV GRIB2 code definitions. */
    if ((ret = g2c_csv_init()))
        return ret;

    /* Section 0 discipline flag. */
    if ((ret = g2c_find_desc("Code table 0.0", msg->discipline, desc)))
        return ret;
    LOG((4, "Discipline: %s", desc));

    /* Section 1 flags. */
    LOG((4, "Identification of originating/generating center: %d", msg->center));
    LOG((4, "Identification of originating/generating subcenter: %d", msg->subcenter));
    if ((ret = g2c_find_desc("Code table 1.0", msg->master_version, desc)))
        return ret;
    LOG((4, "GRIB master tables version number: %s", desc));
    if ((ret = g2c_find_desc("Code table 1.1", msg->local_version, desc)))
        return ret;
    LOG((4, "Version number of GRIB local tables used to augment Master Tables: %s", desc));
    if ((ret = g2c_find_desc("Code table 1.2", msg->sig_ref_time, desc)))
        return ret;
    LOG((4, "Significance of reference time: %s", desc));
    LOG((4, "Reference time: %d/%d/%d %d:%d:%d", msg->year, msg->month, msg->day,
         msg->hour, msg->minute, msg->second));
    if ((ret = g2c_find_desc("Code table 1.3", msg->status, desc)))
        return ret;
    LOG((4, "Production Status of Processed data in the GRIB message: %s", desc));
    if ((ret = g2c_find_desc("Code table 1.4", msg->type, desc)))
        return ret;
    LOG((4, "Type of processed data in this GRIB message: %s", desc));

#endif
    return G2C_NOERROR;
}

/**
 * Log info about a section.
 *
 * @param sec Pointer to ::G2C_SECTION_INFO_T.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett 10/20/22
 */
int
g2c_log_section(G2C_SECTION_INFO_T *sec)
{
#ifdef LOGGING
    int ret;

    LOG((3, "sec_id %d sec_len %d byte_to_sec %ld sec_num %d", sec->sec_id,
         sec->sec_len, sec->bytes_to_sec, sec->sec_num));
    if (sec->sec_num == 4)
    {
        char abbrev[G2C_MAX_NOAA_ABBREV_LEN + 1];

        /* Look up the parameter abbreviation with the discipline,
         * category, and product number. */
        if ((ret = g2c_param_abbrev(sec->msg->discipline, sec->template[0],
                                    sec->template[1], abbrev)))
            return ret;
        LOG((4, "%s", abbrev));
    }

#endif
    return G2C_NOERROR;
}

/**
 * Print a summary of the contents of an open GRIB2 file. If the
 * NCEPLIBS-g2c library is built without the LOGGING option, this
 * function will do nothing.
 *
 * @param g2cid The file ID, returned by g2c_open() or g2c_create().
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID g2cid not found.
 *
 * @author Ed Hartnett 8/22/22
 */
int
g2c_log_file(int g2cid)
{
#ifdef LOGGING
    G2C_MESSAGE_INFO_T *msg;
    int ret;

    /* Read in the CSV GRIB2 code definitions. */
    if ((ret = g2c_csv_init()))
        return ret;

    /* Find the open file struct. */
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    LOG((1, "path %s", g2c_file[g2cid].path));
    LOG((1, "num_messages %ld", g2c_file[g2cid].num_messages));
    for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
    {
        G2C_SECTION_INFO_T *sec;

        LOG((1, "message %ld bytes_to_msg %ld bytes_in_msg %ld num_fields %d num_local %d",
             msg->msg_num, msg->bytes_to_msg, msg->bytes_in_msg, msg->num_fields, msg->num_local));
        LOG((1, "sec1_len %d center %d subcenter %d master_version %d local_version %d",
             msg->sec1_len, msg->center, msg->subcenter, msg->master_version, msg->local_version));
        LOG((1, "sig_ref_time %d %d %d %d %d:%d:%d status %d type %d", msg->sig_ref_time, msg->year,
             msg->month, msg->day, msg->hour, msg->minute, msg->second, msg->status, msg->type));

        /* If we've loaded CSV tables, decode some flags. */
        if (g2c_table)
            if ((ret = g2c_log_section1(msg)))
                return ret;

        /* Section info. */
        for (sec = msg->sec; sec; sec = sec->next)
            if ((ret = g2c_log_section(sec)))
                return ret;
    }

    /* Free CSV code memory. */
    g2c_free_tables();

#endif /* LOGGING */
    return G2C_NOERROR;
}
