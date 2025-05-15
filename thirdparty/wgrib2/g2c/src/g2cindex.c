/**
 * @file
 * @brief Write summary output to an index file, as is done by utility
 * grb2index.
 *
 * @author Ed Hartnett @date 10/12/22
 */
#include "grib2_int.h"
#include <libgen.h>
#include <stdarg.h>
#include <time.h>

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** Length of the two header lines at the top of the index file. */
#define G2C_INDEX_HEADER_LEN 81

/** Length of the basename in header record 2. */
#define G2C_INDEX_BASENAME_LEN 40

/** Length of bitmap section included in the index record. */
#define G2C_INDEX_BITMAP_BYTES 6

/** Length of beginning of index record. */
#define G2C_INDEX_FIXED_LEN 44

/** Length of beginning of index record for large files. */
#define G2C_INDEX_FIXED_LEN_2 48

/** Length of date string in index record. */
#define G2C_INDEX_DATE_STR_LEN 10

/** Length of time string in index record. */
#define G2C_INDEX_TIME_STR_LEN 8

/** Length of str1 string in index record. */
#define G2C_INDEX_STR1_LEN 7

/** Use externally-defined mutex for thread-safety. */
EXTERN_MUTEX(m);

/**
 * Read or write the start of a version 2 index record.
 *
 * @param f FILE * to open index file.
 * @param rw_flag True if function should write, false if it should read.
 * @param reclen Pointer to reclen.
 * @param msg Pointer to msg.
 * @param local Pointer to local.
 * @param gds Pointer to gds.
 * @param pds Pointer to pds.
 * @param drs Pointer to drs.
 * @param bms Pointer to bms.
 * @param data Pointer to data.
 * @param msglen Pointer to msglen.
 * @param version Pointer to version.
 * @param discipline Pointer to discipline.
 * @param fieldnum Pointer to fieldnum, 0- based. (It is 1-based in
 * the index file.)
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_EFILE File I/O error.
 *
 * @author Ed Hartnett 10/26/22
 */
int
g2c_start_index_record(FILE *f, int rw_flag, int *reclen, int *msg, int *local, int *gds,
                       int *pds, int *drs, int *bms, int *data, size_t *msglen,
                       unsigned char *version, unsigned char *discipline, short *fieldnum)
{
    /* size_t size_t_be; */
    short fieldnum1; /* This is for the 1-based fieldnum in the index file. */
    int ret;

    /* All pointers must be provided. */
    if (!f || !reclen || !msg || !local || !gds || !pds || !drs || !bms || !data || !msglen || !version || !discipline || !fieldnum)
        return G2C_EINVAL;

    LOG((4, "g2c_start_index_record rw_flag %d", rw_flag));

    /* When writing, set the fieldnum1 to be a 1-based index, just
     * like in Fortran. */
    if (rw_flag)
        fieldnum1 = *fieldnum + 1;

    /* Read or write the values at the beginning of each index
     * record. */
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)reclen)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)msg)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)local)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)gds)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)pds)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)drs)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)bms)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)data)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)msglen)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, version)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, discipline)))
        return ret;
    if ((ret = g2c_file_io_short(f, rw_flag, &fieldnum1)))
        return ret;

    /* When reading, translate the 1-based fieldnum1 into the 0-based
     * fieldnum that C programmers will expect and love. */
    if (!rw_flag)
        *fieldnum = fieldnum1 - 1;

    return G2C_NOERROR;
}

/**
 * Read or write the start of a version 2 index record.
 *
 * @param f FILE pointer to open index file.
 * @param rw_flag True if function should write, false if it should read.
 * @param reclen Pointer to reclen, the length of the index record in bytes.
 * @param msg Pointer to bytes to skip in file to reach msg.
 * @param local Pointer to bytes to skip in message to reach local.
 * @param gds Pointer to bytes to skip in message to reach gds.
 * @param pds Pointer to bytes to skip in message to reach pds.
 * @param drs Pointer to bytes to skip in message to reach drs.
 * @param bms Pointer to bytes to skip in message to reach bms.
 * @param data Pointer to bytes to skip in message to reach data.
 * @param msglen Pointer to msglen.
 * @param version Pointer to version.
 * @param discipline Pointer to discipline.
 * @param fieldnum Pointer to fieldnum, 0- based. (It is 1-based in
 * the index file.)
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_EFILE File I/O error.
 *
 * @author Ed Hartnett 10/26/22
 */
int
g2c_start_index_record_lf(FILE *f, int rw_flag, int *reclen, size_t *msg, size_t *local, size_t *gds,
                          size_t *pds, size_t *drs, size_t *bms, size_t *data, size_t *msglen,
                          unsigned char *version, unsigned char *discipline, short *fieldnum)
{
    /* size_t size_t_be; */
    short fieldnum1; /* This is for the 1-based fieldnum in the index file. */
    int ret;

    LOG((4, "g2c_start_index_record_lf rw_flag %d", rw_flag));

    /* All pointers must be provided. */
    if (!f || !reclen || !msg || !local || !gds || !pds || !drs || !bms || !data || !msglen || !version || !discipline || !fieldnum)
        return G2C_EINVAL;

    /* When writing, set the fieldnum1 to be a 1-based index, just
     * like in Fortran. */
    if (rw_flag)
        fieldnum1 = *fieldnum + 1;

    /* Read or write the values at the beginning of each index
     * record. */
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)reclen)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)msg)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)local)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)gds)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)pds)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)drs)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)bms)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)data)))
        return ret;
    if ((ret = g2c_file_io_ulonglong(f, rw_flag, (unsigned long long *)msglen)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, version)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, discipline)))
        return ret;
    if ((ret = g2c_file_io_short(f, rw_flag, &fieldnum1)))
        return ret;

    /* When reading, translate the 1-based fieldnum1 into the 0-based
     * fieldnum that C programmers will expect and love. */
    if (!rw_flag)
        *fieldnum = fieldnum1 - 1;

    return G2C_NOERROR;
}

/**
 * Read or write the start of a version 1 index record.
 *
 * For more detail on version 1 of the index format, see the
 * [grbindex](https://noaa-emc.github.io/NCEPLIBS-grib_util/grbindex/grbindex_8f.html)
 * documentation in the
 * [NCEPLIBS-grib_util](https://github.com/NOAA-EMC/NCEPLIBS-grib_util).
 *
 * @param f FILE * to open index file.
 * @param rw_flag True if function should write, false if it should read.
 * @param b2_msg Pointer that gets the bytes to skip in file before msg.
 * @param b2_pds Pointer that gets bytes to skip in message before pds.
 * @param b2_gds Pointer that gets bytes to skip in message before gds (0 if no gds).
 * @param b2_bms Pointer that gets bytes to skip in message before bms (0 if no bms).
 * @param b2_bds Pointer that gets bytes to skip in message before bds.
 * @param msglen Pointer that gets bytes total in the message.
 * @param version Pointer that gets grib version number (always 1 for this function).
 * @param pds_val Pointer that gets an arry of 27 bytes of the product definition section (pds).
 * @param gds_val Pointer that gets an arry of 41 bytes of the gds.
 * @param bms_val Pointer that gets an arry of 5 bytes of the bms.
 * @param bds_val Pointer that gets an arry of 10 bytes, bytes 41-100 of the bds.
 * @param pds_val2 Pointer that gets an arry of 59 bytes 41-100 of the pds. Ignored if null.
 * @param pds_val3 Pointer that gets an arry of 11 bytes 29-40 of the pds. Ignored if null.
 * @param gds_val2 Pointer that gets an arry of 135 bytes 43-178 of the gds. Ignored if null.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_EFILE File I/O error.
 *
 * @author Ed Hartnett 9/11/23
 */
int
g2c_start_index1_record(FILE *f, int rw_flag, unsigned int *b2_msg, unsigned int *b2_pds,
                        unsigned int *b2_gds, unsigned int *b2_bms, unsigned int *b2_bds,
                        unsigned int *msglen, unsigned char *version, unsigned char *pds_val,
                        unsigned char *gds_val, unsigned char *bms_val, unsigned char *bds_val,
                        unsigned char *pds_val2, unsigned char *pds_val3, unsigned char *gds_val2)
{
    size_t bytes_read;
    int ret;

    /* All pointers must be provided. */
    if (!f || !b2_msg || !b2_pds || !b2_gds || !b2_bms || !b2_bds ||
        !msglen || !version)
        return G2C_EINVAL;

    /* Read or write the values at the beginning of each index
     * record. */
    if ((ret = g2c_file_io_uint(f, rw_flag, b2_msg)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, b2_pds)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, b2_gds)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, b2_bms)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, b2_bds)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, msglen)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, version)))
        return ret;

    /* The index record contains some metadata copied directly from
     * the file. These are arrays of unsigned char, of known
     * length. For more detail see
     * https://noaa-emc.github.io/NCEPLIBS-grib_util/grbindex/grbindex_8f.html. */
    if ((bytes_read = fread(pds_val, 1, G2C_INDEX1_PDS_VAL_LEN, f)) != G2C_INDEX1_PDS_VAL_LEN)
        return G2C_EFILE;
    if ((bytes_read = fread(gds_val, 1, G2C_INDEX1_GDS_VAL_LEN, f)) != G2C_INDEX1_GDS_VAL_LEN)
        return G2C_EFILE;
    if ((bytes_read = fread(bms_val, 1, G2C_INDEX1_BMS_VAL_LEN, f)) != G2C_INDEX1_BMS_VAL_LEN)
        return G2C_EFILE;
    if ((bytes_read = fread(bds_val, 1, G2C_INDEX1_BDS_VAL_LEN, f)) != G2C_INDEX1_BDS_VAL_LEN)
        return G2C_EFILE;

    return G2C_NOERROR;
}

/**
 * Given a pointer to a message, and a field number, return pointers
 * to all relevent section structs for that product.
 *
 * Each product is defined in a section 4, and has an associated
 * section 3, 5, 6, and 7.
 *
 * @param msg Pointer to a G2C_MESSAGE_INFO_T with information about
 * the message.
 * @param fieldnum The field number (first field in message is 0).
 * @param sec3 Pointer that gets a pointer to the G2C_SECTION_INFO_T
 * struct for the section 3 associated with this product.
 * @param sec4 Pointer that gets a pointer to the G2C_SECTION_INFO_T
 * struct for the section 4 associated with this product.
 * @param sec5 Pointer that gets a pointer to the G2C_SECTION_INFO_T
 * struct for the section 5 associated with this product.
 * @param sec6 Pointer that gets a pointer to the G2C_SECTION_INFO_T
 * struct for the section 6 associated with this product. NULL is
 * returned if there is no section 6.
 * @param sec7 Pointer that gets a pointer to the G2C_SECTION_INFO_T
 * struct for the section 7 associated with this product.
 *
 * @note This is an internal function and should not be called by users.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOSECTION Section not found.
 *
 * @author Ed Hartnett @date 10/27/22
 */
int
g2c_get_prod_sections(G2C_MESSAGE_INFO_T *msg, int fieldnum, G2C_SECTION_INFO_T **sec3,
                      G2C_SECTION_INFO_T **sec4, G2C_SECTION_INFO_T **sec5,
                      G2C_SECTION_INFO_T **sec6, G2C_SECTION_INFO_T **sec7)
{
    G2C_SECTION_INFO_T *s3, *s4, *s5, *s6, *s7;

    /* Check inputs. */
    if (!msg || fieldnum < 0 || !sec3 || !sec4 || !sec5 || !sec6 || !sec7)
        return G2C_EINVAL;

    /* Find the product with matching fieldnum. */
    for (s4 = msg->sec; s4; s4 = s4->next)
    {
        G2C_SECTION4_INFO_T *s4info = s4->sec_info;
        if (s4->sec_num != 4)
            continue;
        if (s4info->field_num == fieldnum)
            break;
    }
    if (!s4)
        return G2C_ENOSECTION;

    /* Find the section 3, grid definition section, which is
     * associated with this product. */
    for (s3 = s4->prev; s3; s3 = s3->prev)
        if (s3->sec_num == 3)
            break;
    if (!s3)
        return G2C_ENOSECTION;

    /* Find the section 5, data representation section, which
     * is associated with this product. */
    for (s5 = s4->next; s5; s5 = s5->next)
        if (s5->sec_num == 5)
            break;
    if (!s5)
        return G2C_ENOSECTION;

    /* Find the section 6, the bit map section. There may not be
     * one. */
    for (s6 = s5->next; s6; s6 = s6->next)
    {
        if (s6->sec_num == 6)
            break;

        /* If we hit section 7, there's no bitmap. */
        if (s6->sec_num == 7)
        {
            s6 = NULL;
            break;
        }
    }

    /* Find the section 7, data section, which is associated with this
     * product. */
    for (s7 = s5->next; s7; s7 = s7->next)
        if (s7->sec_num == 7)
            break;
    if (!s7)
        return G2C_ENOSECTION;

    /* Return results to caller. */
    *sec3 = s3;
    *sec4 = s4;
    *sec5 = s5;
    *sec6 = s6;
    *sec7 = s7;

    return G2C_NOERROR;
}

/**
 * Create an index file from a GRIB2 file, just like those created by
 * the grb2index utility.
 *
 * The index file starts with two header records:
 * 1. 81-byte header with 'gb2ix1' in columns 42-47.
 * 2. 81-byte header with number of bytes to skip before index records,
 * total length in bytes of the index records, number of index records,
 * and grib file basename written in format ('ix1form:',3i10,2x,a40).
 *
 * Each following index record corresponds to a grib message
 * and has the internal format:
 * - byte 001 - 004 length of index record
 * - byte 005 - 008 bytes to skip in data file before grib message
 * - byte 009 - 012 bytes to skip in message before lus (local use) (0, if none).
 * - byte 013 - 016 bytes to skip in message before gds
 * - byte 017 - 020 bytes to skip in message before pds
 * - byte 021 - 024 bytes to skip in message before drs
 * - byte 025 - 028 bytes to skip in message before bms
 * - byte 029 - 032 bytes to skip in message before data section
 * - byte 033 - 040 bytes total in the message
 * - byte 041 - 041 grib version number (currently 2)
 * - byte 042 - 042 message discipline
 * - byte 043 - 044 field number within grib2 message (1-based)
 * - byte 045 -  ii identification section (ids)
 * - byte ii+1-  jj grid definition section (gds)
 * - byte jj+1-  kk product definition section (pds)
 * - byte kk+1-  ll the data representation section (drs)
 * - byte ll+1-ll+6 first 6 bytes of the bit map section (bms)
 *
 * @param g2cid File it for an open GRIB2 file, as returned by
 * g2c_open().
 * @param mode Mode flags. Set ::G2C_NOCLOBBER to avoid overwriting
 * and existing file.
 * @param index_file The name that will be given to the index file. An
 * existing file will be overwritten.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 10/12/22
 */
int
g2c_write_index(int g2cid, int mode, const char *index_file)
{
    FILE *f;
    char h1[G2C_INDEX_HEADER_LEN * 2 + 1]; /* need extra space to silence GNU warnings */
    char h2[G2C_INDEX_HEADER_LEN + 10];
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    size_t items_written;
    char my_path[G2C_INDEX_BASENAME_LEN + 1];
    G2C_MESSAGE_INFO_T *msg;
    int total_index_size = 0; /* Does not include size of header records. */
    int index_version = 1;    /* 1 for legacy, 2 if indexed file may be > 2 GB. */
    int reclen;
    int ret = G2C_NOERROR;

    /* Is this an open GRIB2 file? */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (!index_file)
        return G2C_EINVAL;

    LOG((2, "g2c_write_index g2cid %d mode %d index_file %s", g2cid, mode,
         index_file));

    /* If NOCLOBBER, check if file exists. */
    if (mode & G2C_NOCLOBBER)
    {
        FILE *f;
        if ((f = fopen(index_file, "r")))
        {
            fclose(f);
            return G2C_EFILE;
        }
    }

    /* If LARGE_INDEX_FILE, set index version. */
    if (mode & G2C_LARGE_FILE_INDEX)
        index_version = 2;

    /* Create the index file. */
    if (!(f = fopen(index_file, "wb+")))
        return G2C_EFILE;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    if (!ret)
    {
        /* Create header 1. */
        snprintf(h1, G2C_INDEX_HEADER_LEN + 1,
                 "!GFHDR!  1   1   162 %4.4u-%2.2u-%2.2u %2.2u:%2.2u:%2.2u %s        hfe08           grb2index\n",
                 (tm.tm_year + 1900), (tm.tm_mon + 1), tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec,
                 "GB2IX1");

        /* Write header 1. */
        if ((items_written = fwrite(h1, G2C_INDEX_HEADER_LEN, 1, f)) != 1)
            ret = G2C_EFILE;
    }

    /* Find the total length of the index we are generating. */
    if (!ret)
    {
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
        {
            short fieldnum;

            /* Find information for each field in the message. */
            for (fieldnum = 0; fieldnum < msg->num_fields; fieldnum++)
            {
                G2C_SECTION_INFO_T *sec3, *sec4, *sec5, *sec6, *sec7;

                if ((ret = g2c_get_prod_sections(msg, fieldnum, &sec3, &sec4, &sec5, &sec6, &sec7)))
                    break;

                /* What will be the length of this index record? */
                reclen = (index_version == 1 ? G2C_INDEX_FIXED_LEN : G2C_INDEX_FIXED_LEN_2) +
                         msg->sec1_len + sec3->sec_len + sec4->sec_len + sec5->sec_len + G2C_INDEX_BITMAP_BYTES;
                total_index_size += reclen;
                LOG((4, "fieldnum %d reclen %d total_index_size %d", fieldnum, reclen, total_index_size));
            } /* next product */

            /* If there was a problem, give up. */
            if (ret)
                break;
        } /* next message */
    }

    /* Create header 2. */
    if (!ret)
    {
        strncpy(my_path, basename(g2c_file[g2cid].path), G2C_INDEX_BASENAME_LEN);
        sprintf(h2, "IX%dFORM:       162    %6d    %6ld  %s    \n", index_version, total_index_size,
                g2c_file[g2cid].num_messages, my_path);
        LOG((5, "header 2: %s", h2));

        /* Write header 2. */
        if ((items_written = fwrite(h2, G2C_INDEX_HEADER_LEN, 1, f)) != 1)
            ret = G2C_EFILE;
    }

    /* Write a record of index file for each message in the file. */
    if (!ret)
    {
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
        {
            short fieldnum;

            /* Find information for each field in the message. */
            for (fieldnum = 0; fieldnum < msg->num_fields; fieldnum++)
            {
                G2C_SECTION_INFO_T *sec3, *sec4, *sec5, *sec6, *sec7;
                int bs3, bs4, bs5, bs6, bs7;              /* bytes to each section, as ints. */
                size_t bs3_8, bs4_8, bs5_8, bs6_8, bs7_8; /* bytes to each section, as size_t. */
                unsigned char sec_num;
                int ret;

                if ((ret = g2c_get_prod_sections(msg, fieldnum, &sec3, &sec4, &sec5, &sec6, &sec7)))
                    return ret;
                bs3 = (int)sec3->bytes_to_sec;
                bs4 = (int)sec4->bytes_to_sec;
                bs5 = (int)sec5->bytes_to_sec;
                bs6 = (int)sec6->bytes_to_sec;
                bs7 = (int)sec7->bytes_to_sec;

                bs3_8 = sec3->bytes_to_sec;
                bs4_8 = sec4->bytes_to_sec;
                bs5_8 = sec5->bytes_to_sec;
                bs6_8 = sec6->bytes_to_sec;
                bs7_8 = sec7->bytes_to_sec;

                /* What will be the length of this index record? */
                reclen = (index_version == 1 ? G2C_INDEX_FIXED_LEN : G2C_INDEX_FIXED_LEN_2) + msg->sec1_len + sec3->sec_len + sec4->sec_len + sec5->sec_len + G2C_INDEX_BITMAP_BYTES;
                LOG((4, "fieldnum %d reclen %d", fieldnum, reclen));

                /* Write the beginning of the index record. */
                if (index_version == 2)
                {
                    if ((ret = g2c_start_index_record_lf(f, G2C_FILE_WRITE, &reclen, &msg->bytes_to_msg, &msg->bytes_to_local,
                                                         &bs3_8, &bs4_8, &bs5_8, &bs6_8, &bs7_8, &msg->bytes_in_msg, &msg->master_version,
                                                         &msg->discipline, &fieldnum)))
                        break;
                }
                else
                {
                    int bytes_to_msg = (int)msg->bytes_to_msg;
                    int b2l;

                    b2l = (int)msg->bytes_to_local;
                    if ((ret = g2c_start_index_record(f, G2C_FILE_WRITE, &reclen, &bytes_to_msg, &b2l,
                                                      &bs3, &bs4, &bs5, &bs6, &bs7, &msg->bytes_in_msg, &msg->master_version,
                                                      &msg->discipline, &fieldnum)))
                        break;
                }

                /* Write the section 1, identification section. */
                if ((ret = g2c_rw_section1_metadata(f, G2C_FILE_WRITE, msg)))
                    break;

                /* Write the section 3, grid definition section. */
                sec_num = 3;
                if ((ret = g2c_file_io_uint(f, G2C_FILE_WRITE, &sec3->sec_len)))
                    return ret;
                if ((ret = g2c_file_io_ubyte(f, G2C_FILE_WRITE, &sec_num)))
                    return ret;
                if ((ret = g2c_rw_section3_metadata(f, G2C_FILE_WRITE, sec3)))
                    break;

                /* Write the section 4, product definition section. */
                sec_num = 4;
                if ((ret = g2c_file_io_uint(f, G2C_FILE_WRITE, &sec4->sec_len)))
                    return ret;
                if ((ret = g2c_file_io_ubyte(f, G2C_FILE_WRITE, &sec_num)))
                    return ret;
                if ((ret = g2c_rw_section4_metadata(f, G2C_FILE_WRITE, sec4)))
                    break;

                /* Write the section 5, data representation section. */
                sec_num = 5;
                if ((ret = g2c_file_io_uint(f, G2C_FILE_WRITE, &sec5->sec_len)))
                    return ret;
                if ((ret = g2c_file_io_ubyte(f, G2C_FILE_WRITE, &sec_num)))
                    return ret;
                if ((ret = g2c_rw_section5_metadata(f, G2C_FILE_WRITE, sec5)))
                    break;

                /* Write the first 6 bytes of the bitmap section, if there
                 * is one. */
                if (!ret && sec6)
                {
                    unsigned char sample[G2C_INDEX_BITMAP_BYTES];
                    int b;

                    /* Read the first 6 bytes of the bitmap section. */
                    if (fseek(msg->file->f, msg->bytes_to_msg + sec6->bytes_to_sec, SEEK_SET))
                    {
                        ret = G2C_EFILE;
                        break;
                    }
                    if ((fread(sample, ONE_BYTE, G2C_INDEX_BITMAP_BYTES, msg->file->f)) != G2C_INDEX_BITMAP_BYTES)
                    {
                        ret = G2C_EFILE;
                        break;
                    }

                    /* Now write these bytes to the end of the index record. */
                    for (b = 0; b < G2C_INDEX_BITMAP_BYTES; b++)
                        if ((ret = g2c_file_io_ubyte(f, G2C_FILE_WRITE, &sample[b])))
                            return ret;
                }
            } /* next product */

            /* If there was a problem, give up. */
            if (ret)
                break;
        } /* next message */
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    /* Close the index file. */
    if (!ret)
        if (fclose(f))
            return G2C_EFILE;

    return ret;
}

/**
 * Read the header record apparently named after Steve Lord.
 *
 * This function reads the first of two 81-byte header records of an
 * index file.
 *
 * @param f Pointer to open FILE.
 * @param ip Pointer that gets i value. Ignored if NULL.
 * @param jp Pointer that gets j value. Ignored if NULL.
 * @param kp Pointer that gets k value. Ignored if NULL.
 * @param date_str Pointer to char array of size
 * ::G2C_INDEX_DATE_STR_LEN + 1 which will get the date string from the
 * header. Ignored if NULL.
 * @param time_str Pointer to char array of size
 * ::G2C_INDEX_TIME_STR_LEN + 1 which will get the time string from the
 * header. Ignored if NULL.
 *
 * @returns 0 for success, error code otherwise.
 *
 * @author Edward Hartnett @date 9/10/23
*/
static int
read_hdr_rec1(FILE *f, int *ip, int *jp, int *kp, char *date_str, char *time_str)
{
    size_t bytes_read;
    char line[G2C_INDEX_HEADER_LEN + 1];
    char str1[G2C_INDEX_STR1_LEN + 1];
    char my_date_str[G2C_INDEX_DATE_STR_LEN + 1];
    char my_time_str[G2C_INDEX_TIME_STR_LEN + 1];
    int i, j, k;

    /* Read the first line of header. */
    if ((bytes_read = fread(line, 1, G2C_INDEX_HEADER_LEN, f)) != G2C_INDEX_HEADER_LEN)
        return G2C_EFILE;
    line[G2C_INDEX_HEADER_LEN] = 0;

    /* Scan the line. */
    {
        char long_date_str[G2C_INDEX_HEADER_LEN + 1], long_time_str[G2C_INDEX_HEADER_LEN + 1];
        char long_str1[G2C_INDEX_HEADER_LEN + 1];

        sscanf(line, "%s %d %d %d %s %s GB2IX1", long_str1, &i, &j, &k, long_date_str, long_time_str);
        memcpy(str1, long_str1, G2C_INDEX_STR1_LEN);
        date_str[G2C_INDEX_STR1_LEN] = 0;
        memcpy(date_str, long_date_str, G2C_INDEX_DATE_STR_LEN);
        date_str[G2C_INDEX_DATE_STR_LEN] = 0;
        memcpy(time_str, long_time_str, G2C_INDEX_TIME_STR_LEN);
        time_str[G2C_INDEX_TIME_STR_LEN] = 0;
    }
    LOG((4, "str1 %s i %d j %d k %d date_str %s time_str %s", str1, i, j, k, date_str,
         time_str));

    /* Return info to caller where desired. */
    if (ip)
        *ip = i;
    if (jp)
        *jp = j;
    if (kp)
        *kp = k;
    if (date_str)
        strncpy(date_str, my_date_str, G2C_INDEX_DATE_STR_LEN + 1);
    if (time_str)
        strncpy(time_str, my_time_str, G2C_INDEX_TIME_STR_LEN + 1);

    return G2C_NOERROR;
}

/**
 * Read the second header record of an index file
 *
 * This function reads the second of two 81-byte header records of an
 * index file.
 *
 * @param f Pointer to open FILE.
 * @param skipp Pointer that gets number of bytes to skip before index
 * records. Ignored if NULL.
 * @param total_lenp Pointer that gets number of bytes in each index
 * record. Ignored if NULL.
 * @param num_recp Pointer that gets number of index records in the
 * file. Ignored if NULL.
 * @param basename Pointer to char array of size
 * ::G2C_INDEX_BASENAME_LEN + 1 which will get the basename string from the
 * second header record. Ignored if NULL.
 * @param index_version The version of the index, 1 for legacy, 2 to
 * allow for > 2 GB GRIB2 files.
 *
 * @returns 0 for success, error code otherwise.
 *
 * @author Edward Hartnett @date 9/10/23
*/
static int
read_hdr_rec2(FILE *f, int *skipp, int *total_lenp, int *num_recp,
              char *basename, int *index_version)
{
    size_t bytes_read;
    char line[G2C_INDEX_HEADER_LEN + 1];
    int skip;
    int total_len, num_rec;
    char my_basename[G2C_INDEX_BASENAME_LEN + 1];

    /* Read the second line of header. */
    if ((bytes_read = fread(line, 1, G2C_INDEX_HEADER_LEN, f)) != G2C_INDEX_HEADER_LEN)
        return G2C_EFILE;
    line[G2C_INDEX_HEADER_LEN] = 0;
    /* Scan the line. Hard! */
    {
        char long_basename[G2C_INDEX_HEADER_LEN + 1];
        sscanf(line, "IX%dFORM: %d %d %d %s", index_version, &skip, &total_len,
               &num_rec, long_basename);
        memcpy(my_basename, long_basename, G2C_INDEX_BASENAME_LEN);
        my_basename[G2C_INDEX_BASENAME_LEN] = 0;
    }

    /* Return info to caller where desired. */
    if (skipp)
        *skipp = skip;
    if (total_lenp)
        *total_lenp = total_len;
    if (num_recp)
        *num_recp = num_rec;
    if (basename)
        strncpy(basename, my_basename, G2C_INDEX_BASENAME_LEN + 1);

    return G2C_NOERROR;
}

/**
 * Open a GRIB1 index file and read the contents.
 *
 * This function opens the GRIB2 index file and reads its metadata,
 * and opens the accompanying GRIB2 file.
 *
 * GRIB2 messages in the file are assigned a message ID, starting with
 * 0 for the first message in the file.
 *
 * Each product within a message is assigned a product ID, starting
 * with 0 for the first product in the message.
 *
 * Files opened with this function should be closed with a call
 * g2c_close() to release resources.
 *
 * @param index_file The name that will be given to the index file. An
 * existing file will be overwritten.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 10/12/22
 */
int
g2c_open_index1(const char *index_file)
{
    FILE *f;
    int i, j, k;
    char date_str[G2C_INDEX_DATE_STR_LEN + 1];
    char time_str[G2C_INDEX_TIME_STR_LEN + 1];
    int skip, total_len, num_rec;
    char basename[G2C_INDEX_BASENAME_LEN + 1];
    size_t file_pos = G2C_INDEX_HEADER_LEN * 2;
    unsigned char pds_val[G2C_INDEX1_PDS_VAL_LEN];
    unsigned char gds_val[G2C_INDEX1_GDS_VAL_LEN];
    unsigned char bms_val[G2C_INDEX1_BMS_VAL_LEN];
    unsigned char bds_val[G2C_INDEX1_BDS_VAL_LEN];
    int index_version;
    int rec;
    int ret = G2C_NOERROR;

    /* Check inputs. */
    if (!index_file)
        return G2C_EINVAL;

    LOG((2, "g2c_open_index1 index_file %s", index_file));

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Open the index file. */
    if (!(f = fopen(index_file, "rb")))
        return G2C_EFILE;

    /* Read header record apparently named after Steve Lord. */
    if ((ret = read_hdr_rec1(f, &i, &j, &k, date_str, time_str)))
        return ret;
    LOG((4, "i %d j %d k %d date_str %s time_str %s", i, j, k, date_str, time_str));

    /* Read second header record. */
    if ((ret = read_hdr_rec2(f, &skip, &total_len, &num_rec, basename, &index_version)))
        return ret;
    LOG((4, "skip %d total_len %d num_rec %d basename %s", skip, total_len, num_rec, basename));

    /* Read each index record. These is one record for each message in
       the original GRIB1 file. */
    for (rec = 0; rec < num_rec; rec++)
    {
        unsigned int b2_msg, b2_gds, b2_pds, b2_bms, b2_bds, msglen;
        unsigned char version;

        /* Move to beginning of index record. */
        if (fseek(f, file_pos, SEEK_SET))
        {
            ret = G2C_EFILE;
            break;
        }

        /* Read the index1 record. */
        LOG((4, "reading index1 record at file position %ld", ftell(f)));
        if ((ret = g2c_start_index1_record(f, G2C_FILE_READ, &b2_msg, &b2_pds, &b2_gds,
                                           &b2_bms, &b2_bds, &msglen, &version, pds_val,
                                           gds_val, bms_val, bds_val, NULL, NULL, NULL)))
            break;

        LOG((4, "b2_msg %d b2_pds %d b2_gds %d b2_bms %d b2_bds %d msglen %d version %d",
             b2_msg, b2_gds, b2_pds, b2_bms, b2_bds, msglen, version));
        printf("b2_msg %d b2_pds %d b2_gds %d b2_bms %d b2_bds %d msglen %d version %d\n",
               b2_msg, b2_gds, b2_pds, b2_bms, b2_bds, msglen, version);

        /* Move the file position to the start of the next index record. */
        file_pos += total_len;
    } /* next rec */

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    /* Close the index file. */
    if (!ret)
        fclose(f);

    return ret;
}

/**
 * Open a GRIB2 file with the help of an index file.
 *
 * The index file, generated by the grb2index utility, of the
 * g2c_write_index() function, contains the byte offsets for the
 * sections of each message in the GRIB2 file. When a GRIB2 file is
 * opened with an index file, the library does not have to scan the
 * file to locate all metadata.
 *
 * @param data_file The name of the data file to which the index applies.
 * @param index_file The name that will be given to the index file. An
 * existing file will be overwritten.
 * @param mode Open mode flags.
 * @param g2cid Pointer that gets the g2cid for this file. Ignored if
 * NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 10/12/22
 */
int
g2c_open_index(const char *data_file, const char *index_file, int mode,
               int *g2cid)
{
    FILE *f;
    size_t bytes_read;
    int ret = G2C_NOERROR;

    /* Check inputs. */
    if (!index_file || !data_file || !g2cid)
        return G2C_EINVAL;
    if (strlen(data_file) > G2C_MAX_NAME)
        return G2C_ENAMETOOLONG;

    LOG((2, "g2c_open_index index_file %s", index_file));

    /* Open the index file. */
    if (!(f = fopen(index_file, "rb")))
        return G2C_EFILE;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Remember file metadata. */
    ret = g2c_add_file(data_file, mode, g2cid);

    /* Read the header. */
    if (!ret)
    {
        char line[G2C_INDEX_HEADER_LEN + 1];
        char str1[G2C_INDEX_STR1_LEN + 1], date_str[G2C_INDEX_DATE_STR_LEN + 1], time_str[G2C_INDEX_TIME_STR_LEN + 1];
        int i, j, k;
        int skip, total_len, num_rec;
        char basename[G2C_INDEX_BASENAME_LEN + 1];
        size_t file_pos = G2C_INDEX_HEADER_LEN * 2;
        int index_version;
        int rec;

        /* Read the first line of header. */
        if ((bytes_read = fread(line, 1, G2C_INDEX_HEADER_LEN, f)) != G2C_INDEX_HEADER_LEN)
            return G2C_EFILE;
        line[G2C_INDEX_HEADER_LEN] = 0;
        /* Scan the line. */
        {
            char long_date_str[G2C_INDEX_HEADER_LEN + 1], long_time_str[G2C_INDEX_HEADER_LEN + 1];
            char long_str1[G2C_INDEX_HEADER_LEN + 1];

            sscanf(line, "%s %d %d %d %s %s GB2IX1", long_str1, &i, &j, &k, long_date_str, long_time_str);
            memcpy(str1, long_str1, G2C_INDEX_STR1_LEN);
            date_str[G2C_INDEX_STR1_LEN] = 0;
            memcpy(date_str, long_date_str, G2C_INDEX_DATE_STR_LEN);
            date_str[G2C_INDEX_DATE_STR_LEN] = 0;
            memcpy(time_str, long_time_str, G2C_INDEX_TIME_STR_LEN);
            time_str[G2C_INDEX_TIME_STR_LEN] = 0;
        }
        LOG((1, "str1 %s i %d j %d k %d date_str %s time_str %s", str1, i, j, k, date_str,
             time_str));

        /* Read the second line of header. */
        if ((bytes_read = fread(line, 1, G2C_INDEX_HEADER_LEN, f)) != G2C_INDEX_HEADER_LEN)
            return G2C_EFILE;
        line[G2C_INDEX_HEADER_LEN] = 0;
        /* Scan the line. Hard! */
        {
            char long_basename[G2C_INDEX_HEADER_LEN + 1];
            sscanf(line, "IX%dFORM: %d %d %d %s", &index_version, &skip, &total_len,
                   &num_rec, long_basename);
            memcpy(basename, long_basename, G2C_INDEX_BASENAME_LEN);
            basename[G2C_INDEX_BASENAME_LEN] = 0;
        }
        LOG((1, "skip %d total_len %d num_rec %d basename %s", skip, total_len, num_rec, basename));

        /* Read each index record. */
        for (rec = 0; rec < num_rec; rec++)
        {
            int reclen, msgint, local, gds, pds, drs, bms, data;
            size_t local8, gds8, pds8, drs8, bms8, data8;
            size_t msglen, msg;
            unsigned char version, discipline;
            short fieldnum;

            /* Move to beginning of index record. */
            if (fseek(f, file_pos, SEEK_SET))
            {
                ret = G2C_EFILE;
                break;
            }

            /* Read the index record. */
            LOG((4, "reading index record at file position %ld, index_version %d",
                 ftell(f), index_version));
            if (index_version == 1)
            {
                if ((ret = g2c_start_index_record(f, G2C_FILE_READ, &reclen, &msgint, &local, &gds, &pds,
                                                  &drs, &bms, &data, &msglen, &version, &discipline, &fieldnum)))
                    break;
                msg = msgint;
                local8 = local;
                gds8 = gds;
                pds8 = pds;
                drs8 = drs;
                bms8 = bms;
                data8 = data;
            }
            else
            {
                if ((ret = g2c_start_index_record_lf(f, G2C_FILE_READ, &reclen, &msg, &local8, &gds8, &pds8,
                                                     &drs8, &bms8, &data8, &msglen, &version, &discipline, &fieldnum)))
                    break;
            }

            LOG((1, "reclen %d msg %ld local8 %d gds8 %d pds8 %d drs8 %d bms8 %d data8 %d "
                    "msglen %ld version %d discipline %d fieldnum %d",
                 reclen, msg, local8, gds8, pds8, drs8, bms8, data8, msglen,
                 version, discipline, fieldnum));

            /* Read the metadata for sections 3, 4, and 5 from
             * the index record. */
            {
                unsigned int sec_len;
                unsigned char sec_num;
                int s;
                G2C_MESSAGE_INFO_T *msgp;
                int sec_id = 0;
                int sec8_found = 0;

                /* Allocate storage for message. */
                if ((ret = add_msg(&g2c_file[*g2cid], rec, msg, msglen, 0, &msgp)))
                    break;
                msgp->discipline = discipline;
                msgp->bytes_to_local = local8;
                msgp->bytes_to_bms = bms8;
                msgp->bytes_to_data = data8;
                msgp->master_version = version;

                /* Read section 1. */
                if ((ret = g2c_rw_section1_metadata(f, G2C_FILE_READ, msgp)))
                    break;
                if ((ret = g2c_log_section1(msgp)))
                    break;

                LOG((4, "reading section info at file position %ld", ftell(f)));

                /* Add a new section to our list of sections. */
                for (s = 3; s < 8 && !sec8_found; s++)
                {
                    size_t bytes_to_sec = gds8; /* Correct value for section 3. */

                    /* For sections 3, 4, 5, read the section length
                     * and number from the index record. */
                    if (s < 7)
                    {
                        if ((ret = g2c_file_io_uint(f, G2C_FILE_READ, &sec_len)))
                            return ret;
                        if ((ret = g2c_file_io_ubyte(f, G2C_FILE_READ, &sec_num)))
                            return ret;
                    }
                    else
                    {
                        /* For section 7, the length of the section is
                         * not in the index file, but is needed for
                         * data read operations. So we will use the
                         * open data file and get the length of this
                         * section. */
                        LOG((4, "seeking to section 7 in data file, position msgp->bytes_to_msg %ld data %ld",
                             msgp->bytes_to_msg, data8));
                        if (fseek(g2c_file[*g2cid].f, msgp->bytes_to_msg + data8, SEEK_SET))
                        {
                            ret = G2C_EFILE;
                            break;
                        }
                        if ((ret = g2c_file_io_uint(g2c_file[*g2cid].f, G2C_FILE_READ, &sec_len)))
                            return ret;
                        if ((ret = g2c_file_io_ubyte(g2c_file[*g2cid].f, G2C_FILE_READ, &sec_num)))
                            return ret;
                        LOG((4, "read section 7 info from data file. sec_len %d sec_num %d",
                             sec_len, sec_num));
                    }

                    /* Select the value from the index record which is
                     * the number of bytes to section s. */
                    if (sec_num == 4)
                        bytes_to_sec = pds8;
                    else if (sec_num == 5)
                        bytes_to_sec = drs8;
                    else if (sec_num == 6)
                        bytes_to_sec = bms8;
                    else if (sec_num == 7)
                        bytes_to_sec = data8;
                    else if (sec_num == 8)
                        sec8_found++;

                    /* Check some stuff. */
                    if (s < 6 && sec_num != s)
                    {
                        ret = G2C_EBADSECTION;
                        break;
                    }
                    if (sec_num == 4)
                        if (fieldnum < 0) /* to silence warning */
                        {
                            ret = G2C_EBADSECTION;
                            break;
                        }

                    /* Read the section info from the index file,
                     * using the same functions that read it from the
                     * GRIB2 data file. */
                    LOG((4, "about to add_section sec_id %d sec_len %d bytes_to_sec %ld sec_num %d",
                         sec_id, sec_len, bytes_to_sec, sec_num));
                    if ((ret = add_section(f, msgp, sec_id++, sec_len, bytes_to_sec, sec_num)))
                        break;
                } /* next section */

                /* If anything went wrong, give up. */
                if (ret)
                    break;
            }

            /* Move the file position to the start of the next index record. */
            file_pos += reclen;
        } /* next rec */
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    /* Close the index file. */
    if (!ret)
        fclose(f);

    return ret;
}
