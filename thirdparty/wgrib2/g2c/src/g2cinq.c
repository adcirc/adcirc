/**
 * @file
 * @brief Inquiry functions.
 *
 * @author Ed Hartnett @date 10/21/22
 */
#include "grib2_int.h"
#include <stdarg.h>

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** If pthreads are enabled, use externally-defined mutex for
 * thread-safety. */
EXTERN_MUTEX(m);

/**
 * Learn about a GRIB2 file.
 *
 * @param g2cid ID of the opened file, as from g2c_open().
 * @param num_msg Pointer that gets the number of messages in the
 * file. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 *
 * @author Ed Hartnett @date 10/21/22
 */
int
g2c_inq(int g2cid, int *num_msg)
{
    int ret = G2C_NOERROR;

    /* Check input parameters. */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Find the open file. */
    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    /* If the caller wants to know the number of messages, tell
     * them. */
    if (!ret)
        if (num_msg)
            *num_msg = g2c_file[g2cid].num_messages;

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

/**
 * Learn about a GRIB2 message.
 *
 * @param g2cid ID of the opened file, as from g2c_open().
 * @param msg_num Number of the message in the file, starting with the
 * first message as 0.
 * @param discipline Pointer that gets the discipline from the
 * message. Ignored if NULL.
 * @param num_fields Pointer that gets the number of fields in the
 * message. Ignored if NULL.
 * @param num_local Pointer that gets the number of local sections in
 * the message. Ignored if NULL.
 * @param center Pointer that gets the code for the producing center from
 * the message. Ignored if NULL.
 * @param subcenter Pointer that gets the code for the producing subcenter from
 * the message. Ignored if NULL.
 * @param master_version Pointer that gets the master version from
 * the message. Ignored if NULL.
 * @param local_version Pointer that gets the local version from
 * the message. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 * - ::G2C_ENOMSG Message not found.
 *
 * @author Ed Hartnett @date 10/21/22
 */
int
g2c_inq_msg(int g2cid, int msg_num, unsigned char *discipline, int *num_fields,
            int *num_local, short *center, short *subcenter, unsigned char *master_version,
            unsigned char *local_version)
{
    G2C_MESSAGE_INFO_T *msg;
    int ret = G2C_NOERROR;

    /* Check input parameters. */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (msg_num < 0)
        return G2C_EINVAL;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Find the open file. */
    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    /* Find the file and message. */
    if (!ret)
    {
        ret = G2C_ENOMSG;
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
        {
            if (msg->msg_num == msg_num)
            {
                if (discipline)
                    *discipline = msg->discipline;
                if (num_fields)
                    *num_fields = msg->num_fields;
                if (num_local)
                    *num_local = msg->num_local;
                if (center)
                    *center = msg->center;
                if (subcenter)
                    *subcenter = msg->subcenter;
                if (master_version)
                    *master_version = msg->master_version;
                if (local_version)
                    *local_version = msg->local_version;
                ret = G2C_NOERROR;
                break;
            }
        }
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

/**
 * Learn about the date/time information in a GRIB2 message.
 *
 * @param g2cid ID of the opened file, as from g2c_open().
 * @param msg_num Number of the message in the file, starting with the
 * first message as 0.
 * @param sig_ref_time Pointer that gets Significane of reference time
 * value from the message. Ignored if NULL.
 * @param year Pointer that gets the year from the message. Ignored if
 * NULL.
 * @param month Pointer that gets the month from the message. Ignored
 * if NULL.
 * @param day Pointer that gets the day from the message. Ignored if
 * NULL.
 * @param hour Pointer that gets the hour from the message. Ignored if
 * NULL.
 * @param minute Pointer that gets the minute from the
 * message. Ignored if NULL.
 * @param second Pointer that gets the seconds from the
 * message. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 * - ::G2C_ENOMSG Message not found.
 *
 * @author Ed Hartnett @date 10/22/22
 */
int
g2c_inq_msg_time(int g2cid, int msg_num, unsigned char *sig_ref_time, short *year,
                 unsigned char *month, unsigned char *day, unsigned char *hour,
                 unsigned char *minute, unsigned char *second)
{
    G2C_MESSAGE_INFO_T *msg;
    int ret = G2C_NOERROR;

    /* Check input parameters. */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (msg_num < 0)
        return G2C_EINVAL;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Find the open file. */
    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    /* Find the file and message. */
    if (!ret)
    {
        ret = G2C_ENOMSG;
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
        {
            if (msg->msg_num == msg_num)
            {
                if (sig_ref_time)
                    *sig_ref_time = msg->sig_ref_time;
                if (year)
                    *year = msg->year;
                if (month)
                    *month = msg->month;
                if (day)
                    *day = msg->day;
                if (hour)
                    *hour = msg->hour;
                if (minute)
                    *minute = msg->minute;
                if (second)
                    *second = msg->second;
                ret = G2C_NOERROR;
                break;
            }
        }
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

/**
 * Inquire about a product.
 *
 * @param g2cid File ID.
 * @param msg_num Message number.
 * @param prod_num Product number.
 * @param pds_template_len PDS template length. Ignored if NULL.
 * @param pds_template Pointer that gets the PDS template. Ignored if NULL.
 * @param gds_template_len GDS template length. Ignored if NULL.
 * @param gds_template Pointer that gets the GDS template. Ignored if NULL.
 * @param drs_template_len The DRS template length. Ignored if NULL.
 * @param drs_template Pointer that gets the DRS template. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 * - ::G2C_ENOMSG Message not found.
 * - ::G2C_ENOPRODUCT Product not found.
 * - ::G2C_ENOSECTION Section not found.
 *
 * @author Ed Hartnett @date 10/21/22
 */
int
g2c_inq_prod(int g2cid, int msg_num, int prod_num, int *pds_template_len,
             long long int *pds_template, int *gds_template_len, long long int *gds_template,
             int *drs_template_len, long long int *drs_template)
{
    G2C_MESSAGE_INFO_T *msg;
    G2C_SECTION_INFO_T *sec4, *sec3, *sec5;
    int t;
    int ret = G2C_NOERROR;

    /* Is this an open GRIB2 file? */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (msg_num < 0 || prod_num < 0)
        return G2C_EINVAL;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    /* Find the message. */
    if (!ret)
    {
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
            if (msg->msg_num == msg_num)
                break;
        if (!msg)
            ret = G2C_ENOMSG;
    }

    /* Find the product. After this, sec4 will point to the
     * appropropriate section 4 G2C_SECTION_INFO_T. */
    if (!ret)
    {
        for (sec4 = msg->sec; sec4; sec4 = sec4->next)
            if (sec4->sec_num == 4 && ((G2C_SECTION4_INFO_T *)sec4->sec_info)->field_num == prod_num)
                break;
        if (!sec4)
            ret = G2C_ENOPRODUCT;
        /* sec4_info = (G2C_SECTION4_INFO_T *)sec4->sec_info; */
    }

    /* Return the info to the caller. */
    if (!ret)
    {
        if (pds_template_len)
            *pds_template_len = sec4->template_len;
        if (pds_template)
            for (t = 0; t < sec4->template_len; t++)
                pds_template[t] = sec4->template[t];
    }

    /* Find the GDS. */
    if (!ret)
    {
        for (sec3 = sec4->prev; sec3; sec3 = sec3->prev)
            if (sec3->sec_num == 3)
                break;
        if (!sec3)
            ret = G2C_ENOSECTION;
    }

    /* Return the info to the caller. */
    if (!ret)
    {
        if (gds_template_len)
            *gds_template_len = sec3->template_len;
        if (gds_template)
            for (t = 0; t < sec3->template_len; t++)
                gds_template[t] = sec3->template[t];
    }

    /* Find the DRS. */
    if (!ret)
    {
        for (sec5 = sec4->next; sec5; sec5 = sec5->next)
            if (sec5->sec_num == 5)
                break;
        if (!sec5)
            ret = G2C_ENOSECTION;
    }

    /* Return the info to the caller. */
    if (!ret)
    {
        if (drs_template_len)
            *drs_template_len = sec5->template_len;
        if (drs_template)
            for (t = 0; t < sec5->template_len; t++)
                drs_template[t] = sec5->template[t];
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

/**
 * Learn about the one of the dimensions of a GRIB2 product. This
 * function will return the size, name, and values along the
 * dimension.
 *
 * @param g2cid ID of the opened file, as from g2c_open().
 * @param msg_num Number of the message in the file, starting with the
 * first message as 0.
 * @param prod_num Product number.
 * @param dim_num Dimension number, with the first dimension as 0.
 * @param len Pointer that gets the length of this dimension. Ignored if NULL.
 * @param name Pointer that gets the name of this dimension. Must have
 * memory of size G2C_MAX_NAME. Ignored if NULL.
 * @param val Pointer that gets array of dimension values, of length
 * len. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 * - ::G2C_ENOMSG Message not found.
 * - ::G2C_ENOPRODUCT Product not found.
 * - ::G2C_ENOSECTION GDS not found.
 *
 * @author Ed Hartnett @date 10/21/22
 */
int
g2c_inq_dim(int g2cid, int msg_num, int prod_num, int dim_num, size_t *len,
            char *name, float *val)
{
    G2C_MESSAGE_INFO_T *msg;
    G2C_SECTION_INFO_T *sec4, *sec3;
    G2C_DIM_INFO_T *dim;
    int d;
    int ret = G2C_NOERROR;

    /* Are these valid IDs? */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (msg_num < 0 || prod_num < 0 || dim_num < 0)
        return G2C_EINVAL;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Find the file. */
    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    /* Find the message. */
    if (!ret)
    {
        for (msg = g2c_file[g2cid].msg; msg; msg = msg->next)
            if (msg->msg_num == msg_num)
                break;
        if (!msg)
            ret = G2C_ENOMSG;
    }

    /* Find the product. After this, sec4 will point to the
     * appropropriate section 4 G2C_SECTION_INFO_T. */
    if (!ret)
    {
        for (sec4 = msg->sec; sec4; sec4 = sec4->next)
            if (sec4->sec_num == 4 && ((G2C_SECTION4_INFO_T *)sec4->sec_info)->field_num == prod_num)
                break;
        if (!sec4)
            ret = G2C_ENOPRODUCT;
        /* sec4_info = (G2C_SECTION4_INFO_T *)sec4->sec_info; */
    }

    /* Find the GDS. */
    if (!ret)
    {
        for (sec3 = sec4->prev; sec3; sec3 = sec3->prev)
            if (sec3->sec_num == 3)
                break;
        if (!sec3)
            ret = G2C_ENOSECTION;
        dim = &((G2C_SECTION3_INFO_T *)sec3->sec_info)->dim[dim_num];
    }

    /* Give the caller the info they want. */
    if (!ret)
    {
        if (len)
            *len = dim->len;
        if (name)
            strncpy(name, dim->name, G2C_MAX_NAME);
        if (val)
            for (d = 0; d < dim->len; d++)
                val[d] = dim->value[d];
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

/**
 * Learn about the one of the dimensions of a GRIB2 product. This
 * function will return the size and name of the dimension.
 *
 * @param g2cid ID of the opened file, as from g2c_open().
 * @param msg_num Number of the message in the file, starting with the
 * first message as 0.
 * @param prod_num Product number.
 * @param dim_num Dimension number, with the first dimension as 0.
 * @param len Pointer that gets the length of this dimension. Ignored if NULL.
 * @param name Pointer that gets the name of this dimension. Must have
 * memory of size G2C_MAX_NAME. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID File ID not found.
 * - ::G2C_ENOMSG Message not found.
 * - ::G2C_ENOPRODUCT Product not found.
 * - ::G2C_ENOSECTION GDS not found.
 *
 * @author Ed Hartnett @date 10/21/22
 */
int
g2c_inq_dim_info(int g2cid, int msg_num, int prod_num, int dim_num, size_t *len,
                 char *name)
{
    return g2c_inq_dim(g2cid, msg_num, prod_num, dim_num, len, name, NULL);
}
