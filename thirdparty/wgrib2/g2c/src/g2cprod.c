/**
 * @file
 * @brief Product functions for the g2c library.
 * @author Ed Hartnett @date Oct 1, 2022
 */

#include "grib2_int.h"

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** If pthreads are enabled, use externally-defined mutex for
 * thread-safety. */
EXTERN_MUTEX(m);

/**
 * Read the data for a product.
 *
 * @param g2cid File ID.
 * @param msg_num Message number in file (first message in file is message 0).
 * @param prod_num Product number in message (first product in message is product 0).
 * @param num_data_points Pointer that gets the number of data points
 * in the product. Ignored if NULL.
 * @param data Pointer that gets the data. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date Sep 28, 2022
 */
int
g2c_get_prod(int g2cid, int msg_num, int prod_num, int *num_data_points, float *data)
{
    G2C_MESSAGE_INFO_T *msg;
    G2C_SECTION_INFO_T *sec3, *sec4, *sec5, *sec7;
    G2C_SECTION3_INFO_T *sec3_info;
    /* G2C_SECTION4_INFO_T *sec4_info; */
    G2C_SECTION5_INFO_T *sec5_info;
    unsigned char *buf;
    size_t bytes_read;
    int ret = G2C_NOERROR;

    /* Check inputs. */
    if (g2cid < 0 || g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (msg_num < 0 || prod_num < 0)
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

    /* Find the grid definiton section, section 3. It will come
     * earlier in the list. */
    if (!ret)
    {
        for (sec3 = sec4; sec3; sec3 = sec3->prev)
            if (sec3->sec_num == 3)
                break;
        if (!sec3)
            ret = G2C_ENOSECTION;
        sec3_info = (G2C_SECTION3_INFO_T *)sec3->sec_info;
    }

    /* Find the section 5, data representation section, to learn how
     * this product is compressed. Section 5 is after section 4 in the
     * list. */
    if (!ret)
    {
        for (sec5 = sec4; sec5; sec5 = sec5->next)
            if (sec5->sec_num == 5)
                break;
        if (!sec5)
            ret = G2C_ENOSECTION;
        sec5_info = (G2C_SECTION5_INFO_T *)sec5->sec_info;
    }

    /* Find the section 7, data section. */
    if (!ret)
    {
        for (sec7 = sec5; sec7; sec7 = sec7->next)
            if (sec7->sec_num == 7)
                break;
        if (!sec7)
            ret = G2C_ENOSECTION;
    }

    /* Give the caller number of data points, if desired. */
    if (!ret)
    {
        if (num_data_points)
            *num_data_points = sec5_info->num_data_points;
    }

    /* Allocate a char buffer to hold the packed data. */
    if (data)
    {
        if (!ret)
            if (!(buf = malloc(sizeof(char) * sec7->sec_len)))
                ret = G2C_ENOMEM;

        /* Jump to this section in the file. */
        if (!ret)
            if (fseek(g2c_file[g2cid].f, sec7->bytes_to_sec + sec7->msg->bytes_to_msg, SEEK_SET))
                ret = G2C_ERROR;

        /* Read the product into a char buffer. */
        if (!ret)
            if ((bytes_read = fread(buf, 1, sec7->sec_len, g2c_file[g2cid].f)) != sec7->sec_len)
                ret = G2C_EFILE;

        /* Unpack the char buffer into a float array, which must be
         * allocated by the caller. */
        if (!ret)
            ret = g2c_unpack7(buf, sec3_info->grid_def, sec3->template_len, sec3->template,
                              sec5_info->data_def, sec5->template_len, sec5->template,
                              sec5_info->num_data_points, data);

        /* Free the char buffer. */
        free(buf);
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}
