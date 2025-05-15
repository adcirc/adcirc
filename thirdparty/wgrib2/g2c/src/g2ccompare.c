/**
 * @file
 * @brief Compare the metadata of two open GRIB2 files.
 *
 * @author Ed Hartnett @date Dec 28, 2022
 */
#include "grib2_int.h"
#include <math.h>
#include <stdarg.h>

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/**
 * Compare the metadata of two open GRIB2 files.
 *
 * @param g2cid1 Indentifier for one file.
 * @param g2cid2 Indentifier for the other file.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid parameters.
 * - ::G2C_EFILE File I/O error.
 *
 * @author Ed Hartnett @date Dec 28, 2022
 */
int
g2c_compare(int g2cid1, int g2cid2)
{
    G2C_MESSAGE_INFO_T *msg1, *msg2;
    int m;
    /* int total_fields = 0; */
    /* int i; */
    /* int ret; */

    /* Check inputs. */
    if (g2cid1 < 0 || g2cid1 > G2C_MAX_FILES || g2c_file[g2cid1].g2cid != g2cid1)
        return G2C_EBADID;
    if (g2cid2 < 0 || g2cid2 > G2C_MAX_FILES || g2c_file[g2cid2].g2cid != g2cid2)
        return G2C_EBADID;

    LOG((2, "g2c_metadata_cmp g2cid1 %d g2cid2 %d", g2cid1, g2cid2));

    /* Same number of messages? */
    if (g2c_file[g2cid1].num_messages != g2c_file[g2cid2].num_messages)
        return G2C_ERROR;

    /* Check metadata of each message. */
    msg1 = g2c_file[g2cid1].msg;
    msg2 = g2c_file[g2cid2].msg;
    for (m = 0; m < g2c_file[g2cid1].num_messages; m++)
    {
        int fld;

        if (msg1->discipline != msg2->discipline)
            return G2C_ERROR;
        if (msg1->center != msg2->center || msg1->subcenter != msg2->subcenter ||
            msg1->master_version != msg2->master_version || msg1->local_version != msg2->local_version ||
            msg1->sig_ref_time != msg2->sig_ref_time || msg1->year != msg2->year ||
            msg1->month != msg2->month || msg1->day != msg2->day || msg1->hour != msg2->hour ||
            msg1->minute != msg2->minute || msg1->second != msg2->second || msg1->status != msg2->status ||
            msg1->type != msg2->type)
            return G2C_ERROR;
        if (msg1->num_local != msg2->num_local || msg1->num_fields != msg2->num_fields)
            return G2C_ERROR;

        /* For each field, print info. */
        for (fld = 0; fld < msg1->num_fields; fld++)
        {
            G2C_SECTION_INFO_T *sec_1, *sec_2;
            G2C_SECTION_INFO_T *sec3_1, *sec3_2;
            G2C_SECTION_INFO_T *sec5_1, *sec5_2;
            G2C_SECTION3_INFO_T *sec3_info_1, *sec3_info_2;
            G2C_SECTION4_INFO_T *sec4_info_1, *sec4_info_2;
            G2C_SECTION5_INFO_T *sec5_info_1, *sec5_info_2;
            int t;

            /* Find this field (a.k.a. product, a.k.a. section 4). */
            for (sec_1 = msg1->sec; sec_1; sec_1 = sec_1->next)
            {
                sec4_info_1 = (G2C_SECTION4_INFO_T *)sec_1->sec_info;
                if (sec_1->sec_num == 4 && sec4_info_1->field_num == fld)
                    break;
            }
            if (!sec_1)
                return G2C_ENOSECTION;
            for (sec_2 = msg2->sec; sec_2; sec_2 = sec_2->next)
            {
                sec4_info_2 = (G2C_SECTION4_INFO_T *)sec_2->sec_info;
                if (sec_2->sec_num == 4 && sec4_info_2->field_num == fld)
                    break;
            }
            if (!sec_2)
                return G2C_ENOSECTION;
            if (sec4_info_1->num_coord != sec4_info_2->num_coord ||
                sec4_info_1->prod_def != sec4_info_2->prod_def)
                return G2C_ERROR;
            if (sec_1->template_len != sec_2->template_len)
                return G2C_ERROR;
            for (t = 0; t < sec_1->template_len; t++)
                if (sec_1->template[t] != sec_2->template[t])
                    return G2C_ERROR;

            /* Find the sec3 that applies to this field. */
            for (sec3_1 = sec_1; sec3_1; sec3_1 = sec3_1->prev)
                if (sec3_1->sec_num == 3)
                    break;
            if (!sec3_1)
                return G2C_ENOSECTION;
            sec3_info_1 = (G2C_SECTION3_INFO_T *)sec3_1->sec_info;
            for (sec3_2 = sec_2; sec3_2; sec3_2 = sec3_2->prev)
                if (sec3_2->sec_num == 3)
                    break;
            if (!sec3_2)
                return G2C_ENOSECTION;
            sec3_info_2 = (G2C_SECTION3_INFO_T *)sec3_2->sec_info;
            if (sec3_info_1->source_grid_def != sec3_info_2->source_grid_def ||
                sec3_info_1->num_data_points != sec3_info_2->num_data_points ||
                sec3_info_1->num_opt != sec3_info_2->num_opt ||
                sec3_info_1->interp_list != sec3_info_2->interp_list ||
                sec3_info_1->grid_def != sec3_info_2->grid_def)
                return G2C_ERROR;
            if (sec3_1->template_len != sec3_2->template_len)
                return G2C_ERROR;
            for (t = 0; t < sec3_1->template_len; t++)
                if (sec3_1->template[t] != sec3_2->template[t])
                    return G2C_ERROR;

            /* Find the sec5 that applies to this field. */
            for (sec5_1 = sec_1; sec5_1; sec5_1 = sec5_1->next)
                if (sec5_1->sec_num == 5)
                    break;
            if (!sec5_1)
                return G2C_ENOSECTION;
            sec5_info_1 = (G2C_SECTION5_INFO_T *)sec5_1->sec_info;
            for (sec5_2 = sec_2; sec5_2; sec5_2 = sec5_2->next)
                if (sec5_2->sec_num == 5)
                    break;
            if (!sec5_2)
                return G2C_ENOSECTION;
            sec5_info_2 = (G2C_SECTION5_INFO_T *)sec5_2->sec_info;
            if (sec5_info_1->num_data_points != sec5_info_2->num_data_points ||
                sec5_info_1->data_def != sec5_info_2->data_def)
                return G2C_ERROR;
            if (sec5_1->template_len != sec5_2->template_len)
                return G2C_ERROR;
            for (t = 0; t < sec5_1->template_len; t++)
                if (sec5_1->template[t] != sec5_2->template[t])
                    return G2C_ERROR;
        }
        msg1 = msg1->next;
        msg2 = msg2->next;
    }

    /* fprintf(f, "  \n  Total Number of Fields Found =  %d\n", total_fields); */

    return G2C_NOERROR;
}
