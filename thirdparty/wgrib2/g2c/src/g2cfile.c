/**
 * @file
 * @brief File functions for the g2c library.
 * @author Ed Hartnett @date Aug 16, 2022
 */

#include "grib2_int.h"

/** Global file information. */
G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/** Next g2cid file ID - used when opening or creating a file. */
int g2c_next_g2cid = 1;

/** Find a minimum. */
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/** Default size of read-buffer. */
#define READ_BUF_SIZE 4092

/** Number of bytes to discipline field in GRIB2 message. */
#define BYTES_TO_DISCIPLINE 6

/** Define mutex for thread-safety. */
MUTEX(m);

#define G2C_SEEKMSG_BUFSIZE 4092 /**< Size of buffer used in g2c_seekmsg(). */

/** Global file information. */
extern G2C_FILE_INFO_T g2c_file[G2C_MAX_FILES + 1];

/**
 * Search a file for the next GRIB2 Message.
 *
 * The search is terminated when a GRIB2 message is found, or the end
 * of the file is reached.
 *
 * @param g2cid ID of an open GRIB2 file, returned from
 * g2c_open()/g2c_create().
 * @param skip The number of bytes in the file to skip before
 * starting the search.
 * @param offset Pointer that gets the number of bytes to skip from
 * the beggining of the file to where the GRIB message starts. Ignored
 * if NULL.
 * @param msglen Pointer that gets the number of bytes in message (set
 * to 0, if no message found). Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID Bad g2cid.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date 2022-09-11
 */
int
g2c_seekmsg(int g2cid, size_t skip, size_t *offset, size_t *msglen)
{
    size_t k4;
    int k, lim;
    int end;
    unsigned char *cbuf;
    size_t bytes_read = G2C_SEEKMSG_BUFSIZE;
    size_t my_msglen = 0, my_offset = 0, ipos;

    /* Find the open file struct. */
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    LOG((3, "g2c_seekgb skip %ld", skip));

    /* Get memory to read in some of the file. */
    if (!(cbuf = malloc(G2C_SEEKMSG_BUFSIZE)))
        return G2C_ENOMEM;

    ipos = skip;

    /* Loop until grib message is found. */
    while (my_msglen == 0 && bytes_read == G2C_SEEKMSG_BUFSIZE)
    {
        /* Read partial section. */
        if (fseek(g2c_file[g2cid].f, ipos, SEEK_SET))
        {
            free(cbuf);
            return G2C_EFILE;
        }
        bytes_read = fread(cbuf, sizeof(unsigned char), G2C_SEEKMSG_BUFSIZE, g2c_file[g2cid].f);
        lim = bytes_read - 8;

        /* Look for 'GRIB...2' in partial section. */
        for (k = 0; k < lim; k++)
        {
            if (!strncmp((char *)&cbuf[k], G2C_MAGIC_HEADER, strlen(G2C_MAGIC_HEADER)) && cbuf[k + 7] == 2)
            {
                /* Find the length of the message. This is stored as
                 * an 8-byte big-endian integer starting at positon
                 * 8 in cbuf. */
                my_msglen = hton64(*(size_t *)&cbuf[k + 8]);

                LOG((4, "my_msglen %ld", my_msglen));

                /* Read the last 4 bytes of the message. */
                if (fseek(g2c_file[g2cid].f, ipos + k + my_msglen - 4, SEEK_SET))
                {
                    free(cbuf);
                    return G2C_EFILE;
                }

                if ((k4 = fread(&end, 4, 1, g2c_file[g2cid].f)) != 1)
                {
                    free(cbuf);
                    return G2C_EFILE;
                }

                /* Look for '7777' at end of grib message. */
                if (k4 == 1 && end == 926365495)
                {
                    /* GRIB message found. */
                    my_offset = ipos + k;
                    LOG((4, "found end of message my_offset %ld", my_offset));
                    break;
                }
            }
        }
        ipos = ipos + lim;
    }

    /* Free resources. */
    free(cbuf);

    /* Return information to the caller. */
    if (offset)
        *offset = my_offset;
    if (msglen)
        *msglen = my_msglen;

    return G2C_NOERROR;
}

/** Search a file for the next GRIB1 or GRIB2 message.
 *
 * A grib message is identified by its indicator section,
 * i.e. an 8-byte sequence with 'GRIB' in bytes 1-4 and a '1' or '2'
 * in byte 8. If found, the length of the message is decoded from
 * bytes 5-7. The search is done over a given section of the file. The
 * search is terminated if an eof or i/o error is encountered.
 *
 * @param g2cid ID of the opened grib file, returned by g2c_open().
 * @param skip_bytes Number of bytes to skip before search.
 * @param max_bytes Maximum number of bytes to search.
 * @param bytes_to_msg Pointer that gets the number of bytes to skip
 * before message.
 * @param bytes_in_msg Pointer that gets the number of bytes in
 * message (or 0 if no message found)
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID g2cid not found.
 * - ::G2C_EFILE File error.
 * - ::G2C_EINVAL Invalid input.
 *
 * @author Ed Hartnett @date 2022-08-19
 */
int
g2c_find_msg2(int g2cid, size_t skip_bytes, size_t max_bytes, size_t *bytes_to_msg,
              size_t *bytes_in_msg)
{
    size_t bytes_to_read = MIN(READ_BUF_SIZE, max_bytes);
    size_t bytes_read;
    unsigned char *buf;
    int grib_version;
    int eof = 0;
    int msg_found = 0;
    size_t num_blocks;
    size_t ftell_pos;
    int i;
    int done = 0;
    int ret = G2C_NOERROR;

    /* Check inputs. */
    if (!bytes_to_msg || !bytes_in_msg)
        return G2C_EINVAL;

    /* Find the open file struct. */
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    /* Skip some bytes if desired. */
    if (fseek(g2c_file[g2cid].f, (off_t)skip_bytes, SEEK_SET))
        return G2C_ERROR;

    /* Allocate storage to read into. */
    if (!(buf = calloc(bytes_to_read, sizeof(char))))
        return G2C_ENOMEM;

    for (num_blocks = 0; !eof && !done; num_blocks++)
    {
        /* Read some bytes. If we don't get the number expected, either a
         * read error occured, or we reached the end of file. */
        if ((ftell_pos = ftell(g2c_file[g2cid].f)) == -1)
            return G2C_EFILE;
        LOG((4, "before read ftell() is %ld (0x%x) reading %ld bytes", ftell_pos,
             ftell_pos, bytes_to_read));
        if ((bytes_read = fread(buf, 1, bytes_to_read, g2c_file[g2cid].f)) != bytes_to_read)
        {
            if (ferror(g2c_file[g2cid].f))
                ret = G2C_EFILE;
            eof++;
        }

        /* Scan for 'GRIB2' in the bytes we have read. */
        if (!ret)
        {
            for (i = 0; i < bytes_read; i++)
            {
#ifdef LOGGING
                /* if (i < 10) LOG((3, "buf[%ld] = %2.2x", i, buf[i])); */
#endif
                /* Find the beginning of a GRIB message. */
                if (buf[i] == 'G' && i < bytes_read - G2C_MAGIC_HEADER_LEN && buf[i + 1] == 'R' && buf[i + 2] == 'I' && buf[i + 3] == 'B')
                {
                    msg_found++;
                    *bytes_to_msg = ftell_pos + i;
                    grib_version = buf[i + 7];
                    LOG((3, "bytes_to_msg %ld grib_version %d", *bytes_to_msg, grib_version));
                    if (grib_version != 1 && grib_version != 2)
                    {
                        ret = G2C_EMSG;
                        done++;
                        break;
                    }
                }

                /* Find the end of a GRIB message. And then we're done. */
                if (msg_found && buf[i] == '7' && i < bytes_read - G2C_MAGIC_HEADER_LEN && buf[i + 1] == '7' && buf[i + 2] == '7' && buf[i + 3] == '7')
                {
                    msg_found--;
                    *bytes_in_msg = ftell_pos + i - *bytes_to_msg + 4;
                    LOG((3, "bytes_in_msg %ld", *bytes_in_msg));
                    ret = G2C_NOERROR;
                    done++;
                    break;
                }
            }
        }

        /* Back up 8 bytes in case the "GRIB" magic header occurred
         * within the last 8 bytes of the previous read. */
        if (!done)
            if (fseek(g2c_file[g2cid].f, (off_t)(ftell(g2c_file[g2cid].f) - G2C_MAGIC_HEADER_LEN),
                      SEEK_SET))
                return G2C_ERROR;
    }

    /* Free storage. */
    free(buf);

    return ret;
}

/** Search a file for the next GRIB1 or GRIB2 message, and read it,
 * allocating space in memory to hold the message.
 *
 * A grib message is identified by its indicator section,
 * i.e. an 8-byte sequence with 'GRIB' in bytes 1-4 and a '1' or '2'
 * in byte 8. If found, the length of the message is decoded from
 * bytes 5-7. The search is done over a given section of the file. The
 * search is terminated if an EOF or I/O error is encountered.
 *
 * @param g2cid ID of the opened grib file, returned by g2c_open().
 * @param skip_bytes The number of bytes to skip before search.
 * @param max_bytes The maximum number of bytes to search. Must be at
 * least 16.
 * @param bytes_to_msg A pointer that gets the number of bytes to skip
 * before message.
 * @param bytes_in_msg A pointer that gets the number of bytes in
 * message (or 0 if no message found)
 * @param cbuf A pointer that gets allocation of memory, into which the
 * message is copied. This memory must be freed by the caller.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EBADID g2cid not found.
 * - ::G2C_EFILE File error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOMSG No GRIB message found.
 *
 * @author Ed Hartnett @date 2022-08-20
 */
int
g2c_get_msg(int g2cid, size_t skip_bytes, size_t max_bytes, size_t *bytes_to_msg,
            size_t *bytes_in_msg, unsigned char **cbuf)
{
    size_t bytes_read;
    int ret = G2C_NOERROR;

    /* Check inputs. */
    if (!bytes_to_msg || !bytes_in_msg || !cbuf || max_bytes < G2C_MIN_MAX_BYTES)
        return G2C_EINVAL;

    LOG((2, "g2c_get_msg g2cid %d skip_bytes %ld max_bytes %ld", g2cid, skip_bytes,
         max_bytes));

    /* Find the open file struct. */
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    /* Find the start and length of the GRIB message. */
    {
        g2int bytes_to_msg_g, bytes_in_msg_g;
        seekgb(g2c_file[g2cid].f, (g2int)skip_bytes, (g2int)max_bytes, &bytes_to_msg_g,
               &bytes_in_msg_g);
        *bytes_to_msg = bytes_to_msg_g;
        *bytes_in_msg = bytes_in_msg_g;
    }
    LOG((4, "*bytes_to_msg %ld *bytes_in_msg %ld", *bytes_to_msg, *bytes_in_msg));

    /* If no message was found, return an error. */
    if (*bytes_in_msg == 0)
        return G2C_ENOMSG;

    /* Allocate storage for the GRIB message. */
    if (!(*cbuf = malloc(*bytes_in_msg)))
        return G2C_ENOMEM;

    /* Position file at start of GRIB message. */
    if (fseek(g2c_file[g2cid].f, (off_t)*bytes_to_msg, SEEK_SET))
    {
#ifdef LOGGING
        int my_errno = errno;
        LOG((0, "fseek error %s", strerror(my_errno)));
#endif
        return G2C_ERROR;
    }

    /* Read the message from the file into the buffer. */
    if ((bytes_read = fread(*cbuf, 1, *bytes_in_msg, g2c_file[g2cid].f)) != *bytes_in_msg)
        return G2C_EFILE;

#ifdef LOGGING
    {
        int i;
        for (i = 0; i < 10; i++)
            LOG((4, "cbuf[%d] = %2x", i, (*cbuf)[i]));
    }
#endif

    return ret;
}

/** Find a g2cid to use for a newly opened or created file.
 *
 * @param g2cid Pointer that gets the next available g2cid.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EINVAL - Invalid input.
 * - ::G2C_ETOOMANYFILES - Trying to open too many files at the same time.
 * '
 * @author Ed Hartnett 8/18/22
 */
static int
find_available_g2cid(int *g2cid)
{
    int i;
    int ret;

    /* Check input. */
    if (!g2cid)
        return G2C_EINVAL;

    /* Find a new g2cid. */
    for (i = 0; i < G2C_MAX_FILES + 1; i++)
    {
        int id = (i + g2c_next_g2cid) % (G2C_MAX_FILES + 1);

        /* Is this ID available? If so, we're done. */
        if (!g2c_file[id].g2cid)
        {
            *g2cid = id;
            g2c_next_g2cid = id + 1;
            ret = G2C_NOERROR;
            break;
        }
    }

    /* If we couldn't find one, all available files are already
     * open. */
    if (i == G2C_MAX_FILES + 1)
        ret = G2C_ETOOMANYFILES;

    return ret;
}

/**
 * Determine the dimension information from the section 3 metadata.
 *
 * See (GRIB2 - SECTION 3 GRID DEFINITION
 * SECTION)[https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml].
 *
 * For a list of grid definitions see [GRIB2 - TABLE 3.1 Grid
 * Definition Template
 * Number](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-1.shtml).
 *
 * @param sec G2C_SECTION3_INFO_T struct.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date Sep 15, 2022
 */
static int
determine_dims(G2C_SECTION_INFO_T *sec)
{
    G2C_DIM_INFO_T *d0, *d1;
    G2C_SECTION3_INFO_T *sec3_info;
    int d;

    sec3_info = (G2C_SECTION3_INFO_T *)(sec->sec_info);
    d0 = &(sec3_info->dim[0]);
    d1 = &(sec3_info->dim[1]);

    /* Based on the grid definition template number, and the contents
     * of the template, decide the len, name, and values of the two
     * dimensions. */
    switch (sec3_info->grid_def)
    {
    case 0:
        LOG((5, "determine_dim allocating storage for lat/lon values"));
        d0->len = sec->template[8];
        strncpy(d0->name, LATITUDE, G2C_MAX_NAME);
        if (!(d0->value = malloc(d0->len * sizeof(float))))
            return G2C_ENOMEM;
        d0->value[0] = sec->template[11];
        for (d = 1; d < d0->len; d++)
            d0->value[d] = d0->value[d - 1] - sec->template[16];

        d1->len = sec->template[7];
        strncpy(d1->name, LONGITUDE, G2C_MAX_NAME);
        if (!(d1->value = malloc(d1->len * sizeof(float))))
            return G2C_ENOMEM;
        d1->value[0] = sec->template[12];
        for (d = 1; d < d1->len; d++)
            d1->value[d] = d1->value[d - 1] - sec->template[17];
        break;
    default:
        break;
    }

    return G2C_NOERROR;
}

/**
 * Read the metadata from section 3 (Grid Definition Section) of a
 * GRIB2 message.
 *
 * When this function is called, the file cursor is positioned just
 * after the section number field in the section. The size of the
 * section, and the section number, have already been read when this
 * function is called.
 *
 * @param f FILE pointer to open GRIB2 file.
 * @param rw_flag ::G2C_FILE_WRITE if function should write,
 * ::G2C_FILE_READ (0) if it should read.
 * @param sec Pointer to the G2C_SECTION_INFO_T struct.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOTEMPLATE Can't find template.
 *
 * @author Ed Hartnett @date Sep 15, 2022
 */
int
g2c_rw_section3_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec)
{
    G2C_SECTION3_INFO_T *sec3_info = NULL;
    int maplen, needsext, map[G2C_MAX_GDS_TEMPLATE_MAPLEN];
    int t;
    int ret;

    /* Check input. */
    if (!f || !sec)
        return G2C_EINVAL;
    if (!rw_flag && sec->sec_num != 3)
        return G2C_EINVAL;

    LOG((6, "g2c_rw_section3_metadata starting to %s section 3 at file position %ld",
         rw_flag ? "write" : "read", ftell(f)));

    /* If reading, allocate storage for a new section 3. */
    if (!rw_flag)
    {
        if (!(sec3_info = calloc(sizeof(G2C_SECTION3_INFO_T), 1)))
            return G2C_ENOMEM;
    }
    else
        sec3_info = sec->sec_info;

    /* Read or write section 3. */
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &sec3_info->source_grid_def)))
        return ret;
    if ((ret = g2c_file_io_uint(f, rw_flag, &sec3_info->num_data_points)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &sec3_info->num_opt)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &sec3_info->interp_list)))
        return ret;
    if ((ret = g2c_file_io_ushort(f, rw_flag, &sec3_info->grid_def)))
        return ret;
    LOG((5, "rw_section3_metadata source_grid_def %d num_data_points %d num_opt %d interp_list %d grid_def %d",
         sec3_info->source_grid_def, sec3_info->num_data_points, sec3_info->num_opt, sec3_info->interp_list,
         sec3_info->grid_def));

    /* Look up the information about this grid. */
    if ((ret = g2c_get_grid_template(sec3_info->grid_def, &maplen, map, &needsext)))
        return ret;

    /* When reading, allocate space to hold the template info. */
    if (!rw_flag)
    {
        sec->template_len = maplen;
        if (!(sec->template = calloc(sizeof(long long int) * maplen, 1)))
            return G2C_ENOMEM;
    }

    /* Read or write the template info. */
    for (t = 0; t < maplen; t++)
    {
        if ((ret = g2c_file_io_template(f, rw_flag, map[t], &sec->template[t])))
            return ret;
        LOG((7, "template[%d] %d", t, sec->template[t]));
    }

    /* Attach sec3_info to our section data. */
    if (!rw_flag)
        sec->sec_info = sec3_info;

    /* Figure out the dimensions, if we can. */
    if (!rw_flag)
        if ((ret = determine_dims(sec)))
            return ret;

    LOG((6, "finished reading or writing section 3 at file position %ld", ftell(f)));
    return G2C_NOERROR;
}

/**
 * Read or write the metadata from section 4 (Product Definition
 * Section) of a GRIB2 message.
 *
 * When this function is called, the file cursor is positioned just
 * after the section number field in the section. The size of the
 * section, and the section number, have already been read/written
 * when this function is called.
 *
 * @param f FILE pointer to open GRIB2 file.
 * @param rw_flag ::G2C_FILE_WRITE if function should write,
 * ::G2C_FILE_READ if it should read.
 * @param sec Pointer to the G2C_SECTION_INFO_T struct.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOTEMPLATE Can't find template.
 *
 * @author Ed Hartnett @date Sep 16, 2022
 */
int
g2c_rw_section4_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec)
{
    G2C_SECTION4_INFO_T *sec4_info;
    int maplen, needsext, map[G2C_MAX_PDS_TEMPLATE_MAPLEN];
    int t;
    int ret;

    /* Check input. */
    if (!f || !sec)
        return G2C_EINVAL;
    if (!rw_flag && sec->sec_num != 4)
        return G2C_EINVAL;

    LOG((3, "read_section4_metadata rw_flag %d", rw_flag));

    /* When reading, allocate storage for a new section 4. */
    if (!rw_flag)
    {
        if (!(sec4_info = calloc(sizeof(G2C_SECTION4_INFO_T), 1)))
            return G2C_ENOMEM;
    }
    else
        sec4_info = sec->sec_info;

    /* When reading, assign a number to this field, and count the
     * number of fields in the message. */
    if (!rw_flag)
        sec4_info->field_num = sec->msg->num_fields++;

    LOG((6, "reading section 4 starting at file position %ld", ftell(f)));

    /* Read section 4. */
    if ((ret = g2c_file_io_ushort(f, rw_flag, &sec4_info->num_coord)))
        return ret;
    if ((ret = g2c_file_io_ushort(f, rw_flag, &sec4_info->prod_def)))
        return ret;
    LOG((5, "read_section4_metadata num_coord %d prod_def %d", sec4_info->num_coord, sec4_info->prod_def));

    /* Look up the information about this grid. */
    if ((ret = g2c_get_pds_template(sec4_info->prod_def, &maplen, map, &needsext)))
        return ret;
    LOG((5, "pds template maplen %d", maplen));

    /* When reading, allocate space to hold the template info. */
    if (!rw_flag)
    {
        sec->template_len = maplen;
        if (!(sec->template = calloc(sizeof(long long int) * maplen, 1)))
            return G2C_ENOMEM;
    }

    /* Read or write the template info. */
    for (t = 0; t < maplen; t++)
    {
        if ((ret = g2c_file_io_template(f, rw_flag, map[t], &sec->template[t])))
            return ret;
        LOG((7, "template[%d] %d", t, sec->template[t]));
    }

    /* When reading, attach sec4_info to our section data. */
    if (!rw_flag)
        sec->sec_info = sec4_info;

    return G2C_NOERROR;
}

/**
 * Read or write the metadata from section 5 (Data Representation
 * Section) of a GRIB2 message.
 *
 * When this function is called, the file cursor is positioned just
 * after the section number field in the section. The size of the
 * section, and the section number, have already been read when this
 * function is called.
 *
 * @param f FILE pointer to open GRIB2 file.
 * @param rw_flag ::G2C_FILE_WRITE if function should write,
 * ::G2C_FILE_READ if it should read.
 * @param sec Pointer to the G2C_SECTION_INFO_T struct.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOTEMPLATE Can't find template.
 *
 * @author Ed Hartnett @date Sep 16, 2022
 */
int
g2c_rw_section5_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec)
{
    G2C_SECTION5_INFO_T *sec5_info;
    int maplen, needsext, map[G2C_MAX_PDS_TEMPLATE_MAPLEN];
    int t;
    int ret;

    /* Check input. */
    if (!f || !sec)
        return G2C_EINVAL;
    LOG((5, "g2c_rw_section5_metadata rw_flag %d at file position %ld", rw_flag,
         ftell(f)));

    /* When reading, allocate storage for a new section 5. When
     * writing, get a pointer to the exitsing sec5_info. */
    if (!rw_flag)
    {
        if (!(sec5_info = calloc(sizeof(G2C_SECTION5_INFO_T), 1)))
            return G2C_ENOMEM;
    }
    else
        sec5_info = sec->sec_info;

    /* Read section 5. */
    if ((ret = g2c_file_io_uint(f, rw_flag, &sec5_info->num_data_points)))
        return ret;
    if ((ret = g2c_file_io_ushort(f, rw_flag, &sec5_info->data_def)))
        return ret;
    LOG((5, "g2c_rw_section5_metadata num_data_points %d data_def %d",
         sec5_info->num_data_points, sec5_info->data_def));

    /* Look up the information about this grid. */
    if ((ret = g2c_get_drs_template(sec5_info->data_def, &maplen, map, &needsext)))
        return ret;
    LOG((5, "grid template maplen %d", maplen));

    /* WHen reading, allocate space to hold the template info. */
    if (!rw_flag)
    {
        sec->template_len = maplen;
        if (!(sec->template = calloc(sizeof(long long int) * maplen, 1)))
            return G2C_ENOMEM;
    }

    /* Read or write the template info. */
    for (t = 0; t < maplen; t++)
    {
        if ((ret = g2c_file_io_template(f, rw_flag, map[t], &sec->template[t])))
            return ret;
        LOG((7, "template[%d] %d", t, sec->template[t]));
    }

    /* When reading, attach sec5_info to our section data. */
    if (!rw_flag)
        sec->sec_info = sec5_info;

    return G2C_NOERROR;
}

/**
 * Read or write the metadata from section 6 (Data Representation
 * Section) of a GRIB2 message.
 *
 * When this function is called, the file cursor is positioned just
 * after the section number field in the section. The size of the
 * section, and the section number, have already been read when this
 * function is called.
 *
 * @param f FILE pointer to open GRIB2 file.
 * @param rw_flag ::G2C_FILE_WRITE if function should write,
 * ::G2C_FILE_READ if it should read.
 * @param sec Pointer to the G2C_SECTION_INFO_T struct.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOTEMPLATE Can't find template.
 *
 * @author Ed Hartnett @date Sep 16, 2022
 */
int
g2c_rw_section6_metadata(FILE *f, int rw_flag, G2C_SECTION_INFO_T *sec)
{
    G2C_SECTION6_INFO_T *sec6_info;
    int ret;

    /* Check input. */
    if (!f || !sec)
        return G2C_EINVAL;
    LOG((4, "g2c_rw_section6_metadata rw_flag %d at file position %ld", rw_flag,
         ftell(f)));

    /* When reading, allocate storage for a new section 6. When
     * writing, get a pointer to the exitsing sec6_info. */
    if (!rw_flag)
    {
        if (!(sec6_info = calloc(sizeof(G2C_SECTION6_INFO_T), 1)))
            return G2C_ENOMEM;
    }
    else
        sec6_info = sec->sec_info;

    /* Read or write section 6. */
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &sec6_info->indicator)))
        return ret;
    LOG((4, "g2c_rw_section6_metadata indicator %d", sec6_info->indicator));

    /* When reading, attach sec6_info to our section data. */
    if (!rw_flag)
        sec->sec_info = sec6_info;

    return G2C_NOERROR;
}

/**
 * Add metadata about a new section 3, 4, 5, 6, or 7.
 *
 * @param f FILE pointer to open GRIB2 file.
 * @param msg Pointer to the G2C_MESSAGE_INFO_T struct.
 * @param sec_id 0-based section ID.
 * @param sec_len Length of section.
 * @param bytes_to_sec Number of bytes from start of message to this section.
 * @param sec_num Section number.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 *
 * @author Ed Hartnett @date Sep 12, 2022
 */
int
add_section(FILE *f, G2C_MESSAGE_INFO_T *msg, int sec_id, unsigned int sec_len,
            size_t bytes_to_sec, unsigned char sec_num)
{
    G2C_SECTION_INFO_T *sec;
    int ret;

    LOG((4, "add_section sec_id %d sec_len %d, bytes_to_sec %ld, sec_num %d",
         sec_id, sec_len, bytes_to_sec, sec_num));

    /* Allocate storage for a new section. */
    if (!(sec = calloc(sizeof(G2C_SECTION_INFO_T), 1)))
        return G2C_ENOMEM;

    /* Add sec to end of linked list. */
    if (!msg->sec)
        msg->sec = sec;
    else
    {
        G2C_SECTION_INFO_T *s;

        for (s = msg->sec; s->next; s = s->next)
            ;
        s->next = sec;
        sec->prev = s;
    }

    /* Remember values. */
    sec->msg = msg;
    sec->sec_id = sec_id;
    sec->sec_len = sec_len;
    sec->sec_num = sec_num;
    sec->bytes_to_sec = bytes_to_sec;

    switch (sec_num)
    {
    case 1:
        break;
    case 2:
        msg->num_local++;
        break;
    case 3:
        if ((ret = g2c_rw_section3_metadata(f, G2C_FILE_READ, sec)))
            return ret;
        break;
    case 4:
        if ((ret = g2c_rw_section4_metadata(f, G2C_FILE_READ, sec)))
            return ret;
        break;
    case 5:
        if ((ret = g2c_rw_section5_metadata(f, G2C_FILE_READ, sec)))
            return ret;
        break;
    case 6:
        if ((ret = g2c_rw_section6_metadata(f, G2C_FILE_READ, sec)))
            return ret;
        break;
    case 7:
        break;
    default:
        return G2C_EBADSECTION;
    }

    return G2C_NOERROR;
}

/**
 * Read Section 1.
 *
 * @param f Pointer to open file.
 * @param rw_flag ::G2C_FILE_WRITE if function should write,
 * ::G2C_FILE_READ if it should read.
 * @param msg Pointer to G2C_MESSAGE_INFO_T which will be populated
 * with the values of section 0.
 *
 * @return
 * -G2C_NOERROR No error.
 *
 * @author Ed Hartnett @date 10/16/22
 */
int
g2c_rw_section1_metadata(FILE *f, int rw_flag, G2C_MESSAGE_INFO_T *msg)
{
    unsigned char sec_num = 1;
    int ret;

    LOG((3, "g2c_rw_section1_metadata rw_flag %d", rw_flag));

    /* Read the section. */
    if ((ret = g2c_file_io_uint(f, rw_flag, (unsigned int *)&msg->sec1_len)))
        return ret;

    if ((ret = g2c_file_io_ubyte(f, rw_flag, &sec_num)))
        return ret;
    if (!rw_flag && sec_num != 1) /* When reading sec num must be 1. */
        return G2C_ENOSECTION;
    if ((ret = g2c_file_io_short(f, rw_flag, &msg->center)))
        return ret;
    if ((ret = g2c_file_io_short(f, rw_flag, &msg->subcenter)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->master_version)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->local_version)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->sig_ref_time)))
        return ret;
    if ((ret = g2c_file_io_short(f, rw_flag, &msg->year)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->month)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->day)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->hour)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->minute)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->second)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->status)))
        return ret;
    if ((ret = g2c_file_io_ubyte(f, rw_flag, &msg->type)))
        return ret;

    /* Section 1 may contain optional numbers at the end of the
     * section. The sec1_len tells us if there are extra values. If
     * so, skip them. */
    if (msg->sec1_len > G2C_SECTION1_BYTES)
        fseek(f, msg->sec1_len - G2C_SECTION1_BYTES, SEEK_CUR);

    return G2C_NOERROR;
}

/**
 * Read the file to get metadata about a message.
 *
 * @param msg Pointer to the G2C_MESSAGE_INFO_T struct for this
 * message.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 *
 * @author Ed Hartnett @date Sep 12, 2022
 */
static int
read_msg_metadata(G2C_MESSAGE_INFO_T *msg)
{
    int total_read = G2C_SECTION0_BYTES;
    int sec_id = 0;
    int ret;

    LOG((6, "read_msg_metadata file position %ld", ftell(msg->file->f)));

    /* Read section 0. */
    if (fseek(msg->file->f, msg->bytes_to_msg + BYTES_TO_DISCIPLINE, SEEK_SET))
        return G2C_EFILE;
    if ((fread(&msg->discipline, ONE_BYTE, 1, msg->file->f)) != 1)
        return G2C_EFILE;

    /* Skip to section 1. */
    if (fseek(msg->file->f, 9, SEEK_CUR))
        return G2C_EFILE;

    /* Read section 1. */
    if ((ret = g2c_rw_section1_metadata(msg->file->f, G2C_FILE_READ, msg)))
        return ret;
    total_read += msg->sec1_len;

    /* Read the sections. */
    while (total_read < msg->bytes_in_msg - FOUR_BYTES)
    {
        unsigned int sec_len;
        unsigned char sec_num;

        LOG((4, "reading new section at file position %ld", ftell(msg->file->f)));

        /* Read section length. */
        if ((ret = g2c_file_io_uint(msg->file->f, G2C_FILE_READ, &sec_len)))
            return ret;

        /* A section length of 926365495 indicates we've reached
         * section 8, the end of the message. */
        if (sec_len != 926365495)
        {
            /* Read section number. */
            if ((ret = g2c_file_io_ubyte(msg->file->f, G2C_FILE_READ, &sec_num)))
                return ret;
            LOG((4, "sec_len %d sec_num %d", sec_len, sec_num));

            /* Add a new section to our list of sections. */
            if ((ret = add_section(msg->file->f, msg, sec_id++, sec_len, total_read, sec_num)))
                return G2C_EBADSECTION;

            /* Skip to next section. */
            total_read += sec_len;
            LOG((4, "total_read %d", total_read));
            if (fseek(msg->file->f, msg->bytes_to_msg + total_read, SEEK_SET))
                return G2C_EFILE;
        }
        else
            break;
    }

    return G2C_NOERROR;
}

/**
 * Add new message to linked list.
 *
 * @param file Pointer to the G2C_FILE_INFO_T for this file.
 * @param msg_num Number of the message in file (0-based).
 * @param bytes_to_msg Number of bytes to the start of the message in
 * the file.
 * @param bytes_in_msg Length of message in bytes.
 * @param read_file Set to true to cause metadata to be read from a
 * GRIB2 data file.
 * @param msg Pointer to a pointer that will get the location of the
 * newly created ::G2C_MESSAGE_INFO_T object. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 *
 * @author Ed Hartnett @date Sep 12, 2022
 */
int
add_msg(G2C_FILE_INFO_T *file, int msg_num, size_t bytes_to_msg, size_t bytes_in_msg,
        int read_file, G2C_MESSAGE_INFO_T **msg)
{
    G2C_MESSAGE_INFO_T *my_msg;
    int ret;

    LOG((4, "add_msg msg_num %d bytes_to_msg %ld bytes_in_msg %ld read_file %d",
         msg_num, bytes_to_msg, bytes_in_msg, read_file));

    /* Allocate storage for a new message. */
    if (!(my_msg = calloc(sizeof(G2C_MESSAGE_INFO_T), 1)))
        return G2C_ENOMEM;

    /* Add my_msg to end of linked list. */
    if (!file->msg)
        file->msg = my_msg;
    else
    {
        G2C_MESSAGE_INFO_T *m;

        for (m = file->msg; m->next; m = m->next)
            ;
        m->next = my_msg;
    }

    /* Remember values. */
    my_msg->msg_num = msg_num;
    my_msg->bytes_to_msg = bytes_to_msg;
    my_msg->bytes_in_msg = bytes_in_msg;
    my_msg->file = file;

    /* Read message metadata. We do this if we are opening a GRIB2
     * data file, but not if we are reading a GRIB2 index file. */
    if (read_file)
        if ((ret = read_msg_metadata(my_msg)))
            return ret;

    /* Increment number of messages in the file. */
    my_msg->file->num_messages++;

    /* Return pointer to caller, if desired. */
    if (msg)
        *msg = my_msg;

    return G2C_NOERROR;
}

/** Read metadata from a GRIB2 file being opened with g2c_open().
 *
 * @param g2cid The indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EBADID g2cid not found.
 * - ::G2C_EFILE File error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOMEM Out of memory.
 * - ::G2C_ENOMSG No GRIB message found.
 *
 * @author Ed Hartnett @date Aug 22, 2022
 */
static int
read_metadata(int g2cid)
{
    size_t msg_num;
    size_t file_pos = 0;
    size_t bytes_to_msg, bytes_in_msg;
    int ret = G2C_NOERROR;

    /* Find the open file struct. */
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    LOG((4, "read_metadata g2cid %d", g2cid));

    /* Read each message in the file. */
    for (msg_num = 0; !ret; msg_num++)
    {
        /* Find the message. */
        if ((ret = g2c_seekmsg(g2cid, file_pos, &bytes_to_msg, &bytes_in_msg)))
            return ret;
        LOG((5, "msg_num %d bytes_to_msg %ld bytes_in_msg %ld", msg_num, bytes_to_msg,
             bytes_in_msg));

        /*  When there are 0 bytes_in_msg, we are done. */
        if (!bytes_in_msg)
            break;

        /* Add new message to our list of messages. */
        if ((ret = add_msg(&g2c_file[g2cid], msg_num, bytes_to_msg, bytes_in_msg,
                           1, NULL)))
            return ret;

        /* Move the file position to the end of this message, ready to
         * scan for the next message. */
        file_pos = bytes_to_msg + bytes_in_msg;
        LOG((6, "file_pos %ld", file_pos));
    }

    /* If we read some messages, then run out, that's success. */
    if (ret == G2C_ENOMSG && msg_num)
        ret = G2C_NOERROR;

#ifdef LOGGING
    /* Print the file contents for library debugging. */
    g2c_log_file(g2cid);
#endif

    return ret;
}

/**
 * Open a GRIB2 file and add it to the list of open files.
 *
 * @param path Path of the file.
 * @param mode Open mode flags.
 * @param g2cid Pointer that gets an indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EINVAL - Invalid input.
 * - ::G2C_ETOOMANYFILES - Trying to open too many files at the same time.
 *
 * @author Ed Hartnett @date Aug 16, 2022
 */
int
g2c_add_file(const char *path, int mode, int *g2cid)
{
    int ret;

    /* Check inputs. */
    if (strlen(path) > G2C_MAX_NAME)
        return G2C_ENAMETOOLONG;
    if (!g2cid)
        return G2C_EINVAL;

    LOG((3, "g2c_add_file path %s mode %d", path, mode));

    /* Find a file ID. */
    if ((ret = find_available_g2cid(g2cid)))
        return ret;

    /* Open the file. */
    if (!(g2c_file[*g2cid].f = fopen(path, (mode & G2C_WRITE ? "rb+" : "rb"))))
        return G2C_EFILE;

    /* Copy the path. */
    strncpy(g2c_file[*g2cid].path, path, G2C_MAX_NAME);

    /* Remember the id. */
    g2c_file[*g2cid].g2cid = *g2cid;

    /* Initialize other values in struct. */
    g2c_file[*g2cid].msg = NULL;
    g2c_file[*g2cid].num_messages = 0;

    return G2C_NOERROR;
}

/** Open an existing GRIB2 file.
 *
 * This function opens the GRIB2 file and reads its metadata.
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
 * @param path Path of the file.
 * @param mode Open mode flags.
 * @param g2cid Pointer that gets an indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EINVAL - Invalid input.
 * - ::G2C_ETOOMANYFILES - Trying to open too many files at the same time.
 *
 * @author Ed Hartnett @date Aug 16, 2022
 */
int
g2c_open(const char *path, int mode, int *g2cid)
{
    int ret;

    LOG((2, "g2c_open path %s mode %d", path, mode));

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Open the file and add it to the list of open files. */
    ret = g2c_add_file(path, mode, g2cid);

    /* Read the metadata. */
    if (!ret)
        ret = read_metadata(*g2cid);

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}

#if 0
/** Create a new GRIB2 file.
 *
 * @param path Path of the file.
 * @param cmode Open mode flags.
 * @param g2cid Pointer that gets an indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EFILE - File exsists and NOCLOBBER, or error opening file.
 *
 * @author Ed Hartnett @date Aug 16, 2022
 */
int
g2c_create(const char *path, int cmode, int *g2cid)
{
    int my_g2cid;
    int ret;

    /* Check inputs. */
    if (strlen(path) > G2C_MAX_NAME)
        return G2C_ENAMETOOLONG;
    if (!g2cid)
        return G2C_EINVAL;

    LOG((2, "g2c_create path %s cmode %d", path, cmode));

    /* If NOCLOBBER, check if file exists. */
    if (cmode & G2C_NOCLOBBER)
    {
        FILE *f;
        if ((f = fopen(path, "r")))
        {
            fclose(f);
            return G2C_EFILE;
        }
    }

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    /* Find a file ID. */
    ret = find_available_g2cid(&my_g2cid);

    /* Create the file. */
    if (!ret)
        if (!(g2c_file[my_g2cid].f = fopen(path, "bw+")))
            ret = G2C_EFILE;

    if (!ret)
    {
        /* Read the metadata. */

        /* Copy the path. */
        strncpy(g2c_file[my_g2cid].path, path, G2C_MAX_NAME);

        /* Remember the id. */
        g2c_file[my_g2cid].g2cid = my_g2cid;

        /* Pass id back to user. */
        *g2cid = my_g2cid;
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return G2C_NOERROR;
}
#endif

/** Free resources holding the file metadata.
 *
 * @param g2cid Indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EBADID - Bad file ID.
 *
 * @author Ed Hartnett @date Aug 16, 2022
 */
static int
free_metadata(int g2cid)
{
    G2C_MESSAGE_INFO_T *msg;

    /* Check input. */
    if (g2cid > G2C_MAX_FILES)
        return G2C_EBADID;
    if (g2c_file[g2cid].g2cid != g2cid)
        return G2C_EBADID;

    /* Free message resources. */
    msg = g2c_file[g2cid].msg;
    while (msg)
    {
        G2C_MESSAGE_INFO_T *mtmp;
        G2C_SECTION_INFO_T *sec;

        /* Free section metadata. */
        sec = msg->sec;
        while (sec)
        {
            G2C_SECTION_INFO_T *stmp;

            stmp = sec->next;
            if (sec->template)
                free(sec->template);
            /* Free dim info in section 3. */
            if (sec->sec_num == 3)
            {
                LOG((5, "free_metadata freeing storage for lat/lon values"));
                float *v0 = ((G2C_SECTION3_INFO_T *)(sec->sec_info))->dim[0].value;
                float *v1 = ((G2C_SECTION3_INFO_T *)(sec->sec_info))->dim[1].value;
                if (v0)
                    free(v0);
                if (v1)
                    free(v1);
            }
            if (sec->sec_info)
                free(sec->sec_info);
            free(sec);
            sec = stmp;
        }

        /* Free message. */
        mtmp = msg->next;
        free(msg);
        msg = mtmp;
    }

    return G2C_NOERROR;
}

/** Close a GRIB2 file, freeing resources.
 *
 * @param g2cid Indentifier for the file.
 *
 * @return
 * - ::G2C_NOERROR - No error.
 * - ::G2C_EBADID - Bad file ID.
 *
 * @author Ed Hartnett @date Aug 16, 2022
 */
int
g2c_close(int g2cid)
{
    int ret = G2C_NOERROR;

    /* Check input. */
    if (g2cid > G2C_MAX_FILES)
        return G2C_EBADID;

    /* If using threading, lock the mutex. */
    MUTEX_LOCK(m);

    if (g2c_file[g2cid].g2cid != g2cid)
        ret = G2C_EBADID;

    LOG((2, "g2c_close %d", g2cid));

    /* Free resources. */
    if (!ret)
        ret = free_metadata(g2cid);

    /* Close the file. */
    if (!ret)
        if (fclose(g2c_file[g2cid].f))
            ret = G2C_EFILE;

    /* Reset the file data. */
    if (!ret)
    {
        g2c_file[g2cid].path[0] = 0;
        g2c_file[g2cid].g2cid = 0;
        g2c_file[g2cid].num_messages = 0;
        g2c_file[g2cid].f = NULL;
    }

    /* If using threading, unlock the mutex. */
    MUTEX_UNLOCK(m);

    return ret;
}
