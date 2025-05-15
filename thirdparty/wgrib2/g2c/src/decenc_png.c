/** @file
 * @brief Decode/encode a PNG stream.
 * @author Alyson Stahl @date 2024-13-08
 */

#include "grib2_int.h"
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Struct for PNG stream.
 */
struct png_stream
{
    unsigned char *stream_ptr; /**< Location to write PNG stream. */
    g2int stream_len;          /**< Number of bytes written. */
};

typedef struct png_stream png_stream; /**< Typedef for PNG stream. */

void user_read_data(png_structp, png_bytep, png_uint_32);
void user_write_data(png_structp, png_bytep, png_uint_32);
void user_flush_data(png_structp);

/**
 * Custom read function used so that libpng will read a PNG stream
 * from memory instead of a file on disk.
 *
 * @param png_ptr Pointer to PNG.
 * @param data Pointer to data.
 * @param length Length.
 *
 * @author Stephen Gilbert
 */
void
user_read_data(png_structp png_ptr, png_bytep data, png_uint_32 length)
{
    char *ptr;
    g2int offset;
    png_stream *mem;

    mem = (png_stream *)png_get_io_ptr(png_ptr);
    ptr = (void *)mem->stream_ptr;
    offset = mem->stream_len;
    memcpy(data, ptr + offset, length);
    mem->stream_len += length;
}

/**
 * Custom write function used to that libpng will write to memory
 * location instead of a file on disk.
 *
 * @param png_ptr pointer
 * @param data data
 * @param length length
 *
 * @author Stephen Gilbert
 */
void
user_write_data(png_structp png_ptr, png_bytep data, png_uint_32 length)
{
    unsigned char *ptr;
    g2int offset;
    png_stream *mem;

    mem = (png_stream *)png_get_io_ptr(png_ptr);
    ptr = mem->stream_ptr;
    offset = mem->stream_len;
    memcpy(ptr + offset, data, length);
    mem->stream_len += length;
}

/**
 * Dummy Custom flush function.
 *
 * @param png_ptr Pointer to PNG struct.
 *
 * @author Stephen Gilbert
 */
void
user_flush_data(png_structp png_ptr)
{
}

/**
 * Decode PNG.
 *
 * @param pngbuf Pointer to PNG buffer.
 * @param width Pointer to width.
 * @param height Pointer to height.
 * @param cout Output buffer.
 *
 * @return 0 for success, error code otherwise.
 *
 * @author Alyson Stahl
 */
int
g2c_dec_png(unsigned char *pngbuf, int *width, int *height,
            unsigned char *cout)
{
    g2int width8 = *width, height8 = *height;
    int ret;

    ret = dec_png(pngbuf, &width8, &height8, cout);

    *width = (g2int)width8;
    *height = (g2int)height8;

    return ret;
}

/**
 * Decode PNG.
 *
 * @param pngbuf Pointer to PNG buffer.
 * @param width Pointer to width.
 * @param height Pointer to height.
 * @param cout Output buffer.
 *
 * @return 0 for success, error code otherwise.
 *
 * @author Stephen Gilbert
 */
int
dec_png(unsigned char *pngbuf, g2int *width, g2int *height,
        unsigned char *cout)
{
    int interlace, color, compres, filter, bit_depth;
    g2int j, k, n, bytes;
    png_structp png_ptr;
    png_infop info_ptr, end_info;
    png_bytepp row_pointers;
    png_stream read_io_ptr;
    png_uint_32 h32, w32;

    /* Check if stream is a valid PNG format. */
    if (png_sig_cmp(pngbuf, 0, 8) != 0)
        return -3;

    /* Create and initialize png_structs. */
    if (!(png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, (png_voidp)NULL,
                                           NULL, NULL)))
        return -1;

    if (!(info_ptr = png_create_info_struct(png_ptr)))
    {
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        return -2;
    }

    if (!(end_info = png_create_info_struct(png_ptr)))
    {
        png_destroy_read_struct(&png_ptr, (png_infopp)info_ptr, (png_infopp)NULL);
        return -2;
    }

    /* Set Error callback. */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        return -3;
    }

    /* Initialize info for reading PNG stream from memory. */
    read_io_ptr.stream_ptr = (png_voidp)pngbuf;
    read_io_ptr.stream_len = 0;

    /* Set new custom read function. */
    png_set_read_fn(png_ptr, (png_voidp)&read_io_ptr, (png_rw_ptr)user_read_data);

    /*     support for larger grids   */
    png_set_user_limits(png_ptr, G2C_PNG_WIDTH_MAX, G2C_PNG_HEIGHT_MAX);

    /* Read and decode PNG stream. */
    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* Get pointer to each row of image data. */
    row_pointers = png_get_rows(png_ptr, info_ptr);

    /* Get image info, such as size, depth, colortype, etc... */
    /*printf("SAGT:png %d %d %d\n", info_ptr->width, info_ptr->height, info_ptr->bit_depth);*/
    // (void)png_get_IHDR(png_ptr,  info_ptr,  (png_uint_32 *)width,  (png_uint_32 *)height,
    (void)png_get_IHDR(png_ptr, info_ptr, &w32, &h32,
                       &bit_depth, &color, &interlace, &compres, &filter);

    *height = h32;
    *width = w32;

    /* Check if image was grayscale. */
    /*
      if (color != PNG_COLOR_TYPE_GRAY) {
      fprintf(stderr, "dec_png: Grayscale image was expected. \n");
      }
    */
    if (color == PNG_COLOR_TYPE_RGB)
        bit_depth = 24;
    else if (color == PNG_COLOR_TYPE_RGB_ALPHA)
        bit_depth = 32;

    /* Copy image data to output string   */
    n = 0;
    bytes = (*width * bit_depth) / 8;
    if ((*width * bit_depth) % 8 != 0)
    {
        bytes++;
    }
    for (j = 0; j < *height; j++)
    {
        for (k = 0; k < bytes; k++)
        {
            cout[n] = *(row_pointers[j] + k);
            n++;
        }
    }

    /* Clean up. */
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    return 0;
}

/**
 * Encode PNG.
 *
 * @param data data.
 * @param width width.
 * @param height height.
 * @param nbits number of bits.
 * @param pngbuf PNG buffer.
 *
 * @return PNG length, or negative number for error.
 *
 * @author Alyson Stahl
 */
int
g2c_enc_png(unsigned char *data, int width, int height, int nbits,
            unsigned char *pngbuf)
{
    g2int width8 = width, height8 = height, nbits8 = nbits;

    return enc_png(data, width8, height8, nbits8, pngbuf);
}

/**
 * Encode PNG.
 *
 * @param data data.
 * @param width width.
 * @param height height.
 * @param nbits number of bits.
 * @param pngbuf PNG buffer.
 *
 * @return PNG length, or negative number for error.
 *
 * @author Stephen Gilbert
 */
int
enc_png(unsigned char *data, g2int width, g2int height, g2int nbits,
        unsigned char *pngbuf)
{
    int color_type;
    g2int j, bytes, pnglen, bit_depth;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep **row_pointers;
    png_stream write_io_ptr;

    /* Create and initialize png_structs. */
    if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL)))
        return -1;

    if (!(info_ptr = png_create_info_struct(png_ptr)))
    {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        return -2;
    }

    /* Set Error callback. */
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return -3;
    }

    /* Initialize info for writing PNG stream to memory. */
    write_io_ptr.stream_ptr = (png_voidp)pngbuf;
    write_io_ptr.stream_len = 0;

    /* Set new custom write functions. */
    png_set_write_fn(png_ptr, (png_voidp)&write_io_ptr, (png_rw_ptr)user_write_data,
                     (png_flush_ptr)user_flush_data);

    /* Set the image size, colortype, filter type, etc... */
    bit_depth = nbits;
    color_type = PNG_COLOR_TYPE_GRAY;
    if (nbits == 24)
    {
        bit_depth = 8;
        color_type = PNG_COLOR_TYPE_RGB;
    }
    else if (nbits == 32)
    {
        bit_depth = 8;
        color_type = PNG_COLOR_TYPE_RGB_ALPHA;
    }
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    /* Put image data into the PNG info structure. */
    bytes = (width * nbits) / 8;
    if ((width * nbits) % 8 != 0)
    {
        bytes++;
    }

    row_pointers = malloc(height * sizeof(png_bytep));
    for (j = 0; j < height; j++)
        row_pointers[j] = (png_bytep *)(data + (j * bytes));

    png_set_rows(png_ptr, info_ptr, (png_bytepp)row_pointers);

    /* Do the PNG encoding, and write out PNG stream. */
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* Clean up. */
    png_destroy_write_struct(&png_ptr, &info_ptr);
    free(row_pointers);
    pnglen = write_io_ptr.stream_len;
    return pnglen;
}
