/** @file
 * @brief Decode/encode a JPEG2000 code stream.
 * @author Alyson Stahl @date 2024-14-08
 */

#include "grib2_int.h"
#include "jasper/jasper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXOPTSSIZE 1024 /**< Maximum size of options. */

/**
 * Encode a grayscale image into a JPEG2000 code stream
 * specified in the JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1)
 * using [JasPer Software](https://github.com/jasper-software/jasper).
 *
 * @param cin Packed matrix of Grayscale image values to encode.
 * @param width width of image.
 * @param height height of image.
 * @param nbits depth (in bits) of image.  i.e number of bits used to
 * hold each data value.
 * @param ltype indicator of lossless or lossy compression.
 * - 1, for lossy compression
 * - != 1, for lossless compression
 * @param ratio target compression ratio. (ratio:1) Used only when
 * ltype == 1.
 * @param retry If 1 try increasing number of guard bits.
 * @param outjpc Output encoded JPEG2000 code stream.
 * @param jpclen Number of bytes allocated for the output JPEG2000
 * code stream in outjpc.
 *
 * @return
 * - > 0 = Length in bytes of encoded JPEG2000 code stream
 * - ::G2_JASPER_INIT Error initializing jasper library.
 * - ::G2_JASPER_ENCODE Error encode jpeg2000 code stream.
 *
 * @note Requires JasPer Software version 1.500.4 or 1.700.2 or later.
 *
 * @author Stephen Gilbert @date 2002-12-02
 * @author Ed Hartnett
 */
int
g2c_enc_jpeg2000(unsigned char *cin, int width, int height, int nbits,
                 int ltype, int ratio, int retry, char *outjpc,
                 size_t jpclen)
{
    g2int width8 = width, height8 = height, nbits8 = nbits, ltype8 = ltype;
    g2int ratio8 = ratio, retry8 = retry, jpclen8 = jpclen;

    return enc_jpeg2000(cin, width8, height8, nbits8, ltype8, ratio8, retry8,
                        outjpc, jpclen8);
}

/**
 * Encode a grayscale image into a JPEG2000 code stream
 * specified in the JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1)
 * using [JasPer Software](https://github.com/jasper-software/jasper).
 *
 *  * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-12-02 | Gilbert | Initial
 * 2004-12-16 | Gilbert | Added retry argument allowing increased guard bits.
 * 2022-04-15 | Hartnett | Converted to use jas_ instead of jpc_ functions.
 *
 * @param cin Packed matrix of Grayscale image values to encode.
 * @param width width of image.
 * @param height height of image.
 * @param nbits depth (in bits) of image.  i.e number of bits used to
 * hold each data value.
 * @param ltype indicator of lossless or lossy compression.
 * - 1, for lossy compression
 * - != 1, for lossless compression
 * @param ratio target compression ratio. (ratio:1) Used only when
 * ltype == 1.
 * @param retry If 1 try increasing number of guard bits.
 * @param outjpc Output encoded JPEG2000 code stream.
 * @param jpclen Number of bytes allocated for the output JPEG2000
 * code stream in outjpc.
 *
 * @return
 * - > 0 = Length in bytes of encoded JPEG2000 code stream
 * - ::G2_JASPER_INIT Error initializing jasper library.
 * - ::G2_JASPER_ENCODE Error encode jpeg2000 code stream.
 *
 * @note Requires JasPer Software version 1.500.4 or 1.700.2 or later.
 *
 * @author Stephen Gilbert @date 2002-12-02
 * @author Ed Hartnett
 */
int
enc_jpeg2000(unsigned char *cin, g2int width, g2int height, g2int nbits,
             g2int ltype, g2int ratio, g2int retry, char *outjpc,
             g2int jpclen)
{
    int ier, rwcnt;
    jas_image_t image;
    jas_stream_t *jpcstream, *istream;
    jas_image_cmpt_t cmpt, *pcmpt;
    char opts[MAXOPTSSIZE];
    int fmt;

    LOG((3, "enc_jpeg2000 width %ld height %ld nbits %ld ltype %ld ratio %ld retry %ld jpclen %d",
         width, height, nbits, ltype, ratio, retry, jpclen));

    /* Set lossy compression options, if requested. */
    if (ltype != 1)
        opts[0] = (char)0;
    else
        snprintf(opts, MAXOPTSSIZE, "mode=real\nrate=%f", 1.0 / (float)ratio);

    if (retry == 1) /* option to increase number of guard bits */
        strcat(opts, "\nnumgbits=4");

    /* Initialize the JasPer image structure describing the grayscale
     * image to encode into the JPEG2000 code stream. */
    image.tlx_ = 0;
    image.tly_ = 0;
    image.brx_ = (jas_image_coord_t)width;
    image.bry_ = (jas_image_coord_t)height;
    image.numcmpts_ = 1;
    image.maxcmpts_ = 1;
    image.clrspc_ = JAS_CLRSPC_SGRAY; /* grayscale Image */
    image.cmprof_ = 0;

    cmpt.tlx_ = 0;
    cmpt.tly_ = 0;
    cmpt.hstep_ = 1;
    cmpt.vstep_ = 1;
    cmpt.width_ = (jas_image_coord_t)width;
    cmpt.height_ = (jas_image_coord_t)height;
    cmpt.type_ = JAS_IMAGE_CT_COLOR(JAS_CLRSPC_CHANIND_GRAY_Y);
    cmpt.prec_ = nbits;
    cmpt.sgnd_ = 0;
    cmpt.cps_ = (nbits + 7) / 8;

    pcmpt = &cmpt;
    image.cmpts_ = &pcmpt;

    /* Initialize Jasper. */
#ifdef JASPER3
#define HUNDRED_MB 100000000
    jas_conf_clear();
    /* static jas_std_allocator_t allocator; */
    /* jas_std_allocator_init(&allocator); */
    /* jas_conf_set_allocator(JAS_CAST(jas_std_allocator_t *, &allocator)); */
    jas_conf_set_max_mem_usage(HUNDRED_MB);
    jas_conf_set_multithread(true);
    if (jas_init_library())
        return G2_JASPER_INIT;
    if (jas_init_thread())
        return G2_JASPER_INIT;
#else
    if (jas_init())
        return G2_JASPER_INIT;
#endif /* JASPER3 */

    /* Open a JasPer stream containing the input grayscale values. */
    istream = jas_stream_memopen((char *)cin, height * width * cmpt.cps_);
    cmpt.stream_ = istream;

    /* Open an output stream that will contain the encoded jpeg2000
     * code stream. */
    jpcstream = jas_stream_memopen(outjpc, (int)jpclen);

    /* Get jasper ID of JPEG encoder. */
    fmt = jas_image_strtofmt(G2C_JASPER_JPEG_FORMAT_NAME);

    /* Encode image. */
    if ((ier = jas_image_encode(&image, jpcstream, fmt, opts)))
        return G2_JASPER_ENCODE;

    /* Rememeber the length in bytes of the encoded JPEG code
     * stream. */
    rwcnt = jpcstream->rwcnt_;

    /* Clean up JasPer work structures. */
    ier = jas_stream_close(istream);
    ier = jas_stream_close(jpcstream);

    /* Finalize jasper. */
#ifdef JASPER3
    jas_cleanup_thread();
    jas_cleanup_library();
#else
    jas_cleanup();
#endif /* JASPER3 */

    /* Return size of jpeg2000 code stream. */
    return rwcnt;
}

/**
 * Decode a JPEG2000 code stream specified in the
 * JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) using [JasPer
 * Software](https://github.com/jasper-software/jasper).
 *
 * @param injpc Pointer to buffer that holds the input JPEG2000 code
 * stream.
 * @param bufsize Length (in bytes) of the buffer that holds the input
 * JPEG2000 code stream.
 * @param outfld Pointer to either int or g2int array, already
 * allocated, that gets the unpacked data.
 * @param out_is_g2int Non-zero if the output array is of type g2int
 * (i.e. 64-bit ints), zero if output is an int array (32-bits).
 *
 * @return
 * - 0 Successful decode
 * - ::G2_JASPER_DECODE Error decode jpeg2000 code stream.
 * - ::G2_JASPER_DECODE_COLOR decoded image had multiple color
 *     components. Only grayscale is expected.
 * - ::G2_JASPER_INIT Error inializing Jasper library.
 *
 * @author Stephen Gilbert @date 2002-12-02
 * @author Ed Hartnett
 * @author Eric Engle
 */
static int
int_dec_jpeg2000(char *injpc, g2int bufsize, void *outfld, int out_is_g2int)
{
    g2int i, j, k;
    jas_image_t *image = NULL;
    jas_stream_t *jpcstream;
    jas_image_cmpt_t *pcmpt;
    char *opts = NULL;
    jas_matrix_t *data;
    int fmt;

    LOG((3, "int_dec_jpeg2000 bufsize %ld out_is_g2int %d", bufsize, out_is_g2int));

    /* Initialize Jasper. */
#ifdef JASPER3
    jas_conf_clear();
    /* static jas_std_allocator_t allocator; */
    /* jas_std_allocator_init(&allocator); */
    /* jas_conf_set_allocator(JAS_CAST(jas_std_allocator_t *, &allocator)); */
    jas_conf_set_max_mem_usage(G2C_JASPER_MAX_MEM);
    jas_conf_set_multithread(true);
    if (jas_init_library())
        return G2_JASPER_INIT;
    if (jas_init_thread())
        return G2_JASPER_INIT;
#else
    if (jas_init())
        return G2_JASPER_INIT;
#endif /* JASPER3 */

    /* Create jas_stream_t containing input JPEG200 codestream in
     * memory. */
    jpcstream = jas_stream_memopen(injpc, bufsize);

    /* Get jasper ID of JPEG encoder. */
    fmt = jas_image_strtofmt(G2C_JASPER_JPEG_FORMAT_NAME);

    /* Decode JPEG200 codestream into jas_image_t structure. */
    if (!(image = jas_image_decode(jpcstream, fmt, opts)))
        return G2_JASPER_DECODE;

    pcmpt = image->cmpts_[0];
    /*
      printf(" SAGOUT DECODE:\n");
      printf(" tlx %d \n", image->tlx_);
      printf(" tly %d \n", image->tly_);
      printf(" brx %d \n", image->brx_);
      printf(" bry %d \n", image->bry_);
      printf(" numcmpts %d \n", image->numcmpts_);
      printf(" maxcmpts %d \n", image->maxcmpts_);
      printf(" colorspace %d \n", image->clrspc_);
      printf(" inmem %d \n", image->inmem_);
      printf(" COMPONENT:\n");
      printf(" tlx %d \n", pcmpt->tlx_);
      printf(" tly %d \n", pcmpt->tly_);
      printf(" hstep %d \n", pcmpt->hstep_);
      printf(" vstep %d \n", pcmpt->vstep_);
      printf(" width %d \n", pcmpt->width_);
      printf(" height %d \n", pcmpt->height_);
      printf(" prec %d \n", pcmpt->prec_);
      printf(" sgnd %d \n", pcmpt->sgnd_);
      printf(" cps %d \n", pcmpt->cps_);
      printf(" type %d \n", pcmpt->type_);
    */

    /* Expecting jpeg2000 image to be grayscale only. No color components. */
    if (image->numcmpts_ != 1)
        return G2_JASPER_DECODE_COLOR;

    /* Create a data matrix of grayscale image values decoded from the
     * jpeg2000 codestream. */
    data = jas_matrix_create(jas_image_height(image), jas_image_width(image));
    jas_image_readcmpt(image, 0, 0, 0, jas_image_width(image),
                       jas_image_height(image), data);

    LOG((3, "pcmpt->height_ %d pcmpt->width_ %d", pcmpt->height_, pcmpt->width_));

    /* Copy data matrix to output integer array. */
    k = 0;
    if (out_is_g2int)
    {
        for (i = 0; i < pcmpt->height_; i++)
            for (j = 0; j < pcmpt->width_; j++)
                ((g2int *)outfld)[k++] = data->rows_[i][j];
    }
    else
    {
        for (i = 0; i < pcmpt->height_; i++)
            for (j = 0; j < pcmpt->width_; j++)
                ((int *)outfld)[k++] = data->rows_[i][j];
    }

    /* Clean up JasPer work structures. */
    jas_matrix_destroy(data);
    jas_stream_close(jpcstream);
    jas_image_destroy(image);

    /* Finalize jasper. */
#ifdef JASPER3
    jas_cleanup_thread();
    jas_cleanup_library();
#else
    jas_cleanup();
#endif /* JASPER3 */

    return 0;
}

/**
 * Decode a JPEG2000 code stream specified in the
 * JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) using [JasPer
 * Software](https://github.com/jasper-software/jasper).
 *
 * @param injpc Pointer to buffer that holds the input JPEG2000 code
 * stream.
 * @param bufsize Length (in bytes) of the buffer that holds the input
 * JPEG2000 code stream.
 * @param outfld Pointer to int array, already allocated, that gets
 * the unpacked data.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2_JASPER_DECODE Error decode jpeg2000 code stream.
 * - ::G2_JASPER_DECODE_COLOR decoded image had multiple color
 *     components. Only grayscale is expected.
 * - ::G2_JASPER_INIT Error inializing Jasper library.
 *
 * @author Ed Hartnett @date 9/7/22
 */
int
g2c_dec_jpeg2000(char *injpc, size_t bufsize, int *outfld)
{
    return int_dec_jpeg2000(injpc, bufsize, outfld, 0);
}

/**
 * Decode a JPEG2000 code stream specified in the
 * JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) using [JasPer
 * Software](https://github.com/jasper-software/jasper).
 *
 * @param injpc Pointer to buffer that holds the input JPEG2000 code
 * stream.
 * @param bufsize Length (in bytes) of the buffer that holds the input
 * JPEG2000 code stream.
 * @param outfld Pointer to g2int array, already allocated, that gets
 * the unpacked data.
 *
 * @return
 * - 0 Successful decode
 * - ::G2_JASPER_DECODE Error decode jpeg2000 code stream.
 * - ::G2_JASPER_DECODE_COLOR decoded image had multiple color
 *     components. Only grayscale is expected.
 * - ::G2_JASPER_INIT Error inializing Jasper library.
 *
 * @author Stephen Gilbert, Ed Hartnett
 */
int
dec_jpeg2000(char *injpc, g2int bufsize, g2int *outfld)
{
    return int_dec_jpeg2000(injpc, bufsize, outfld, 1);
}
