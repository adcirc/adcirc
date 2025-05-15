/** @file
 * @brief Decode/encode an AEC code stream.
 * @author Eric Engle @date 2023-10-16
 */

#include "grib2_int.h"
#include <libaec.h>
#include <stdint.h>
#include <stdio.h>

/**
 * Decode an AEC code stream specified in the [CCSDS 121.0-B-3 Blue
 * Book](https://public.ccsds.org/Pubs/121x0b3.pdf).
 *
 * @param cpack Pointer to buffer that holds the input AEC code
 * stream.
 * @param len Length (in bytes) of the buffer that holds the input
 * AEC code stream.
 * @param nbits CCSDS bits per sample.
 * @param flags CCSDS compression options mask.
 * @param block_size CCSDS block size.
 * @param rsi CCSDS reference sample interval.
 * @param cfld Pointer to output buffer from the AEC decoder.
 * @param cfldlen length of output buffer.
 *
 * @return
 * - >0 Length of data from AEC decoder
 * - 0 Successful decode (AEC_OK)
 * - -1 AEC_CONF_ERROR
 * - -2 AEC_STREAM_ERROR
 * - -3 AEC_DATA_ERROR
 * - -4 AEC_MEM_ERROR
 * - -5 AEC_RSI_OFFSETS_ERROR
 *
 * @author Eric Engle @date 2023-10-16
 */
int
dec_aec(unsigned char *cpack, g2int len, g2int nbits, g2int flags,
        g2int block_size, g2int rsi, unsigned char *cfld, g2int cfldlen)
{
    struct aec_stream strm;
    int ret = 0;

    /* Define bits per sample */
    strm.bits_per_sample = nbits;
    LOG((3, "dec_aec: bits_per_sample = %d", strm.bits_per_sample));

    /* Define a block size */
    LOG((3, "dec_aec: block_size = %d", block_size));
    strm.block_size = block_size;

    /* Define the reference sample interval */
    LOG((3, "dec_aec: rsi = %d", rsi));
    strm.rsi = rsi;

    /* Define the AEC compression flags. */
    strm.flags = flags;
    LOG((3, "dec_aec: flags = %d", strm.flags));

    /* Pointer to input */
    strm.next_in = cpack;

    /* Length of input in bytes */
    strm.avail_in = len;
    LOG((3, "dec_aec: avail_in = %d", strm.avail_in));

    /* Pointer to output buffer */
    strm.next_out = cfld;

    /* Length of output buffer in bytes */
    strm.avail_out = (size_t)cfldlen;

    /* Decode. */
    ret = aec_buffer_decode(&strm);
    LOG((3, "dec_aec: return from aec_buffer_encode = %d", ret));
    if (ret == AEC_OK)
        ret = strm.total_out;

    return ret;
}

/**
 * Encode data into an AEC code stream specified in the
 * [CCSDS 121.0-B-3 Blue Book](https://public.ccsds.org/Pubs/121x0b3.pdf).
 *
 * @param data Pointer to buffer that holds the input data.
 * @param ctemplen Length (in bytes) of the buffer that holds
 * the input data..
 * @param nbits CCSDS bits per sample.
 * @param flags CCSDS compression options mask.
 * @param block_size CCSDS block size.
 * @param rsi CCSDS reference sample interval.
 * @param aecbuf Pointer to buffer holding the AEC encoded stream.
 * @param aecbuflen Length of AEC code stream.
 *
 * @return
 * - >0 Exact length of AEC encoded data.
 * - 0 Successful decode (AEC_OK)
 * - -1 AEC_CONF_ERROR
 * - -2 AEC_STREAM_ERROR
 * - -3 AEC_DATA_ERROR
 * - -4 AEC_MEM_ERROR
 * - -5 AEC_RSI_OFFSETS_ERROR
 *
 * @author Eric Engle @date 2023-10-16
 */
int
enc_aec(unsigned char *data, g2int ctemplen, g2int nbits, g2int flags,
        g2int block_size, g2int rsi, unsigned char *aecbuf, g2int *aecbuflen)
{
    struct aec_stream strm;
    int ret = 0;

    /* Define bits per sample */
    strm.bits_per_sample = nbits;
    LOG((3, "enc_aec: bits_per_sample = %d", strm.bits_per_sample));

    /* Define a block size */
    LOG((3, "enc_aec: block_size = %d", block_size));
    strm.block_size = block_size;

    /* Define the reference sample interval */
    LOG((3, "enc_aec: rsi = %d", rsi));
    strm.rsi = rsi;

    /* Define the AEC compression flags. */
    strm.flags = flags;
    LOG((3, "enc_aec: flags = %d", strm.flags));

    /* Pointer to input */
    strm.next_in = data;

    /* Length of input in bytes */
    strm.avail_in = ctemplen;
    LOG((3, "enc_aec: avail_in = %d", strm.avail_in));

    /* Pointer to output buffer */
    strm.next_out = aecbuf;

    /* Length of output buffer in bytes */
    strm.avail_out = (size_t)aecbuflen;

    /* Encode into AEC. */
    ret = aec_buffer_encode(&strm);
    LOG((3, "enc_aec: return from aec_buffer_encode = %d", ret));
    if (ret == AEC_OK)
        ret = strm.total_out;

    return ret;
}
