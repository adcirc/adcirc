/* unpack_grib 
 * 3/2008 public domain Wesley Ebisuzaki
 * 5/2016 public domain DWD
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>

#include "wgrib2.h"
#include "grb2.h"

#ifdef USE_G2CLIB_LOW
    #include <grib2.h>
#endif

#ifdef USE_AEC
   #include <libaec.h>
#endif

/*
 * unpack grib -- only some formats (code table 5.0) are supported
 *
 * supported: 0 (simple), 4 (ieee), 40 (jpeg), 41(png), 42(aec)
 *
 * input:  sec[]
 *         float data[npnts]
 *
 */

int unpk_grib(unsigned char **sec, float *data) {

    int packing, bitmap_flag, nbits;
    unsigned int ndata, ii;
    unsigned char *mask_pointer, mask;
    unsigned char *ieee, *p;
    float tmp;
    // float reference, tmp;
    double reference;
    double bin_scale, dec_scale, b;

#if G2_PNG_ENABLED == 1
    int i, width, height;
#endif

#if G2_JPEG2000_ENABLED == 1
    int *ifld, err;
    unsigned int kk;
#endif

#if (G2_PNG_ENABLED == 1 || defined USE_AEC)
    unsigned char *c;
#endif

#ifdef USE_AEC
    struct aec_stream strm;
    int status;
    int numBitsNeeded;
    size_t size;
#endif

    packing = code_table_5_0(sec);
    ndata = GB2_Sec3_npts(sec);
    bitmap_flag = code_table_6_0(sec);

    if (bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255)
         fatal_error("unknown bitmap", "");

    if (packing == 4) {			// ieee
        if (sec[5][11] != 1) fatal_error_i("unpk ieee grib file precision %d not supported", 
		(int) sec[5][11]);

        // ieee depacking -- simple no bitmap
        if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
            for (ii = 0; ii < ndata; ii++) {
                data[ii] = ieee2flt_nan(sec[7]+5+ii*4);
            }
	    return 0;
        }
        if (bitmap_flag == 0 || bitmap_flag == 254) {
	    mask_pointer = sec[6] + 6;
	    ieee = sec[7]+5;
	    mask = 0;
	    for (ii = 0; ii < ndata; ii++) {
		if ((ii & 7) == 0) mask = *mask_pointer++;
		if (mask & 128) {
		    data[ii] = ieee2flt_nan(ieee);
		    ieee += 4;
		}
		else {
		    data[ii] = UNDEFINED;
		}
		mask <<= 1;
	    }
	    return 0;
	}
        fatal_error("unknown bitmap", "");
    }
    else if (packing == 0 || packing == 61) {			// simple grib1 packing  61 -- log preprocessing

	p = sec[5];
	reference = ieee2flt(p+11);
	bin_scale = Int_Power(2.0, int2(p+15));
	dec_scale = Int_Power(10.0, -int2(p+17));
	nbits = p[19];
	b = 0.0;
        if (packing == 61) b = ieee2flt(p+20);

	if (bitmap_flag != 0 && bitmap_flag != 254 && bitmap_flag != 255)
            fatal_error("unknown bitmap", "");

        if (nbits == 0) {
            tmp = reference*dec_scale;
	    if (packing == 61) tmp = exp(tmp) - b;		// remove log prescaling
            if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
		mask = 0;
                for (ii = 0; ii < ndata; ii++) {
		    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ?  tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
        }

	mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

	unpk_0(data, sec[7]+5, mask_pointer, nbits, ndata, reference, 
		bin_scale,dec_scale);

	if (packing == 61) {		// remove log prescaling
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
            for (ii = 0; ii < ndata; ii++) {
                if (DEFINED_VAL(data[ii])) data[ii] = exp(data[ii]) - b;
            }
	}
	return 0;
    }

    else if (packing == 2 || packing == 3) {		// complex
	return unpk_complex(sec, data, ndata);
    }
    else if (packing == 200) {				// run length
	return unpk_run_length(sec, data, ndata);
    }

#if G2_JPEG2000_ENABLED == 1
    else if (packing == 40 ||  packing == 40000) {		// jpeg2000
	p = sec[5];
	reference = ieee2flt(p+11);
	bin_scale = Int_Power(2.0, int2(p+15));
	dec_scale = Int_Power(10.0, -int2(p+17));
	nbits = p[19];
	
	if (nbits == 0) {
	    tmp = reference*dec_scale;
            if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
		for (ii = 0; ii < ndata; ii++) {
		    data[ii] = tmp;
		}
		return 0;
	    }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                ieee = sec[7]+5;
                mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ? tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
	    }
            fatal_error("unknown bitmap", "");
        }

	// decode jpeg2000

    ifld = (int *) malloc(ndata * sizeof(int));
	if (ifld == 0) fatal_error("unpk: memory allocation error","");
	err = g2c_dec_jpeg2000((char *) sec[7]+5, (size_t) GB2_Sec7_size(sec)-5, ifld);
	if (err != 0) fatal_error_i("dec_jpeg2000, error %d",err);

    if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii)
#endif
        for (ii = 0; ii < ndata; ii++) {
            data[ii] = ((ifld[ii]*bin_scale)+reference)*dec_scale;
        }
            }
    else if (bitmap_flag == 0 || bitmap_flag == 254) {
        mask_pointer = sec[6] + 6;
        mask = 0;
	    kk = 0;
            for (ii = 0; ii < ndata; ii++) {
                if ((ii & 7) == 0) mask = *mask_pointer++;
                data[ii] = (mask & 128) ? ((ifld[kk++]*bin_scale)+reference)*dec_scale : UNDEFINED;
                mask <<= 1;
            }
    } else {
        fatal_error_i("unknown bitmap: %d", bitmap_flag);
    }

    free(ifld);
	return 0;
    }
#endif

#if G2_PNG_ENABLED == 1
    else if (packing == 41) {		// png
	p = sec[5];
	reference = ieee2flt(p+11);
	bin_scale = Int_Power(2.0, int2(p+15));
	dec_scale = Int_Power(10.0, -int2(p+17));
	nbits = p[19];

        if (nbits == 0) {
            tmp = reference*dec_scale;
            if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
                for (ii = 0; ii < ndata; ii++) {
                    data[ii] = tmp;
                }
                return 0;
            }
            if (bitmap_flag == 0 || bitmap_flag == 254) {
                mask_pointer = sec[6] + 6;
                ieee = sec[7]+5;
		mask = 0;
                for (ii = 0; ii < ndata; ii++) {
                    if ((ii & 7) == 0) mask = *mask_pointer++;
                    data[ii] = (mask & 128) ?  tmp : UNDEFINED;
                    mask <<= 1;
                }
                return 0;
            }
            fatal_error("unknown bitmap", "");
        }

	/* allocate max size of buffer needed */
	if ((c = (unsigned char *) malloc(4*sizeof(char) * (size_t) ndata)) == NULL)
            fatal_error("unpk: png decode allocation error", "");

    i = g2c_dec_png(sec[7]+5, &width, &height, (unsigned char *) c);
	if (i) fatal_error_i("unpk: png decode error %d",i);
	mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

//	check sizes

	if (mask_pointer == NULL) {
	    if (ndata != width*height) 
    		fatal_error_i("png size mismatch w*h=%d", width*height);
	}
	else {
	    if (ndata != width*height + missing_points(mask_pointer, GB2_Sec3_npts(sec)) )
                fatal_error("png size mismatch", "");
	}

	unpk_0(data, c, mask_pointer, nbits, ndata, reference, 
		bin_scale,dec_scale);
	free(c);
	return 0;
    }
#endif
#ifdef USE_AEC
    else if (packing == 42) {		// aec

    	p = sec[5];
    	reference = ieee2flt(p+11);
    	bin_scale = Int_Power(2.0, int2(p+15));
    	dec_scale = Int_Power(10.0, -int2(p+17));
    	nbits = p[19];

    	if (nbits == 0) {
    	    tmp = reference*dec_scale;
    	    if (bitmap_flag == 255) {
#ifdef USE_OPENMP
#pragma omp parallel for private(ii) schedule(static)
#endif
    		for (ii = 0; ii < ndata; ii++) {
    		    data[ii] = tmp;
    		}
    		return 0;
    	    }
    	    if (bitmap_flag == 0 || bitmap_flag == 254) {
    		mask_pointer = sec[6] + 6;
    		mask = 0;
    		for (ii = 0; ii < ndata; ii++) {
    		    if ((ii & 7) == 0) mask = *mask_pointer++;
    		    data[ii] = (mask & 128) ?  tmp : UNDEFINED;
    		    mask <<= 1;
    		}
    		return 0;
    	    }
    	    fatal_error("unknown bitmap", "");
    	}


        strm.flags = (int) sec[5][21];
        strm.bits_per_sample = (int) sec[5][19];
        strm.block_size = (int) sec[5][22];
        strm.rsi = uint2(sec[5]+23);

	strm.next_in = sec[7]+5;
	strm.avail_in = uint4(sec[7]) - 5;

	numBitsNeeded = (int) sec[5][19];
	size = ((numBitsNeeded + 7)/8) * (size_t) ndata;

        if ((c = (unsigned char *) malloc(size)) == NULL) fatal_error("unpk: allocation error", "");

	strm.next_out = c;
	strm.avail_out = size;

        status = aec_buffer_decode(&strm);

        if (status != AEC_OK) fatal_error_i("unpk: aec decode error %d",status);

    	mask_pointer = (bitmap_flag == 255) ? NULL : sec[6] + 6;

    	unpk_0(data, c, mask_pointer, ((nbits+7)/8)*8, ndata, reference, bin_scale, dec_scale);

    	free(c);
    	return 0;
    }
#endif
    fatal_error_i("packing type %d not supported", packing);
    return 1;
}
