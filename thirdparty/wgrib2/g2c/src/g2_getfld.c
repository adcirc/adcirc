/** @file
 * @brief Return all the metadata, template values, Bit-map (if
 * applicable), and the unpacked data for a given data field.
 * @author Stephen Gilbert @date 2002-10-28
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Return all the metadata, template values, bit-map
 * (if applicable), and the unpacked data for a data field.
 *
 * All of the information returned is stored in a ::gribfield
 * structure, which is defined in file grib2.h. Users of this routine
 * will need to include grib2.h in their source code that calls this
 * routine.
 *
 * Since there can be multiple data fields packed into a GRIB2
 * message, the calling routine indicates which field is being
 * requested with the ifldnum argument.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-28 | Gilbert | Initial
 * 2013-08-08 | Vuong | Free up memory in array igds - free(igds)
 *
 * @param cgrib Character pointer to the GRIB2 message.
 * @param ifldnum Specifies which field in the GRIB2 message to
 * return. The first field is number 1, Fortran style.
 * @param unpack Boolean value indicating whether to unpack
 * bitmap/data field.
 * - 1 unpack bitmap (if present) and data values.
 * - 0 do not unpack bitmap and data values.
 * @param expand Boolean value indicating whether the data points
 * should be expanded to the correspond grid, if a bit-map is present.
 * This argument is ignored if unpack == 0 OR if the returned field
 * does not contain a bit-map.
 * - 1 if possible, expand data field to grid, inserting zero values
 * at gridpoints that are bitmapped out.  (SEE REMARKS2)
 * - 0 do not expand data field, leaving it an array of consecutive
 * data points for each "1" in the bitmap.
 * @param gfld pointer to structure ::gribfield containing all decoded
 * data for the data field.
 *
 * @return
 * - ::G2_NO_ERROR no error
 * - ::G2_GETFLD_NO_GRIB Beginning characters "GRIB" not found.
 * - ::G2_GETFLD_GRIB_VERSION GRIB message is not Edition 2.
 * - ::G2_GETFLD_INVAL The data field request number was not positive.
 * - ::G2_GETFLD_WRONG_END End string "7777" found, but not where expected.
 * - ::G2_GETFLD_WRONG_NFLDS GRIB message did not contain the requested number of data fields.
 * - ::G2_GETFLD_BAD_END End string "7777" not found at end of message.
 * - ::G2_GETFLD_INVAL_SEC Unrecognized Section encountered.
 * - ::G2_GETFLD_NO_DRT Data Representation Template 5.NN not yet implemented.
 * - ::G2_GETFLD_BAD_SEC1 Error unpacking Section 1.
 * - ::G2_GETFLD_BAD_SEC2 Error unpacking Section 2.
 * - ::G2_GETFLD_BAD_SEC3 Error unpacking Section 3.
 * - ::G2_GETFLD_BAD_SEC4 Error unpacking Section 4.
 * - ::G2_GETFLD_BAD_SEC5 Error unpacking Section 5.
 * - ::G2_GETFLD_BAD_SEC6 Error unpacking Section 6.
 * - ::G2_GETFLD_BAD_SEC7 Error unpacking Section 7.
 * - ::G2_GETFLD_NO_BITMAP Previous bitmap specified, yet none exists.
 *
 * @note Struct ::gribfield is allocated by this routine and it also
 * contains pointers to many arrays of data that were allocated during
 * decoding. Users are encouraged to free up this memory, when it is
 * no longer needed, by an explicit call to routine g2_free().
 *
 * ## Example:
 * @code
 *              #include "grib2.h"
 *              gribfield *gfld;
 *              ret=g2_getfld(cgrib,1,1,1,&gfld);
 *                ...
 *              g2_free(gfld);
 * @endcode
 *
 * Routine g2_info() can be used to first determine how many data
 * fields exist in a given GRIB message.
 *
 * @note It may not always be possible to expand a bit-mapped data
 * field. If a pre-defined bit-map is used and not included in the
 * GRIB2 message itself, this routine would not have the necessary
 * information to expand the data. In this case, gfld->expanded would
 * would be set to 0 (false), regardless of the value of input
 * argument expand.
 *
 * @author Stephen Gilbert @date 2002-10-28
 */
g2int
g2_getfld(unsigned char *cgrib, g2int ifldnum, g2int unpack, g2int expand,
          gribfield **gfld)
{

    g2int have3 = 0, have4 = 0, have5 = 0, have6 = 0, have7 = 0;
    g2int numfld = 0, j, n, istart, iofst, ipos;
    g2int disc, ver, lensec0, lengrib, lensec, isecnum;
    g2int *igds;
    g2int *bmpsave;
    float *newfld;
    gribfield *lgfld;
    g2int jerr;

    /* Allocate storage for the field struct. The caller must free
     * this memory. */
    lgfld = calloc(1, sizeof(gribfield));
    *gfld = lgfld;

    /* Check for valid request number. */
    if (ifldnum <= 0)
    {
        printf("g2_getfld: Request for field number must be positive.\n");
        free(lgfld);
        return G2_GETFLD_INVAL;
    }

    /*  Check for beginning of GRIB message in the first 100 bytes. */
    istart = -1;
    for (j = 0; j < 100; j++)
    {
        if (cgrib[j] == 'G' && cgrib[j + 1] == 'R' && cgrib[j + 2] == 'I' &&
            cgrib[j + 3] == 'B')
        {
            istart = j;
            break;
        }
    }
    if (istart == -1)
    {
        printf("g2_getfld:  Beginning characters GRIB not found.\n");
        free(lgfld);
        return G2_GETFLD_NO_GRIB;
    }

    /*  Unpack Section 0 - Indicator Section. */
    iofst = 8 * (istart + 6);
    gbit(cgrib, &disc, iofst, 8); /* Discipline */
    iofst = iofst + 8;
    gbit(cgrib, &ver, iofst, 8); /* GRIB edition number */
    iofst = iofst + 8;
    iofst = iofst + 32;
    gbit(cgrib, &lengrib, iofst, 32); /* Length of GRIB message */
    iofst = iofst + 32;
    lensec0 = 16;
    ipos = istart + lensec0;

    /*  Currently handles only GRIB Edition 2. */
    if (ver != 2)
    {
        printf("g2_getfld: can only decode GRIB edition 2.\n");
        free(lgfld);
        return G2_GETFLD_GRIB_VERSION;
    }

    /*  Loop through the remaining sections keeping track of the
     *  length of each.  Also keep the latest Grid Definition Section
     *  info. Unpack the requested field number. */
    for (;;)
    {
        /*    Check to see if we are at end of GRIB message */
        if (cgrib[ipos] == '7' && cgrib[ipos + 1] == '7' && cgrib[ipos + 2] == '7' &&
            cgrib[ipos + 3] == '7')
        {
            ipos = ipos + 4;

            /*    If end of GRIB message not where expected, issue
             *    error. */
            if (ipos != istart + lengrib)
            {
                printf("g2_getfld: '7777' found, but not where expected.\n");
                free(lgfld);
                return G2_GETFLD_WRONG_END;
            }
            break;
        }

        /* Get length of Section and Section number. */
        iofst = ipos * 8;
        gbit(cgrib, &lensec, iofst, 32); /* Get Length of Section */
        iofst = iofst + 32;
        gbit(cgrib, &isecnum, iofst, 8); /* Get Section number */
        iofst = iofst + 8;

        /* Check to see if section number is valid. */
        if (isecnum < 1 || isecnum > 7)
        {
            printf("g2_getfld: Unrecognized Section Encountered=%ld\n", isecnum);
            free(lgfld);
            return G2_GETFLD_INVAL_SEC;
        }

        /*   If found Section 1, decode elements in Identification Section. */
        if (isecnum == 1)
        {
            iofst = iofst - 40; /* reset offset to beginning of section */
            if (g2_unpack1(cgrib, &iofst, &lgfld->idsect, &lgfld->idsectlen))
            {
                g2_free(lgfld);
                return G2_GETFLD_BAD_SEC1;
            }
        }

        /*   If found Section 2, Grab local section. Save in case this
         *   is the latest one before the requested field. */
        if (isecnum == 2)
        {
            iofst = iofst - 40; /* reset offset to beginning of section */
            if (lgfld->local)
                free(lgfld->local);
            if (g2_unpack2(cgrib, &iofst, &lgfld->locallen, &lgfld->local))
            {
                g2_free(lgfld);
                return G2_GETFLD_BAD_SEC2;
            }
        }

        /* If found Section 3, unpack the GDS info using the
         * appropriate template.  Save in case this is the latest grid
         * before the requested field. */
        if (isecnum == 3)
        {
            iofst = iofst - 40; /* reset offset to beginning of section */
            if (lgfld->igdtmpl)
                free(lgfld->igdtmpl);
            if (lgfld->list_opt)
                free(lgfld->list_opt);
            if (g2_unpack3(cgrib, &iofst, &igds, &lgfld->igdtmpl,
                           &lgfld->igdtlen, &lgfld->list_opt, &lgfld->num_opt))
            {
                g2_free(lgfld);
                return G2_GETFLD_BAD_SEC3;
            }

            have3 = 1;
            lgfld->griddef = igds[0];
            lgfld->ngrdpts = igds[1];
            lgfld->numoct_opt = igds[2];
            lgfld->interp_opt = igds[3];
            lgfld->igdtnum = igds[4];
            free(igds);
        }

        /* If found Section 4, check to see if this field is the one
         * requested. */
        if (isecnum == 4)
        {
            numfld = numfld + 1;
            if (numfld == ifldnum)
            {
                lgfld->discipline = disc;
                lgfld->version = ver;
                lgfld->ifldnum = ifldnum;
                lgfld->unpacked = unpack;
                lgfld->expanded = 0;
                iofst = iofst - 40; /* reset offset to beginning of section */
                if (g2_unpack4(cgrib, &iofst, &lgfld->ipdtnum, &lgfld->ipdtmpl,
                               &lgfld->ipdtlen, &lgfld->coord_list, &lgfld->num_coord))
                {
                    g2_free(lgfld);
                    return G2_GETFLD_BAD_SEC4;
                }

                have4 = 1;
            }
        }

        /*   If found Section 5, check to see if this field is the one
         *   requested. */
        if (isecnum == 5 && numfld == ifldnum)
        {
            iofst = iofst - 40; /* reset offset to beginning of section */
            if (g2_unpack5(cgrib, &iofst, &lgfld->ndpts, &lgfld->idrtnum,
                           &lgfld->idrtmpl, &lgfld->idrtlen))
            {
                g2_free(lgfld);
                return G2_GETFLD_BAD_SEC5;
            }

            have5 = 1;
        }

        /*   If found Section 6, Unpack bitmap. Save in case this is
         *   the latest bitmap before the requested field. */
        if (isecnum == 6)
        {
            if (unpack)
            {                          /* unpack bitmap */
                iofst = iofst - 40;    /* reset offset to beginning of section */
                bmpsave = lgfld->bmap; /* save pointer to previous bitmap */
                if (g2_unpack6(cgrib, &iofst, lgfld->ngrdpts, &lgfld->ibmap,
                               &lgfld->bmap))
                {
                    g2_free(lgfld);
                    return G2_GETFLD_BAD_SEC6;
                }

                have6 = 1;
                if (lgfld->ibmap == 254) /* use previously specified bitmap */
                {
                    if (bmpsave)
                        lgfld->bmap = bmpsave;
                    else
                    {
                        printf("g2_getfld: Prev bit-map specified, but none exist.\n");
                        return G2_GETFLD_NO_BITMAP;
                    }
                }
                else /* get rid of it */
                    if (bmpsave)
                        free(bmpsave);
            }
            else
            {                                         /* do not unpack bitmap */
                gbit(cgrib, &lgfld->ibmap, iofst, 8); /* Get BitMap Indicator */
                have6 = 1;
            }
        }

        /*   If found Section 7, check to see if this field is the one
         *   requested. */
        if (isecnum == 7 && numfld == ifldnum && unpack)
        {
            iofst = iofst - 40; /* reset offset to beginning of section */
            if ((jerr = g2_unpack7(cgrib, &iofst, lgfld->igdtnum, lgfld->igdtmpl,
                                   lgfld->idrtnum, lgfld->idrtmpl, lgfld->ndpts,
                                   &lgfld->fld)))
            {
                printf("g2_getfld: return from g2_unpack7 = %d \n", (int)jerr);
                g2_free(lgfld);
                return G2_GETFLD_BAD_SEC7;
            }

            have7 = 1;

            /*  If bitmap is used with this field,  expand data field */
            /*  to grid, if possible. */
            if (lgfld->ibmap != 255 && lgfld->bmap != 0)
            {
                if (expand == 1)
                {
                    n = 0;
                    newfld = calloc(lgfld->ngrdpts, sizeof(float));
                    for (j = 0; j < lgfld->ngrdpts; j++)
                    {
                        if (lgfld->bmap[j] == 1)
                            newfld[j] = lgfld->fld[n++];
                    }
                    free(lgfld->fld);
                    lgfld->fld = newfld;
                    lgfld->expanded = 1;
                }
                else
                {
                    lgfld->expanded = 0;
                }
            }
            else
            {
                lgfld->expanded = 1;
            }
        }

        /* Check to see if we read pass the end of the GRIB message
         * and missed the terminator string '7777'. */
        ipos = ipos + lensec; /* Update beginning of section pointer */
        if (ipos > istart + lengrib)
        {
            printf("g2_getfld: '7777'  not found at end of GRIB message.\n");
            g2_free(lgfld);
            return G2_GETFLD_BAD_END;
        }

        /*  If unpacking requested, return when all sections have been
         *  processed. */
        if (unpack && have3 && have4 && have5 && have6 && have7)
            return G2_NO_ERROR;

        /*  If unpacking is not requested, return when sections 3
         *  through 6 have been processed. */
        if (!unpack && have3 && have4 && have5 && have6)
            return G2_NO_ERROR;
    }

    /*  If exited from above loop, the end of the GRIB message was
     *  reached before the requested field was found. */
    printf("g2_getfld: GRIB message contained %ld different fields.\n", numfld);
    printf("g2_getfld: The request was for field %ld.\n", ifldnum);
    g2_free(lgfld);
    return G2_GETFLD_WRONG_NFLDS;
}
