/**
 * @file
 * @brief Functions, Structures, and Data for Data Representation
 * Templates (DRT).
 *
 * Each Template has three parts:
 * 1. The number of entries in the template (mapdrslen).
 * 2. A map of the template (mapdrs), which contains the number of
 * octets in which to pack each of the template values.
 * 3. A logical value (needext) that indicates whether the
 * Template needs to be extended. In some cases the number of entries
 * in a template can vary depending upon values specified in the
 * "static" part of the template. (See Template 5.1 as an example.)
 *
 * @note Array mapdrs contains the number of octets in which the
 * corresponding template values will be stored. A negative value in
 * mapdrs is used to indicate that the corresponding template entry
 * can contain negative values. This information is used later when
 * packing (or unpacking) the template data values. Negative data
 * values in GRIB are stored with the left most bit set to one, and a
 * negative number of octets value in mapdrs indicates that this
 * possibility should be considered. The number of octets used to
 * store the data value in this case would be the absolute value of
 * the negative value in mapdrs.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2001-06-28 | Gilbert | Initial
 * 2009-01-14 | Vuong | Changed structure name template to gtemplate
 * 2022-10-18 | Hartnett | Added g2c_get_drs_template().
 * 2023-09-08 | Engle | Added template 5.42.
 * 2025-01-06 | Stahl | Added template 5.53
 *
 * @author Stephen Gilbert @date 2001-06-28
 */
#include "grib2_int.h"

/**
 * Stuct for GRIB2 Data Representation Section (DRS) template.
 */
struct drstemplate
{
    g2int template_num;                        /**< The number of entries in the template. */
    g2int mapdrslen;                           /**< Length of map of the template. */
    g2int needext;                             /**< Whether the Template needs to be extended. */
    g2int mapdrs[G2C_MAX_DRS_TEMPLATE_MAPLEN]; /**< A map of the template. */
};

/**
 * Stuct holding data for GRIB2 Data Representation Section (DRS)
 * template.
 */
static const struct drstemplate templatesdrs[G2C_MAX_DRS_TEMPLATE] =
    {
        /** [5.0: Grid point data - Simple
     * Packing](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-0.shtml) */
        {0, 5, 0, {4, -2, -2, 1, 1}},

        /** [5.2: Grid point data - Complex
     * Packing](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-2.shtml) */
        {2, 16, 0, {4, -2, -2, 1, 1, 1, 1, 4, 4, 4, 1, 1, 4, 1, 4, 1}},

        /** [5.3: Grid point data - Complex Packing and spatial
     * differencing](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-3.shtml) */
        {3, 18, 0, {4, -2, -2, 1, 1, 1, 1, 4, 4, 4, 1, 1, 4, 1, 4, 1, 1, 1}},

        /** [5.50: Spectral Data - Simple
     * Packing](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-50.shtml) */
        {50, 5, 0, {4, -2, -2, 1, 4}},

        /** [5.51: Spherical Harmonics data - Complex
     * packing](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-51.shtml) */
        {51, 10, 0, {4, -2, -2, 1, -4, 2, 2, 2, 4, 1}},

        /* 5.1: Matrix values at gridpoint - Simple packing.
     * Comment from Stephen Gilbert in 2021:
     *
     * This encoder/decoder was written in the early days of GRIB2
     * adoption as a standard. It was used to help WMO validate the
     * templates in the specification by sharing GRIB2 encoded message
     * with other organizations to verify that the data could be
     * transmitted and processed successfully.
     *
     * We did not have a use case for DRS template 5.1 at that time
     * and did not produce any GRIB2 messages using that template. It
     * appears that other organizations did not work on it as
     * well. The latest GRIB2 specification still includes the DRS
     * Template 5.1 definition, but there is a disclaimer to use it
     * with caution, since it has not yet been validated. I assume we
     * commented it out because it was not validated, which means it's
     * definition could possibly change during any validation attempts
     * in the future.
     */

        /* {1, 15, 1, {4, -2, -2, 1, 1, 1, 4, 2, 2, 1, 1, 1, 1, 1, 1}}, */

        /** [5.40: Grid point data - JPEG2000
     * encoding](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-40.shtml) */
        {40, 7, 0, {4, -2, -2, 1, 1, 1, 1}},

        /** [5.41: Grid point data - PNG
     * encoding](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-41.shtml) */
        {41, 5, 0, {4, -2, -2, 1, 1}},

        /** [5.42: Grid point data - CCSDS recommended lossless compress (libaec)
     * encoding](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-42.shtml) */
        {42, 8, 0, {4, -2, -2, 1, 1, 1, 1, 2}},

        /** [5.53: Spectral data for limited area models - Complex packing
     * encoding](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp5-53.shtml) */
        {53, 11, 0, {4, -2, -2, 1, 1, 1, -4, 2, 2, 4, 1}},

        /** 5.40000: Grid point data - JPEG2000 encoding
     *
     * This is a local template number, from a time before WMO standardized use of JPEG2000
     * with 5.40. This should not be used in new data files. Use 5.40 instead. */
        {40000, 7, 0, {4, -2, -2, 1, 1, 1, 1}},

        /** 5.40010: Grid point data - PNG encoding
     *
     * This is a local template number, from a time before WMO standardized use of PNG
     * with 5.41. This should not be used in new data files. Use 5.41 instead. */
        {40010, 5, 0, {4, -2, -2, 1, 1}},
};

/**
 * This function returns the index of specified Data Representation
 * Template.
 *
 * @param number The number of the Data Representation Template that
 * is being requested.
 *
 * @return Index of the DRT in array gtemplates, if it exists. -1,
 * otherwise.
 *
 * @author Stephen Gilbert @date 2001-06-28
 */
static g2int
getdrsindex(g2int number)
{
    g2int j, getdrsindex = -1;

    for (j = 0; j < G2C_MAX_DRS_TEMPLATE; j++)
    {
        if (number == templatesdrs[j].template_num)
        {
            getdrsindex = j;
            return (getdrsindex);
        }
    }

    return (getdrsindex);
}

/**
 * This subroutine returns DRS template information for a specified
 * Data Representation Template. The number of entries in the
 * template is returned along with a map of the number of octets
 * occupied by each entry. Also, a flag is returned to indicate
 * whether the template would need to be extended.
 *
 * @param number The number of the Data Representation
 * Template that is being requested.
 *
 * @return Pointer to the returned template struct. Returns NULL
 * if template not found.
 *
 * @author Stephen Gilbert @date 2000-05-11
 */
gtemplate *
getdrstemplate(g2int number)
{
    g2int index;
    gtemplate *new;

    index = getdrsindex(number);

    if (index != -1)
    {
        new = malloc(sizeof(gtemplate));
        new->type = 5;
        new->num = templatesdrs[index].template_num;
        new->maplen = templatesdrs[index].mapdrslen;
        new->needext = templatesdrs[index].needext;
        new->map = (g2int *)templatesdrs[index].mapdrs;
        new->extlen = 0;
        new->ext = NULL;
        return (new);
    }
    else
    {
        printf("getdrstemplate: DRS Template 5.%d not defined.\n", (int)number);
        return (NULL);
    }

    return (NULL);
}

/**
 * This subroutine generates the remaining octet map for a given Data
 * Representation Template, if required. Some Templates can vary
 * depending on data values given in an earlier part of the Template,
 * and it is necessary to know some of the earlier entry values to
 * generate the full octet map of the Template.
 *
 * @param number The number of the Data Representation Template that
 * is being requested.
 * @param list The list of values for each entry in the the Data
 * Representation Template.
 *
 * @return Pointer to the returned template struct. Returns NULL
 * pointer, if template not found.
 *
 * @author Stephen Gilbert @date 2000-05-11
 */
gtemplate *
extdrstemplate(g2int number, g2int *list)
{
    gtemplate *new;

    if (getdrsindex(number) == -1)
        return NULL;

    new = getdrstemplate(number);

    if (!new->needext)
        return (new);

    /* This template is commented out (see comment in struct
     * drstemplate for explanation). */
    /* if (number == 1) */
    /* { */
    /*     new->extlen = list[10] + list[12]; */
    /*     new->ext = malloc(sizeof(g2int) * new->extlen); */
    /*     for (i = 0; i < new->extlen; i++) */
    /*     { */
    /*         new->ext[i] = 4; */
    /*     } */
    /* } */

    return new;
}

/**
 * Get DRS template information.
 *
 * The DRS template consists of a template map, and its length. There
 * are no supported DRS templates with extensions.
 *
 * @param drs_template_num The DRS template number.
 * @param maplen Pointer that gets the length of the map. Ignored if
 * NULL.
 * @param map Pointer that gets the map as an array of int. Memory
 * must be allocated by caller. Ignored if NULL.
 * @param needext Pointer that a non-zero value if an extension to
 * this template is needed. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOTEMPLATE Template not found.
 *
 * @author Ed Hartnett @date 10/18/22
 */
int
g2c_get_drs_template(int drs_template_num, int *maplen, int *map, int *needext)
{
    int j, m;

    /* Look through the array of templates to find a matching one. */
    for (j = 0; j < G2C_MAX_DRS_TEMPLATE; j++)
    {
        if (drs_template_num == templatesdrs[j].template_num)
        {
            /* Copy maplen and map if the caller wants them. */
            if (maplen)
                *maplen = templatesdrs[j].mapdrslen;
            if (map)
                for (m = 0; m < templatesdrs[j].mapdrslen; m++)
                    map[m] = templatesdrs[j].mapdrs[m];
            if (needext)
                *needext = templatesdrs[j].needext;

            /* Done. */
            return G2C_NOERROR;
        }
    }

    /* If we didn't find a template, return an error. */
    return G2C_ENOTEMPLATE;
}
