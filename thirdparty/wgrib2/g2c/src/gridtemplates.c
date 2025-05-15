/** @file
 * @brief Returns grid template information for a specified Grid
 * Definition Template for Section 3 - the Grid Definition Section
 * (GDS).
 *
 * Each Template has three parts:
 * 1. The number of entries in the template (mapgridlen).
 * 2. A map of the template (mapgrid), which contains the number of
 * octets in which to pack each of the template values.
 * 3. A logical value (needext) that indicates whether the Template
 * needs to be extended.
 *
 * In some cases the number of entries in a template can vary
 * depending upon values specified in the "static" part of the
 * template. (See Template 3.120 as an example).
 *
 * @note Array mapgrid contains the number of octets in which the
 * corresponding template values will be stored. A negative value in
 * mapgrid is used to indicate that the corresponding template entry
 * can contain negative values. This information is used later when
 * packing (or unpacking) the template data values. Negative data
 * values in GRIB are stored with the left most bit set to one, and a
 * negative number of octets value in mapgrid indicates that this
 * possibility should be considered. The number of octets used to
 * store the data value in this case would be the absolute value of
 * the negative value in mapgrid.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2001-06-28 | Gilbert | Initial
 * 2007-08-16 | Vuong | Added GDT 3.204  Curvilinear Orthogonal Grid
 * 2008-07-08 | Vuong | Added GDT 3.32768 Rotate Lat/Lon E-grid (Arakawa)
 * 2009-01-14 | Vuong | Changed structure name template to gtemplate
 * 2010-05-11 | Vuong | Added GDT 3.32769 Rotate Lat/Lon Non-E Staggered grid (Arakawa)
 * 2013-08-06 | Vuong | Added GDT 3.4, 3.5, 3.12, 3.101, 3.140
 * 2022-10-16 | Hartnett | Added g2c_get_grid_template().
 * 2025-01-06 | Stahl| Added GDT 3.13, 3.23, 3.33, 3.61, 3.62, 3.63, 3.150
 *
 * @author Stephen Gilbert @date 2001-06-28
 */
#include "grib2_int.h"

/**
 * Struct for grid template.
 */
struct gridtemplate
{
    g2int template_num;                         /**< Template number. */
    g2int mapgridlen;                           /**< The number of entries in the template. */
    g2int needext;                              /**< Does template need extension? */
    g2int mapgrid[G2C_MAX_GDS_TEMPLATE_MAPLEN]; /**< Number of bytes for each template value. */
};

/**
 * Templates grid.
 */
static const struct gridtemplate templatesgrid[G2C_MAX_GDS_TEMPLATE] =
    {
        /* 3.0: Lat/Lon grid */
        {0, 19, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1}},
        /* 3.1: Rotated Lat/Lon grid */
        {1, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, 4}},
        /* 3.2: Stretched Lat/Lon grid */
        {2, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, -4}},
        /* 3.3: Stretched & Rotated Lat/Lon grid */
        {3, 25, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, 4, -4, 4, -4}},
        /* Added GDT 3.4,3.5    (08/05/2013) */
        /* 3.4: Variable resolution Latitude/Longitude */
        {4, 13, 1, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, 1, 1}},
        /* 3.5: Variable resolution rotate Latitude/Longitude */
        {5, 16, 1, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, 1, 1, -4, 4, 4}},
        /* 3.12: Transverse Mercator */
        {12, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 1, 4, 4, -4, -4, -4, -4}},
        /* 3.101: General unstructured grid */
        {101, 4, 0, {1, 4, 1, -4}},
        /* 3.140: Lambert Azimuthal Equal Area Projection */
        {140, 17, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 4, 4, 1, 4, 4, 1}},

        /* 3.10: Mercator */
        {10, 19, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, -4, 4, 1, 4, 4, 4}},
        /* 3.20: Polar Stereographic Projection */
        {20, 18, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 1}},
        /* 3.30: Lambert Conformal */
        {30, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 1, -4, -4, -4, 4}},
        /* 3.31: Albers equal area */
        {31, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 1, -4, -4, -4, 4}},
        /* 3.40: Guassian Lat/Lon */
        {40, 19, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1}},
        /* 3.41: Rotated Gaussian Lat/Lon */
        {41, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, 4}},
        /* 3.42: Stretched Gaussian Lat/Lon */
        {42, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, -4}},
        /* 3.43: Stretched and Rotated Gaussian Lat/Lon */
        {43, 25, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, -4, 4, 4, -4, 4, -4}},
        /* 3.50: Spherical Harmonic Coefficients */
        {50, 5, 0, {4, 4, 4, 1, 1}},
        /* 3.51: Rotated Spherical Harmonic Coefficients */
        {51, 8, 0, {4, 4, 4, 1, 1, -4, 4, 4}},
        /* 3.52: Stretched Spherical Harmonic Coefficients */
        {52, 8, 0, {4, 4, 4, 1, 1, -4, 4, -4}},
        /* 3.53: Stretched and Rotated Spherical Harmonic Coefficients */
        {53, 11, 0, {4, 4, 4, 1, 1, -4, 4, 4, -4, 4, -4}},
        /* 3.90: Space View Perspective or orthographic */
        {90, 21, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, 4, 4, 4, 4, 1, 4, 4, 4, 4}},
        /* 3.100: Triangular grid based on an icosahedron */
        {100, 11, 0, {1, 1, 2, 1, -4, 4, 4, 1, 1, 1, 4}},
        /* 3.110: Equatorial Azimuthal equidistant */
        {110, 16, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, 4, 4, 1, 1}},
        /* 3.120: Azimuth-range projection */
        {120, 7, 1, {4, 4, -4, 4, 4, 4, 1}},
        /* 3.204: Curvilinear Orthogonal Grid */
        {204, 19, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1}},
        /* 3.32768: Rot Lat/Lon E-grid (Arakawa) */
        {32768, 19, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1}},
        /* 3.32769: Rot Lat/Lon Non-E Staggered grid (Arakawa) */
        {32769, 21, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 4, 4}},
        /* 3.1000: Cross Section Grid */
        {1000, 20, 1, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, -4, 4, 1, 4, 4, 1, 2, 1, 1, 2}},
        /* 3.1100: Hovmoller Diagram Grid */
        {1100, 28, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, 4, -4, 4, 1, -4, 4, 1, 4, 1, -4, 1, 1, -4, 2, 1, 1, 1, 1, 1}},
        /* 3.1200: Time Section Grid */
        {1200, 16, 1, {4, 1, -4, 1, 1, -4, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2}},

        /* 3.13: Mercator with modelling subdomains definition */
        {13, 23, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, -4, 4, 1, 4, 4, 4, 4, 4, 4, 4}},
        /* 3.23: Polar stereographic projection with modelling subdomains definition */
        {23, 22, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 1, 4, 4, 4, 4}},
        /* 3.30: Lambert conformal with modelling subdomains definition */
        {33, 26, 0, {1, 1, 4, 1, 4, 1, 4, 4, 4, -4, 4, 1, -4, 4, 4, 4, 1, 1, -4, -4, -4, 4, 4, 4, 4, 4}},
        /* 3.61: Spectral Mercator with modelling subdomains definition */
        {61, 23, 0, {1, 4, 4, 1, 8, 8, 8, 8, 8, 8, 1, 1, 4, 1, 4, 1, 4, -4, 4, -4, -4, 4, 4}},
        /* 3.62: Spectral Polar Stereographic with modelling subdomains definition */
        {62, 23, 0, {1, 4, 4, 1, 8, 8, 8, 8, 8, 8, 1, 1, 4, 1, 4, 1, 4, -4, 4, 1, -4, 4, 1}},
        /* 3.63: Spectral Lambert Conformal with modelling subdomains definition */
        {63, 26, 0, {1, 4, 4, 1, 8, 8, 8, 8, 8, 8, 1, 1, 4, 1, 4, 1, 4, -4, 4, -4, 4, 1, -4, -4, -4, 4}},
        /* 3.150: Hierarchical Equal Area isoLatitude Pixelization grid (HEALPix) */
        {150, 13, 0, {1, 1, 4, 1, 4, 1, 4, 1, 4, 4, 1, 1, 1}}};

/**
 * This function returns the index of specified Grid Definition
 * Template in array templates for [Section 3 - the Grid Definition
 * Section
 * (GDS)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml).
 *
 * @param number The number of the Grid Definition Template being
 * requested.
 *
 * @return Index of the Grid Definition Template(GDT) in array
 * templates, if template exists. -1, otherwise.
 *
 * @author Stephen Gilbert @date 2001-06-28
 */
static g2int
getgridindex(g2int number)
{
    g2int j, getgridindex = -1;

    for (j = 0; j < G2C_MAX_GDS_TEMPLATE; j++)
    {
        if (number == templatesgrid[j].template_num)
        {
            getgridindex = j;
            return (getgridindex);
        }
    }

    return (getgridindex);
}

/**
 * This subroutine returns grid template information for a specified
 * Grid Definition Template for [Section 3 - the Grid Definition
 * Section
 * (GDS)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml). The
 * number of entries in the template is returned along with a map of
 * the number of octets occupied by each entry. Also, a flag is
 * returned to indicate whether the template would need to be
 * extended.
 *
 * This function allocates storage for the template. The returned
 * pointer must be freed by the caller.
 *
 * @param number The number of the Grid Definition Template that is
 * being requested.
 *
 * @return Pointer to the returned template struct (must be freed by
 * caller). Returns NULL pointer, if template not found.
 *
 * @author Stephen Gilbert @date 2000-05-09
 */
gtemplate *
getgridtemplate(g2int number)
{
    g2int index;
    gtemplate *new;

    index = getgridindex(number);

    if (index != -1)
    {
        new = malloc(sizeof(gtemplate));
        new->type = 3;
        new->num = templatesgrid[index].template_num;
        new->maplen = templatesgrid[index].mapgridlen;
        new->needext = templatesgrid[index].needext;
        new->map = (g2int *)templatesgrid[index].mapgrid;
        new->extlen = 0;
        new->ext = NULL;
        return (new);
    }
    else
    {
        printf("getgridtemplate: GDT Template 3.%d not defined.\n", (int)number);
        return (NULL);
    }

    return (NULL);
}

/**
 * This subroutine generates the remaining octet map for a given Grid
 * Definition Template, if required. Some Templates can vary
 * depending on data values given in an earlier part of the Template,
 * and it is necessary to know some of the earlier entry values to
 * generate the full octet map of the Template.
 *
 * This function allocates memory for the extension. The pointer ext
 * in the gtemplate struct must be freed to prevent memory leaks.
 *
 * @param number The number of the Grid Definition
 * Template that is being requested.
 * @param template The grid definition template array.
 *
 * @return Pointer to the returned template struct. Returns NULL
 * pointer, if template not found.
 *
 * @author Stephen Gilbert @date 2000-05-09
 */
gtemplate *
extgridtemplate(g2int number, g2int *template)
{
    gtemplate *new;
    g2int index, i;

    index = getgridindex(number);
    if (index == -1)
        return NULL;

    new = getgridtemplate(number);

    if (!new->needext)
        return (new);

    if (number == 120)
    {
        new->extlen = template[1] * 2;
        new->ext = malloc(sizeof(g2int) * new->extlen);
        for (i = 0; i < new->extlen; i++)
        {
            if (i % 2 == 0)
                new->ext[i] = 2;
            else
                new->ext[i] = -2;
        }
    }
    else if (number == 4 || number == 5)
    {
        /* The extension is of length template[7] + template[8]. The first
         * template[7] values are 4, the next template[8] values are -4. */
        new->extlen = template[7] + template[8];
        new->ext = malloc(sizeof(g2int) * new->extlen);
        for (i = 0; i < new->extlen; i++)
            new->ext[i] = i < template[7] ? 4 : -4;
    }
    else if (number == 1000)
    {
        new->extlen = template[19];
        new->ext = malloc(sizeof(g2int) * new->extlen);
        for (i = 0; i < new->extlen; i++)
        {
            new->ext[i] = 4;
        }
    }
    else if (number == 1200)
    {
        new->extlen = template[15];
        new->ext = malloc(sizeof(g2int) * new->extlen);
        for (i = 0; i < new->extlen; i++)
        {
            new->ext[i] = 4;
        }
    }

    return (new);
}

/**
 * Get grid template extension information.
 *
 * @param grid_template_num The grid template number.
 * @param template Pointer to array that contains the template values.
 * @param extlen Pointer that gets the length of the extension. Ignored if NULL.
 * @param ext Pointer that gets template extension array, if there is
 * one. Memory must be allocated by the caller. Ignored if NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_EINVAL Invalid input.
 * - ::G2C_ENOTEMPLATE Template not found.
 * - ::G2C_ENOMEM Out of memory.
 *
 * @author Ed Hartnett @date 10/16/22
 */
int
g2c_get_grid_template_extension(int grid_template_num, int *template,
                                int *extlen, int *ext)
{
    int j, t;

    /* Check input. */
    if (!template)
        return G2C_EINVAL;

    /* Look through the array of templates to find a matching one. */
    for (j = 0; j < G2C_MAX_GDS_TEMPLATE; j++)
    {
        if (grid_template_num == templatesgrid[j].template_num)
        {
            /* Is there an extension to this template? */
            if (templatesgrid[j].needext)
            {
                gtemplate *gt;
                g2int *template8;
                int e;

                /* Copy templage to g2int for extgridtemplate() function. */
                if (!(template8 = malloc(sizeof(g2int) * templatesgrid[j].mapgridlen)))
                    return G2C_ENOMEM;
                for (t = 0; t < templatesgrid[j].mapgridlen; t++)
                    template8[t] = template[t];
                if (!(gt = extgridtemplate(grid_template_num, template8)))
                    return G2C_ENOTEMPLATE;
                free(template8);
                if (extlen)
                    *extlen = gt->extlen;
                if (ext)
                    for (e = 0; e < gt->extlen; e++)
                        ext[e] = gt->ext[e];
                free(gt->ext);
                free(gt);
            }
            else
            {
                if (extlen)
                    *extlen = 0;
            }

            /* Done. */
            return G2C_NOERROR;
        }
    }

    /* If we didn't find a template, return an error. */
    return G2C_ENOTEMPLATE;
}

/**
 * Get grid template information.
 *
 * The grid template consists of a template map, its length, and, for
 * some templates, an extra extension map, and its length. If an
 * extension is needed, use g2c_get_grid_template_extension() to get
 * it.
 *
 * @param grid_template_num The grid template number.
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
 * @author Ed Hartnett @date 10/16/22
 */
int
g2c_get_grid_template(int grid_template_num, int *maplen, int *map, int *needext)
{
    int j, m;

    /* Look through the array of templates to find a matching one. */
    for (j = 0; j < G2C_MAX_GDS_TEMPLATE; j++)
    {
        if (grid_template_num == templatesgrid[j].template_num)
        {
            /* Copy maplen and map if the caller wants them. */
            if (maplen)
                *maplen = templatesgrid[j].mapgridlen;
            if (map)
                for (m = 0; m < templatesgrid[j].mapgridlen; m++)
                    map[m] = templatesgrid[j].mapgrid[m];
            if (needext)
                *needext = templatesgrid[j].needext;

            /* Done. */
            return G2C_NOERROR;
        }
    }

    /* If we didn't find a template, return an error. */
    return G2C_ENOTEMPLATE;
}

/**
 * Get initial length (number of entries) in static part of 
 * Grid Definition Template.
 *
 * @param grid_template_num The Grid template number.
 * @param maplen Pointer that gets the length of the map. Ignored if
 * NULL.
 *
 * @return
 * - ::G2C_NOERROR No error.
 * - ::G2C_ENOTEMPLATE Template not found.
 *
 * @author Alyson Stahl @date 01/21/25
 */
int
g2c_get_gdt_len(int grid_template_num, int *maplen)
{
    int j;

    /* Look through the array of templates to find a matching one. */
    for (j = 0; j < G2C_MAX_GDS_TEMPLATE; j++)
    {
        if (grid_template_num == templatesgrid[j].template_num)
        {
            if (maplen)
                *maplen = templatesgrid[j].mapgridlen;

            /* Done. */
            return G2C_NOERROR;
        }
    }

    /* If we didn't find a template, return an error. */
    return G2C_ENOTEMPLATE;
}