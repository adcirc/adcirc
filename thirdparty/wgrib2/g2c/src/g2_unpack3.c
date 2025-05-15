/** @file
 * @brief Unpack Section 3 (Grid Definition Section) of a GRIB2
 * message.
 * @author Stephen Gilbert @date 2002-10-31
 */
#include "grib2_int.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Unpack [Section 3 (Grid Definition
 * Section)](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_sect3.shtml)
 * of a GRIB2 message.
 *
 * ### Program History Log
 * Date | Programmer | Comments
 * -----|------------|---------
 * 2002-10-31 | Gilbert | Initial
 * 2009-01-14 | Vuong | Changed structure name template to gtemplate
 *
 * @param cgrib Char array ontaining Section 3 of the GRIB2 message.
 * @param iofst Pointer to g2int which contains the bit offset for the
 * beginning of Section 3 in cgrib.
 * @param igds Pointer to a pointer which will get a pointer to memory
 * allocated for the GDS array, of length 5. The array will contain
 * information read from the appropriate GRIB Grid Definition Section
 * 3 for the field being returned.
 * - igds[0] Source of grid definition (see [Table
 * 3.0](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-0.shtml)).
 * - igds[1] Number of grid points in the defined grid.
 * - igds[2] Number of octets needed for each additional grid points
 *   definition. Used to define number of points in each row (or
 *   column) for non-regular grids. = 0, if using regular grid.
 * - igds[3] Interpretation of list for optional points
 *   definition. (See [Table
 *   3.11](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-11.shtml))
 * - igds[4] Grid Definition Template Number (see [Table
 3.1](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table3-1.shtml)).
 * @param igdstmpl Pointer a pointer to g2int, which will get a
 * pointer to an allocated array that contians containing the data
 * values from the Grid Definition Template specified by igds[4].
 * @param mapgridlen A pointer tat gets the number of elements in
 * igdstmpl. i.e. number of entries in Grid Defintion Template
 * specified by igds[4].
 * @param ideflist (Used if igds[2] .ne. 0) Pointer to integer array
 * containing the number of grid points contained in each row (or
 * column).
 * @param idefnum (Used if igds[2] != 0.) A pointer that gets the
 * number of entries in array ideflist - i.e. number of rows (or
 * columns) for which optional grid points are defined.
 *
 * @return
 * - ::G2_NO_ERROR No error.
 * - ::G2_UNPACK_BAD_SEC Array passed had incorrect section number.
 * - ::G2_UNPACK3_BAD_GDT message contains an undefined Grid Definition Template.
 * - ::G2_UNPACK_NO_MEM Memory allocation error.
 *
 * @author Stephen Gilbert @date 2002-10-31
 */

g2int
g2_unpack3(unsigned char *cgrib, g2int *iofst, g2int **igds, g2int **igdstmpl,
           g2int *mapgridlen, g2int **ideflist, g2int *idefnum)
{
    g2int i, j, nbits, isecnum;
    g2int lensec, ibyttem = 0, isign, newlen;
    g2int *ligds, *ligdstmpl = NULL, *lideflist = NULL;
    gtemplate *mapgrid;

    *igds = NULL;
    *igdstmpl = NULL;
    *ideflist = NULL;

    gbit(cgrib, &lensec, *iofst, 32); /* Get Length of Section */
    *iofst = *iofst + 32;
    gbit(cgrib, &isecnum, *iofst, 8); /* Get Section Number */
    *iofst = *iofst + 8;

    if (isecnum != 3)
    {
        *idefnum = 0;
        *mapgridlen = 0;
        return G2_UNPACK_BAD_SEC;
    }

    ligds = calloc(5, sizeof(g2int));
    *igds = ligds;

    gbit(cgrib, &ligds[0], *iofst, 8); /* Get source of Grid def. */
    *iofst = *iofst + 8;
    gbit(cgrib, &ligds[1], *iofst, 32); /* Get number of grid pts. */
    *iofst = *iofst + 32;
    gbit(cgrib, &ligds[2], *iofst, 8); /* Get num octets for opt. list */
    *iofst = *iofst + 8;
    gbit(cgrib, &ligds[3], *iofst, 8); /* Get interpret. for opt. list */
    *iofst = *iofst + 8;
    gbit(cgrib, &ligds[4], *iofst, 16); /* Get Grid Def Template num. */
    *iofst = *iofst + 16;

    if (ligds[4] != 65535)
    {
        /*   Get Grid Definition Template */
        if (!(mapgrid = getgridtemplate(ligds[4])))
        { /* undefined template */
            free(ligds);
            return G2_UNPACK3_BAD_GDT;
        }
        *mapgridlen = mapgrid->maplen;

        /*   Unpack each value into array igdstmpl from the the
         *   appropriate number of octets, which are specified in
         *   corresponding entries in array mapgrid. */
        if (*mapgridlen > 0)
        {
            if (!(ligdstmpl = calloc(*mapgridlen, sizeof(g2int))))
            {
                *mapgridlen = 0;
                *igdstmpl = NULL;
                if (mapgrid)
                    free(mapgrid);
                return G2_UNPACK_NO_MEM;
            }
            *igdstmpl = ligdstmpl;
        }
        ibyttem = 0;
        for (i = 0; i < *mapgridlen; i++)
        {
            nbits = abs(mapgrid->map[i]) * 8;
            if (mapgrid->map[i] >= 0)
            {
                gbit(cgrib, ligdstmpl + i, *iofst, nbits);
            }
            else
            {
                gbit(cgrib, &isign, *iofst, 1);
                gbit(cgrib, ligdstmpl + i, *iofst + 1, nbits - 1);
                if (isign == 1)
                    ligdstmpl[i] = -1 * ligdstmpl[i];
            }
            *iofst = *iofst + nbits;
            ibyttem = ibyttem + abs(mapgrid->map[i]);
        }

        /*   Check to see if the Grid Definition Template needs to be
         *   extended. The number of values in a specific template may
         *   vary depending on data specified in the "static" part of
         *   the gtemplate. */
        if (mapgrid->needext == 1)
        {
            free(mapgrid);
            mapgrid = extgridtemplate(ligds[4], ligdstmpl);
            /*   Unpack the rest of the Grid Definition Template */
            newlen = mapgrid->maplen + mapgrid->extlen;
            ligdstmpl = realloc(ligdstmpl, newlen * sizeof(g2int));
            *igdstmpl = ligdstmpl;
            j = 0;
            for (i = *mapgridlen; i < newlen; i++)
            {
                nbits = abs(mapgrid->ext[j]) * 8;
                if (mapgrid->ext[j] >= 0)
                {
                    gbit(cgrib, ligdstmpl + i, *iofst, nbits);
                }
                else
                {
                    gbit(cgrib, &isign, *iofst, 1);
                    gbit(cgrib, ligdstmpl + i, *iofst + 1, nbits - 1);
                    if (isign == 1)
                        ligdstmpl[i] = -1 * ligdstmpl[i];
                }
                *iofst = *iofst + nbits;
                ibyttem = ibyttem + abs(mapgrid->ext[j]);
                j++;
            }
            *mapgridlen = newlen;
        }
        if (mapgrid->ext)
            free(mapgrid->ext);
        if (mapgrid)
            free(mapgrid);
    }
    else
    { /* No Grid Definition Template */
        *mapgridlen = 0;
        *igdstmpl = 0;
    }

    /* Unpack optional list of numbers defining number of points in
     * each row or column, if included.  This is used for non
     * regular grids. */
    if (ligds[2] != 0)
    {
        nbits = ligds[2] * 8;
        *idefnum = (lensec - 14 - ibyttem) / ligds[2];
        if (*idefnum > 0)
            lideflist = calloc(*idefnum, sizeof(g2int));
        if (!lideflist)
        {
            *idefnum = 0;
            *ideflist = NULL;
            return G2_UNPACK_NO_MEM;
        }
        *ideflist = lideflist;
        gbits(cgrib, lideflist, *iofst, nbits, 0, *idefnum);
        *iofst = *iofst + (nbits * *idefnum);
    }
    else
    {
        *idefnum = 0;
        *ideflist = NULL;
    }

    return G2_NO_ERROR; /* End of Section 3 processing */
}
