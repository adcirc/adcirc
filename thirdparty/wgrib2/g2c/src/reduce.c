/**
 * @file
 * @brief Determines whether the number of groups should be increased
 * in order to reduce the size of the large groups, and to make that
 * adjustment.
 *
 * @author Harry Glahn, The Meteorological Development Lab (MDL)
 * @date November 2001
 */

#include "grib2_int.h"
#include <stdlib.h>

/**
 * Determines whether the number of groups should be increased in
 * order to reduce the size of the large groups, and to make that
 * adjustment. By reducing the size of the large groups, less bits
 * may be necessary for packing the group sizes and all the
 * information about the groups.
 *
 * The reference for nov was removed in the calling routine so that
 * kbit could be determined. This furnishes a starting point for the
 * iterations in reduce. However, the reference must be considered.
 *
 * PROGRAM HISTORY LOG
 * - November 2001 Harry Glahn, The Meteorological Development Lab
 *    (MDL) tdl grib2
 * - March 2002 Harry Glahn, The Meteorological Development Lab (MDL)
 *   comment ier = 715
 * - March 2002 Harry Glahn, The Meteorological Development Lab (MDL)
 *   modified to accommodate lx=1 on entry
 *
 * DATA SET USE
 * - kfildo - unit number for output (print) file. (output)
 *
 * @param kfildo unit number for output (print) file. (input)
 * @param jmin the minimum of each group (j=1,lx). it is possible
 * after splitting the groups, jmin() will not be the minimum of the
 * new group.  this doesn't matter; jmin() is really the group
 * reference and doesn't have to be the smallest value. (input/output)
 * @param jmax the maximum of each group (j=1,lx). (input/output)
 * @param lbit the number of bits necessary to pack each group
 * (j=1,lx). (input/output)
 * @param nov the number of values in each group (j=1,lx).
 * (input/output)
 * @param lx the number of groups. this will be increased if groups
 * are split. (input/output)
 * @param ndg the dimension of jmin(), jmax(), lbit(), and nov(
 *). (input)
 * @param ibit the number of bits necessary to pack the jmin(j)
 * values, j=1,lx. (input)
 * @param jbit the number of bits necessary to pack the lbit(j)
 * values, j=1,lx. (input)
 * @param kbit the number of bits necessary to pack the nov(j) values,
 * j=1,lx. if the groups are split, kbit is reduced. (input/output)
 * @param novref reference value for nov(). (input)
 * @param ibxx2 2**j (j=0,30). (input)
 * @param ier error return. (output)
 * - 0 = good return.
 * - 714 = problem in algorithm. reduce aborted.
 * - 715 = ngp not large enough. reduce aborted.
 *
 * @return always returns 0, see parameter ier for error code.
 *
 * <pre>
 * ntotbt(j) = the total bits used for the packing bits j
 *                       (j=1,30). (internal)
 *  nboxj(j) = new boxes needed for the packing bits j
 *                       (j=1,30). (internal)
 *           newbox(l) = number of new boxes (groups) for each original
 *                       group (l=1,lx) for the current j. (automatic)
 *                       (internal)
 *          newboxp(l) = same as newbox() but for the previous j.
 *                       this eliminates recomputation. (automatic)
 *                       (internal)
 *               cfeed = contains the character representation
 *                       of a printer form feed. (character) (internal)
 *               ifeed = contains the g2int value of a printer
 *                       form feed. (internal)
 *              iorigb = the original number of bits necessary
 *                       for the group values. (internal)
 * </pre>
 *
 * @author Harry Glahn, The Meteorological Development Lab (MDL)
 * @date November 2001
 */
int
reduce(g2int *kfildo, g2int *jmin, g2int *jmax,
       g2int *lbit, g2int *nov, g2int *lx, g2int *ndg, g2int *ibit,
       g2int *jbit, g2int *kbit, g2int *novref, g2int *ibxx2,
       g2int *ier)
{
    /* Initialized data */

    static g2int ifeed = 12;

    /* System generated locals */
    g2int i__1, i__2;

    /* Local variables */
    static g2int newboxtp, j, l, m, jj, lxn, left;
    static float pimp;
    static g2int move, novl;
    static char cfeed[1];
    /* static g2int nboxj[31]; */
    static g2int lxnkp, iorigb, ibxx2m1, movmin,
        ntotbt[31], ntotpr, newboxt;
    g2int *newbox, *newboxp;

    /*        NEWBOX() AND NEWBOXP() were AUTOMATIC ARRAYS. */
    newbox = (g2int *)calloc(*ndg, sizeof(g2int));
    newboxp = (g2int *)calloc(*ndg, sizeof(g2int));

    /* Parameter adjustments */
    --nov;
    --lbit;
    --jmax;
    --jmin;

    /* Function Body */

    *ier = 0;
    if (*lx == 1)
    {
        goto L410;
    }
    /*        IF THERE IS ONLY ONE GROUP, RETURN. */

    *(unsigned char *)cfeed = (char)ifeed;

    /*        INITIALIZE NUMBER OF NEW BOXES PER GROUP TO ZERO. */

    i__1 = *lx;
    for (l = 1; l <= i__1; ++l)
    {
        newbox[l - 1] = 0;
        /* L110: */
    }

    /*        INITIALIZE NUMBER OF TOTAL NEW BOXES PER J TO ZERO. */

    for (j = 1; j <= 31; ++j)
    {
        ntotbt[j - 1] = 999999999;
        /* nboxj[j - 1] = 0; */
        /* L112: */
    }

    iorigb = (*ibit + *jbit + *kbit) * *lx;
    /*        IBIT = BITS TO PACK THE JMIN(). */
    /*        JBIT = BITS TO PACK THE LBIT(). */
    /*        KBIT = BITS TO PACK THE NOV(). */
    /*        LX = NUMBER OF GROUPS. */
    ntotbt[*kbit - 1] = iorigb;
    /*           THIS IS THE VALUE OF TOTAL BITS FOR THE ORIGINAL LX */
    /*           GROUPS, WHICH REQUIRES KBITS TO PACK THE GROUP */
    /*           LENGHTS.  SETTING THIS HERE MAKES ONE LESS LOOPS */
    /*           NECESSARY BELOW. */

    /*        COMPUTE BITS NOW USED FOR THE PARAMETERS DEFINED. */

    /*        DETERMINE OTHER POSSIBILITES BY INCREASING LX AND DECREASING */
    /*        NOV() WITH VALUES GREATER THAN THRESHOLDS.  ASSUME A GROUP IS */
    /*        SPLIT INTO 2 OR MORE GROUPS SO THAT KBIT IS REDUCED WITHOUT */
    /*        CHANGING IBIT OR JBIT. */

    jj = 0;

    /* Computing MIN */
    i__1 = 30, i__2 = *kbit - 1;
    /*for (j = min(i__1,i__2); j >= 2; --j) {*/
    for (j = (i__1 < i__2) ? i__1 : i__2; j >= 2; --j)
    {
        /*           VALUES GE KBIT WILL NOT REQUIRE SPLITS.  ONCE THE TOTAL */
        /*           BITS START INCREASING WITH DECREASING J, STOP.  ALSO, THE */
        /*           NUMBER OF BITS REQUIRED IS KNOWN FOR KBITS = NTOTBT(KBIT). */

        newboxt = 0;

        i__1 = *lx;
        for (l = 1; l <= i__1; ++l)
        {

            if (nov[l] < ibxx2[j])
            {
                newbox[l - 1] = 0;
                /*                 NO SPLITS OR NEW BOXES. */
                goto L190;
            }
            else
            {
                novl = nov[l];

                m = (nov[l] - 1) / (ibxx2[j] - 1) + 1;
                /*                 M IS FOUND BY SOLVING THE EQUATION BELOW FOR M: */
                /*                 (NOV(L)+M-1)/M LT IBXX2(J) */
                /*                 M GT (NOV(L)-1)/(IBXX2(J)-1) */
                /*                 SET M = (NOV(L)-1)/(IBXX2(J)-1)+1 */
            L130:
                novl = (nov[l] + m - 1) / m;
                /*                 THE +M-1 IS NECESSARY.  FOR INSTANCE, 15 WILL FIT */
                /*                 INTO A BOX 4 BITS WIDE, BUT WON'T DIVIDE INTO */
                /*                 TWO BOXES 3 BITS WIDE EACH. */

                if (novl < ibxx2[j])
                {
                    goto L185;
                }
                else
                {
                    ++m;
                    /* ***                  WRITE(KFILDO,135)L,NOV(L),NOVL,M,J,IBXX2(J) */
                    /* *** 135              FORMAT(/' AT 135--L,NOV(L),NOVL,M,J,IBXX2(J)',6I10) */
                    goto L130;
                }

                /*                 THE ABOVE DO LOOP WILL NEVER COMPLETE. */
            }

        L185:
            newbox[l - 1] = m - 1;
            newboxt = newboxt + m - 1;
        L190:;
        }

        /* nboxj[j - 1] = newboxt; */
        ntotpr = ntotbt[j];
        ntotbt[j - 1] = (*ibit + *jbit) * (*lx + newboxt) + j * (*lx +
                                                                 newboxt);

        if (ntotbt[j - 1] >= ntotpr)
        {
            jj = j + 1;
            /*              THE PLUS IS USED BECAUSE J DECREASES PER ITERATION. */
            goto L250;
        }
        else
        {

            /*              SAVE THE TOTAL NEW BOXES AND NEWBOX() IN CASE THIS */
            /*              IS THE J TO USE. */

            newboxtp = newboxt;

            i__1 = *lx;
            for (l = 1; l <= i__1; ++l)
            {
                newboxp[l - 1] = newbox[l - 1];
                /* L195: */
            }

            /*           WRITE(KFILDO,197)NEWBOXT,IBXX2(J) */
            /* 197        FORMAT(/' *****************************************' */
            /*    1             /' THE NUMBER OF NEWBOXES PER GROUP OF THE TOTAL', */
            /*    2              I10,' FOR GROUP MAXSIZE PLUS 1 ='I10 */
            /*    3             /' *****************************************') */
            /*           WRITE(KFILDO,198) (NEWBOX(L),L=1,LX) */
            /* 198        FORMAT(/' '20I6/(' '20I6)) */
        }

        /* 205     WRITE(KFILDO,209)KBIT,IORIGB */
        /* 209     FORMAT(/' ORIGINAL BITS WITH KBIT OF',I5,' =',I10) */
        /*        WRITE(KFILDO,210)(N,N=2,10),(IBXX2(N),N=2,10), */
        /*    1                    (NTOTBT(N),N=2,10),(NBOXJ(N),N=2,10), */
        /*    2                    (N,N=11,20),(IBXX2(N),N=11,20), */
        /*    3                    (NTOTBT(N),N=11,20),(NBOXJ(N),N=11,20), */
        /*    4                    (N,N=21,30),(IBXX2(N),N=11,20), */
        /*    5                    (NTOTBT(N),N=21,30),(NBOXJ(N),N=21,30) */
        /* 210     FORMAT(/' THE TOTAL BYTES FOR MAXIMUM GROUP LENGTHS BY ROW'// */
        /*    1      '   J         = THE NUMBER OF BITS PER GROUP LENGTH'/ */
        /*    2      '   IBXX2(J)  = THE MAXIMUM GROUP LENGTH PLUS 1 FOR THIS J'/ */
        /*    3      '   NTOTBT(J) = THE TOTAL BITS FOR THIS J'/ */
        /*    4      '   NBOXJ(J)  = THE NEW GROUPS FOR THIS J'/ */
        /*    5      4(/10X,9I10)/4(/10I10)/4(/10I10)) */

        /* L200: */
    }

L250:
    pimp = (iorigb - ntotbt[jj - 1]) / (float)iorigb * 100.f;
    /*     WRITE(KFILDO,252)PIMP,KBIT,JJ */
    /* 252  FORMAT(/' PERCENT IMPROVEMENT =',F6.1, */
    /*    1        ' BY DECREASING GROUP LENGTHS FROM',I4,' TO',I4,' BITS') */
    if (pimp >= 2.f)
    {

        /*        WRITE(KFILDO,255)CFEED,NEWBOXTP,IBXX2(JJ) */
        /* 255     FORMAT(A1,/' *****************************************' */
        /*    1             /' THE NUMBER OF NEWBOXES PER GROUP OF THE TOTAL', */
        /*    2             I10,' FOR GROUP MAXSIZE PLUS 1 ='I10 */
        /*    2             /' *****************************************') */
        /*        WRITE(KFILDO,256) (NEWBOXP(L),L=1,LX) */
        /* 256     FORMAT(/' '20I6) */

        /*           ADJUST GROUP LENGTHS FOR MAXIMUM LENGTH OF JJ BITS. */
        /*           THE MIN PER GROUP AND THE NUMBER OF BITS REQUIRED */
        /*           PER GROUP ARE NOT CHANGED.  THIS MAY MEAN THAT A */
        /*           GROUP HAS A MIN (OR REFERENCE) THAT IS NOT ZERO. */
        /*           THIS SHOULD NOT MATTER TO THE UNPACKER. */

        lxnkp = *lx + newboxtp;
        /*           LXNKP = THE NEW NUMBER OF BOXES */

        if (lxnkp > *ndg)
        {
            /*              DIMENSIONS NOT LARGE ENOUGH.  PROBABLY AN ERROR */
            /*              OF SOME SORT.  ABORT. */
            /*           WRITE(KFILDO,257)NDG,LXNPK */
            /*        1         2         3         4         5         6         7 X */
            /* 257        FORMAT(/' DIMENSIONS OF JMIN, ETC. IN REDUCE =',I8, */
            /*    1              ' NOT LARGE ENOUGH FOR THE EXPANDED NUMBER OF', */
            /*    2              ' GROUPS =',I8,'.  ABORT REDUCE.') */
            *ier = 715;
            goto L410;
            /*              AN ABORT CAUSES THE CALLING PROGRAM TO REEXECUTE */
            /*              WITHOUT CALLING REDUCE. */
        }

        lxn = lxnkp;
        /*           LXN IS THE NUMBER OF THE BOX IN THE NEW SERIES BEING */
        /*           FILLED.  IT DECREASES PER ITERATION. */
        ibxx2m1 = ibxx2[jj] - 1;
        /*           IBXX2M1 IS THE MAXIMUM NUMBER OF VALUES PER GROUP. */

        for (l = *lx; l >= 1; --l)
        {

            /*              THE VALUES IS NOV() REPRESENT THOSE VALUES + NOVREF. */
            /*              WHEN VALUES ARE MOVED TO ANOTHER BOX, EACH VALUE */
            /*              MOVED TO A NEW BOX REPRESENTS THAT VALUE + NOVREF. */
            /*              THIS HAS TO BE CONSIDERED IN MOVING VALUES. */

            if (newboxp[l - 1] * (ibxx2m1 + *novref) + *novref > nov[l] + *novref)
            {
                /*                 IF THE ABOVE TEST IS MET, THEN MOVING IBXX2M1 VALUES */
                /*                 FOR ALL NEW BOXES WILL LEAVE A NEGATIVE NUMBER FOR */
                /*                 THE LAST BOX.  NOT A TOLERABLE SITUATION. */
                movmin = (nov[l] - newboxp[l - 1] * *novref) / newboxp[l - 1];
                left = nov[l];
                /*                 LEFT = THE NUMBER OF VALUES TO MOVE FROM THE ORIGINAL */
                /*                 BOX TO EACH NEW BOX EXCEPT THE LAST.  LEFT IS THE */
                /*                 NUMBER LEFT TO MOVE. */
            }
            else
            {
                movmin = ibxx2m1;
                /*                 MOVMIN VALUES CAN BE MOVED FOR EACH NEW BOX. */
                left = nov[l];
                /*                 LEFT IS THE NUMBER OF VALUES LEFT TO MOVE. */
            }

            if (newboxp[l - 1] > 0)
            {
                if ((movmin + *novref) * newboxp[l - 1] + *novref <= nov[l] +
                                                                         *novref &&
                    (movmin + *novref) * (newboxp[l - 1] + 1) >= nov[l] + *novref)
                {
                    goto L288;
                }
                else
                {
                    /* ***D                 WRITE(KFILDO,287)L,MOVMIN,NOVREF,NEWBOXP(L),NOV(L) */
                    /* ***D287              FORMAT(/' AT 287 IN REDUCE--L,MOVMIN,NOVREF,', */
                    /* ***D    1                    'NEWBOXP(L),NOV(L)',5I12 */
                    /* ***D    2                    ' REDUCE ABORTED.') */
                    /*              WRITE(KFILDO,2870) */
                    /* 2870          FORMAT(/' AN ERROR IN REDUCE ALGORITHM.  ABORT REDUCE.') */
                    *ier = 714;
                    goto L410;
                    /*                 AN ABORT CAUSES THE CALLING PROGRAM TO REEXECUTE */
                    /*                 WITHOUT CALLING REDUCE. */
                }
            }

        L288:
            i__1 = newboxp[l - 1] + 1;
            for (j = 1; j <= i__1; ++j)
            {
                /*move = min(movmin,left);*/
                move = (movmin < left) ? movmin : left;
                jmin[lxn] = jmin[l];
                jmax[lxn] = jmax[l];
                lbit[lxn] = lbit[l];
                nov[lxn] = move;
                --lxn;
                left -= move + *novref;
                /*                 THE MOVE OF MOVE VALUES REALLY REPRESENTS A MOVE OF */
                /*                 MOVE + NOVREF VALUES. */
                /* L290: */
            }

            if (left != -(*novref))
            {
                /* ***               WRITE(KFILDO,292)L,LXN,MOVE,LXNKP,IBXX2(JJ),LEFT,NOV(L), */
                /* ***     1                          MOVMIN */
                /* *** 292           FORMAT(' AT 292 IN REDUCE--L,LXN,MOVE,LXNKP,', */
                /* ***     1                'IBXX2(JJ),LEFT,NOV(L),MOVMIN'/8I12) */
            }

            /* L300: */
        }

        *lx = lxnkp;
        /*           LX IS NOW THE NEW NUMBER OF GROUPS. */
        *kbit = jj;
        /*           KBIT IS NOW THE NEW NUMBER OF BITS REQUIRED FOR PACKING */
        /*           GROUP LENGHTS. */
    }

    /*     WRITE(KFILDO,406)CFEED,LX */
    /* 406  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP SIZES NOV() AFTER REDUCTION IN SIZE', */
    /*    2           ' FOR'I10,' GROUPS', */
    /*    3          /' *****************************************') */
    /*     WRITE(KFILDO,407) (NOV(J),J=1,LX) */
    /* 407  FORMAT(/' '20I6) */
    /*     WRITE(KFILDO,408)CFEED,LX */
    /* 408  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP MINIMA JMIN() AFTER REDUCTION IN SIZE', */
    /*    2           ' FOR'I10,' GROUPS', */
    /*    3          /' *****************************************') */
    /*     WRITE(KFILDO,409) (JMIN(J),J=1,LX) */
    /* 409  FORMAT(/' '20I6) */

L410:
    if (newbox)
        free(newbox);
    if (newboxp)
        free(newboxp);
    return 0;
} /* reduce_ */
