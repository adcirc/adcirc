/** @file
 * @author Glahn @date february 1994
 */

/*#include "f2c.h"*/
#include "grib2_int.h"
#include <stdlib.h>

#define TRUE_ (1)  /**< True. */
#define FALSE_ (0) /**< False. */

/**
 * Determines groups of variable size, but at least of size minpk, the
 * associated max (jmax( )) and min (jmin( )), the number of bits
 * necessary to hold the values in each group (lbit( )), the number of
 * values in each group (nov( )), the number of bits necessary to pack
 * the jmin( ) values (ibit), the number of bits necessary to pack the
 * lbit( ) values (jbit), and the number of bits necessary to pack the
 * nov( ) values (kbit). The routine is designed to determine the
 * groups such that a small number of bits is necessary to pack the
 * data without excessive computations. If all values in the group are
 * zero, the number of bits to use in packing is defined as zero when
 * there can be no missing values; when there can be missing values,
 * the number of bits must be at least 1 to have the capability to
 * recognize the missing value. However, if all values in a group are
 * missing, the number of bits needed is 0, and the unpacker
 * recognizes this. All variables are g2int. Even though the groups
 * are initially of size minpk or larger, an adjustment between two
 * groups (the lookback procedure) may make a group smaller than
 * minpk. The control on group size is that the sum of the sizes of
 * the two consecutive groups, each of size minpk or larger, is not
 * decreased. When determining the number of bits necessary for
 * packing, the largest value that can be accommodated in, say, mbits,
 * is 2**mbits-1; this largest value (and the next smallest value) is
 * reserved for the missing value indicator (only) when is523 ne 0. If
 * the dimension ndg is not large enough to hold all the groups, the
 * local value of minpk is increased by 50 percent. This is repeated
 * until ndg will suffice. A diagnostic is printed whenever this
 * happens, which should be very rarely. If it happens often, ndg in
 * subroutine pack should be increased and a corresponding increase in
 * subroutine unpack made. Considerable code is provided so that no
 * more checking for missing values within loops is done than
 * necessary; the added efficiency of this is relatively minor, but
 * does no harm. For grib2, the reference value for the length of
 * groups in nov( ) and for the number of bits necessary to pack group
 * values are determined, and subtracted before jbit and kbit are
 * determined.
 *
 * When 1 or more groups are large compared to the others, the width
 * of all groups must be as large as the largest. A subroutine reduce
 * breaks up large groups into 2 or more to reduce total bits
 * required. If reduce should abort, pack_gp will be executed again
 * without the call to reduce.
 *
 * PROGRAM HISTORY LOG:
 * - February 1994 Glahn tdl mos-2000
 * - June 1995 Glahn modified for lmiss error.
 * - July 1996 Glahn added misss
 * - February 1997 Glahn removed 4 redundant tests for missp.eq.0;
 * inserted a test to better handle a string of 9999's
 * - February 1997 Glahn added loops to eliminate test for misss when
 * misss = 0
 * - March 1997 Glahn corrected for secondary missing value
 * - March 1997 Glahn corrected for use of local value of minpk
 * - March 1997 Glahn corrected for secondary missing value
 * - March 1997 Glahn changed calculating number of bits through
 * exponents to an array (improved overall packing performance by
 * about 35 percent!). allowed 0 bits for packing jmin( ), lbit( ),
 * and nov( ).
 * - May 1997 Glahn a number of changes for efficiency.  mod functions
 * eliminated and one ifthen added. jount removed.  recomputation of
 * bits not made unless necessary after moving points from one group
 * to another. nendb adjusted to eliminate possibility of very small
 * group at the end.  about 8 percent improvement in overall
 * packing. iskipa removed; there is always a group b that can become
 * group a. control on size of group b (statement below 150)
 * added. added adda, and use of ge and le instead of gt and lt in
 * loops between 150 and 160.  ibitbs added to shorten trips through
 * loop.
 * - March 2000 Glahn modified for grib2; changed name from packgp
 * - january 2001 Glahn comments; ier = 706 substituted for stops;
 * added return1; removed statement number 110; added ier and * return
 * - November 2001 Glahn changed some diagnostic formats to allow
 * printing larger numbers
 * - November 2001 Glahn added misslx( ) to put maximum value into
 * jmin( ) when all values missing to agree with grib standard.
 * - November 2001 Glahn changed two tests on missp and misss eq 0 to
 * tests on is523. however, missp and misss cannot in general be = 0.
 * - November 2001 Glahn added call to reduce; defined itest before
 * loops to reduce computation; started large group when all same
 * value
 * - December 2001 Glahn modified and added a few comments
 * - January 2002 Glahn removed loop before 150 to determine a group
 * of all same value
 * - January 2002 Glahn changed mallow from 9999999 to 2**30+1, and
 * made it a parameter
 * - March 2002 Glahn added non fatal ier = 716, 717; removed
 * nendb=nxy above 150; added iersav=0; comments
 *
 * DATA SET USE
 * - kfildo - unit number for output (print) file. (output)
 *
 * @param kfildo unit number for output (print) file. (input)
 * @param ic array to hold data for packing. the values do not have to
 * be positive at this point, but must be in the range -2**30 to
 * +2**30 (the the value of mallow). these g2int values will be
 * retained exactly through packing and unpacking. (input)
 * @param nxy number of values in ic( ). also treated as its
 * dimension. (input)
 * @param is523 missing value management 0=data contains no missing
 * values 1=data contains primary missing values 2=data contains
 * primary and secondary missing values (input)
 * @param minpk the minimum size of each group, except possibly the
 * last one. (input)
 * @param inc the number of values to add to an already existing group
 * in determining whether or not to start a new group. ideally, this
 * would be 1, but each time inc values are attempted, the max and min
 * of the next minpk values must be found. this is "a loop within a
 * loop," and a slightly larger value may give about as good results
 * with slightly less computational time.  if inc is le 0, 1 is used,
 * and a diagnostic is output. note: it is expected that inc will
 * equal 1. the code uses inc primarily in the loops starting at
 * statement 180. if inc were 1, there would not need to be loops as
 * such. however, kinc (the local value of inc) is set ge 1 when near
 * the end of the data to forestall a very small group at the end.
 * (input)
 * @param missp when missing points can be present in the data, they
 * will have the value missp or misss.  missp is the primary missing
 * value and misss is the secondary missing value . these must not be
 * values that would occur with subtracting the minimum (reference)
 * value or scaling.  for example, missp = 0 would not be advisable.
 * (input)
 * @param misss secondary missing value indicator (see missp).
 * (input)
 * @param jmin the minimum of each group (j=1,lx). (output)
 * @param jmax the maximum of each group (j=1,lx). this is not really
 * needed, but since the max of each group must be found, saving it
 * here is cheap in case the user wants it. (output)
 * @param lbit the number of bits necessary to pack each group
 * (j=1,lx). it is assumed the minimum of each group will be removed
 * before packing, and the values to pack will, therefore, all be
 * positive.  however, ic( ) does not necessarily contain all positive
 * values. if the overall minimum has been removed (the usual case),
 * then ic( ) will contain only positive values. (output)
 * @param nov the number of values in each group (j=1,lx). (output)
 * @param ndg the dimension of jmin, jmax, lbit, and nov. (input)
 * @param lx the number of groups determined. (output)
 * @param ibit the number of bits necessary to pack the jmin(j)
 * values, j=1,lx. (output)
 * @param jbit the number of bits necessary to pack the lbit(j)
 * values, j=1,lx. (output)
 * @param kbit the number of bits necessary to pack the nov(j) values,
 * j=1,lx. (output)
 * @param novref reference value for nov( ). (output)
 * @param lbitref reference value for lbit( ). (output)
 * @param ier Error code
 * - 0 No error.
 * - 706 value will not pack in 30 bits--fatal
 * - 714 error in reduce--non-fatal
 * - 715 ngp not large enough in reduce--non-fatal
 * - 716 minpk inceased--non-fatal
 * - 717 inc set
 * - 1--non-fatal
 * * alternate return when ier ne 0 and fatal error.
 *
 * @return 0 - check ier for error code.
 *
 *        INTERNAL VARIABLES
 * <pre>
 *               cfeed = contains the character representation
 *                       of a printer form feed.
 *               ifeed = contains the g2int value of a printer
 *                       form feed.
 *                kinc = working copy of inc. may be modified.
 *                mina = minimum value in group a.
 *                maxa = maximum value in group a.
 *               nenda = the place in ic( ) where group a ends.
 *              kstart = the place in ic( ) where group a starts.
 *               ibita = number of bits needed to hold values in group a.
 *                minb = minimum value in group b.
 *                maxb = maximum value in group b.
 *               nendb = the place in ic( ) where group b ends.
 *               ibitb = number of bits needed to hold values in group b.
 *                minc = minimum value in group c.
 *                maxc = maximum value in group c.
 *              ktotal = count of number of values in ic( ) processed.
 *               nount = number of values added to group a.
 *               lmiss = 0 when is523 = 0. when packing into a
 *                       specific number of bits, say mbits,
 *                       the maximum value that can be handled is
 *                       2**mbits-1. when is523 = 1, indicating
 *                       primary missing values, this maximum value
 *                       is reserved to hold the primary missing value
 *                       indicator and lmiss = 1. when is523 = 2,
 *                       the value just below the maximum i.e.,
 *                       2**mbits-2 is reserved to hold the secondary
 *                       missing value indicator and lmiss = 2.
 *              lminpk = local value of minpk. this will be adjusted
 *                       upward whenever ndg is not large enough to hold
 *                       all the groups.
 *              mallow = the largest allowable value for packing.
 *              mislla = set to 1 when all values in group a are missing.
 *                       this is used to distinguish between a real
 *                       minimum when all values are not missing
 *                       and a minimum that has been set to zero when
 *                       all values are missing. 0 otherwise.
 *                       note that this does not distinguish between
 *                       primary and secondary missings when secondary
 *                       missings are present. this means that
 *                       lbit( ) will not be zero with the resulting
 *                       compression efficiency when secondary missings
 *                       are present. also note that a check has been
 *                       made earlier to determine that secondary
 *                       missings are really there.
 *              misllb = set to 1 when all values in group b are missing.
 *                       this is used to distinguish between a real
 *                       minimum when all values are not missing
 *                       and a minimum that has been set to zero when
 *                       all values are missing. 0 otherwise.
 *              misllc = performs the same function for group c that
 *                       mislla and misllb do for groups b and c,
 *                       respectively.
 *            ibxx2(j) = an array that when this routine is first entered
 *                       is set to 2**j, j=0,30. ibxx2(30) = 2**30, which
 *                       is the largest value packable, because 2**31
 *                       is larger than the g2int word size.
 *              ifirst = set by data statement to 0. changed to 1 on
 *                       first
 *                       entry when ibxx2( ) is filled.
 *               minak = keeps track of the location in ic( ) where the
 *                       minimum value in group a is located.
 *               maxak = does the same as minak, except for the maximum.
 *               minbk = the same as minak for group b.
 *               maxbk = the same as maxak for group b.
 *               minck = the same as minak for group c.
 *               maxck = the same as maxak for group c.
 *                adda = keeps track whether or not an attempt to add
 *                       points to group a was made. if so, then adda
 *                       keeps from trying to put one back into b.
 *                       (g2int)
 *              ibitbs = keeps current value if ibitb so that loop
 *                       ending at 166 doesn't have to start at
 *                       ibitb = 0 every time.
 *           misslx(j) = mallow except when a group is all one value (and
 *                       lbit(j) = 0) and that value is missing. in
 *                       that case, misslx(j) is missp or misss. this
 *                       gets inserted into jmin(j) later as the
 *                       missing indicator; it can't be put in until
 *                       the end, because jmin( ) is used to calculate
 *                       the maximum number of bits (ibits) needed to
 *                       pack jmin( ).
 * </pre>
 */
int
pack_gp(g2int *kfildo, g2int *ic, g2int *nxy,
        g2int *is523, g2int *minpk, g2int *inc, g2int *missp, g2int *misss,
        g2int *jmin, g2int *jmax, g2int *lbit, g2int *nov,
        g2int *ndg, g2int *lx, g2int *ibit, g2int *jbit, g2int *kbit,
        g2int *novref, g2int *lbitref, g2int *ier)
{
    /* Initialized data */

    const g2int mallow = 1073741825; /*  MALLOW=2**30+1  */
    static g2int ifeed = 12;
    static g2int ifirst = 0;

    /* System generated locals */
    g2int i__1, i__2, i__3;

    /* Local variables */
    static g2int j, k, l;
    static g2int adda;
    static g2int ired, kinc, mina, maxa, minb, maxb, minc, maxc, ibxx2[31];
    static char cfeed[1];
    static g2int nenda, nendb, ibita, ibitb, minak, minbk, maxak, maxbk,
        minck, maxck, nouta, lmiss, itest, nount;
    extern /* Subroutine */ int reduce(g2int *, g2int *, g2int *,
                                       g2int *, g2int *, g2int *, g2int *, g2int *, g2int *,
                                       g2int *, g2int *, g2int *, g2int *);
    static g2int ibitbs, mislla, misllb, misllc, iersav, lminpk, ktotal,
        kounta, kountb, kstart, mstart, mintst, maxtst,
        kounts, mintstk, maxtstk;
    g2int *misslx;

    /*        NON SYSTEM SUBROUTINES CALLED */
    /*           NONE */

    /*        MISSLX( ) was AN AUTOMATIC ARRAY. */
    misslx = (g2int *)calloc(*ndg, sizeof(g2int));

    /* Parameter adjustments */
    --ic;
    --nov;
    --lbit;
    --jmax;
    --jmin;

    /* Function Body */

    *ier = 0;
    iersav = 0;
    /*     CALL TIMPR(KFILDO,KFILDO,'START PACK_GP        ') */
    *(unsigned char *)cfeed = (char)ifeed;

    ired = 0;
    /*        IRED IS A FLAG.  WHEN ZERO, REDUCE WILL BE CALLED. */
    /*        IF REDUCE ABORTS, IRED = 1 AND IS NOT CALLED.  IN */
    /*        THIS CASE PACK_GP EXECUTES AGAIN EXCEPT FOR REDUCE. */

    if (*inc <= 0)
    {
        iersav = 717;
        /*        WRITE(KFILDO,101)INC */
        /* 101     FORMAT(/' ****INC ='I8,' NOT CORRECT IN PACK_GP.  1 IS USED.') */
    }

    /*        THERE WILL BE A RESTART OF PACK_GP IF SUBROUTINE REDUCE */
    /*        ABORTS.  THIS SHOULD NOT HAPPEN, BUT IF IT DOES, PACK_GP */
    /*        WILL COMPLETE WITHOUT SUBROUTINE REDUCE.  A NON FATAL */
    /*        DIAGNOSTIC RETURN IS PROVIDED. */

L102:
    /*kinc = max(*inc,1);*/
    kinc = (*inc > 1) ? *inc : 1;
    lminpk = *minpk;

    /*         CALCULATE THE POWERS OF 2 THE FIRST TIME ENTERED. */

    if (ifirst == 0)
    {
        ifirst = 1;
        ibxx2[0] = 1;

        for (j = 1; j <= 30; ++j)
        {
            ibxx2[j] = ibxx2[j - 1] << 1;
            /* L104: */
        }
    }

    /*        THERE WILL BE A RESTART AT 105 IS NDG IS NOT LARGE ENOUGH. */
    /*        A NON FATAL DIAGNOSTIC RETURN IS PROVIDED. */

L105:
    kstart = 1;
    ktotal = 0;
    *lx = 0;
    adda = FALSE_;
    lmiss = 0;
    if (*is523 == 1)
    {
        lmiss = 1;
    }
    if (*is523 == 2)
    {
        lmiss = 2;
    }

    /*        ************************************* */

    /*        THIS SECTION COMPUTES STATISTICS FOR GROUP A.  GROUP A IS */
    /*        A GROUP OF SIZE LMINPK. */

    /*        ************************************* */

    ibita = 0;
    mina = mallow;
    maxa = -mallow;
    minak = mallow;
    maxak = -mallow;

    /*        FIND THE MIN AND MAX OF GROUP A.  THIS WILL INITIALLY BE OF */
    /*        SIZE LMINPK (IF THERE ARE STILL LMINPK VALUES IN IC( )), BUT */
    /*        WILL INCREASE IN SIZE IN INCREMENTS OF INC UNTIL A NEW */
    /*        GROUP IS STARTED.  THE DEFINITION OF GROUP A IS DONE HERE */
    /*        ONLY ONCE (UPON INITIAL ENTRY), BECAUSE A GROUP B CAN ALWAYS */
    /*        BECOME A NEW GROUP A AFTER A IS PACKED, EXCEPT IF LMINPK */
    /*        HAS TO BE INCREASED BECAUSE NDG IS TOO SMALL.  THEREFORE, */
    /*        THE SEPARATE LOOPS FOR MISSING AND NON-MISSING HERE BUYS */
    /*        ALMOST NOTHING. */

    /* Computing MIN */
    i__1 = kstart + lminpk - 1;
    /*nenda = min(i__1,*nxy);*/
    nenda = (i__1 < *nxy) ? i__1 : *nxy;
    if (*nxy - nenda <= lminpk / 2)
    {
        nenda = *nxy;
    }
    /*        ABOVE STATEMENT GUARANTEES THE LAST GROUP IS GT LMINPK/2 BY */
    /*        MAKING THE ACTUAL GROUP LARGER.  IF A PROVISION LIKE THIS IS */
    /*        NOT INCLUDED, THERE WILL MANY TIMES BE A VERY SMALL GROUP */
    /*        AT THE END.  USE SEPARATE LOOPS FOR MISSING AND NO MISSING */
    /*        VALUES FOR EFFICIENCY. */

    /*        DETERMINE WHETHER THERE IS A LONG STRING OF THE SAME VALUE */
    /*        UNLESS NENDA = NXY.  THIS MAY ALLOW A LARGE GROUP A TO */
    /*        START WITH, AS WITH MISSING VALUES.   SEPARATE LOOPS FOR */
    /*        MISSING OPTIONS.  THIS SECTION IS ONLY EXECUTED ONCE, */
    /*        IN DETERMINING THE FIRST GROUP.  IT HELPS FOR AN ARRAY */
    /*        OF MOSTLY MISSING VALUES OR OF ONE VALUE, SUCH AS */
    /*        RADAR OR PRECIP DATA. */

    if (nenda != *nxy && ic[kstart] == ic[kstart + 1])
    {
        /*           NO NEED TO EXECUTE IF FIRST TWO VALUES ARE NOT EQUAL. */

        if (*is523 == 0)
        {
            /*              THIS LOOP IS FOR NO MISSING VALUES. */

            i__1 = *nxy;
            for (k = kstart + 1; k <= i__1; ++k)
            {

                if (ic[k] != ic[kstart])
                {
                    /* Computing MAX */
                    i__2 = nenda, i__3 = k - 1;
                    /*nenda = max(i__2,i__3);*/
                    nenda = (i__2 > i__3) ? i__2 : i__3;
                    goto L114;
                }

                /* L111: */
            }

            nenda = *nxy;
            /*              FALL THROUGH THE LOOP MEANS ALL VALUES ARE THE SAME. */
        }
        else if (*is523 == 1)
        {
            /*              THIS LOOP IS FOR PRIMARY MISSING VALUES ONLY. */

            i__1 = *nxy;
            for (k = kstart + 1; k <= i__1; ++k)
            {

                if (ic[k] != *missp)
                {

                    if (ic[k] != ic[kstart])
                    {
                        /* Computing MAX */
                        i__2 = nenda, i__3 = k - 1;
                        /*nenda = max(i__2,i__3);*/
                        nenda = (i__2 > i__3) ? i__2 : i__3;
                        goto L114;
                    }
                }

                /* L112: */
            }

            nenda = *nxy;
            /*              FALL THROUGH THE LOOP MEANS ALL VALUES ARE THE SAME. */
        }
        else
        {
            /*              THIS LOOP IS FOR PRIMARY AND SECONDARY MISSING VALUES. */

            i__1 = *nxy;
            for (k = kstart + 1; k <= i__1; ++k)
            {

                if (ic[k] != *missp && ic[k] != *misss)
                {

                    if (ic[k] != ic[kstart])
                    {
                        /* Computing MAX */
                        i__2 = nenda, i__3 = k - 1;
                        /*nenda = max(i__2,i__3);*/
                        nenda = (i__2 > i__3) ? i__2 : i__3;
                        goto L114;
                    }
                }

                /* L113: */
            }

            nenda = *nxy;
            /*              FALL THROUGH THE LOOP MEANS ALL VALUES ARE THE SAME. */
        }
    }

L114:
    if (*is523 == 0)
    {

        i__1 = nenda;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] < mina)
            {
                mina = ic[k];
                minak = k;
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
                maxak = k;
            }
            /* L115: */
        }
    }
    else if (*is523 == 1)
    {

        i__1 = nenda;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp)
            {
                goto L117;
            }
            if (ic[k] < mina)
            {
                mina = ic[k];
                minak = k;
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
                maxak = k;
            }
        L117:;
        }
    }
    else
    {

        i__1 = nenda;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp || ic[k] == *misss)
            {
                goto L120;
            }
            if (ic[k] < mina)
            {
                mina = ic[k];
                minak = k;
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
                maxak = k;
            }
        L120:;
        }
    }

    kounta = nenda - kstart + 1;

    /*        INCREMENT KTOTAL AND FIND THE BITS NEEDED TO PACK THE A GROUP. */

    ktotal += kounta;
    mislla = 0;
    if (mina != mallow)
    {
        goto L125;
    }
    /*        ALL MISSING VALUES MUST BE ACCOMMODATED. */
    mina = 0;
    maxa = 0;
    mislla = 1;
    ibitb = 0;
    if (*is523 != 2)
    {
        goto L130;
    }
    /*        WHEN ALL VALUES ARE MISSING AND THERE ARE NO */
    /*        SECONDARY MISSING VALUES, IBITA = 0. */
    /*        OTHERWISE, IBITA MUST BE CALCULATED. */

L125:
    itest = maxa - mina + lmiss;

    for (ibita = 0; ibita <= 30; ++ibita)
    {
        if (itest < ibxx2[ibita])
        {
            goto L130;
        }
        /* ***        THIS TEST IS THE SAME AS: */
        /* ***     IF(MAXA-MINA.LT.IBXX2(IBITA)-LMISS)GO TO 130 */
        /* L126: */
    }

    /*     WRITE(KFILDO,127)MAXA,MINA */
    /* 127  FORMAT(' ****ERROR IN PACK_GP.  VALUE WILL NOT PACK IN 30 BITS.', */
    /*    1       '  MAXA ='I13,'  MINA ='I13,'.  ERROR AT 127.') */
    *ier = 706;
    goto L900;

L130:

    /* ***D     WRITE(KFILDO,131)KOUNTA,KTOTAL,MINA,MAXA,IBITA,MISLLA */
    /* ***D131  FORMAT(' AT 130, KOUNTA ='I8,'  KTOTAL ='I8,'  MINA ='I8, */
    /* ***D    1       '  MAXA ='I8,'  IBITA ='I3,'  MISLLA ='I3) */

L133:
    if (ktotal >= *nxy)
    {
        goto L200;
    }

    /*        ************************************* */

    /*        THIS SECTION COMPUTES STATISTICS FOR GROUP B.  GROUP B IS A */
    /*        GROUP OF SIZE LMINPK IMMEDIATELY FOLLOWING GROUP A. */

    /*        ************************************* */

L140:
    minb = mallow;
    maxb = -mallow;
    minbk = mallow;
    maxbk = -mallow;
    ibitbs = 0;
    mstart = ktotal + 1;

    /*        DETERMINE WHETHER THERE IS A LONG STRING OF THE SAME VALUE. */
    /*        THIS WORKS WHEN THERE ARE NO MISSING VALUES. */

    nendb = 1;

    if (mstart < *nxy)
    {

        if (*is523 == 0)
        {
            /*              THIS LOOP IS FOR NO MISSING VALUES. */

            i__1 = *nxy;
            for (k = mstart + 1; k <= i__1; ++k)
            {

                if (ic[k] != ic[mstart])
                {
                    nendb = k - 1;
                    goto L150;
                }

                /* L145: */
            }

            nendb = *nxy;
            /*              FALL THROUGH THE LOOP MEANS ALL REMAINING VALUES */
            /*              ARE THE SAME. */
        }
    }

L150:
    /* Computing MAX */
    /* Computing MIN */
    i__3 = ktotal + lminpk;
    /*i__1 = nendb, i__2 = min(i__3,*nxy);*/
    i__1 = nendb, i__2 = (i__3 < *nxy) ? i__3 : *nxy;
    /*nendb = max(i__1,i__2);*/
    nendb = (i__1 > i__2) ? i__1 : i__2;
    /* **** 150  NENDB=MIN(KTOTAL+LMINPK,NXY) */

    if (*nxy - nendb <= lminpk / 2)
    {
        nendb = *nxy;
    }
    /*        ABOVE STATEMENT GUARANTEES THE LAST GROUP IS GT LMINPK/2 BY */
    /*        MAKING THE ACTUAL GROUP LARGER.  IF A PROVISION LIKE THIS IS */
    /*        NOT INCLUDED, THERE WILL MANY TIMES BE A VERY SMALL GROUP */
    /*        AT THE END.  USE SEPARATE LOOPS FOR MISSING AND NO MISSING */

    /*        USE SEPARATE LOOPS FOR MISSING AND NO MISSING VALUES */
    /*        FOR EFFICIENCY. */

    if (*is523 == 0)
    {

        i__1 = nendb;
        for (k = mstart; k <= i__1; ++k)
        {
            if (ic[k] <= minb)
            {
                minb = ic[k];
                /*              NOTE LE, NOT LT.  LT COULD BE USED BUT THEN A */
                /*              RECOMPUTE OVER THE WHOLE GROUP WOULD BE NEEDED */
                /*              MORE OFTEN.  SAME REASONING FOR GE AND OTHER */
                /*              LOOPS BELOW. */
                minbk = k;
            }
            if (ic[k] >= maxb)
            {
                maxb = ic[k];
                maxbk = k;
            }
            /* L155: */
        }
    }
    else if (*is523 == 1)
    {

        i__1 = nendb;
        for (k = mstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp)
            {
                goto L157;
            }
            if (ic[k] <= minb)
            {
                minb = ic[k];
                minbk = k;
            }
            if (ic[k] >= maxb)
            {
                maxb = ic[k];
                maxbk = k;
            }
        L157:;
        }
    }
    else
    {

        i__1 = nendb;
        for (k = mstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp || ic[k] == *misss)
            {
                goto L160;
            }
            if (ic[k] <= minb)
            {
                minb = ic[k];
                minbk = k;
            }
            if (ic[k] >= maxb)
            {
                maxb = ic[k];
                maxbk = k;
            }
        L160:;
        }
    }

    kountb = nendb - ktotal;
    misllb = 0;
    if (minb != mallow)
    {
        goto L165;
    }
    /*        ALL MISSING VALUES MUST BE ACCOMMODATED. */
    minb = 0;
    maxb = 0;
    misllb = 1;
    ibitb = 0;

    if (*is523 != 2)
    {
        goto L170;
    }
    /*        WHEN ALL VALUES ARE MISSING AND THERE ARE NO SECONDARY */
    /*        MISSING VALUES, IBITB = 0.  OTHERWISE, IBITB MUST BE */
    /*        CALCULATED. */

L165:
    for (ibitb = ibitbs; ibitb <= 30; ++ibitb)
    {
        if (maxb - minb < ibxx2[ibitb] - lmiss)
        {
            goto L170;
        }
        /* L166: */
    }

    /*     WRITE(KFILDO,167)MAXB,MINB */
    /* 167  FORMAT(' ****ERROR IN PACK_GP.  VALUE WILL NOT PACK IN 30 BITS.', */
    /*    1       '  MAXB ='I13,'  MINB ='I13,'.  ERROR AT 167.') */
    *ier = 706;
    goto L900;

    /*        COMPARE THE BITS NEEDED TO PACK GROUP B WITH THOSE NEEDED */
    /*        TO PACK GROUP A.  IF IBITB GE IBITA, TRY TO ADD TO GROUP A. */
    /*        IF NOT, TRY TO ADD A'S POINTS TO B, UNLESS ADDITION TO A */
    /*        HAS BEEN DONE.  THIS LATTER IS CONTROLLED WITH ADDA. */

L170:

    /* ***D     WRITE(KFILDO,171)KOUNTA,KTOTAL,MINA,MAXA,IBITA,MISLLA, */
    /* ***D    1                               MINB,MAXB,IBITB,MISLLB */
    /* ***D171  FORMAT(' AT 171, KOUNTA ='I8,'  KTOTAL ='I8,'  MINA ='I8, */
    /* ***D    1       '  MAXA ='I8,'  IBITA ='I3,'  MISLLA ='I3, */
    /* ***D    2       '  MINB ='I8,'  MAXB ='I8,'  IBITB ='I3,'  MISLLB ='I3) */

    if (ibitb >= ibita)
    {
        goto L180;
    }
    if (adda)
    {
        goto L200;
    }

    /*        ************************************* */

    /*        GROUP B REQUIRES LESS BITS THAN GROUP A.  PUT AS MANY OF A'S */
    /*        POINTS INTO B AS POSSIBLE WITHOUT EXCEEDING THE NUMBER OF */
    /*        BITS NECESSARY TO PACK GROUP B. */

    /*        ************************************* */

    kounts = kounta;
    /*        KOUNTA REFERS TO THE PRESENT GROUP A. */
    mintst = minb;
    maxtst = maxb;
    mintstk = minbk;
    maxtstk = maxbk;

    /*        USE SEPARATE LOOPS FOR MISSING AND NO MISSING VALUES */
    /*        FOR EFFICIENCY. */

    if (*is523 == 0)
    {

        i__1 = kstart;
        for (k = ktotal; k >= i__1; --k)
        {
            /*           START WITH THE END OF THE GROUP AND WORK BACKWARDS. */
            if (ic[k] < minb)
            {
                mintst = ic[k];
                mintstk = k;
            }
            else if (ic[k] > maxb)
            {
                maxtst = ic[k];
                maxtstk = k;
            }
            if (maxtst - mintst >= ibxx2[ibitb])
            {
                goto L174;
            }
            /*           NOTE THAT FOR THIS LOOP, LMISS = 0. */
            minb = mintst;
            maxb = maxtst;
            minbk = mintstk;
            maxbk = maxtstk;
            --kounta;
            /*           THERE IS ONE LESS POINT NOW IN A. */
            /* L1715: */
        }
    }
    else if (*is523 == 1)
    {

        i__1 = kstart;
        for (k = ktotal; k >= i__1; --k)
        {
            /*           START WITH THE END OF THE GROUP AND WORK BACKWARDS. */
            if (ic[k] == *missp)
            {
                goto L1718;
            }
            if (ic[k] < minb)
            {
                mintst = ic[k];
                mintstk = k;
            }
            else if (ic[k] > maxb)
            {
                maxtst = ic[k];
                maxtstk = k;
            }
            if (maxtst - mintst >= ibxx2[ibitb] - lmiss)
            {
                goto L174;
            }
            /*           FOR THIS LOOP, LMISS = 1. */
            minb = mintst;
            maxb = maxtst;
            minbk = mintstk;
            maxbk = maxtstk;
            misllb = 0;
            /*           WHEN THE POINT IS NON MISSING, MISLLB SET = 0. */
        L1718:
            --kounta;
            /*           THERE IS ONE LESS POINT NOW IN A. */
            /* L1719: */
        }
    }
    else
    {

        i__1 = kstart;
        for (k = ktotal; k >= i__1; --k)
        {
            /*           START WITH THE END OF THE GROUP AND WORK BACKWARDS. */
            if (ic[k] == *missp || ic[k] == *misss)
            {
                goto L1729;
            }
            if (ic[k] < minb)
            {
                mintst = ic[k];
                mintstk = k;
            }
            else if (ic[k] > maxb)
            {
                maxtst = ic[k];
                maxtstk = k;
            }
            if (maxtst - mintst >= ibxx2[ibitb] - lmiss)
            {
                goto L174;
            }
            /*           FOR THIS LOOP, LMISS = 2. */
            minb = mintst;
            maxb = maxtst;
            minbk = mintstk;
            maxbk = maxtstk;
            misllb = 0;
            /*           WHEN THE POINT IS NON MISSING, MISLLB SET = 0. */
        L1729:
            --kounta;
            /*           THERE IS ONE LESS POINT NOW IN A. */
            /* L173: */
        }
    }

    /*        AT THIS POINT, KOUNTA CONTAINS THE NUMBER OF POINTS TO CLOSE */
    /*        OUT GROUP A WITH.  GROUP B NOW STARTS WITH KSTART+KOUNTA AND */
    /*        ENDS WITH NENDB.  MINB AND MAXB HAVE BEEN ADJUSTED AS */
    /*        NECESSARY TO REFLECT GROUP B (EVEN THOUGH THE NUMBER OF BITS */
    /*        NEEDED TO PACK GROUP B HAVE NOT INCREASED, THE END POINTS */
    /*        OF THE RANGE MAY HAVE). */

L174:
    if (kounta == kounts)
    {
        goto L200;
    }
    /*        ON TRANSFER, GROUP A WAS NOT CHANGED.  CLOSE IT OUT. */

    /*        ONE OR MORE POINTS WERE TAKEN OUT OF A.  RANGE AND IBITA */
    /*        MAY HAVE TO BE RECOMPUTED; IBITA COULD BE LESS THAN */
    /*        ORIGINALLY COMPUTED.  IN FACT, GROUP A CAN NOW CONTAIN */
    /*        ONLY ONE POINT AND BE PACKED WITH ZERO BITS */
    /*        (UNLESS MISSS NE 0). */

    nouta = kounts - kounta;
    ktotal -= nouta;
    kountb += nouta;
    if (nenda - nouta > minak && nenda - nouta > maxak)
    {
        goto L200;
    }
    /*        WHEN THE ABOVE TEST IS MET, THE MIN AND MAX OF THE */
    /*        CURRENT GROUP A WERE WITHIN THE OLD GROUP A, SO THE */
    /*        RANGE AND IBITA DO NOT NEED TO BE RECOMPUTED. */
    /*        NOTE THAT MINAK AND MAXAK ARE NO LONGER NEEDED. */
    ibita = 0;
    mina = mallow;
    maxa = -mallow;

    /*        USE SEPARATE LOOPS FOR MISSING AND NO MISSING VALUES */
    /*        FOR EFFICIENCY. */

    if (*is523 == 0)
    {

        i__1 = nenda - nouta;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] < mina)
            {
                mina = ic[k];
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
            }
            /* L1742: */
        }
    }
    else if (*is523 == 1)
    {

        i__1 = nenda - nouta;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp)
            {
                goto L1744;
            }
            if (ic[k] < mina)
            {
                mina = ic[k];
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
            }
        L1744:;
        }
    }
    else
    {

        i__1 = nenda - nouta;
        for (k = kstart; k <= i__1; ++k)
        {
            if (ic[k] == *missp || ic[k] == *misss)
            {
                goto L175;
            }
            if (ic[k] < mina)
            {
                mina = ic[k];
            }
            if (ic[k] > maxa)
            {
                maxa = ic[k];
            }
        L175:;
        }
    }

    mislla = 0;
    if (mina != mallow)
    {
        goto L1750;
    }
    /*        ALL MISSING VALUES MUST BE ACCOMMODATED. */
    mina = 0;
    maxa = 0;
    mislla = 1;
    if (*is523 != 2)
    {
        goto L177;
    }
    /*        WHEN ALL VALUES ARE MISSING AND THERE ARE NO SECONDARY */
    /*        MISSING VALUES IBITA = 0 AS ORIGINALLY SET.  OTHERWISE, */
    /*        IBITA MUST BE CALCULATED. */

L1750:
    itest = maxa - mina + lmiss;

    for (ibita = 0; ibita <= 30; ++ibita)
    {
        if (itest < ibxx2[ibita])
        {
            goto L177;
        }
        /* ***        THIS TEST IS THE SAME AS: */
        /* ***         IF(MAXA-MINA.LT.IBXX2(IBITA)-LMISS)GO TO 177 */
        /* L176: */
    }

    /*     WRITE(KFILDO,1760)MAXA,MINA */
    /* 1760 FORMAT(' ****ERROR IN PACK_GP.  VALUE WILL NOT PACK IN 30 BITS.', */
    /*    1       '  MAXA ='I13,'  MINA ='I13,'.  ERROR AT 1760.') */
    *ier = 706;
    goto L900;

L177:
    goto L200;

    /*        ************************************* */

    /*        AT THIS POINT, GROUP B REQUIRES AS MANY BITS TO PACK AS GROUPA. */
    /*        THEREFORE, TRY TO ADD INC POINTS TO GROUP A WITHOUT INCREASING */
    /*        IBITA.  THIS AUGMENTED GROUP IS CALLED GROUP C. */

    /*        ************************************* */

L180:
    if (mislla == 1)
    {
        minc = mallow;
        minck = mallow;
        maxc = -mallow;
        maxck = -mallow;
    }
    else
    {
        minc = mina;
        maxc = maxa;
        minck = minak;
        maxck = minak;
    }

    nount = 0;
    if (*nxy - (ktotal + kinc) <= lminpk / 2)
    {
        kinc = *nxy - ktotal;
    }
    /*        ABOVE STATEMENT CONSTRAINS THE LAST GROUP TO BE NOT LESS THAN */
    /*        LMINPK/2 IN SIZE.  IF A PROVISION LIKE THIS IS NOT INCLUDED, */
    /*        THERE WILL MANY TIMES BE A VERY SMALL GROUP AT THE END. */

    /*        USE SEPARATE LOOPS FOR MISSING AND NO MISSING VALUES */
    /*        FOR EFFICIENCY.  SINCE KINC IS USUALLY 1, USING SEPARATE */
    /*        LOOPS HERE DOESN'T BUY MUCH.  A MISSING VALUE WILL ALWAYS */
    /*        TRANSFER BACK TO GROUP A. */

    if (*is523 == 0)
    {

        /* Computing MIN */
        i__2 = ktotal + kinc;
        /*i__1 = min(i__2,*nxy);*/
        i__1 = (i__2 < *nxy) ? i__2 : *nxy;
        for (k = ktotal + 1; k <= i__1; ++k)
        {
            if (ic[k] < minc)
            {
                minc = ic[k];
                minck = k;
            }
            if (ic[k] > maxc)
            {
                maxc = ic[k];
                maxck = k;
            }
            ++nount;
            /* L185: */
        }
    }
    else if (*is523 == 1)
    {

        /* Computing MIN */
        i__2 = ktotal + kinc;
        /*i__1 = min(i__2,*nxy);*/
        i__1 = (i__2 < *nxy) ? i__2 : *nxy;
        for (k = ktotal + 1; k <= i__1; ++k)
        {
            if (ic[k] == *missp)
            {
                goto L186;
            }
            if (ic[k] < minc)
            {
                minc = ic[k];
                minck = k;
            }
            if (ic[k] > maxc)
            {
                maxc = ic[k];
                maxck = k;
            }
        L186:
            ++nount;
            /* L187: */
        }
    }
    else
    {

        /* Computing MIN */
        i__2 = ktotal + kinc;
        /*i__1 = min(i__2,*nxy);*/
        i__1 = (i__2 < *nxy) ? i__2 : *nxy;
        for (k = ktotal + 1; k <= i__1; ++k)
        {
            if (ic[k] == *missp || ic[k] == *misss)
            {
                goto L189;
            }
            if (ic[k] < minc)
            {
                minc = ic[k];
                minck = k;
            }
            if (ic[k] > maxc)
            {
                maxc = ic[k];
                maxck = k;
            }
        L189:
            ++nount;
            /* L190: */
        }
    }

    /* ***D     WRITE(KFILDO,191)KOUNTA,KTOTAL,MINA,MAXA,IBITA,MISLLA, */
    /* ***D    1   MINC,MAXC,NOUNT,IC(KTOTAL),IC(KTOTAL+1) */
    /* ***D191  FORMAT(' AT 191, KOUNTA ='I8,'  KTOTAL ='I8,'  MINA ='I8, */
    /* ***D    1       '  MAXA ='I8,'  IBITA ='I3,'  MISLLA ='I3, */
    /* ***D    2       '  MINC ='I8,'  MAXC ='I8, */
    /* ***D    3       '  NOUNT ='I5,'  IC(KTOTAL) ='I9,'  IC(KTOTAL+1) =',I9) */

    /*        IF THE NUMBER OF BITS NEEDED FOR GROUP C IS GT IBITA, */
    /*        THEN THIS GROUP A IS A GROUP TO PACK. */

    if (minc == mallow)
    {
        minc = mina;
        maxc = maxa;
        minck = minak;
        maxck = maxak;
        misllc = 1;
        goto L195;
        /*           WHEN THE NEW VALUE(S) ARE MISSING, THEY CAN ALWAYS */
        /*           BE ADDED. */
    }
    else
    {
        misllc = 0;
    }

    if (maxc - minc >= ibxx2[ibita] - lmiss)
    {
        goto L200;
    }

    /*        THE BITS NECESSARY FOR GROUP C HAS NOT INCREASED FROM THE */
    /*        BITS NECESSARY FOR GROUP A.  ADD THIS POINT(S) TO GROUP A. */
    /*        COMPUTE THE NEXT GROUP B, ETC., UNLESS ALL POINTS HAVE BEEN */
    /*        USED. */

L195:
    ktotal += nount;
    kounta += nount;
    mina = minc;
    maxa = maxc;
    minak = minck;
    maxak = maxck;
    mislla = misllc;
    adda = TRUE_;
    if (ktotal >= *nxy)
    {
        goto L200;
    }

    if (minbk > ktotal && maxbk > ktotal)
    {
        mstart = nendb + 1;
        /*           THE MAX AND MIN OF GROUP B WERE NOT FROM THE POINTS */
        /*           REMOVED, SO THE WHOLE GROUP DOES NOT HAVE TO BE LOOKED */
        /*           AT TO DETERMINE THE NEW MAX AND MIN.  RATHER START */
        /*           JUST BEYOND THE OLD NENDB. */
        ibitbs = ibitb;
        nendb = 1;
        goto L150;
    }
    else
    {
        goto L140;
    }

    /*        ************************************* */

    /*        GROUP A IS TO BE PACKED.  STORE VALUES IN JMIN( ), JMAX( ), */
    /*        LBIT( ), AND NOV( ). */

    /*        ************************************* */

L200:
    ++(*lx);
    if (*lx <= *ndg)
    {
        goto L205;
    }
    lminpk += lminpk / 2;
    /*     WRITE(KFILDO,201)NDG,LMINPK,LX */
    /* 201  FORMAT(' ****NDG ='I5,' NOT LARGE ENOUGH.', */
    /*    1       '  LMINPK IS INCREASED TO 'I3,' FOR THIS FIELD.'/ */
    /*    2       '  LX = 'I10) */
    iersav = 716;
    goto L105;

L205:
    jmin[*lx] = mina;
    jmax[*lx] = maxa;
    lbit[*lx] = ibita;
    nov[*lx] = kounta;
    kstart = ktotal + 1;

    if (mislla == 0)
    {
        misslx[*lx - 1] = mallow;
    }
    else
    {
        misslx[*lx - 1] = ic[ktotal];
        /*           IC(KTOTAL) WAS THE LAST VALUE PROCESSED.  IF MISLLA NE 0, */
        /*           THIS MUST BE THE MISSING VALUE FOR THIS GROUP. */
    }

    /* ***D     WRITE(KFILDO,206)MISLLA,IC(KTOTAL),KTOTAL,LX,JMIN(LX),JMAX(LX), */
    /* ***D    1                 LBIT(LX),NOV(LX),MISSLX(LX) */
    /* ***D206  FORMAT(' AT 206,  MISLLA ='I2,'  IC(KTOTAL) ='I5,'  KTOTAL ='I8, */
    /* ***D    1       '  LX ='I6,'  JMIN(LX) ='I8,'  JMAX(LX) ='I8, */
    /* ***D    2       '  LBIT(LX) ='I5,'  NOV(LX) ='I8,'  MISSLX(LX) =',I7) */

    if (ktotal >= *nxy)
    {
        goto L209;
    }

    /*        THE NEW GROUP A WILL BE THE PREVIOUS GROUP B.  SET LIMITS, ETC. */

    ibita = ibitb;
    mina = minb;
    maxa = maxb;
    minak = minbk;
    maxak = maxbk;
    mislla = misllb;
    nenda = nendb;
    kounta = kountb;
    ktotal += kounta;
    adda = FALSE_;
    goto L133;

    /*        ************************************* */

    /*        CALCULATE IBIT, THE NUMBER OF BITS NEEDED TO HOLD THE GROUP */
    /*        MINIMUM VALUES. */

    /*        ************************************* */

L209:
    *ibit = 0;

    i__1 = *lx;
    for (l = 1; l <= i__1; ++l)
    {
    L210:
        if (jmin[l] < ibxx2[*ibit])
        {
            goto L220;
        }
        ++(*ibit);
        goto L210;
    L220:;
    }

    /*        INSERT THE VALUE IN JMIN( ) TO BE USED FOR ALL MISSING */
    /*        VALUES WHEN LBIT( ) = 0.  WHEN SECONDARY MISSING */
    /*        VALUES CAN BE PRESENT, LBIT(L) WILL NOT = 0. */

    if (*is523 == 1)
    {

        i__1 = *lx;
        for (l = 1; l <= i__1; ++l)
        {

            if (lbit[l] == 0)
            {

                if (misslx[l - 1] == *missp)
                {
                    jmin[l] = ibxx2[*ibit] - 1;
                }
            }

            /* L226: */
        }
    }

    /*        ************************************* */

    /*        CALCULATE JBIT, THE NUMBER OF BITS NEEDED TO HOLD THE BITS */
    /*        NEEDED TO PACK THE VALUES IN THE GROUPS.  BUT FIND AND */
    /*        REMOVE THE REFERENCE VALUE FIRST. */

    /*        ************************************* */

    /*     WRITE(KFILDO,228)CFEED,LX */
    /* 228  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP WIDTHS LBIT( ) FOR ',I8,' GROUPS' */
    /*    2          /' *****************************************') */
    /*     WRITE(KFILDO,229) (LBIT(J),J=1,MIN(LX,100)) */
    /* 229  FORMAT(/' '20I6) */

    *lbitref = lbit[1];

    i__1 = *lx;
    for (k = 1; k <= i__1; ++k)
    {
        if (lbit[k] < *lbitref)
        {
            *lbitref = lbit[k];
        }
        /* L230: */
    }

    if (*lbitref != 0)
    {

        i__1 = *lx;
        for (k = 1; k <= i__1; ++k)
        {
            lbit[k] -= *lbitref;
            /* L240: */
        }
    }

    /*     WRITE(KFILDO,241)CFEED,LBITREF */
    /* 241  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP WIDTHS LBIT( ) AFTER REMOVING REFERENCE ', */
    /*    2             I8, */
    /*    3          /' *****************************************') */
    /*     WRITE(KFILDO,242) (LBIT(J),J=1,MIN(LX,100)) */
    /* 242  FORMAT(/' '20I6) */

    *jbit = 0;

    i__1 = *lx;
    for (k = 1; k <= i__1; ++k)
    {
    L310:
        if (lbit[k] < ibxx2[*jbit])
        {
            goto L320;
        }
        ++(*jbit);
        goto L310;
    L320:;
    }

    /*        ************************************* */

    /*        CALCULATE KBIT, THE NUMBER OF BITS NEEDED TO HOLD THE NUMBER */
    /*        OF VALUES IN THE GROUPS.  BUT FIND AND REMOVE THE */
    /*        REFERENCE FIRST. */

    /*        ************************************* */

    /*     WRITE(KFILDO,321)CFEED,LX */
    /* 321  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP SIZES NOV( ) FOR ',I8,' GROUPS' */
    /*    2          /' *****************************************') */
    /*     WRITE(KFILDO,322) (NOV(J),J=1,MIN(LX,100)) */
    /* 322  FORMAT(/' '20I6) */

    *novref = nov[1];

    i__1 = *lx;
    for (k = 1; k <= i__1; ++k)
    {
        if (nov[k] < *novref)
        {
            *novref = nov[k];
        }
        /* L400: */
    }

    if (*novref > 0)
    {

        i__1 = *lx;
        for (k = 1; k <= i__1; ++k)
        {
            nov[k] -= *novref;
            /* L405: */
        }
    }

    /*     WRITE(KFILDO,406)CFEED,NOVREF */
    /* 406  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP SIZES NOV( ) AFTER REMOVING REFERENCE ',I8, */
    /*    2          /' *****************************************') */
    /*     WRITE(KFILDO,407) (NOV(J),J=1,MIN(LX,100)) */
    /* 407  FORMAT(/' '20I6) */
    /*     WRITE(KFILDO,408)CFEED */
    /* 408  FORMAT(A1,/' *****************************************' */
    /*    1          /' THE GROUP REFERENCES JMIN( )' */
    /*    2          /' *****************************************') */
    /*     WRITE(KFILDO,409) (JMIN(J),J=1,MIN(LX,100)) */
    /* 409  FORMAT(/' '20I6) */

    *kbit = 0;

    i__1 = *lx;
    for (k = 1; k <= i__1; ++k)
    {
    L410:
        if (nov[k] < ibxx2[*kbit])
        {
            goto L420;
        }
        ++(*kbit);
        goto L410;
    L420:;
    }

    /*        DETERMINE WHETHER THE GROUP SIZES SHOULD BE REDUCED */
    /*        FOR SPACE EFFICIENCY. */

    if (ired == 0)
    {
        reduce(kfildo, &jmin[1], &jmax[1], &lbit[1], &nov[1], lx, ndg, ibit,
               jbit, kbit, novref, ibxx2, ier);

        if (*ier == 714 || *ier == 715)
        {
            /*              REDUCE HAS ABORTED.  REEXECUTE PACK_GP WITHOUT REDUCE. */
            /*              PROVIDE FOR A NON FATAL RETURN FROM REDUCE. */
            iersav = *ier;
            ired = 1;
            *ier = 0;
            goto L102;
        }
    }

    if (misslx != 0)
    {
        free(misslx);
        misslx = 0;
    }
    /*     CALL TIMPR(KFILDO,KFILDO,'END   PACK_GP        ') */
    if (iersav != 0)
    {
        *ier = iersav;
        return 0;
    }

    /* 900  IF(IER.NE.0)RETURN1 */

L900:
    if (misslx != 0)
        free(misslx);
    return 0;
} /* pack_gp__ */
