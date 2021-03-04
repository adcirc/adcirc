/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1.1.1 2003/01/23 17:21:05 estrade Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif

/*TCM added windows dependency*/
#ifndef WINDOWS
#include <strings.h>
#endif

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "rename.h"
#include "proto.h"

