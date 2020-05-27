/*******************************************************************************
NAME                            TRANSVERSE MERCATOR

PURPOSE:	Transforms input longitude and latitude to Easting and
		Northing for the Transverse Mercator projection.  The
		longitude and latitude must be in radians.  The Easting
		and Northing values will be returned in meters.

PROGRAMMER              DATE		REASON
----------              ----		------
D. Steinwand, EROS      Nov, 1991
T. Mittan		Mar, 1993
S. Nelson		Feb, 1995	Divided tmfor.c into two files, one
                                        for UTM (utmfor.c) and one for 
					TM (tmfor.c).  This was a
                                        necessary change to run forward
                                        projection conversions for both UTM
					and TM in the same process.

ALGORITHM REFERENCES

1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
    Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
    State Government Printing Office, Washington D.C., 1987.

2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
    U.S. Geological Survey Professional Paper 1453 , United State Government
    Printing Office, Washington D.C., 1989.
*******************************************************************************/
#include "cproj.h"

/* Variables common to all subroutines in this code file
  -----------------------------------------------------*/
static double r_major;		/* major axis 				*/
static double r_minor;		/* minor axis 				*/
static double scale_factor;	/* scale factor				*/
static double lon_center;	/* Center longitude (projection center) */
static double lat_origin;	/* center latitude			*/
static double e0,e1,e2,e3;	/* eccentricity constants		*/
static double e,es,esp;		/* eccentricity constants		*/
static double ml0;		/* small value m			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/
static double ind;		/* spherical flag			*/


/* Initialize the Transverse Mercator (TM) projection
  -------------------------------------------------*/
long tmforint(r_maj,r_min,scale_fact,center_lon,center_lat,false_east,
	false_north) 

double r_maj;			/* major axis			*/
double r_min;			/* minor axis			*/
double scale_fact;		/* scale factor			*/
double center_lon;		/* center longitude		*/
double center_lat;		/* center latitude		*/
double false_east;		/* x offset in meters		*/
double false_north;		/* y offset in meters		*/
{
double temp;			/* temporary variable		*/

/* Place parameters in static storage for common use
  -------------------------------------------------*/
r_major = r_maj;
r_minor = r_min;
scale_factor = scale_fact;
lon_center = center_lon;
lat_origin = center_lat;
false_northing = false_north;
false_easting = false_east;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = e0fn(es);
e1 = e1fn(es);
e2 = e2fn(es);
e3 = e3fn(es);
ml0 = r_major * mlfn(e0, e1, e2, e3, lat_origin);
esp = es / (1.0 - es);

if (es < .00001)
   ind = 1;
else 
   ind = 0;

/* Report parameters to the user
  -----------------------------*/
ptitle("TRANSVERSE MERCATOR (TM)"); 
radius2(r_major, r_minor);
genrpt(scale_factor,"Scale Factor at C. Meridian:    ");
cenlonmer(lon_center);
origin(lat_origin);
offsetp(false_easting,false_northing);
return(OK);
}

/* Transverse Mercator forward equations--mapping lat,long to x,y
   Note:  The algorithm for TM is exactly the same as UTM and therefore
	  if a change is implemented, also make the change to UTMFOR.c
  --------------------------------------------------------------*/
long tmfor(lon, lat, x, y)
double lon;			/* (I) Longitude 		*/
double lat;			/* (I) Latitude 		*/
double *x;			/* (O) X projection coordinate 	*/
double *y;			/* (O) Y projection coordinate 	*/
{
double delta_lon;	/* Delta longitude (Given longitude - center 	*/
double theta;		/* angle					*/
double delta_theta;	/* adjusted longitude				*/
double sin_phi, cos_phi;/* sin and cos value				*/
double al, als;		/* temporary values				*/
double b;		/* temporary values				*/
double c, t, tq;	/* temporary values				*/
double con, n, ml;	/* cone constant, small m			*/

/* Forward equations
  -----------------*/
delta_lon = adjust_lon(lon - lon_center);
sincos(lat, &sin_phi, &cos_phi);

/* This part was in the fortran code and is for the spherical form 
----------------------------------------------------------------*/
if (ind != 0)
  {
  b = cos_phi * sin(delta_lon);
  if ((fabs(fabs(b) - 1.0)) < .0000000001)
     {
     p_error("Point projects into infinity","tm-for");
     return(93);
     }
  else
     {
     *x = .5 * r_major * scale_factor * log((1.0 + b)/(1.0 - b));
     con = acos(cos_phi * cos(delta_lon)/sqrt(1.0 - b*b));
     if (lat < 0)
        con = - con;
     *y = r_major * scale_factor * (con - lat_origin); 
     return(OK);
     }
  }

al  = cos_phi * delta_lon;
als = SQUARE(al);
c   = esp * SQUARE(cos_phi);
tq  = tan(lat);
t   = SQUARE(tq);
con = 1.0 - es * SQUARE(sin_phi);
n   = r_major / sqrt(con);
ml  = r_major * mlfn(e0, e1, e2, e3, lat);

*x  = scale_factor * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 *
      (5.0 - 18.0 * t + SQUARE(t) + 72.0 * c - 58.0 * esp))) + false_easting;

*y  = scale_factor * (ml - ml0 + n * tq * (als * (0.5 + als / 24.0 *
      (5.0 - t + 9.0 * c + 4.0 * SQUARE(c) + als / 30.0 * (61.0 - 58.0 * t
      + SQUARE(t) + 600.0 * c - 330.0 * esp))))) + false_northing;

return(OK);
}
