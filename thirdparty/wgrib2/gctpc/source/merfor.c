/*******************************************************************************
NAME                            MERCATOR

PURPOSE:	Transforms input longitude and latitude to Easting and
		Northing for the Mercator projection.  The
		longitude and latitude must be in radians.  The Easting
		and Northing values will be returned in meters.

PROGRAMMER              DATE
----------              ----
D. Steinwand, EROS      Nov, 1991
T. Mittan		Mar, 1993

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
static double lon_center;	/* Center longitude (projection center) */
static double lat_origin;	/* center latitude			*/
static double e,es;		/* eccentricity constants		*/
static double m1;		/* small value m			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/


/* Initialize the Mercator projection
  -------------------------------------------------*/
long merforint(r_maj,r_min,center_lon,center_lat,false_east,false_north) 

double r_maj;			/* major axis			*/
double r_min;			/* minor axis			*/
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
lon_center = center_lon;
lat_origin = center_lat;
false_northing = false_north;
false_easting = false_east;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
m1 = cos(center_lat)/(sqrt(1.0 - es * sin(center_lat) * sin(center_lat)));

/* Report parameters to the user
  -----------------------------*/
ptitle("MERCATOR"); 
radius2(r_major, r_minor);
cenlonmer(lon_center);
origin(lat_origin);
offsetp(false_easting,false_northing);
return(OK);
}


/* Mercator forward equations--mapping lat,long to x,y
  --------------------------------------------------*/
long merfor(lon, lat, x, y)
double lon;			/* (I) Longitude 		*/
double lat;			/* (I) Latitude 		*/
double *x;			/* (O) X projection coordinate 	*/
double *y;			/* (O) Y projection coordinate 	*/
{
double ts;		/* small t value				*/
double sinphi;		/* sin value					*/

/* Forward equations
  -----------------*/
if (fabs(fabs(lat) - HALF_PI)  <= EPSLN)
   {
   p_error("Transformation cannot be computed at the poles","mer-forward");
   return(53);
   }
else
   {
   sinphi = sin(lat);
   ts = tsfnz(e,lat,sinphi);
   *x = false_easting + r_major * m1 * adjust_lon(lon - lon_center);
   *y = false_northing - r_major * m1 * log(ts);
   }
return(OK);
}
