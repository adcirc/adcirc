/*******************************************************************************
NAME                           AZIMUTHAL EQUIDISTANT 

PURPOSE:	Transforms input Easting and Northing to longitude and
		latitude for the Azimuthal Equidistant projection.  The
		Easting and Northing must be in meters.  The longitude
		and latitude values will be returned in radians.

PROGRAMMER              DATE
----------              ----
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
static double lon_center;	/* Center longitude (projection center) */
static double lat_origin;	/* center latitude			*/
static double false_northing;	/* y offset in meters			*/
static double false_easting;	/* x offset in meters			*/
static double sin_p12;		/* sin of center latitude		*/
static double cos_p12;		/* cos of center latitude		*/

/* Initialize the Azimuthal projection
  ----------------------------------*/
long aziminvint(r_maj,center_lon,center_lat,false_east,false_north) 

double r_maj;			/* major axis			*/
double center_lon;		/* center longitude		*/
double center_lat;		/* center latitude		*/
double false_east;		/* x offset in meters		*/
double false_north;		/* y offset in meters		*/
{

/* Place parameters in static storage for common use
  -------------------------------------------------*/
r_major = r_maj;
lon_center = center_lon;
lat_origin = center_lat;
false_northing = false_north;
false_easting = false_east;

sincos(center_lat,&sin_p12,&cos_p12);

/* Report parameters to the user
  -----------------------------*/
ptitle("AZIMUTHAL EQUIDISTANT"); 
radius(r_major);
cenlonmer(lon_center);
origin(lat_origin);
offsetp(false_easting,false_northing);
return(OK);
}


/* Azimuthal inverse equations--mapping x,y to lat/long
  ---------------------------------------------------*/
long aziminv(x, y, lon, lat)
double x;			/* (O) X projection coordinate 	*/
double y;			/* (O) Y projection coordinate 	*/
double *lon;			/* (I) Longitude 		*/
double *lat;			/* (I) Latitude 		*/
{
double rh;		/* height above ellipsoid			*/
double z;		/* angle					*/
double sinz,cosz;	/* sin of z and cos of z			*/
double temp;
double con;

/* Inverse equations
  -----------------*/
x -= false_easting;
y -= false_northing;
rh = sqrt(x * x + y * y);
if (rh > (2.0 * HALF_PI * r_major))
   {
   p_error("Input data error","azim-inv");
   return(125);
   }
z = rh / r_major;
sincos(z,&sinz,&cosz);
*lon = lon_center;
if (fabs(rh) <= EPSLN)
   {
   *lat = lat_origin;
   return(OK);
   }
*lat = asinz(cosz * sin_p12 + (y * sinz * cos_p12) / rh);
con = fabs(lat_origin) - HALF_PI;
if (fabs(con) <= EPSLN)
   {
   if (lat_origin >= 0.0)
      {
      *lon = adjust_lon(lon_center + atan2(x , -y));
      return(OK);
      }
   else
      {
      *lon = adjust_lon(lon_center - atan2(-x , y));
      return(OK);
      }
   }
con = cosz - sin_p12 * sin(*lat);
if ((fabs(con) < EPSLN) && (fabs(x) < EPSLN))
   return(OK);
temp = atan2((x * sinz * cos_p12), (con * rh));
*lon = adjust_lon(lon_center + atan2((x * sinz * cos_p12), (con * rh)));

return(OK);
} 
