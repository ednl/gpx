/*****************************************************************************
 * GREAT CIRCLE DISTANCE
 * Calculates the surface distance between two lat/lon points on Earth. Use as
 * a command line tool with exactly four arguments:
 *
 *     greatcircledist lat1 lon1 lat2 lon2
 *
 * where all arguments are decimal degrees. Outputs distance in metres between
 * (lat1,lon1) and (lat2,lon2) with precision in centimetres.
 *
 * Author: E. Dronkert https://github.com/ednl
 * Licence: MIT (free to use as you like, with attribution)
 *****************************************************************************/

#include <stdio.h>    // printf, fprintf
#include <stdlib.h>   // strtod
#include <string.h>   // strlen
#include <math.h>     // sin, cos, asin, sqrt
#include <errno.h>    // errno, ERANGE

// Error exit codes
#define ERR_NUMARG  1  // number of arguments must be 4
#define ERR_INVALID 2  // argument must be a floating point number in the C locale
#define ERR_RANGE   3  // argument must be a valid double
#define ERR_LAT90   4  // latitude must be between -90 and +90
#define ERR_LON180  5  // longitude must be between -180 and +180

// Distance calculation constants
#define DEG2RAD   1.74532925199432953e-2  // pi/180 (degrees to radians)
#define RADIUS_P  6.357e+6                // earth polar radius in metres
#define RADIUS_E  6.378e+6                // earth equatorial radius in metres
#define RADIUS_P2 (RADIUS_P * RADIUS_P)   // square polar radius
#define RADIUS_E2 (RADIUS_E * RADIUS_E)   // square equatorial radius

// Line segment on the Earth surface defined by 2 lat/lon points
// = four doubles addressable as either coor[0..3] or lat1,lon1,lat2,lon2
typedef union {
    double coor[4];
    struct { double lat1, lon1, lat2, lon2; };
} LineSegment;

int main(int argc, char *argv[])
{
    // Number of arguments = 1 (program name) + 4 (command line arguments)
    if (argc != 5) {
        fprintf(stderr, "Provide 4 arguments: lat1 lon1 lat2 lon2.\n");
        exit(ERR_NUMARG);
    }

    LineSegment a = {0};
    for (int i = 0, j = 1; i < 4; ++i, ++j) {  // value index i, argument index j
        // Parse string argument to double
        char *end;
        a.coor[i] = strtod(argv[j], &end);
        // Whole string used?
        if ((size_t)(end - argv[j]) != strlen(argv[j])) {
            fprintf(stderr, "Not a number: %s.\n", argv[j]);
            exit(ERR_INVALID);
        }
        // Representable as double?
        if (errno == ERANGE) {
            fprintf(stderr, "Out of range: %s.\n", argv[j]);
            exit(ERR_RANGE);
        }
        // Arg 1 and 3 (values 0 and 2) are latitudes in degrees [-90,+90]
        if ((i == 0 || i == 2) && (a.coor[i] < -90 || a.coor[i] > 90)) {
            fprintf(stderr, "Latitude must be between -90 and +90: %s.\n", argv[j]);
            exit(ERR_LAT90);
        }
        // Arg 2 and 4 (values 1 and 3) are longitudes in degrees [-180,+180]
        if ((i == 1 || i == 3) && (a.coor[i] < -180 || a.coor[i] > 180)) {
            fprintf(stderr, "Longitude must be between -180 and +180: %s.\n", argv[j]);
            exit(ERR_LON180);
        }
        // From degrees to radians
        a.coor[i] *= DEG2RAD;
    }

    // Intermediate values for local earth radius
    // Ref.: https://en.wikipedia.org/wiki/Earth_radius#Location-dependent_radii
    double avglat = (a.lat1 + a.lat2) / 2;
    double s = sin(avglat);
    double c = cos(avglat);
    double rs = RADIUS_P2 * s * s;
    double rc = RADIUS_E2 * c * c;

    // Intermediate values for inverse haversine function
    // Ref.: https://en.wikipedia.org/wiki/Haversine_formula#Formulation
    double slat = sin((a.lat2 - a.lat1) / 2);
    double slon = sin((a.lon2 - a.lon1) / 2);

    // Great circle distance in m, precision to cm
    printf("%.2f\n", 2 * sqrt((RADIUS_P2 * rs + RADIUS_E2 * rc) / (rs + rc))
        * asin(sqrt(slat * slat + cos(a.lat1) * cos(a.lat2) * slon * slon)));

    return 0;
}
