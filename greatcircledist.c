/*****************************************************************************
 * GREAT CIRCLE DISTANCE
 * Calculates the surface distance between two lat/lon points on Earth. Use as
 * a command line tool with exactly four arguments:
 *
 *     greatcircledist lat1 lon1 lat2 lon2
 *
 * where all arguments are floating point numbers in degrees. Outputs distance
 * in metres between (lat1,lon1) and (lat2,lon2) with precision in centimetres.
 *
 * Author: E. Dronkert https://github.com/ednl
 * Licence: MIT (free to use as you like, with attribution)
 *****************************************************************************/

#include <stdio.h>    // printf, fprintf
#include <stdlib.h>   // strtod
#include <string.h>   // strlen
#include <math.h>     // sin, cos, asin, sqrt
#include <errno.h>    // errno, ERANGE
#include <stdbool.h>  // bool

#define E_NUMARG  1
#define E_INVALID 2
#define E_RANGE   3
#define E_LAT90   4
#define E_LON180  5

#define EPS  1e-8  // epsilon for double comparison
#define D2R  1.74532925199432953e-2  // pi/180 (degrees to radians)
#define D2RH 8.72664625997164656e-3  // pi/360 (degrees to radians, half)

#define R_POL   6.357e+6         // earth polar radius in metres
#define R_EQT   6.378e+6         // earth equatorial radius in metres
#define R_POL2  (R_POL * R_POL)  // square polar radius
#define R_EQT2  (R_EQT * R_EQT)  // square equatorial radius

// Two lat/lon points in degrees
// = four doubles addressable as either deg[0..3] or lat1,lon1,lat2,lon2
typedef union {
    double deg[4];
    struct { double lat1, lon1, lat2, lon2; };
} Segment;

static bool equal(const double a, const double b)
{
    return fabs(a - b) <= EPS;
}

// Local earth diameter in metres, from latitude in degrees.
// Ref.: https://en.wikipedia.org/wiki/Earth_radius#Location-dependent_radii
static double local_earth_diameter(const double latdeg)
{
    if (equal(latdeg, 0))
        return 2 * R_EQT;
    if (equal(latdeg, 90) || equal(latdeg, -90))
        return 2 * R_POL;
    double lat = latdeg * D2R;
    double s = sin(lat), c = cos(lat);
    double rs2 = R_POL2 * s * s;
    double rc2 = R_EQT2 * c * c;
    return 2 * sqrt((R_POL2 * rs2 + R_EQT2 * rc2) / (rs2 + rc2));
}

// Great circle distance in metres, using local earth radius and
// inverse haversine function, for 2 lat/lon points in degrees.
// Ref.: https://en.wikipedia.org/wiki/Haversine_formula#Formulation
static double great_circle_distance(const Segment a)
{
    double dlat = a.lat2 - a.lat1;
    double dlon = a.lon2 - a.lon1;
    bool eqlat = equal(dlat, 0);
    bool eqlon = equal(dlon, 0);

    if (eqlat && eqlon)
        return 0;

    if (eqlat)
        return local_earth_diameter(a.lat1) * asin(fabs(cos(a.lat1 * D2R) * sin(dlon * D2RH)));

    double D = local_earth_diameter((a.lat1 + a.lat2) / 2);
    if (eqlon)
        return D * asin(fabs(sin(dlat * D2RH)));

    double slat = sin(dlat * D2RH);  // sin lat diff
    double slon = sin(dlon * D2RH);  // sin lon diff
    return D * asin(sqrt(slat * slat + cos(a.lat1 * D2R) * cos(a.lat2 * D2R) * slon * slon));
}

int main(int argc, char *argv[])
{
    if (argc != 5) {
        fprintf(stderr, "Provide 4 arguments: lat1 lon1 lat2 lon2.\n");
        exit(E_NUMARG);
    }

    Segment a = {0};
    for (int i = 0, j = 1; i < 4; ++i, ++j) {
        char *end;
        a.deg[i] = strtod(argv[j], &end);  // parse to double
        if ((size_t)(end - argv[j]) != strlen(argv[j])) {
            fprintf(stderr, "Not a number: %s.\n", argv[j]);
            exit(E_INVALID);
        }
        if (errno == ERANGE) {
            fprintf(stderr, "Out of range: %s.\n", argv[j]);
            exit(E_RANGE);
        }
        if ((i == 0 || i == 2) && (a.deg[i] < -90 || a.deg[i] > 90)) {
            fprintf(stderr, "Latitude must be between -90 and +90: %s.\n", argv[j]);
            exit(E_LAT90);
        }
        if ((i == 1 || i == 3) && (a.deg[i] < -180 || a.deg[i] > 180)) {
            fprintf(stderr, "Longitude must be between -180 and +180: %s.\n", argv[j]);
            exit(E_LON180);
        }
    }

    // Distance in m, precision to cm
    printf("%.2f\n", great_circle_distance(a));
    return 0;
}
