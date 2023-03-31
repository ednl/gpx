/*****************************************************************************
 * GREAT-CIRCLE DISTANCE
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
#include <stdbool.h>  // bool

// Error exit codes
#define ERR_NUMARG  1  // number of arguments must be 4
#define ERR_INVALID 2  // argument must be a floating point number in the C locale
#define ERR_RANGE   3  // argument must be a valid double
#define ERR_LAT90   4  // latitude must be between -90 and +90
#define ERR_LON180  5  // longitude must be between -180 and +180

// Mathematical constants
#define DEG2RAD (M_PI / 180.0)  // pi/180 (degrees to radians)
#define EPSILON 1e-12

// WGS-84 constants
#define RA   6.378137e+6      // earth equatorial radius in metres
#define FINV 298.257223563    // 1/f = inverse flattening of the ellipsoid

// WGS-84 derived constants
#define F   (1.0 / FINV)      // f    ~= 0.00335
#define F1  (1.0 - F)         // 1-f  ~= 0.9966
#define F16 (F / 16.0)        // f/16 ~= 0.00021
#define RB  (RA * F1)         // earth polar radius in metres
#define A2  (RA * RA)         // square equatorial radius
#define B2  (RB * RB)         // square polar radius
#define RF  ((A2 - B2) / B2)  // reduced a/b fraction

// Line segment on the Earth surface defined by 2 lat/lon points
// = four doubles addressable as either coor[0..3] or lat1,lon1,lat2,lon2
typedef union {
    double coor[4];
    struct { double lat1, lon1, lat2, lon2; };
} LineSegment;

static bool equal(const double a, const double b)
{
    return fabs(a - b) <= EPSILON;
}

// Ref.: https://en.wikipedia.org/wiki/Vincenty%27s_formulae
static double vincenty(const LineSegment a)
{
    if (equal(a.lat1, a.lat2) && equal(a.lon1, a.lon2))
        return 0;

    double U1 = atan(F1 * tan(a.lat1));
    double U2 = atan(F1 * tan(a.lat2));
    double L = a.lon2 - a.lon1;
    double l = L, l0, ss, cs, s, c2a, c2sm;
    do {
        l0 = l;
        double sU1 = sin(U1), cU1 = cos(U1);
        double sU2 = sin(U2), cU2 = cos(U2);
        double sl  = sin(l),  cl  = cos(l);
        double sU12 = sU1 * sU2, cU12 = cU1 * cU2;

        double p = cU2 * sl;
        double q = cU1 * sU2 - sU1 * cU2 * cl;
        ss = sqrt(p * p + q * q);
        cs = sU12 + cU12 * cl;
        s = atan2(ss, cs);
        double sa = cU12 * sl / ss;
        c2a = 1 - sa * sa;
        c2sm = cos(s) - 2 * sU12 / c2a;
        double C = F16 * c2a * (4 + F * (4 - 3 * c2a));
        l = L + (1 - C) * F * sa * (s + C * ss * (c2sm + C * cs * (-1 + 2 * c2sm * c2sm)));
    } while (!equal(l, l0));
    double u2 = c2a * RF;
    double t = sqrt(1 + u2);
    double k1 = (t - 1) / (t + 1);
    double k24 = 0.25 * k1 * k1;
    // double A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
    // double B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
    double A = (1 + k24) / (1 - k1);
    double B = k1 * (1 - 1.5 * k24);
    double ds = B * ss * (c2sm + (B / 4) * (cs * (-1 + 2 * c2sm * c2sm) - (B / 6) * c2sm * (-3 + 4 * ss * ss) * (-3 + 4 * c2sm * c2sm)));
    return RB * A * (s - ds);
}

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
        if (i == 1 || i == 3) {
            if (a.coor[i] < -180 || a.coor[i] > 180) {
                fprintf(stderr, "Longitude must be between -180 and +180: %s.\n", argv[j]);
                exit(ERR_LON180);
            } else if (equal(a.coor[i], -180))
                a.coor[i] = 180;
        }
        // From degrees to radians
        a.coor[i] *= DEG2RAD;
    }

    // Intermediate values for local earth radius
    // Ref.: https://en.wikipedia.org/wiki/Earth_radius#Location-dependent_radii
    double avglat = (a.lat1 + a.lat2) / 2;
    double s = sin(avglat);
    double c = cos(avglat);
    double rs = B2 * s * s;
    double rc = A2 * c * c;

    // Intermediate values for inverse haversine function
    // Ref.: https://en.wikipedia.org/wiki/Haversine_formula#Formulation
    double slat = sin((a.lat2 - a.lat1) / 2);
    double slon = sin((a.lon2 - a.lon1) / 2);

    // Great-circle distance in m, precision to cm
    // Ref.: https://en.wikipedia.org/wiki/Great-circle_distance#Computational_formulas
    double dist = 2 * sqrt((B2 * rs + A2 * rc) / (rs + rc))
                    * asin(sqrt(slat * slat + cos(a.lat1) * cos(a.lat2) * slon * slon));
    printf("%.2f\n", dist);
    printf("%.3f\n", vincenty(a));

    return 0;
}
