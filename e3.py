from math import fabs, sin, cos, asin, sqrt
import xml.etree.ElementTree as ET

EPSILON  = 1e-10
DEG2RAD  = 1.74532925199432953e-2
DEG2RAD_HALF = 8.72664625997164656e-3

RADIUS_POLAR  = 6.357e+6
RADIUS_EQUATORIAL  = 6.378e+6
RADIUS_POLAR_SQUARED = RADIUS_POLAR * RADIUS_POLAR
RADIUS_EQUATORIAL_SQUARED = RADIUS_EQUATORIAL * RADIUS_EQUATORIAL

def isequal(a, b):
    return fabs(a - b) < EPSILON

def earth_diameter(latdeg):
    if isequal(latdeg, 0):
        return 2 * RADIUS_EQUATORIAL
    if isequal(latdeg, 90) or isequal(latdeg, -90):
        return 2 * RADIUS_POLAR
    lat = latdeg * DEG2RAD
    slat = sin(lat)
    clat = cos(lat)
    radius_sin_sqr = RADIUS_POLAR_SQUARED * slat * slat
    radius_cos_sqr = RADIUS_EQUATORIAL_SQUARED * clat * clat
    numer = RADIUS_POLAR_SQUARED * radius_sin_sqr + RADIUS_EQUATORIAL_SQUARED * radius_cos_sqr
    denom = radius_sin_sqr + radius_cos_sqr
    return 2 * sqrt(numer / denom)

def distance(p1, p2):
    dlat = p2['lat'] - p1['lat']
    dlon = p2['lon'] - p1['lon']
    eqlat = isequal(dlat, 0)
    eqlon = isequal(dlon, 0)

    if eqlat and eqlon:
        return 0

    if eqlat:
        return earth_diameter(p1['lat']) * asin(fabs(cos(p1['lat'] * DEG2RAD) * sin(dlon * DEG2RAD_HALF)))

    D = earth_diameter((p1['lat'] + p2['lat']) / 2)
    if eqlon:
        return D * asin(fabs(sin(dlat * DEG2RAD_HALF)))

    slat = sin(dlat * DEG2RAD_HALF)
    slon = sin(dlon * DEG2RAD_HALF)
    return D * asin(sqrt(slat * slat + cos(p1['lat'] * DEG2RAD) * cos(p2['lat'] * DEG2RAD) * slon * slon))

tree = ET.parse('e3.gpx')
root = tree.getroot()
ns = {'': 'http://www.topografix.com/GPX/1/1'}

p1 = None
time = 0
speed = 40 / 3.6
first = True
for trkpt in root.findall('.//trk/trkseg/trkpt', ns):
    p2 = dict((k, float(v)) for k, v in trkpt.attrib.items())
    if first:
        first = False
    else:
        time += distance(p1, p2) / speed
    print(p2['lat'], p2['lon'], time)
    p1 = p2.copy()
