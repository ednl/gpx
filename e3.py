from math import fabs, sin, cos, asin, sqrt
import xml.etree.ElementTree as ET

EPS  = 1e-10
D2R  = 1.74532925199432953e-2
D2RH = 8.72664625997164656e-3

R_POL  = 6.357e+6
R_EQT  = 6.378e+6
R_POL2 = R_POL * R_POL
R_EQT2 = R_EQT * R_EQT

def eq(a, b):
    return fabs(a - b) < EPS

def diam(latdeg):
    if eq(latdeg, 0):
        return 2 * R_EQT
    if eq(latdeg, 90) or eq(latdeg, -90):
        return 2 * R_POL
    lat = latdeg * D2R
    s = sin(lat)
    c = cos(lat)
    rs2 = R_POL2 * s * s
    rc2 = R_EQT2 * c * c
    return 2 * sqrt((R_POL2 * rs2 + R_EQT2 * rc2) / (rs2 + rc2))

def dist(p1, p2):
    dlat = p2['lat'] - p1['lat']
    dlon = p2['lon'] - p1['lon']
    eqlat = eq(dlat, 0)
    eqlon = eq(dlon, 0)

    if eqlat and eqlon:
        return 0

    if eqlat:
        return diam(p1['lat']) * asin(fabs(cos(p1['lat'] * D2R) * sin(dlon * D2RH)))

    D = diam((p1['lat'] + p2['lat']) / 2)
    if eqlon:
        return D * asin(fabs(sin(dlat * D2RH)))

    slat = sin(dlat * D2RH)
    slon = sin(dlon * D2RH)
    return D * asin(sqrt(slat * slat + cos(p1['lat'] * D2R) * cos(p2['lat'] * D2R) * slon * slon))

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
        time += dist(p1, p2) / speed
    print(p2['lat'], p2['lon'], time)
    p1 = p2.copy()
