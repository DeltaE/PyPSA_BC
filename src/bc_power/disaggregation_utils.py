import math as m

'''
    FUNCTIONS TO HANDLE CONVERSION BETWEEN LAT/LON AND CARTESIAN COORDINATES

    Source: https://stackoverflow.com/questions/1185408/converting-from-longitude-latitude-to-cartesian-coordinates
'''

#Radius of Earth in km
r = 6371

#Convert a pair of latitude and longitude values to a point in a 3D cartesian coordinate system. Return a list to represent that point.
def to_cartesian(lat, lon):
    x = r * m.cos(m.radians(lat)) * m.cos(m.radians(lon))
    y = r * m.cos(m.radians(lat)) * m.sin(m.radians(lon))
    z = r * m.sin(m.radians(lat))

    return [x, y, z]

#Inverse of the above
def to_lat_lon(x, y, z):
    lat = m.asin(z / r)
    lon = m.atan2(y, x)

    return m.degrees(lat), m.degrees(lon)



#Shapely's built-in distance function only works on a 2D-plane (x, y). It did not account for z values when testing it out.
#p1 and p2 are lists that represent a point in 3D space like [x, y, z]
def dist_3d(p1, p2):
    x = m.pow((p1[0] - p2[0]), 2)
    y = m.pow((p1[1] - p2[1]), 2)
    z = m.pow((p1[2] - p2[2]), 2)

    return m.sqrt(x + y + z)