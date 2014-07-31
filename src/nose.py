"""Support routines for generating nose contours."""

import math
import poly

tau = math.pi * 2.0

def tangent_point(sa, r):
    x = math.cos(sa/2.0) * r
    y = math.sin(sa/2.0) * r

    return [x, y]

def nose_poly(r, h, sa, ra, xfloor, erode=0.0, depth=0.0):
    """
    Generates an open contour polygon to match the described nose.
    r, h, sa and ra are standard radius, height of intercept, splay angle and ridge angle parameters (in radians)
    yfloor is the end point in the y direction of the nose side lines, for getting the right size of contour
    optional erode parameter will erode the contour by the given amount, usually the radius of the cutting tool
    optional depth parameter will generate a contour that matches the nose profile at the given depth
    """

    nose_radius = r
    r -= erode
    h -= erode
    half_sa = sa / 2.0
    xfloor = float(xfloor)
    depth = float(depth)
    depth_displacement = -depth/math.tan(ra);
    translation = h - depth_displacement - nose_radius
    xfloor -= translation

    p2 = tangent_point(sa, r)
    p1 = [p2[0] + math.sin(half_sa) * abs(xfloor) / math.sin(tau/4.0 - half_sa), xfloor]
    rpoly = poly.arc(r, half_sa, tau/2.0 - half_sa)
    p3 = [-p2[0], p2[1]]
    p4 = [-p1[0], p1[1]]

    base = [p1] + [p2] + rpoly + [p3] + [p4]

    ret = poly.translate(base, 0.0, translation)

    return poly.rotate_90(ret)
