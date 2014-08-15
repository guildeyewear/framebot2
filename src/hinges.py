"""Manufacturing information for hinge features."""

import json
import poly

def get_hinge(pn, left=True):
    """
    Get the dictonary describing a hinge.
    By convention, you should only create descriptions for the left hinge.
    Right hinge descriptions can be derived using the flip_hinge function.
    Face coordinates are looking at the back of the face and have an origin with x axis on the vertical center of the hinge, and y axis on the inside edge of the left temple.
    Temple coordinates are looking at the inside of the left temple and have an origin with x axis on the vertical center of the hinge, and y axis on the back edge of the face.
    Temple and face coordinates are specified as if there is no pocketing.  Pocketing adjustments can be performed according to pocket depths.
    face_con and temple_con describe the contour of the face and temple pockets

    left says whether you want the left version or not, setting to false will return Right version

    'description': str(description),
    'drill_dia': float(drill_dia),
    'cb_dia': float(cb_dia),
    'cb_depth': float(cb_depth),
    'face_holes': make_points(face_holes),
    'temple_holes': make_points(temple_holes),
    """

    h = json.load(open("../hinges/%d.json" % pn))
    # hinges are modelled with y in reverse orientation from the mill
#    to_rotate = ['face_holes', 'face_contour']
#    to_flip = ['face_holes', 'face_contour'] #, 'temple_contour']

#    for key in to_flip:
#        h[key] = flip_x(h[key])

#    for key in to_rotate:
#        h[key] = poly.rotate(h[key], 90)

    #h['temple_contour'] = poly.rotate(h['temple_contour'], 90)
    #h['temple_contour'] = flip_y(h['temple_contour'])

    # Special processing to correct for modelling errors.
# NOTE: This is probably a result of the front plane of 6 degrees.  We must understand this better!!!
#    if pn == 1:
        # Rotational angle of temple hinge is slightly off
#        h['temple_holes'] = poly.rotate(h['temple_holes'], -11)
#        h['temple_contour'] = poly.rotate(h['temple_contour'], -11)
#        h['temple_holes'] = poly.translate(h['temple_holes'], -0.5, 0)
#        h['temple_contour'] = poly.translate(h['temple_contour'], -0.5, 0)
#    elif pn == 0:
#        h['temple_holes'] = poly.translate(h['temple_holes'], -1, -1)
#        h['temple_contour'] = poly.translate(h['temple_contour'], -1, -1)
    if left:
        return h
    else:
        return get_right_version(h)

def get_right_version(h):
    """Make a right hinge description from a left hinge description."""
    right_version = h.copy()

    to_flip = ['face_holes', 'temple_holes', 'face_contour', 'temple_contour']

    for key in to_flip:
        right_version[key] = flip_x(h[key])

    return right_version

def rotate(pl):
    return [[p[1], p[0]] for p in pl]

def flip_x(pl):
    """Flip the X coordinates of the provided point list across the Y axis."""
    return [[-1.0 * p[0], p[1]] for p in pl]

def flip_y(pl):
    """Flip the Y coordinates of the provided point list across the X axis."""
    return [[p[0], -p[1]] for p in pl]



