"""Creates and performs operations on polygons."""

import math
import pyclipper
import sys
from cStringIO import StringIO
from decimal import Decimal

tau = math.pi * 2.0  # That's right, I said it

"""
Today we're adding and switching to length-along-polyline functions.  This should make deriving more
complicated toolpaths from polylines easier (e.g. segments forming tabs) because we will not have to
modify a polyline through a sequence of steps, but instead extract information about a line and then
extract the final required features (e.g. gather length-along-line points, then extract segments with
gaps at these lengths).

It will also enable some interesting functionality, like tab-every-inch-of-polyline segments.

"""
def intersection(poly1, poly2):
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    scale_factor = 1000.0
    scaled1 = scale(poly1, scale_factor);
    scaled2 = scale(poly2, scale_factor);
    c = pyclipper.Pyclipper()
    c.add_polygon(scaled1)
    c.sub_polygon(scaled2)
    result = c.execute(0)
    sys.stdout = old_stdout
    if result and len(result) > 0:
        return [[p[0]/scale_factor, p[1]/scale_factor] for p in result[0]]
    else:
        return result
    #return [[p[0]/scale_factor, p[1]/scale_factor] for p in result[0]];

def difference(poly1, poly2):
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    scale_factor = 1000.0
    scaled1 = scale(poly1, scale_factor);
    scaled2 = scale(poly2, scale_factor);
    c = pyclipper.Pyclipper()
    c.add_polygon(scaled1)
    c.sub_polygon(scaled2)
    result = c.execute(3)
    sys.stdout = old_stdout
    if result and len(result) > 0:
        return [[p[0]/scale_factor, p[1]/scale_factor] for p in result[0]]
    else:
        return result
    #return [[p[0]/scale_factor, p[1]/scale_factor] for p in result[0]];

def union(poly1, poly2):
    old_stdout = sys.stdout
    sys.stdout = mystdout = StringIO()
    scale_factor = 1000.0
    scaled1 = scale(poly1, scale_factor);
    scaled2 = scale(poly2, scale_factor);
    c = pyclipper.Pyclipper()
    c.add_polygon(scaled1)
    c.sub_polygon(scaled2)
    result = c.execute(2)
    sys.stdout = old_stdout
    if result and len(result) > 0:
        return [[p[0]/scale_factor, p[1]/scale_factor] for p in result[0]]
    else:
        return result

def line_length(line):
    """Returns length of a line segment [[x0, y0], [x1, y1]]"""
    dx = line[0][0] - line[1][0]
    dy = line[0][1] - line[1][1]

    l = math.sqrt(dx**2.0 + dy**2.0)

    return l

def polyline_length(polyline, closed = True):
    """Returns total length of a polyline.  Includes a line segment from the last vertex to the first if closed=True (the default)."""
    lengths = [line_length(pair) for pair in pairs(polyline, closed)]

    l = sum(lengths)

    return l

def x_intercept(x, line):
    """
    Finds intersection point between given line segment and vertical line at x.

    The possible return values are:
    [], no points of intersection.
    [[x, y]], a single point of intersection (note that this point can be the same point as the start or end of the line segment).
    list(line), co-incident with the provided line segment, i.e. a continuous intersection with the provided line segment.  In this case a copy of the provided line segment is returned.
    """

    (x1, y1) = line[0]
    (x2, y2) = line[1]

    # Continuous intersection with provided line segment
    if x == x1 and x1 == x2:
        return list(line) # A copy of the provided line segment

    # Intersect directly with only the first point
    elif x == x1:
        return [[x1, y1]]

    # Intersect directly with only the second point
    elif x == x2:
        return [[x2, y2]]

    # Intersect somewhere along the line
    elif (x1 < x and x2 > x) or (x2 < x and x1 > x):
        xi = x
        yi = y1 + ((x - x1)/(x2 - x1)) * (y2 - y1)
        return [[xi, yi]]

    # No intersection
    else:
        return []

def y_intercept(y, line):
    """
    Finds intersection point between given line segment and horizontal line at y.
    Please see x_intercept.
    """

    r = x_intercept(y, rotate_90(line, True)) # CW rotation

    return rotate_90(r, False) # CCW rotation (undo previous CW)

def intercepts(polyline, closed, x=[], y=[]):
    """
    Finds x and/or y intercept points on a polyline, expressed in lengths from the start of the polyline.

    Returns a list of lengths along the polyline where the intercepts are found.

    If intercept is coincident with a line on the polyline, the two points that make up that line segment are treated as separate intersections.
    """

    intercept_lengths = []

    l = 0.0 # Length accumulator

    for p in pairs(polyline, closed):
        intercept_points = []
        for yi in y:
            intercept_points += y_intercept(yi, p)
        for xi in x:
            intercept_points += x_intercept(xi, p)
        for point in intercept_points:
            intercept_lengths += [l + line_length([p[0], point])]

        l += line_length(p)

    return intercept_lengths

def rotate(poly, angle):
    if len(poly) == 0:
        return []
    rad = math.radians(angle)
    cosx = math.cos(rad)
    sinx = math.sin(rad)
    if len(poly[0]) > 2:
        return [[p[0]*cosx - p[1]*sinx, p[1]*cosx + p[0]*sinx, p[2]] for p in poly]
    else:
        return [[p[0]*cosx - p[1]*sinx, p[1]*cosx + p[0]*sinx] for p in poly]

def rotate_90(poly, ccw=True):
    """Rotates points in a polyline by 90 degrees counterclockwise (CCW) or clockwise (CW) depending on ccw flag (default CCW, or True)."""

    if ccw:
        return [[p[1], -p[0]] for p in poly]
    else:
        return [[-p[1], p[0]] for p in poly]

def scale(polygon, scale):
    return [[round(p[0]*scale), round(p[1]*scale)] for p in polygon]

def dilate(r, polygon):
    #old_stdout = sys.stdout
    #sys.stdout = mystdout = StringIO()
    scale_factor = 100.0


    scaled = scale(polygon, scale_factor);

    if(not is_ccw(scaled)):
        scaled = reverse(scaled)
    offset = pyclipper.offset([scaled], r*scale_factor, jointype=2);

    #sys.stdout = old_stdout
    return [[p[0]/scale_factor, p[1]/scale_factor] for p in offset[0]];
    """Return the provided polygon, dilated by r."""
    # return m2d.run(circle(r), polygon, mode="dilate")


def erode(r, polygon, jointype=2):
    #old_stdout = sys.stdout
    #sys.stdout = mystdout = StringIO()
    scale_factor = 1000.0
    scaled = scale(polygon, scale_factor);
    if(not is_ccw(scaled)):
        scaled = reverse(scaled)

    offset = pyclipper.offset([scaled], -r*scale_factor, jointype);

    #sys.stdout = old_stdout

    if len(polygon[0]) > 2:
        return  [[[p[0]/scale_factor, p[1]/scale_factor] for p in off] for off in offset]
        #return  [[[p[0]/scale_factor, p[1]/scale_factor, p[2]] for p in off] for off in offset]
    #Holes are messing us up, only return the first eroded path.

    return  [[[p[0]/scale_factor, p[1]/scale_factor] for p in off] for off in offset]
    #return  [[[p[0]/scale_factor, p[1]/scale_factor] for p in off] for off in offset]

#    return m2d.run(circle(r), polygon, mode="erode")

def circle(r, n=32):
    """
    Generate a counter-clockwise (CCW) polygon roughly representing a circle with radius r.
    Argument n specifies the number of vertices in the result.
    """
    step = tau / float(n)

    angles = [step * float(i) for i in range(0, n)]

    points = [[math.cos(a) * float(r), math.sin(a) * float(r)] for a in angles]

    return points

def arc(r, starta, enda, n=32):
    """
    Generate a counter-clockwise (CCW) polygon roughly representing a circular arc with radius r, from angle starta to angle enda.
    Argument n specifies the number of vertices in the result.
    """
    step = (enda - starta) / float(n)

    angles = [starta + step * float(i) for i in range(0, n)]

    points = [[math.cos(a) * float(r), math.sin(a) * float(r)] for a in angles]

    return points

def shorter_segment(line, length):
    """Returns a segment of a segment that is 'length' long."""
    ratio = abs(length / line_length(line))

    #dx = line[0][0] - line[1][0]
    #dy = line[0][1] - line[1][1]

    dx = line[1][0] - line[0][0]
    dy = line[1][1] - line[0][1]

    r = [line[0][0] + (ratio * dx), line[0][1] + (ratio * dy)]

    return r

def point_at_length(polyline, length, closed):
    """Returns the coordinates of a point at length distance along the polyline."""

    # degenerate case: requested length 0
    if length == 0:
        return polyline[0]


    if length > polyline_length(polyline, closed):
        raise Exception("Requested length %f is longer than the polyline's length of %f" % (length, polyline_length(polyline, closed)))

    l = 0.0

    for p in pairs(polyline, closed):
        l += line_length(p)
        if l == length:
            return p[1]
        elif l > length:
            l -= line_length(p) # Roll it back to length at start point
            r = shorter_segment(p, length-l) # Now get the point that gives us the requested length
            return r

def start_point(contour):
    return contour[0]

def new_start(poly, n):
    """Returns a new polygon with poly[n] as its start point."""
    return poly[n:] + poly[:n]

def segment(polyline, start, end, closed):
    """Returns a new polyline which is a segment of polyline from length 'start' to length 'end'."""

    l = 0.0

    # Grab our start and end points
    start_p = point_at_length(polyline, start, closed)
    end_p = point_at_length(polyline, end, closed)
    middle = []

    # Get all our inbetween points
    for p in pairs(polyline, closed):
        l += line_length(p)

        # If one of the endpoints is the same as start_p or end_p, it's automatically
# disqualified from being in the middle
        disqualified = False
        for pt in p:
            if points_are_equal(pt, start_p) or points_are_equal(pt, end_p):
                disqualified = True
        if disqualified:
            continue

        if l > start and l < end:
            middle += [p[0]]
        elif l > end:
            break
    return [start_p] + middle + [end_p]

# Make sure the polygon has a point at length 'length'.  Simplifies
# degenerate cases when modifying polygons later to add tabs
def insert_point_at_length(polyline, length, closed):
    r = [polyline[0]]
    inserted = False
    l = 0.0
    for p in pairs(polyline, closed):
        l += line_length(p)
        if l < length:
            r += [p[1]]
        elif not inserted:
            r += [point_at_length(polyline, length, closed)]
            r += [p[1]]
            inserted = True
        else:
            r += [p[1]]
    return r

def points_are_equal(p1, p2):
    if abs(p1[0] - p2[0]) < 0.001 and abs(p1[1] - p2[1]) < 0.001:
        return True
    return False

def make_ramped_tabs(polyline, lengths, tab_width, bottom, top, closed):
    r = []
    start = 0.0

    # Insert points at the relevant transition points.
    for l in sorted(lengths):
        tab_start = l - tab_width/2.0
        tab_end = l + tab_width/2.0
        polyline = insert_point_at_length(polyline, tab_start, closed)
        polyline = insert_point_at_length(polyline, l, closed)
        polyline = insert_point_at_length(polyline, tab_end, closed)

    for l in sorted(lengths):
        tab_start = l - tab_width/2.0
        tab_end = l + tab_width/2.0
        free = set_z(segment(polyline, start, tab_start, False), bottom)
        tab_sg1 = segment(polyline, tab_start, l, False)
        tab_sg2 = segment(polyline, l, tab_end, False)
        ramp1 = ramp(tab_sg1, bottom, top, False)
        ramp2 = ramp(tab_sg2, top, bottom, False)
        r += free
        r += ramp1
        r += ramp2
        start = tab_end

    # Tail end
    tail = set_z(segment(polyline, start, polyline_length(polyline, closed), closed), bottom)
    r += tail
    return r



def make_square_tabs(polyline, lengths, tab_width, bottom, top, closed):
    r = []
    start = 0.0

# Insert points at the relevant transition points.
    for l in sorted(lengths):
        tab_start = l - tab_width/2.0
        tab_end = l + tab_width/2.0
        polyline = insert_point_at_length(polyline, tab_start, closed)
        polyline = insert_point_at_length(polyline, tab_end, closed)

    for l in sorted(lengths):
        tab_start = l - tab_width/2.0
        tab_end = l + tab_width/2.0
        free = set_z(segment(polyline, start, tab_start, False), bottom)
        tab_sg = segment(polyline, tab_start, tab_end, False)
        tab = set_z(tab_sg, top)
        r += free
        r += tab

        start = tab_end

    # Tail end
    tail = set_z(segment(polyline, start, polyline_length(polyline, closed), closed), bottom)
    r += tail
    return r


def make_tabs(polyline, lengths, tab_width, bottom, top, closed):
    """
    Adds steps to the 'polyline' at the given 'lengths' with a width of 'gap'.
    Non-gap pieces have a height of 'bottom', gap pieces have a height of 'top'.
    """

    r = []
    start = 0.0

    for l in sorted(lengths):
        tab_start = l - tab_width/2.0
        tab_end = l + tab_width/2.0
#        ramp_up_end = tab_start + gap/4.0
#        ramp_down_start = tab_end - gap/4.0

        free = set_z(segment(polyline, start, tab_start, False), bottom)
#        ramp_up = ramp2(segment(polyline, tab_start, ramp_up_end, False), bottom, top)
        tab1 = segment(polyline, tab_start, l, False)
        tab2 = segment(polyline, l, tab_end, False)

        tab1 = ramp(tab1[1:], bottom, top, False)
        tab2 = ramp(tab2[1:], top, bottom, False)

#        ramp_down = ramp2(segment(polyline, ramp_down_start, tab_end, False), top, bottom)

#        tab = set_z(segment(polyline, tab_start, tab_end, False), top)
#        ramp_points = len(tab)/4
#        z_interval = (top-bottom) / ramp_points
#        z_levels = [i*z_interval for i in range(ramp_points)]
#        for i, z_level in enumerate(z_levels):
#            tab[i][2] = bottom + z_level
#            tab[-i][2] = bottom + z_level

#        ramp2(segment(polyline, tab_start, ramp_up_end, False), bottom, top)

        r += free
#        r += ramp_up
        r += tab1
        r += tab2
#        r += ramp_down

        start = tab_end

    # Tail end
    tail = set_z(segment(polyline, start, polyline_length(polyline, closed), closed), bottom)
    r += tail
    return r

def ramp(polyline, start_height, end_height, closed):
    """
    Adds a third coordinate, z, to the polyline points, from 'top' to 'bottom', proportional to the length along the line.
    Will add an extra point that is a duplicate of the first point if closed=True (the default).
    """
    #print "ramp(", start_height, end_height, closed, ")"
    max = polyline_length(polyline, closed)

    r = []
    for (p, l) in lengths(polyline, closed):
        r.append(p + [start_height + (l/max) * (end_height - start_height)])

    return r

def ramp2(polyline, start_z, end_z):
    step= (end_z - start_z)/len(polyline)
    z = [start_z + i*step for i in range(0, len(polyline))]
    ramp = [[p[0], p[1], z] for p, z in zip(polyline, z)]
    return ramp

def set_z(polyline, z):
    """Sets the z coordinate to the polyline of 'z'."""
    return [p + [z] for p in polyline]

def mirror_y(contour, closed):
    """Returns a polygon with mirrored y coordinates, i.e. mirrored across x axis."""
    if len(contour[0]) > 2:
        r = [[p[0], -p[1], p[2]] for p in contour]
    else:
        r = [[p[0], -p[1]] for p in contour]

    if closed:
        r = reverse(r)

    return r


def mirror_x(contour, closed):
    """Returns a polygon with mirrored x coordinates, i.e. mirrored across y axis."""
    if len(contour[0]) > 2:
        r = [[p[0] * -1.0, p[1], p[2]] for p in contour]
    else:
        r = [[p[0] * -1.0, p[1]] for p in contour]

    if closed:
        r = reverse(r)

    return r

def reverse(contour, closed=True):
    if closed:
        return [contour[0]] + contour[1:][::-1]  # Contour[0] is added to the front to maintain the start point
    return contour[::-1]

def translate(poly, dx, dy):
    if len(poly) == 0:
        return []
    if len(poly[0]) > 2:
        return [[p[0] + dx, p[1] + dy, p[2]] for p in poly]
    else:
        return [[p[0] + dx, p[1] + dy] for p in poly]

def bottom(poly):
    """Returns the lowest Y value of the polygon points."""
    return min([p[1] for p in poly])

def top(poly):
    """Returns the highest Y value of the polygon points."""
    return max([p[1] for p in poly])

def right(poly):
    """Returns the highest X value of the polygon points."""
    return max([p[0] for p in poly])


def left(poly):
    """Returns the lowest X value of the polygon points."""
    return min([p[0] for p in poly])

# returns True if p1 is completely contained in p2
def is_contained(poly1, poly2):
    contained = top(poly2) < top(poly1) and \
                bottom(poly2) > bottom(poly1) and \
                left(poly2) > left(poly1) and \
                right(poly2) < right(poly1)
    return contained


def leftmost_index(poly):
    values = [p[0] for p in poly]
    return values.index(min(values))

def bottommost_index(poly):
    values = [p[1] for p in poly]
    return values.index(min(values))

def closest_index(poly, to_point):
    values = [abs(line_length([p, to_point])) for p in poly]
    return values.index(min(values))

def pairs(l, closed):
    """
    Generates pairs of items from a list.  Can be used to generate pairs of points around a polyline.
    Will pair the last item with the first if 'closed' is True (the default).
    """

    for i in range(1, len(l)):
        yield [l[i-1], l[i]]

    if closed:
        yield [l[-1], l[0]]

def area(poly):
    """Returns the area of the polygon.  Mostly clockwise polygons have positive areas, mostly counterclockwise polygons have negative areas."""

    a = 0.0

    for pair in pairs(poly, True):
        (x1, y1) = pair[0]
        (x2, y2) = pair[1]

        a += (x2-x1) * (y2+y1)

    return a/2.0

def is_cw(poly):
    return True if area(poly) >= 0.0 else False

def is_ccw(poly):
    return True if area(poly) <= 0.0 else False

def lengths(polyline, closed):
    """
    Returns the points of the provided line along with the distance along the line at that point in the format (point, length at that point).
    """

    l = 0.0

    yield [polyline[0], l] # First point is always 0.0 along the line

    for p in pairs(polyline, closed):
        l += line_length(p)
        yield [p[1], l]


