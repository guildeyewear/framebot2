import math
from collections import namedtuple
import time

Point = namedtuple('Point', 'x y')
Vector = namedtuple('Vector', 'x y')
PolarPoint = namedtuple('PolarPoint', 'theta r')
LineSegment = namedtuple('Line', 'p1 p2')
Rectangle = namedtuple('Rectangle', 'p1 p2') # p1 is top left, p2 is bottom right

tau = math.pi * 2.0  # Pi is wrong

class Polygon:
    def __init__(self, pts):
        pts = pts or []
        # Accept arrays of Points, lists or tuples
        self.points = [Point(pt[0], pt[1]) for pt in pts]
        if self.points[0] != self.points[-1]:
            self.points.append(Point(self.points[0].x, self.points[0].y))
        if len(self.points) >= 3 and ccw(self.points[0], self.points[1], self.points[2]) <= 0:
            self.points = self.points[::-1]
    def line_segments(self):
        return [LineSegment(a, b) for (a, b) in zip(self.points[0:-1], self.points[1:])]
    def left_extent(self):
        return min([pt.x for pt in self.points])
    def right_extent(self):
        return max([pt.x for pt in self.points])
    def top_extent(self):
        return max([pt.y for pt in self.points])
    def bottom_extent(self):
        return min([pt.y for pt in self.points])
    def bounding_box(self):
        return Rectangle(Point(self.left_extent(), self.top_extent()), Point(self.right_extent(), self.bottom_extent()))

def point_to_polar(pt):
    return PolarPoint(math.atan2(pt[1], pt[0]), math.sqrt(pt[0]*pt[0] + pt[1]*pt[1]))

def polar_to_point(polar):
    return Point(polar.r*math.sin(polar.theta), polar.r*math.cos(polar.theta))

def translate_polygon(polygon, vec):
    return Polygon([(pt.x+vec.x, pt.y+vec.y) for pt in polygon.points])

def triangle_area(a, b, c):
    """ Return the area of a triangle.  Area is positive if the
        points are counterclockwise, negative if the points
        are clockwise, and 0 if the points are colinear."""
    return ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y))/2

def ccw(a, b, c):
    """ Return true if the <abc forms a counterclockwise turn, false otherwise."""
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)

def colinear(a, b, c):
    """ Returns true if a, b, and c are colinear, false otherwise"""
    return triangle_area(a, b, c) == 0

def lines_intersect(seg1, seg2):
    """ Returns true if seg1 intersects seg2, false otherwise.
        Currently ignores the following degeneracies:
        1) one or both line segments are single points, i.e. p1 == p2
        2) one of the endpoints of a segment is colinear with the other segment"""
    a, b = seg1
    c, d = seg2
    if (( a.x > c.x and a.x > d.x and b.x > c.x and b.x > d.x) or
        ( a.x < c.x and a.x < d.x and b.x < c.x and b.x < d.x) or
        ( a.y > c.y and a.y > d.y and b.y > c.y and b.y > d.y) or
        ( a.y < c.y and a.y < d.y and b.y < c.y and b.y < d.y) ):
        return False
   # or
   #      (min(seg1.p1.x, seg1.p2.x) > max(seg2.p1.x, seg2.p2.x)) or
   #      (max(seg1.p1.y, seg1.p2.y) < min(seg2.p1.y, seg2.p2.y)) or
   #     (min(seg1.p1.y, seg1.p2.y) < max(seg2.p1.y, seg2.p2.y)) ):
   #     return False
    if ccw(seg1.p1, seg1.p2, seg2.p1) * ccw(seg1.p1, seg1.p2, seg2.p2) > 0:
        return False
    elif ccw(seg2.p1, seg2.p2, seg1.p1) * ccw(seg2.p1, seg2.p2, seg1.p2) > 0:
        return False
    return True

def line_intersection(seg1, seg2):
    """ Returns a point giving the interection point of two line segments.
        Returns None if the segments do not intersect."""
    if not lines_intersect(seg1, seg2):
        return None
    a, b = seg1
    c, d = seg2

    r = (((a.y-c.y)*(d.x-c.x) - (a.x-c.x)*(d.y-c.y)) /
        ((b.x-a.x)*(d.y-c.y) - (b.y-a.y)*(d.x-c.x)))
    return Point(a.x + r*(b.x-a.x), a.y + r*(b.y-a.y))

def line_polygon_intersection(segment, polygon):
    """ Returns a list of the points where segment intersects polygon,
        or an empty tuple if there are no intersections."""
    intersections =  [line_intersection(segment, pseg) for pseg in polygon.line_segments()]
    return [pt for pt in intersections if pt is not None]

def polygon_contains_point(polygon, pt):
    """ Returns true if the point lies inside the polygon.
        First, returns False if the pt is outside the bounding box
        of the polygon.  Then, creates a horizontal line from the point
        to just outside the right bound of the polygon.  The point
        is inside the polygon if the line intersects the polygon an
        odd number of times."""
    bound = polygon.bounding_box()
    if pt.x < bound.p1.x or pt.x > bound.p2.x or pt.y > bound.p1.y or pt.y < bound.p2.y:
        return False
    seg = LineSegment(pt, Point(bound.p2.x + 10, pt.y))
    return (len(line_polygon_intersection(seg, polygon))%2 == 1)

def polygon_to_uniform_polar(polygon, point_count):
    """ Takes a polygon and converts it to a list of point_count points in polar
        form, with a uniform angle between each point.  A polygon drawn from the
        result will not be exactly the same as the input polygon since this function
        doesn't necessarily hit all the vertices.
        Assumes that the origin is contained in the polygon."""
    time1 = time.time()
    if not polygon_contains_point(polygon, Point(0.0, 0.0)):
        return []
    time2 = time.time()
    bound = polygon.bounding_box()
    r = max(abs(bound.p1.x), abs(bound.p1.y), abs(bound.p2.x), abs(bound.p2.y)) + 10.0

    interval = tau/point_count
    point_circle = [polar_to_point(PolarPoint(interval*i, r)) for i in range(point_count)]
    intersections = [line_polygon_intersection(LineSegment(Point(0,0), pt), polygon) for pt in point_circle]
    polar_form = [point_to_polar(pt) for pts in intersections for pt in pts]
    polar_form.sort(key=lambda pt: pt.theta)
    time3 = time.time()
    return [PolarPoint(pt.theta+tau/2, pt.r) for pt in polar_form]


def _flatten(l):
    '''Flattens a tree structure to a list of strings.  Nodes with a value of None are removed.'''
    if type(l) is type(""):
        return [l]
    elif type(l) is type(None):
        return []
    else:
        r = []
        for e in l:
            r += _flatten(e)
        return r



