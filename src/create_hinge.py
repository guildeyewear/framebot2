'''DXF to JSON utility.'''

import dxfgrabber
import math
import json
import poly
from optparse import OptionParser
import sys

def main():
    """Run as a command line program."""

    parser = OptionParser()
    parser.add_option("-f", "--face", dest="face",
                      help="name of face dxf file", metavar="FACEDXF")
    parser.add_option("-t", "--temple", dest="temple",
                      help="name of temple dxf file", metavar="TEMPLEDXF")
    parser.add_option("-n", "--name", dest="name",
                      help="name of hinge", metavar="HINGENAME")
    parser.add_option("-o", "--outfile", dest="outfile",
                      help="name of output file", metavar="OUTFILE")

    (options, args) = parser.parse_args()

    if options.face == None or options.temple == None or options.name == None or options.outfile == None:
        sys.stderr.write("Insufficient arguments.")
        parser.print_help()
        sys.exit(0)

    facedxf = dxfgrabber.readfile(options.face)
    templedxf = dxfgrabber.readfile(options.temple)
    spec = {'description': options.name, 'drill_dia': 1.0, 'cb_dia': 2.0, 'cb_depth': 0.5, 'pocket_depth': 1.0}
    process_contour('face', facedxf, spec)
    process_contour('temple', templedxf, spec)
    with open(options.outfile, 'w') as outfile:
          json.dump(spec, outfile)

def is_same_point(pt1, pt2):
    if(len(pt1) < 2 or len(pt2) < 2): return False
    return abs(pt1[0]-pt2[0]) < .01 and abs(pt1[1]-pt2[1]) < .01

def get_points(e, startpoint):
    points = []
    if isinstance(e, dxfgrabber.entities.Line):
        points += line_to_points(e)
    elif isinstance(e, dxfgrabber.entities.Polyline):
        points += polyline_to_points(e)
    elif isinstance(e, dxfgrabber.entities.Arc):
        points += arc_to_points(e)
    else:
        raise Exception("Unrecognized DXF entity %s." % str(e))
    if not is_same_point(startpoint, points[0]):
        print 'points different', startpoint, points[0]
        print 'reversing'
        points.reverse()
    return points
def get_start_point(e):
    if isinstance(e, dxfgrabber.entities.Line):
        return e.start
    elif isinstance(e, dxfgrabber.entities.Polyline):
        return e.start
    elif isinstance(e, dxfgrabber.entities.Arc):
        ang = math.radians(e.startangle)
        return [e.center[0]+math.cos(ang)*e.radius, e.center[1]+math.sin(ang)*e.radius]
    return []
def get_end_point(e):
    if isinstance(e, dxfgrabber.entities.Line):
        return e.end
    elif isinstance(e, dxfgrabber.entities.Polyline):
        return e.end
    elif isinstance(e, dxfgrabber.entities.Arc):
        ang = math.radians(e.endangle)
        return [e.center[0]+math.cos(ang)*e.radius, e.center[1]+math.sin(ang)*e.radius]
    return []

def process_contour(name, dxf, spec):

    # need to flip because we've modeled the right-side hinges but display left-side
#    def flip(points):
#        return[[-p[0], p[1]] for p in points]

    points = []
    # find the first entity that isn't a circle
    e = next(ent for ent in dxf.entities if (isinstance(ent, dxfgrabber.entities.Line) or isinstance(ent, dxfgrabber.entities.Polyline)))
    firstpoint = e.start
    lastpoint = e.end
    entities = dxf.entities.get_entities()

    entities.remove(e)
    points.extend(get_points(e, firstpoint))
    print "start:", points

    def get_connected_entity(entity_list, point):
        print 'looking for connection to ', point
        print entity_list
        if not entity_list:
            return None
        for ent in entity_list:
            startpoint = get_start_point(ent)
            endpoint = get_end_point(ent)
            print 'entity:', ent, startpoint, endpoint
            if is_same_point(get_start_point(ent), point) or is_same_point(get_end_point(ent), point):
                return ent
        return None

    while e and not is_same_point(lastpoint, firstpoint):
        try:
            e = get_connected_entity(entities, lastpoint)
            print 'got entity', e
            print get_points(e, lastpoint)
            points.extend(get_points(e, lastpoint))
            lastpoint = points[-1]
            print 'last poitn is ', lastpoint
            entities.remove(e)
        except StopIteration:
            e = None
            break
    points = filter_duplicate_points(points);
    points.append(points[0]);
#    contour = flip(points)
#    if name == 'face':
        # We messed up the origin in the flip, have to translate back
#        points = [[p[0]+4, p[1]] for p in points]

    contour = points
    if not poly.is_ccw(contour):
        contour.reverse()
    contour_name = name + "_contour"
    spec[contour_name] = contour

    entities = dxf.entities.get_entities()
    holes = [flatten_point(c.center) for c in [e for e in entities if isinstance(e, dxfgrabber.entities.Circle)]]
    holes_name = name + "_holes"
    spec[holes_name] = holes



def flatten_point(point):
    return [point[0], point[1]]

def arc_to_points(arc):
    def get_coords(angle):
        ang = math.radians(angle)
        return [arc.center[0] + math.cos(ang)*arc.radius, arc.center[1] + math.sin(ang)*arc.radius]
    def get_range():
        increment = 2
        start = arc.startangle
        end = arc.endangle
        while(start > end):
            end += 360
        return range(int(math.ceil(start)), int(math.ceil(end)), increment)
    print 'getting points for arc', arc.startangle, arc.endangle, arc.center
    print get_range()
    points = [get_start_point(arc)] + [get_coords(angle) for angle in get_range()] + [get_end_point(arc)]
    return points

def polyline_to_points(entity):

    points = [flatten_point(p) for p in entity.points()]

    return points

def line_to_points(entity):
    points = [
        flatten_point(entity.start),
        flatten_point(entity.end)
    ]

    return points

def filter_duplicate_points(points):
    filtered = [points[i] for i in range(1, len(points)) if points[i-1] != points[i]]

    return filtered

def to_polygon(filename):
    """Returns list of [x,y] points from DXF."""

    dxf = dxfgrabber.readfile(filename)

    points = []

    for e in dxf.entities:
        if isinstance(e, dxfgrabber.entities.Line):
            points += line_to_points(e)
        elif isinstance(e, dxfgrabber.entities.Polyline):
            points += polyline_to_points(e)
        elif isinstance(e, dxfgrabber.entities.Arc):
            points += arc_to_points(e)
        elif isinstance(e, dxfgrabber.entities.Circle):
            pass
        else:
            raise Exception("Unrecognized DXF entity %s." % str(e))

    return filter_duplicate_points(points)

if __name__ == '__main__':
    main()
