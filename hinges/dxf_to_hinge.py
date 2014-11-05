import dxfgrabber
import json
import math
from optparse import OptionParser
import sys

def main():
    parser = OptionParser()
    parser.add_option("-n", "--name", dest="name",
            help="Name of the hinge.  Program will load file called <name>.json and output <name>-face.dxf and <name>-temple.dxf.", metavar = "HINGENAME")

    (options, args) = parser.parse_args()

    if options.name == None:
        sys.stderr.write("Not enough arguments")
        parser.print_help()
        sys.exit(0)

    name = options.name
    facename = name + "-face.dxf"
    templename = name + "-temple.dxf"

    facedxf = dxfgrabber.readfile(facename)
    templedxf = dxfgrabber.readfile(templename)

    hinge = { 'description': name }
    process_contour(facedxf, "face", hinge)
    process_contour(templedxf, "temple", hinge)

    filename = name + ".json"
    with open(filename, 'w') as outfile:
        json.dump(hinge, outfile)

def close_poly(points):
    if abs(points[0][0] - points[-1][0]) < 0.001:
        points[-1] = [points[0][0], points[-1][1]]
    if abs(points[0][1] - points[-1][1]) < 0.001:
        points[-1] = [points[-1][0], points[0][1]]

def distance(p1, p2):
    return ((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2) ** 0.5

def process_contour(dxf, name, hinge):
    """ Take a DXF and produce an object that can be used to create the JSON spec for the hinge.
    ASSUMPTION: the outline of the hinge is in lines or polylines, and holes for the rivets or
    screws are represented as circles.
    ASSUMPTION: If a closed polyline is found, it is assumed to be the entire hinge outline."""
    polylines = [entity for entity in dxf.entities if (entity.dxftype == 'POLYLINE')]
    lwpolylines = [entity for entity in dxf.entities if (entity.dxftype == 'LWPOLYLINE')]
    lines = [entity for entity in dxf.entities if entity.dxftype == 'LINE']
    arcs = [entity for entity in dxf.entities if entity.dxftype == 'ARC']
    circles = [entity for entity in dxf.entities if entity.dxftype == 'CIRCLE']

    if len(polylines) == 0 and len(lines) == 0 and len(lwpolylines) == 0:
        sys.stderr.write("No outline found.  Hinge outline should consist of lines or polylines.")
        sys.exit(0)

    holes = [[c.center[0], c.center[1]] for c in circles]
    print 'holes', holes
    hinge[name + "_holes"] = holes


    for polyline in polylines:
        points = []
        bulges = polyline.bulge
        for i, p in enumerate(polyline.points):
            if(bulges[i] == 0):
                points.append([p[0], p[1]])
            else:
                points.append([p[0], p[1]])
                arc_points= bulge_to_points(p, polyline.points[i+1], bulges[i])
                for pt in arc_points:
                    points.append([pt[0], pt[1]])
        print points
        close_poly(points)
        hinge[name + '_contour'] = points
#       points = [[p[0], p[1]] for p in polyline.points()]
#        close_poly(points)
#        startpoint = points[0]
#        endpoint = points[-1]
#        if startpoint == endpoint: # closed polyline
#            hinge[name + '_contour'] = points
        return

    #TODO: Implement lines and multiple polylines

    for polyline in lwpolylines:
        points = []
        bulges = polyline.bulge
        for i, p in enumerate(polyline.points):
            if(bulges[i] == 0):
                points.append([p[0], p[1]])
            else:
                points.append([p[0], p[1]])
                arc_points= bulge_to_points(p, polyline.points[i+1], bulges[i])
                for pt in arc_points:
                    points.append([pt[0], pt[1]])
        print points
        close_poly(points)
        hinge[name + '_contour'] = points
        return

def bulge_to_points(p1, p2, bulge):
    theta = 4*math.atan(bulge) # Included angle of arc
    theta2 = math.pi/2 - abs(theta/2)
    angle_to_next = math.atan2(p2[1]-p1[1], p2[0]-p1[0]) # Angle of vector to next point
    angle_to_cp = angle_to_next - theta2
    chord = distance(p1, p2)
    sagatta = chord/2*abs(bulge)
    radius = ((chord/2)**2 + sagatta**2)/2*sagatta
    cp = (p1[0] + radius*math.cos(angle_to_cp), p1[1] + radius*math.sin(angle_to_cp))

    # Generate 20 points along arc
    start = math.atan2(p1[1]-cp[1], p1[0]-cp[0])
    arc_theta = [start + (theta/20)*i for i in range(20)][1:]
    return [(cp[0]+radius*math.cos(t), cp[1]+radius*math.sin(t)) for t in arc_theta]






if __name__ == '__main__':
    main()


