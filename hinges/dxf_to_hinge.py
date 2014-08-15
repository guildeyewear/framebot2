import dxfgrabber
import json
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


def process_contour(dxf, name, hinge):
    """ Take a DXF and produce an object that can be used to create the JSON spec for the hinge.
    ASSUMPTION: the outline of the hinge is in lines or polylines, and holes for the rivets or
    screws are represented as circles.
    ASSUMPTION: If a closed polyline is found, it is assumed to be the entire hinge outline."""
    polylines = [entity for entity in dxf.entities if (entity.dxftype == 'POLYLINE')]
    lwpolylines = [entity for entity in dxf.entities if (entity.dxftype == 'LWPOLYLINE')]
    lines = [entity for entity in dxf.entities if entity.dxftype == 'LINE']
    circles = [entity for entity in dxf.entities if entity.dxftype == 'CIRCLE']

    if len(polylines) == 0 and len(lines) == 0 and len(lwpolylines) == 0:
        sys.stderr.write("No outline found.  Hinge outline should consist of lines or polylines.")
        sys.exit(0)

    holes = [[c.center[0], c.center[1]] for c in circles]
    print 'holes', holes
    hinge[name + "_holes"] = holes


    for polyline in polylines:
        points = [[p[0], p[1]] for p in polyline.points()]
        startpoint = points[0]
        endpoint = points[-1]
        if startpoint == endpoint: # closed polyline
            hinge[name + '_contour'] = points
    #TODO: Implement lines and multiple polylines

    for polyline in lwpolylines:
        points = [[p[0], p[1]] for p in polyline.points]
        startpoint = points[0]
        endpoint = points[-1]
        if startpoint == endpoint: # closed polyline
            hinge[name + '_contour'] = points


if __name__ == '__main__':
    main()


