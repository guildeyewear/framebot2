from dxfwrite import DXFEngine as dxf
import json
from optparse import OptionParser
import os
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
    infile = open(name + ".json")
    hinge = json.loads(infile.read())

    face = hinge["face_contour"]
    temple = hinge["temple_contour"]
    face_hole_diam = hinge["drill_dia"]
    temple_hole_diam = hinge["drill_dia"]
    face_holes = hinge["face_holes"]
    temple_holes = hinge["temple_holes"]

    face_drawing = dxf.drawing("./%s-face.dxf" % name)
    temple_drawing = dxf.drawing("./%s-temple.dxf" % name)
    face_drawing.add_layer("OUTLINE", color=1)
    face_drawing.add_layer("HOLES", color=2)
    temple_drawing.add_layer("OUTLINE", color=1)
    temple_drawing.add_layer("HOLES", color=2)

    face_outline = dxf.polyline(layer="OUTLINE", thickness = 0.1)
    face_outline.add_vertices(face)
    face_drawing.add(face_outline)
    for hole in face_holes:
        circ = dxf.circle(face_hole_diam/2.0, hole, layer="HOLES", thickness=0.1)
        face_drawing.add(circ)
    face_drawing.save()

    temple_outline = dxf.polyline(layer="OUTLINE", thickness = 0.1)
    temple_outline.add_vertices(temple)
    temple_drawing.add(temple_outline)
    for hole in temple_holes:
        circ = dxf.circle(temple_hole_diam/2.0, hole, layer="HOLES", thickness=0.1)
        temple_drawing.add(circ)
    temple_drawing.save()

if __name__ == '__main__':
    main()


