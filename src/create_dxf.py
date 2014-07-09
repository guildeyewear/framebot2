from dxfwrite import DXFEngine as dxf
import json
from optparse import OptionParser
import os
import sys
import poly

def log(message):
    sys.stdout.write(message + "\n")

def main():
    """Run as a command line program."""

    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="infile",
                      help="read specification from INFILE", metavar="INFILE")
    parser.add_option("-u", "--url", dest="url",
                      help="read specification from URL", metavar="URL")
    parser.add_option("-o", "--outdir", dest="outdir",
                      help="write manufacturing instruction files to OUTDIR", metavar="OUTDIR")

    (options, args) = parser.parse_args()

    if (options.infile == None and options.url == None):
        log("Insufficient arguments.")
        parser.print_help()
        sys.exit(0)

    infile = open(options.infile) if options.infile else urllib.urlopen(options.url)

    o = json.loads(infile.read())

    drawing = dxf.drawing('test.dxf')
    drawing.add_layer('OUTLINE', color=1)

    lhole = o['lhole_con']
    polyline = dxf.polyline(layer="OUTLINE")
    polyline.add_vertices(lhole)
    drawing.add(polyline)

    p2 = poly.erode(3.175/2.0, lhole)[0]
    p2 = poly.erode(2.0, p2);
    print p2[0][0]
    poly2 = dxf.polyline(layer="OUTLINE")
    poly2.add_vertices(p2[0])
    drawing.add(poly2)
    drawing.save()



def create_dxf(p):
    drawing = dxf.drawing('test.dxf')
    drawing.add_layer('OUTLINE', color=1)
    polyline = dxf.polyline(layer="OUTLINE")
    polyline.add_vertices(p)
    drawing.add(polyline)
    drawing.save()

if __name__ == '__main__':
    main()
