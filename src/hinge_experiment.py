from dxfwrite import DXFEngine as dxf
import json
from optparse import OptionParser
import os
import sys
import poly
import cam
import hinges
import math
import nose
import geometry as g

def log(message):
    sys.stdout.write(message + "\n")

def main():
    """Run as a command line program."""
    # Create the output director
    order_dir = "/Users/rodfrey/Dropbox/guild_orders/hinge_experiment"
    if not os.path.exists(order_dir):
        os.makedirs(order_dir)
    mill_hinges(order_dir)

    # Arrange the temples to fit on the stock

def mill_hinges(outdir):
    """Creates the g-code for the first milling operation.  The
    first milling operation is done on an unregistered plastic blank,
    so includes creating registration holes for subsequent operations."""
    y_offsets = [y*15 for y in range(-5, 6)]
    x_offsets = [x*7 for x in range(-3, 4)]

    hinge = hinges.get_hinge(2);
    contour = hinge['face_contour']

    contours = []
    erode = poly.erode(1.5875/2, contour)
    while len(erode) > 0:
        if len(erode) >= 1:
            contours.append(erode[0])
        else:
            break
        erode = poly.erode(1.5875/2, contours[-1])



    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.activate_pin("stock_clamp"),
        cam.change_tool("1/16in endmill"),
        cam.rapid([0,0]),
        cam.start_spindle(20000),
        cam.dwell(3),
#        cam.pocket(contours, -1, -1),
#        cam.rapid([None, None, 20.0]),
#        cam.change_tool("1mm drill"),
#        cam.start_spindle(4500),
#        [cam.rmp(p + [-2.5], retract=10.0) for p in hinge['face_holes']],
        ]



    for y_offset in y_offsets:
        for x_offset in x_offsets:
            program += [
                cam.rapid([x_offset, y_offset]),
                cam.temporary_offset((0,0)),
                cam.pocket(contours, -1.2, -1.2),
                cam.rapid([None, None, 20.0]),
                cam.remove_temporary_offset(),
                    ]
#    program += [
#            cam.change_tool("1mm drill"),
#            cam.start_spindle(4500),
#            ]
#    for y_offset in y_offsets:
#        for x_offset in x_offsets:
#            program += [
#                cam.rapid([x_offset, y_offset]),
#                cam.temporary_offset((0,0)),
#                [cam.rmp(p + [-2.5], retract=10.0) for p in hinge['face_holes']],
#                cam.rapid([None, None, 20.0]),
#                cam.remove_temporary_offset(),
#                    ]

    program += [
        cam.deactivate_pin("stock_clamp"),
        cam.end_program(),
    ]


    open(outdir + "/face_stage1.ngc", "w").write(to_string(program))

def flatten(l):
    '''Flattens a tree structure to a list of strings.  Nodes with a value of None are removed.'''
    if type(l) is type(""):
        return [l]
    elif type(l) is type(None):
        return []
    else:
        r = []
        for e in l:
            r += flatten(e)
        return r

def to_string(l):
    '''Converts a list of strings or a tree structure of strings to a single string for output.'''
    return "\n".join(flatten(l))


if __name__ == '__main__':
    main()
