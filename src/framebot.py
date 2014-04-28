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

    if (options.infile == None and options.url == None) or options.outdir == None:
        log("Insufficient arguments.")
        parser.print_help()
        sys.exit(0)

    infile = open(options.infile) if options.infile else urllib.urlopen(options.url)
    o = json.loads(infile.read())
    outdir = options.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create a dxf for cutting the outline on the laser cutter
    create_dxf('/Users/rodfrey/Dropbox/outer_contour.dxf', [o['face_con']])

    # Create the milling program for the lens holes/groove, hinge pockets, etc.
    mill_fronts(outdir, o)

    temples =  arrange_temple_curves(o['ltemple_con'], o['rtemple_con'], o['lhinge'], o['lhinge_y'], o['rhinge_y'])
    create_dxf("/Users/rodfrey/Dropbox/temple_contour.dxf", [temples['left_temple_contour'], temples['right_temple_contour']])
    mill_temples(outdir, temples, o['temple_length'])


def mill_fronts(outdir, order):
    """Creates the g-code for the first milling operation.  The
    first milling operation is done on an unregistered plastic blank,
    so includes creating registration holes for subsequent operations."""

    #TODO: Replace with information in materials database
    front_surface_thickness = 0
    back_surface_thickness = 4
    final_front_thickness = 0
    final_back_thickness = 4
    thickness = final_front_thickness + final_back_thickness

    front_surface_removal = final_front_thickness - front_surface_thickness
    back_surface_removal = final_back_thickness - back_surface_thickness


   # The machine has the stock clamp oriented 90 degrees to the way the
    # software creates the contours.
    face_c = poly.rotate_90(order["face_con"])
    left_lens_c = poly.rotate_90(order["lhole_con"])
    right_lens_c = poly.rotate_90(order["rhole_con"])

    msg = check_frame_size(face_c)
    if msg:
        print msg
        sys.exit(0)

    print 'milling front with hinge', order['lhinge'], order['rhinge']
    offset = frame_offset(face_c)
    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.activate_pin("stock_clamp"),
        surface_front(front_surface_removal),
        surface_back(back_surface_removal),

        cam.change_tool("1/16in endmill"),
        cam.rapid([0, 0]),
        cam.temporary_offset(offset),
       # Note that X and Y parameters in the order are switched from our system
        #TODO: replace the thickness offset with the thickness of the TEMPLE, not the fronts.
        face_hinge_pockets(order["lhinge"], order["lhinge_y"], order["ltemple_x"]),
        index_holes([face_c], thickness),
        lens_holes(left_lens_c, right_lens_c, thickness),
        nose_pads(order, thickness),
        nose_contour(order["nose_rad"], order["nose_h"], order["nose_sa"], order["nose_ra"], face_c, thickness),


        cam.contour(face_c, True),
        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        cam.end_program(),
    ]
    open(outdir + "/face_stage1.ngc", "w").write(to_string(program))

def mill_temples(outdir, temples, temple_length):
    #TODO: Replace with information in materials database
    front_surface_thickness = 0
    back_surface_thickness = 4
    final_front_thickness = 0
    final_back_thickness = 4
    thickness = final_front_thickness + final_back_thickness

    front_surface_removal = final_front_thickness - front_surface_thickness
    back_surface_removal = final_back_thickness - back_surface_thickness

    r_temple = poly.rotate_90(temples['right_temple_contour'])
    l_temple = poly.rotate_90(temples['left_temple_contour'])

    offset = frame_offset(l_temple)
    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.activate_pin("stock_clamp"),
        surface_front(front_surface_removal),
        surface_back(back_surface_removal),
        cam.change_tool("1/16in endmill"),
        cam.rapid([0,0]),
        cam.temporary_offset(offset),
        temple_hinge_pockets(temples),
        index_holes([l_temple, r_temple], thickness),
        #thin_temples([l_temple, r_temple], temple_length),
    ]
    open(outdir + "/temples_milling.ngc", "w").write(to_string(program))

def thin_temples(temples, temple_length):
    left = temples[0]
    right = temples[1]
    left_side = poly.bottom(left)
    right_side = poly.top(left)
    halfway = left_side + (right_side - left_side)/2

    taperpath = []
    for pt in left:
        if pt[1] == left_side:
            break
        elif pt[1] > halfway:
            taperpath.append(pt)

    print taperpath
    # offset the path so our 1/2 mill cuts the whole thing
    # TODO: gauge width and make sure we're cutting the whole thing.
    flat_begin = left_side + temple_length -10 # 20 mm flat
    flat_end = flat_begin + 20
    front_slope = -2.0 / (halfway-flat_begin)

    print "flat", flat_begin, flat_end
    print left_side, right_side, halfway
    def calc_thinning_z(pt):
        if pt[1] > flat_begin:
            print 'Over flat begin', pt
            return (abs(pt[1]-halfway) * front_slope)
        elif pt[1] > flat_end:
            print 'over flat end'
            return -4
        else:
            return -(pt[0]-flat_end)/4 - 4



    shiftedpath = []
    for idx, pt in enumerate(taperpath):
        if idx == 0:
            shiftedpath.append([pt[0], pt[1]-3])
        else:
            lastpt = taperpath[idx-1]
            line=[pt[0]-lastpt[0], pt[1]-lastpt[1]]
            normal=[-line[1], line[0]]
            length = math.sqrt(normal[0]*normal[0]+normal[1]*normal[1])
            normal = [6*(x/length) for x in normal]
            shiftedpath.append([pt[0]+normal[0], pt[1]+normal[1]])
    thinning_contour_left = [[pt[0], pt[1], calc_thinning_z(pt)] for pt in shiftedpath]

    #thinning_contour_right = poly.mirror_x(thinning_contour_left)
    return [
        cam.rmp(thinning_contour_left[0]),
        cam.contour(thinning_contour_left, False),
    ]


def surface_front(amount):
    log("Surfacing front with amount %f" % amount)
    return [
        cam.flip_stock(),
        cam.change_tool("1/4in endmill"),
        cam.spindle_speed(20000),
        cam.feedrate(2000),
        cam.start_spindle(),
        cam.surface_along_y(-80, -100, -5, 100, 3.175, amount),
        cam.stop_spindle(),
        cam.retract_spindle(),
        cam.flip_stock(),
        ] if amount < 0 else None

def surface_back(amount):
    return [
        cam.change_tool("1/4in endmill"),
        cam.spindle_speed("20000"),
        cam.start_spindle(),
        cam.surface_along_y(-80, -100, -5, 100, 0.25/2, amount),
        cam.stop_spindle(),
        cam.retract_spindle(),
        cam.dactivate_pin("stock_clamp"),
        ] if amount > 0 else None




def face_hinge_pockets(hinge_num, xposition, yposition):
    left_hinge = hinges.get_hinge(hinge_num)
    right_hinge = hinges.get_hinge(hinge_num, False)
    left_translate = [xposition, -yposition]
    #left_translate = [xposition, 0]
    right_translate = [xposition, yposition]
    #right_translate = [xposition, 0]
    # Adjust by pocket depth of hinge
    pocket_depth = left_hinge['pocket_depth']

    left_contour = poly.translate(left_hinge["face_contour"], left_translate[0], left_translate[1])
    right_contour = poly.translate(right_hinge["face_contour"], right_translate[0], right_translate[1])
    left_holes = poly.translate(left_hinge["face_holes"], left_translate[0], left_translate[1])
    right_holes = poly.translate(right_hinge["face_holes"], right_translate[0], right_translate[1])

    if not poly.is_ccw(left_contour):
        left_contour = poly.reverse(left_contour)
    if not poly.is_ccw(right_contour):
        right_contour = poly.reverse(right_contour)

    left_hinge_pocket_contours = [];
    while len(left_contour) > 0:
        left_contour = poly.erode(1.5875/2, left_contour)
        if len(left_contour) > 0:
            left_contour = left_contour[0]
            left_hinge_pocket_contours.append(left_contour)

    right_hinge_pocket_contours = [];
    while len(right_contour) > 0:
            right_contour = poly.erode(1.5875/2, right_contour)
            if len(right_contour) > 0:
                right_contour = right_contour[0]
                right_hinge_pocket_contours.append(right_contour)
    r = [
        cam.comment("Hinge Pockets"),
        cam.feedrate(750),
        cam.change_tool("1/16in endmill"),
        cam.start_spindle(15000),
        cam.dwell(3),
        cam.comment("Right Hinge Pocket"),
        cam.pocket(right_hinge_pocket_contours, -abs(right_hinge['pocket_depth']), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Left Hinge Pocket"),
        cam.pocket(left_hinge_pocket_contours, -abs(left_hinge['pocket_depth']), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Hinge Holes"),
        cam.change_tool("1mm drill"),
        cam.start_spindle(4500),
        cam.dwell(2),
        [cam.rmp(p + [-8.0], retract=10.0) for p in right_holes],
        [cam.rmp(p + [-8.0], retract=10.0) for p in left_holes],
        cam.rapid([None, None, 20.0]),
    ]
    return r

def temple_hinge_pockets(temples):
    # We're operating in a 90 degree rotated fixture
    #l_hinge = poly.rotate_90(temples["left_hinge_contour"])
    #r_hinge = poly.rotate_90(temples["right_hinge_contour"])

    l_hinge = temples["left_hinge_contour"]
    r_hinge = temples["right_hinge_contour"]
    if not poly.is_ccw(l_hinge):
        l_hinge = poly.reverse(l_hinge)
    if not poly.is_ccw(r_hinge):
        r_hinge = poly.reverse(r_hinge)

    left_hinge_pocket_contours = [];
    while len(l_hinge) > 0:
        l_hinge = poly.erode(1.5875/2, l_hinge)
        if len(l_hinge) > 0:
            l_hinge = l_hinge[0]
            left_hinge_pocket_contours.append(l_hinge)

    right_hinge_pocket_contours = [];
    while len(r_hinge) > 0:
            r_hinge = poly.erode(1.5875/2, r_hinge)
            if len(r_hinge) > 0:
                r_hinge = r_hinge[0]
                right_hinge_pocket_contours.append(r_hinge)
    r = [
        cam.comment("Hinge Pockets"),
        cam.feedrate(750),
        cam.change_tool("1/16in endmill"),
        cam.start_spindle(15000),
        cam.dwell(3),
        cam.comment("Right Hinge Pocket"),
        cam.pocket(right_hinge_pocket_contours, -abs(temples['pocket_depth']), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Left Hinge Pocket"),
        cam.pocket(left_hinge_pocket_contours, -abs(temples['pocket_depth']), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Hinge Holes"),
        cam.change_tool("1mm drill"),
        cam.start_spindle(4500),
        cam.dwell(2),
        [cam.rmp(p + [-8.0], retract=10.0) for p in temples['right_hinge_holes']],
        [cam.rmp(p + [-8.0], retract=10.0) for p in temples['left_hinge_holes']],
        cam.rapid([None, None, 20.0]),

        cam.move([None, None, 0]),
        cam.contour(poly.rotate_90(temples['left_temple_contour']), True),
        cam.contour(poly.rotate_90(temples['right_temple_contour']), True),
    ]
    return r





def lens_holes(left_c, right_c, thickness):
    """Generates the toolpath for the lens holes (holes, groove and tabs)."""
    if not poly.is_ccw(left_c):
        left_c = poly.reverse(left_c)
    if not poly.is_ccw(right_c):
        right_c = poly.reverse(right_c)

    lhole = poly.erode(3.175/2.0, left_c)[0]
    rhole = poly.erode(3.175/2.0, right_c)[0]

    right_rough = poly.erode(0.1, rhole)[0]
    left_rough = poly.erode(0.1, lhole)[0]

    lgroove = poly.erode(0.8, left_c)[0]
    rgroove = poly.erode(0.8, right_c)[0]

    left_entry = poly.erode(2.0, lhole)[0][0];
    right_entry = poly.erode(2.0, rhole)[0][0];

    lhole = poly.reverse(lhole)
    rhole = poly.reverse(rhole)

    print "Groove will be at ", (-thickness/2)

    r = [
        "(Lens Holes)",
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(20000),
        cam.feedrate(2000),
        cam.rmh(right_entry + [-thickness - 1.0], 1.5, 0.5, 1.0),
        cam.contour(right_rough, True),
        cam.contour(rhole, True),
        cam.rmh(left_entry + [-thickness - 1.0], 1.5, 0.5, 1.0),
        cam.contour(left_rough, True),
        cam.contour(lhole, True),

        "(Lens Grooves)",
        cam.change_tool("vgroove"),
        cam.start_spindle(20000),
        cam.feedrate(2000),
        cam.rmp(right_entry + [-thickness/2]),
        cam.contour(rgroove, True),
        cam.move(right_entry), # Get out from under the overhang
        cam.rmp(left_entry + [-thickness/2]),
        cam.contour(lgroove, True),
        cam.move(left_entry), # Get out from under the overhang
    ]
    return r

def index_holes(contours, thickness):
# We put the index holes 1/2 between top and bottom, 160mm apart
#    log(str(poly.right(face_c)))
#    log(str(poly.left(face_c)))
#    log(str(poly.right(face_c) - poly.left(face_c)))
#    log(str((poly.right(face_c) - poly.left(face_c))/2))
#    log(str(poly.right(face_c) - (poly.right(face_c) - poly.left(face_c))/2))

    rightmost = -1000000
    leftmost = 1000000
    for contour in contours:
        right = poly.right(contour)
        left = poly.left(contour)
        rightmost = max(right, rightmost)
        leftmost = min(left, leftmost)

#    x_offset = poly.right(face_c) - (poly.right(face_c) - poly.left(face_c))/2
    x_offset = rightmost - (rightmost - leftmost)/2

    hole_radius = 4.85/2 # Measured from dowel pin
    tool_radius = 3.175/2
    helix_radius = hole_radius - tool_radius
    r_hole = [x_offset, 90]
    l_hole = [x_offset, -90]
    r = [
        cam.comment("Index Holes for secondary operations"),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(15000),
        cam.feedrate(1000),
        cam.dwell(3),
        cam.rmh(r_hole + [-thickness - 1.0], helix_radius, 0.5, 1),
        cam.rmh(l_hole + [-thickness - 1.0], helix_radius, 0.5, 1),
    ]
    return r


def nose_pads(order, thickness):
    return []

def nose_contour(nose_rad, nose_h, nose_sa, nose_ra, face_con, thickness):
    """Creates the nose contour feature toolpath.  Angular arguments are in degrees."""

    nr = nose_rad
    h = nose_h
    sa = math.radians(nose_sa)
    ra = math.radians(nose_ra)
    xfloor = poly.left(face_con) - 3.175/2.0  # bottom most point minus tool radius
    xfloor = max(xfloor, -27.0) # miminum safe distance without hitting clamp
    nose_tool_radius = 3.175/2.0

    nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, nose_tool_radius, 0.0)

    r = [
        "(Nose Contour)",
        cam.change_tool("1/8in ballmill"),
        cam.start_spindle(20000),
        cam.feedrate(2000),
        cam.rmp(nextpoly[0] + [2.0])  # Start near our first contour
    ]

    direction = 1
    for i in range(-20, (thickness+2)*10):
        z = -i/10.0
#        r += cam.move(nextpoly[0])
        if(direction < 0):
            nextpoly.reverse()
        direction = direction * -1
        r += cam.contour(nextpoly, False)
        r += cam.move([None, None, z])
        nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, nose_tool_radius, z)

    return r

def arrange_temple_curves(left_temple_contour, right_temple_contour, hinge, lhinge_y, rhinge_y):
    left_hinge = hinges.get_hinge(hinge)
    right_hinge = hinges.get_hinge(hinge, False)

    print 'first and last point of left_temple_contour'
    print left_temple_contour[0], left_temple_contour[-1]
    left_holes = poly.rotate_90(left_hinge['temple_holes'])
    right_holes = poly.rotate(right_hinge['temple_holes'], 90) # opposite direction.  FIX someday
    left_hinge_contour = poly.rotate_90(left_hinge['temple_contour'])
    right_hinge_contour = poly.rotate(right_hinge['temple_contour'], 90)

    #right_hinge_contour = right_hinge['temple_contour']
    #left_holes = left_hinge['temple_holes']
    #right_holes = right_hinge['temple_holes']

    # Get the thing as horizontal as possible
    # When top-bottom distance is minimized, it's horizontal
    height = poly.top(left_temple_contour) - poly.bottom(left_temple_contour)
    opt_angle = 0

    for angle in range(2, 41):
        candidate_contour = poly.rotate(left_temple_contour, -angle)
        candidate_height = poly.top(candidate_contour)-poly.bottom(candidate_contour)
        if candidate_height < height:
            height = candidate_height
            opt_angle = angle
        else:
            break

    # The temple is raised or lowered as compared to the Y axis.  We need
    # to bring it to the Y origin for rotation to avoid offsetting it, then
    # send it back to its original spot.
    original_y_offset = left_temple_contour[0][1] - (left_temple_contour[0][1] - left_temple_contour[-1][1])/2

    left_temple_contour = poly.translate(left_temple_contour, 0, -original_y_offset)
    left_temple_contour = poly.rotate(left_temple_contour, -opt_angle);
    left_temple_contour = poly.translate(left_temple_contour, 0, original_y_offset)

    right_temple_contour = poly.translate(right_temple_contour, 0, -original_y_offset)
    right_temple_contour = poly.rotate(right_temple_contour, opt_angle);
    right_temple_contour = poly.translate(right_temple_contour, 0, original_y_offset)

    #left_holes = poly.rotate(left_hinge['temple_holes'], -opt_angle)
    #right_holes = poly.rotate(right_hinge['temple_holes'], opt_angle)
    #left_hinge_contour = poly.rotate(left_hinge['temple_contour'], -opt_angle)
    #right_hinge_contour = poly.rotate(right_hinge['temple_contour'], opt_angle)

    left_holes = poly.rotate(left_holes, -opt_angle)
    right_holes = poly.rotate(right_holes, opt_angle)
    left_hinge_contour = poly.rotate(left_hinge_contour, -opt_angle)
    right_hinge_contour = poly.rotate(right_hinge_contour, opt_angle)
    left_hinge_contour = poly.translate(left_hinge_contour,  lhinge_y, -left_hinge['pocket_depth'])
    right_hinge_contour = poly.translate(right_hinge_contour, rhinge_y, right_hinge['pocket_depth'])
    # Left and right translate are reversed because we've rotated everything
    left_holes = poly.translate(left_holes, lhinge_y, -left_hinge['pocket_depth'])
    right_holes = poly.translate(right_holes, rhinge_y, right_hinge['pocket_depth'])

    temple_width = poly.right(left_temple_contour) - poly.left(left_temple_contour)
    print "temple width", temple_width
    lt_trans = [temple_width/2, 0]
    rt_trans = [-lt_trans[0], 0]

    #lt_trans = [180-poly.right(left_temple_contour), -poly.top(left_temple_contour) - 30]
    #lt_trans = [-poly.right(left_temple_contour)-30, 180 -poly.top(left_temple_contour)]
    #rt_trans = [-poly.left(right_temple_contour)+25, lt_trans[1]-height]

   # Translate them
    left_temple_contour = poly.translate(left_temple_contour, lt_trans[0], lt_trans[1])
    right_temple_contour = poly.translate(right_temple_contour, rt_trans[0], rt_trans[1])
    left_hinge_contour = poly.translate(left_hinge_contour, lt_trans[1], -lt_trans[0])
    right_hinge_contour = poly.translate(right_hinge_contour, rt_trans[1], -rt_trans[0])
    left_holes = poly.translate(left_holes, lt_trans[1], -lt_trans[0])
    right_holes = poly.translate(right_holes, rt_trans[1], -rt_trans[0])


    # Translate bottom one upward as much as we can
    trans_amount = poly.bottom(left_temple_contour) - poly.top(right_temple_contour) - 8
    for i in range(1, 100):
        candidate = poly.translate(right_temple_contour, 0, trans_amount + i)
        intersection = poly.intersection(candidate, left_temple_contour)
        if len(intersection) > 0:
            trans_amount = trans_amount + i -8
            right_temple_contour = poly.translate(right_temple_contour, 0, trans_amount)
            break
    print 'trans amount:', trans_amount

    right_hinge_contour = poly.translate(right_hinge_contour, trans_amount, 0)
    right_holes = poly.translate(right_holes, trans_amount, 0)

    return {
        "pocket_depth": left_hinge['pocket_depth'],
        "left_hinge_contour": left_hinge_contour,
        "right_hinge_contour": right_hinge_contour,
        "left_hinge_holes": left_holes,
        "right_hinge_holes": right_holes,
        "left_temple_contour": left_temple_contour,
        "right_temple_contour": right_temple_contour
            }



def create_dxf(filename, polys):
    drawing = dxf.drawing(filename)
    drawing.add_layer('OUTLINE', color=1)
    for p in polys:
        polyline = dxf.polyline(layer="OUTLINE")
        p = p + [p[0], p[1]]  # Close the polygon to avoid a cusp
        polyline.add_vertices(p)
        drawing.add(polyline)

    drawing.save()

def check_frame_size(contour):
    # We're limited by our stock size and clamp clearances
    # We can be about 160mm wide max, and about 65mm high max.
    if abs(poly.top(contour)) > 85:
        return "Frame is too wide: %f mm" % (poly.top(contour) * 2)
    if poly.right(contour) - poly.left(contour) > 70:
        return "Frame is too tall: %f mm" % (poly.right(contour) - poly.left(contour))
    return None


def frame_offset(contour):
    """
    The origin of the glasses is in the middle of the pair (i.e. centered at Y) with the
    X origin at the line between the pupils. The fixture on the milling machine has its origin
    centered on the Y axis but with the X axis at the edge of the clamp.  We need to shift the
    machine origin toward the middle of the stock.
    """
    xoffset = poly.right(contour) + 10 # Offset by the furthest point, plus some extra for the tool
    return [xoffset, 0]

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
