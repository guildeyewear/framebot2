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
import datetime

def log(message):
    sys.stdout.write(message + "\n")

def main():
    """Run as a command line program."""
    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="infile",
                      help="read specification from INFILE", metavar="INFILE")
    parser.add_option("-o", "--out", dest="outname",
            help="Name of output director", metavar="OUTDIR")
    parser.add_option("-u", "--url", dest="url",
                      help="read specification from URL", metavar="URL")
    (options, args) = parser.parse_args()
    if options.infile == None and options.url == None:
        log("Insufficient arguments.")
        parser.print_help()
        sys.exit(0)
    infile = open(options.infile) if options.infile else urllib.urlopen(options.url)
    o = json.loads(infile.read())

    order_id = o.get("order_id") or "testorder"
    outname = datetime.date.today().isoformat() + "/" + order_id
#    if options.outname == None:
#        print "extracting out name from ", options.infile
#        outname = options.infile.split("/")[-1]
#        outname = outname.split(".")[0]
#        outname = datetime.date.today().isoformat() + "/" + outname
#        print "out file is ", outname
#    else:
#        outname = options.outname
    # Create the output director
    #order_dir = "/Users/rodfrey/Dropbox/guild_orders/" + outname
    #order_dir = "/Volumes/Untitled/" + outname
    order_dir = "/Users/rodfrey/Development/framebot2/" + outname
    order_dir = options.outname
    if not os.path.exists(order_dir):
        os.makedirs(order_dir)
        print 'Created dir', order_dir
    else:
        for i in range(1, 101):
            candidate = order_dir + "_" + str(i)
            print 'checking directory', candidate
            if i == 100:
                print 'Too many directories!  Could not create order directory'
                sys.exit(1)
            if not os.path.exists(candidate):
                order_dir = candidate
                os.makedirs(order_dir)
                print "Created dir", order_dir
                break



    # Create the milling program for the lens holes/groove, hinge pockets, etc.
    # and the dxf for the laser for the fronts
    mill_fronts(order_dir, o)
    laser_face = o['face_con']
    laser_face = laser_face + poly.reverse(poly.mirror_x(laser_face, False))[1:]

    # Rearrange the face curve a bit so the laser starts outside the curve.
    # Otherwise it leaves a little scar where it starts.
    laser_face = poly.new_start(laser_face, poly.leftmost_index(laser_face))
    laser_face.append(laser_face[0]) # Close the curve

    first_vector = [laser_face[1][0]-laser_face[0][0], laser_face[1][1]-laser_face[0][1]]
    unit = math.sqrt(first_vector[0]**2 + first_vector[1]**2)
    leadin_vector = [(4*first_vector[0])/unit,(4*first_vector[1])/unit]
    leadin_start = [laser_face[0][0] - leadin_vector[0], laser_face[0][1] - leadin_vector[1]]

    last_vector = [laser_face[-2][0]-laser_face[-1][0], laser_face[-2][1]-laser_face[-1][1]]
    unit = math.sqrt(last_vector[0]**2 + last_vector[1]**2)
    leadout_vector = [(4*last_vector[0])/unit,(4*last_vector[1])/unit]
    leadout_end = [laser_face[-1][0] + leadout_vector[0], laser_face[-1][1] + leadout_vector[1]]

    laser_face.insert(0, leadin_start)
    create_dxf(order_dir + "/face_contour.dxf", [laser_face], close=False)

#    temple_dxf(o, order_dir + "/left_temple.dxf", 1, False)
#    temple_dxf(o, order_dir + "/right_temple.dxf", 1, True)


    # Arrange the temples to fit on the stock
    temples =  arrange_temple_curves(o['ltemple_con'], o.get('lhinge') or 1)
    left_hinge = poly.rotate_90(temples['left_hinge_contour'])
    right_hinge = poly.rotate_90(temples['right_hinge_contour'])
#
    l_temple = poly.rotate_90(temples['left_temple_contour']);
    r_temple = poly.rotate_90(temples['right_temple_contour']);
#
    create_dxf(order_dir + "/temple_contour.dxf",
            [poly.rotate_90(temples['left_temple_contour']),
             poly.rotate_90(temples['right_temple_contour']),
#             left_hinge,
#             right_hinge,
            ],
            close_with_arc=False,
            close=False)

    mill_temples(order_dir, temples, o)

#    mill_lenses(order_dir, o)


def mill_lenses(outdir, order):
    def to_polar(polyline):
        return [[-math.sqrt(p[0]**2 + p[1]**2), math.degrees(math.atan2(p[1],p[0])) + 180] for p in polyline]

    """ Creates g-code for milling the left and right lenses for the frames."""
    lens = g.Polygon(order['lhole_con'])
    x = lens.bounding_box()
    box = lens.bounding_box()
    center = g.Point((box.p2.x+box.p1.x)/2, (box.p2.y+box.p1.y)/2)
    shift_to_origin = g.Vector(-center.x, -center.y)
    lens = g.translate_polygon(lens, shift_to_origin)
    polar = g.polygon_to_uniform_polar(lens, 1500)

    # Inflate the polar form by 1/2 the diameter of the cutter, and convert
    # the angles to degrees
    tau = math.pi * 2
    conv = 360/tau
    roughing = [(pt.theta*conv, pt.r+4.77) for pt in polar]

    # Expand for the roughing cut
#    left_lens_rough = poly.dilate(4.77, left_lens)
    # Convert to polar coordinates
#    left_lens_roughing_polar = to_polar(left_lens_rough)
    # Start the cut from 0 degrees
#    angle_list = [p[1] for p in left_lens_roughing_polar]
#    zero_idx = angle_list.index(min(angle_list))
#    left_lens_roughing_polar = left_lens_roughing_polar[zero_idx:] + left_lens_roughing_polar[:zero_idx]
    # Cut it down to every 0.2 degrees
#    coarse = []
#    for idx, pt in enumerate(left_lens_roughing_polar):
#        if idx == 0 or (pt[1]-coarse[-1][1]) >= 0.2:
#            coarse.append(pt)

    # The polar distances aren't correct quite.  That's because the conversion assumed a flat
    # surface, but we'll actually be bending that contour around a sphere.  The lens is already
    # spherical.  So we need to adjust inwards a bit to account for the x distance actually being an arc
    # distance.

    # The radius of the sphere is 88mm (base 6 lens).  ArcLength = Radius * angle, and chord
    # length is sin(angle) * radius.
    roughing = [(p[0], math.sin(-p[1]/88) * 88) for p in roughing]
    roughing.sort(key=lambda pt: pt[0])
    closest = max([p[1] for p in roughing])
    if abs(closest) < 22.5:
        print "Error!  Cannot cut lens, it's too small.", closest

    roughing_reversed = [[-1*(math.degrees(-math.radians(p[1]))+360), p[1]] for p in roughing]

    program = [
        cam.setup(),
        cam.select_fixture("lens_clamp"),
        cam.retract_spindle(),
        cam.change_tool("1/4in endmill"),
        cam.start_spindle(20000),
        cam.rapid([-50, 0]),
        ["G1 Z0 F500"],
        cam.polar_contour(roughing),
        ["G1 X-50"],
        ["G1 A0"],
        cam.stop_spindle(),
        cam.retract_spindle(),
        cam.end_program(),
    ]
    open(outdir + "/left_lens.ngc", "w").write(to_string(program))
    program = [
        cam.setup(),
        cam.select_fixture("lens_clamp"),
        cam.retract_spindle(),
        cam.change_tool("1/4in endmill"),
        cam.start_spindle(20000),
        cam.rapid([-50, 0]),
        ["G1 Z0 F500"],
        cam.polar_contour(roughing_reversed),
        ["G1 X-50"],
        ["G1 A0"],
        cam.stop_spindle(),
        cam.retract_spindle(),
        cam.end_program(),
    ]
    open(outdir + "/right_lens.ngc", "w").write(to_string(program))





def mill_fronts(outdir, order):
    """Creates the g-code for the first milling operation.  The
    first milling operation is done on an unregistered plastic blank,
    so includes creating registration holes for subsequent operations."""

# Initial thickness of the forward-facing lamination.  0 if no lamination.
    print 'Creating milling program for frame fronts'
    top_thickness = 6
    bottom_thickness = 0
    top_raw = 6
    bottom_raw = 0
    if order.get("face_material"):
        top_thickness = order["face_material"].get("top_thickness")
        bottom_thickness = order["face_material"].get("bottom_thickness") or 0
        top_raw = order["face_material"].get("top_raw_thickness") or top_thickness
        bottom_raw = order["face_material"].get("bottom_raw_thickness") or bottom_thickness
    frame_thickness = top_thickness + bottom_thickness
    machining_z_offset = bottom_raw - bottom_thickness

    print 'top raw', top_raw
    print 'bottom raw', bottom_raw
#    top_thickness = 6
#    bottom_thickness = 0
#    total_thickness = 6
#    top_raw = 0


    # Automatically determine some surfacing.  The bottom thickness
    # is a lamination and should be thin.  It should be thinned to 1mm.
    # The top should be thinned to bring total finished thickness to 6mm or less.
    thin_back = 0
    if bottom_raw > bottom_thickness:
        thin_back = bottom_raw - bottom_thickness
    thin_front = 0
    if top_raw > top_thickness:
        thin_front = top_raw - top_thickness

    print 'Calculated surfacing: thinning back by', thin_back, ', thinning front by', thin_front, ', total frame thickess is', frame_thickness

# The machine has the stock clamp oriented 90 degrees to the way the
# software creates the contours.
    face_c = poly.rotate_90(order["face_con"])
    left_lens_c = poly.rotate_90(order["lhole_con"])
    right_lens_c = poly.mirror_y(left_lens_c, True)
    face_c = face_c + poly.reverse(poly.mirror_y(face_c, False))[1:]

    print 'Got order contours and rotated them for the router'
    temple_height = abs(order["ltemple_con"][0][1] - order["ltemple_con"][-1][1])
    print 'Calculated temple height:', temple_height

    msg = check_frame_size(left_lens_c)
    print 'Checked frame size'
    if msg:
        print msg
        sys.exit(1)
    print 'Frame size OK'
#   Instead of offset, we'll make sure this thing is centered on the stock
    offset = frame_offset(face_c)
    top = poly.top(face_c)
    bottom = poly.bottom(face_c)
    left = poly.left(face_c)
    right = poly.right(face_c)
    y_shift = (top + bottom)/2
    x_shift = (left + right)/2

    face_c = poly.translate(face_c, -x_shift, -y_shift)
    left_lens_c = poly.translate(left_lens_c, -x_shift, -y_shift)
    right_lens_c = poly.translate(right_lens_c, -x_shift, -y_shift)
    print 'Moved all contours into place'

# Groove for lenses is 2/3 of the distance from back to front.
# Here we're calculating the actual cutting depth so we need to add
# back the material that we surfaced away from the back.
    groove_depth = -(float(machining_z_offset) + (2.0/3)*float(frame_thickness))
    print 'groove depth is', groove_depth
    print 'ltemple endpoints', order["ltemple_con"][0][1], order["ltemple_con"][-1][1]
    hinge_loc = order["ltemple_con"][0][1] - (order["ltemple_con"][0][1] - order["ltemple_con"][-1][1])/2

    size_info = order.get('size_info') or order.get('usersizes')

    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.activate_pin("stock_clamp"),
        cam.rapid([0,0]),
        surface_front(thin_front),
        surface_back(thin_back),

        cam.change_tool("1/8in endmill"),
        cam.rapid([0, 0]),

        # Note that X and Y parameters in the order are switched from our system
        #TODO: replace the thickness offset with the thickness of the TEMPLE, not the fronts.
        index_holes(frame_thickness+machining_z_offset),
        lens_holes(left_lens_c, right_lens_c, frame_thickness + machining_z_offset),
        lens_groove(left_lens_c, right_lens_c, groove_depth),
        face_hinge_pockets(order.get("lhinge") or 1, hinge_loc, order["ltemple_x"], (-x_shift, -y_shift), machining_z_offset),

        #nose_pads(order, thickness),
        nose_contour(
            float(size_info["noseradius"]),
            float(size_info["noseheight"]),
            float(size_info["splayangle"]),
            float(size_info["ridgeangle"]),
            face_c, frame_thickness, machining_z_offset, -x_shift ),
        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        cam.change_tool("1/8in endmill"),
        #cam.rapid([face_c[0][0], face_c[0][1], -thin_back] ),
        #cam.contour(face_c, True),
        cam.end_program(),
    ]
    print 'Writing face milling program to ', outdir + "/face_stange1.ngc"
    open(outdir + "/face_stage1.ngc", "w").write(to_string(program))

def nose_contour(nose_rad, nose_h, nose_sa, nose_ra, face_con, thickness, thin_back, centering_shift):
    """Creates the nose contour feature toolpath.  Angular arguments are in degrees."""#
    print 'Generating nose contour'
    nr = nose_rad
    nose_tool_radius = 3.175

    # Set the nose height to the bottom of the bridge, but no
    # higher than 10mm.
    #h = min(face_con[len(face_con)/2][0], 10)

    # Now the nose height will be set to the measured height, but
    # no less than 10, or the bridge if that's less than 10
    #h = max(h, nose_h)

    #bridge_top = face_con[0][0]


    # We're cutting with a ball-end mill.  Where it actually cuts is dependent on the
    # ridge angle.  If you draw a line at the ridge angle and put it tangent to the ball mill,
    # that is the cutting line.  The angle between the center of the ball mill and the intersection
    # of the cutting line and the surface is 1/2 of the ridge angle.  From that and the radius
    # of the ball mill we can figure out the offset.
    cutter_offset = (nose_tool_radius)*math.tan(math.radians(nose_ra/2))
    print "Cutter offset", cutter_offset

    sa = math.radians(nose_sa)
    ra = math.radians(nose_ra)
    h = nose_h + centering_shift

    xfloor = poly.left(face_con) - 3.175  # bottom most point minus tool radius
    xfloor = max(xfloor, -27.0) # miminum safe distance without hitting clamp

    z_depth = 0 # Start a bit above the surface of the glasses
    nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, cutter_offset, z_depth)

    r = [
        "(Nose Contour)",
        cam.change_tool("1/4in ballmill"),
        cam.start_spindle(20000),
        cam.feedrate(4000),
        cam.rmp(nextpoly[0] + [2.0]),  # Start near our first contour
    ]

    direction = 1

    z_start = int((z_depth)*10) # We have to use integers for the range, also step in 1/10 mm steps

    for i in range(-z_start, int(((thickness)+3)*10)):
        z = -i/10.0
#        r += cam.move(nextpoly[0])
        if(direction < 0):
            nextpoly.reverse()
        direction = direction * -1
        r += cam.move([None, None, z-thin_back]) # Z adjusted for any surfacing that happened
        r += cam.contour(nextpoly, False)
        nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, cutter_offset, z)
    return r

def mill_temples(outdir, temples, order):
#TODO: Replace with information in materials database
    print 'Creating milling program for temples'
    top_thickness = 4 # Assume 4mm temple
    top_raw = 4
    if order.get("temple_material"):
        top_raw = order["temple_material"].get("top_raw_thickness") or order["temple_material"].get("top_thickness") or top_thickness

    print 'top raw', top_raw

    # Automatically determine some surfacing.  The bottom thickness
    # is a lamination and should be thin.  It should be thinned to 1mm.
    # The top should be thinned to bring total finished thickness to 6mm or less.
    thin_back = 0
    if top_raw > top_thickness:
        thin_back = top_raw - top_thickness

    l_temple = temples['left_temple_contour']
    r_temple = temples['right_temple_contour']
    offset = frame_offset(l_temple)

# Calculate entry points for bevelling operation
    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.rapid([0,0]),
        cam.activate_pin("stock_clamp"),
        surface_back(thin_back),
        index_holes(top_raw),
        #cam.rapid(l_temple[0]),
        #cam.contour(l_temple, True),
        #cam.rapid(r_temple[0]),
        #cam.contour(r_temple, True),
        rough_temple_bevel(l_temple, thin_back),
        rough_temple_bevel(r_temple, thin_back),
        cam.change_tool("1/16in endmill"),
        cam.rapid([0,0]),
        temple_hinge_pockets(temples, thin_back),
        cam.change_tool("dovetail"),
        bevel_temple(l_temple, thin_back),
        bevel_temple(r_temple, thin_back),
        temple_hinge_clearance(l_temple, thin_back),
        temple_hinge_clearance(r_temple, thin_back),
        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        cam.change_tool("1/8in endmill"),
        cam.end_program()
    ]
    open(outdir + "/temples_milling.ngc", "w").write(to_string(program))

def extendLine(p1, p2, distance):
    lineLen = math.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
    return [
        p1[0] + distance * ((p1[0]-p2[0])/lineLen),
        p1[1] + distance * ((p1[1]-p2[1])/lineLen)
        ]

def bevel_temple(temple, thinning):
    # Assume a 20 degree dovetail cutter, cutting 5mm from bottom
    dovetail_offset = 9.52/2 - 1.82

    entry = extendLine(temple[-1], temple[-2], 7.5)
    # Endpoints for the undercut
    p1 = extendLine(temple[-1], temple[-2], dovetail_offset)
    p2 = extendLine(temple[0], temple[1], dovetail_offset)
    p2 = extendLine(p2, p1, 2)

    return [
        cam.change_tool("dovetail"),
        cam.feedrate(750),
        cam.start_spindle(20000),
        cam.rmp(entry + [-5-thinning]),
        cam.move(p1),
        cam.move(p2),
        cam.move(entry),
        cam.move(entry + [10]),
            ]

def temple_hinge_clearance(temple, thinning):
    # Endpoints for the top bevel
    entry = extendLine(temple[-1], temple[-2], 7.5)
    p1 = extendLine(temple[0], temple[-1], 1)
    p2 = extendLine(temple[-1], temple[0], 1)
    return [
       cam.change_tool("engraver"),
       cam.feedrate(750),
       cam.start_spindle(20000),
       cam.rmp(entry + [-2 - thinning]),
       cam.move(p1),
       cam.move(p2),
       cam.move(entry),
       cam.move(entry + [10]),
            ]




def rough_temple_bevel(temple,  thinning):
    p1 = extendLine(temple[-1], temple[-2], 3.175/2 - 0.5)
    p2 = extendLine(temple[0], temple[1], 3.175/2 - 0.5)
    p3 = extendLine(temple[0], temple[1], 15)
    p4 = extendLine(temple[-1], temple[-2], 15) # room for dovetail
# p1 and p2 are just extensions of the temple - move them to the side a bit to
# clearance for when the dovetail cutter comes through
    p1 = extendLine(p1, p2, 3)
    p2 = extendLine(p2, p1, 3)
    p3 = extendLine(p3, p4, 3)
    p4 = extendLine(p4, p3, 3)


# Move to the dovetail cutter entry point, helix through stock.
# Cut a circle big enough to admit the dovetail cutter.
# Rough cut the end of the temple
# Clear a return path for the dovetail cutter.
    return [
        cam.change_tool("1/8in endmill"),
        cam.feedrate(1000),
        cam.rmh(p4 + [-5-thinning], 1, pitch=1),
        cam.move(p1),
        cam.move(p2),
        cam.move(p3),
        cam.move(p4),
        cam.rapid(p4 + [10])
        ]


def surface_front(amount):
    if amount < 0.1:
        print 'Not surfacing front, returning'
        return None
    print "Surfacing front with amount", amount
    surface_amount = min(amount, 4)
    surface_heights = []
    while surface_amount <= amount:
        surface_heights.append(surface_amount)
        surface_amount = max(surface_amount+2, amount)
    print "surface amounts", surface_heights

    program =  [
        cam.comment("Surface front by %f mm" % amount),
        cam.flip_stock(),
        cam.change_tool("3/4in surfacer"),
        cam.spindle_speed(15000),
        cam.feedrate(1500),
        cam.start_spindle(),]
    for height in surface_heights:
        program = program + cam.surface_along_y(-40, -110, 40, 110, 9.525, -height)
    program = program + [
        cam.stop_spindle(),
        cam.retract_spindle(),
        cam.flip_stock(),
        ]
    return program

def surface_back(amount):
    if amount < 0.1:
        return None

    print 'surfacing back with', amount
    surface_amount = min(amount, 4)
    surface_heights = []
    while surface_amount <= amount:
        surface_heights.append(surface_amount);
        print surface_heights
        surface_amount = max(surface_amount+2, amount)
    print "surface amounts on back", surface_heights


    return [
        cam.comment("Surface back by %f mm" % amount),
        cam.change_tool("3/4in surfacer"),
        cam.spindle_speed(15000),
        cam.feedrate(1500),
        cam.start_spindle(),
        cam.dwell(5),
        cam.surface_along_y(-40, -110, 40, 110, 9.525, -amount),
        cam.rapid([None, None, 20]),
        ] if amount > 0 else None
    program =  [
        cam.comment("Surface back by %f mm" % amount),
        cam.change_tool("3/4in surfacer"),
        cam.spindle_speed(15000),
        cam.feedrate(1500),
        cam.start_spindle(),
        cam.dwell(5),
        ]
    for height in surface_heights:
        program = program + cam.surface_along_y(-40, -110, 40, 110, 9.525, -height),
    program = program + [
        cam.stop_spindle(),
        cam.retract_spindle(),
        ]
    return program


def face_hinge_pockets(hinge_num, hinge_height, temple_position, centering_shift, thin_back):

    print 'Generating face hinge pockets', hinge_num
    xposition = hinge_height;
    yposition = temple_position+3.5; # 4mm for the temple, but less 0.5mm for temple hinge pocket
    print 'Position is ', xposition, yposition
    left_hinge = hinges.get_hinge(hinge_num)
    print 'Got left hinge'
    right_hinge = hinges.get_hinge(hinge_num, False)
    print 'Retrieved hinge contours'
    left_translate = [xposition, yposition]
    right_translate = [xposition, -yposition]
    print 'Calculated hinge translations'
    #right_translate = [xposition, 0]
    # Adjust by pocket depth of hinge
    #pocket_depth = left_hinge['pocket_depth']+thin_back

    pocket_depth = 1 + thin_back
    drill_depth = -thin_back - 2.0

    left_contour = poly.mirror_x(poly.rotate_90(left_hinge["face_contour"]), False)
    right_contour = poly.mirror_x(poly.rotate_90(right_hinge["face_contour"]), False)
    left_holes = poly.mirror_x(poly.rotate_90(left_hinge["face_holes"]), False)
    right_holes = poly.mirror_x(poly.rotate_90(right_hinge["face_holes"]), False)
    left_contour = poly.translate(left_contour, left_translate[0], -left_translate[1])
    right_contour = poly.translate(right_contour, right_translate[0], -right_translate[1])
    left_holes = poly.translate(left_holes, left_translate[0], -left_translate[1])
    right_holes = poly.translate(right_holes, right_translate[0], -right_translate[1])

    # Now center everything on the stock
    left_contour = poly.translate(left_contour, centering_shift[0], -centering_shift[1])
    right_contour = poly.translate(right_contour, centering_shift[0], -centering_shift[1])
    left_holes = poly.translate(left_holes, centering_shift[0], -centering_shift[1])
    right_holes = poly.translate(right_holes, centering_shift[0], -centering_shift[1])

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
        cam.start_spindle(20000),
        cam.dwell(3),
        cam.comment("Right Hinge Pocket"),
        cam.pocket(right_hinge_pocket_contours, -abs(pocket_depth), retract=-abs(pocket_depth)),
        cam.rapid([None, None, 20.0]),
        cam.comment("Left Hinge Pocket"),
        cam.pocket(left_hinge_pocket_contours, -abs(pocket_depth), retract=-abs(pocket_depth)),
        cam.rapid([None, None, 20.0]),
        cam.comment("Hinge Holes"),
        cam.change_tool("1mm drill"),
        cam.start_spindle(5000),
        cam.dwell(2),
        [cam.rmp(p + [drill_depth], retract=10.0) for p in right_holes],
        [cam.rmp(p + [drill_depth], retract=10.0) for p in left_holes],
        cam.rapid([None, None, 20.0]),
    ]
    return r

def bevel_temple_ends(path1, path2, thinning):
    depth = -4 - thinning
# Assume a 14 degree dovetail cutter, cutting 5mm from bottom
    dovetail_offset = 12.35/2 - 1.25
    return [
            cam.comment("Bevel temple ends"),
            cam.change_tool("dovetail"),
            cam.start_spindle(20000),
            cam.feedrate(750),
            cam.rapid(path1[0]),
            cam.rapid([None, None, depth]),
            cam.move(path1[1]),
            cam.move(path1[2]),
            cam.move(path1[0]),
            cam.rapid([None, None, 10]),
            cam.rapid(path2[0]),
            cam.rapid([None, None, depth]),
            cam.move(path2[1]),
            cam.move(path2[2]),
            cam.move(path2[0]),
            cam.move([None, None, 20])
            ]

def temple_hinge_pockets(temples, thinned):
    # We're operating in a 90 degree rotated fixture
    #l_hinge = poly.rotate_90(temples["left_hinge_contour"])
    #r_hinge = poly.rotate_90(temples["right_hinge_contour"])
    print 'Generating temple hinge pockets'

    l_hinge = temples["left_hinge_contour"]
    r_hinge = temples["right_hinge_contour"]
    if not poly.is_ccw(l_hinge):
        l_hinge = poly.reverse(l_hinge)
    if not poly.is_ccw(r_hinge):
        r_hinge = poly.reverse(r_hinge)

    #pocket_depth = temples['pocket_depth'] + thinned
    pocket_depth = 1 + thinned;

    def pocket_contours(contour):
        contours = []
        erode = poly.erode(1.5875/2, contour)
        making_progress = True
        while len(erode) > 0 and making_progress:
            making_progress = False
            for path in erode:
                if len(path) > 5:
                    making_progress = True
                    contours.append(path)
            erode = poly.erode(1.5875/2, contours[-1])
        return contours

    left_hinge_pocket_contours = pocket_contours(l_hinge)
    right_hinge_pocket_contours = pocket_contours(r_hinge)
#    left_hinge_pocket_contours = [];
#    while len(l_hinge) > 0:
#        l_hinge = poly.erode(1.5875/2, l_hinge)
#        if len(l_hinge) > 0:
#            l_hinge = l_hinge[0]
#            left_hinge_pocket_contours.append(l_hinge)
#    right_hinge_pocket_contours = [];
#    while len(r_hinge) > 0:
#            r_hinge = poly.erode(1.5875/2, r_hinge)
#            if len(rhinge_) == 1:
#                right_hinge_pocket
#            if len(r_hinge) > 0:
#                r_hinge = r_hinge[0]
#                right_hinge_pocket_contours.append(r_hinge)
    r = [
        cam.comment("Hinge Pockets"),
        cam.feedrate(750),
        cam.change_tool("1/16in endmill"),
        cam.start_spindle(22000),
        cam.dwell(5),
        cam.comment("Right Hinge Pocket"),
        cam.pocket(right_hinge_pocket_contours, -abs(pocket_depth), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Left Hinge Pocket"),
        cam.pocket(left_hinge_pocket_contours, -abs(pocket_depth), retract=0),
        cam.rapid([None, None, 20.0]),
        cam.comment("Hinge Holes"),
        cam.change_tool("1mm drill"),
        cam.start_spindle(5000),
        cam.dwell(2),
        [cam.rmp(p + [-8.0], retract=10.0) for p in temples['right_hinge_holes']],
        [cam.rmp(p + [-8.0], retract=10.0) for p in temples['left_hinge_holes']],
        cam.rapid([None, None, 20.0]),
    ]
    return r



def lens_groove(left_c, right_c, height):
    """Generates the toolpath for the lens holes (holes, groove and tabs)."""
    print 'Generating lens grooves'
    if not poly.is_ccw(left_c):
        left_c = poly.reverse(left_c)
    if not poly.is_ccw(right_c):
        right_c = poly.reverse(right_c)

    lgroove = poly.erode(0.85, left_c)[0]
    rgroove = poly.erode(0.85, right_c)[0]

    left_entry = poly.erode(7.0, left_c)[0][0];
    right_entry = poly.erode(7.0, right_c)[0][0];
    r = [
        "(Lens Grooves)",
        cam.change_tool("vgroove"),
        cam.start_spindle(20000),
        cam.dwell(5),
        cam.feedrate(2000),
        cam.rmp(right_entry + [height]),
        cam.contour(rgroove, True),
        cam.move(right_entry), # Get out from under the overhang
        cam.rmp(left_entry + [height]),
        cam.contour(lgroove, True),
        cam.move(left_entry), # Get out from under the overhang
    ]
    return r


def lens_holes(left_c, right_c, thickness):
    """Generates the toolpath for the lens holes (holes, groove and tabs)."""
    print 'Calculating the lens holes'
    tool_radius = 3.175
    if not poly.is_ccw(left_c):
        left_c = poly.reverse(left_c)
    if not poly.is_ccw(right_c):
        right_c = poly.reverse(right_c)

#    drawing = dxf.drawing('test.dxf')
#    drawing.add_layer('OUTLINE', color=1)
#    polyline = dxf.polyline(layer="OUTLINE")
#    polyline.add_vertices(left_c)
#    drawing.add(polyline)


    lhole = poly.erode(tool_radius/2.0, left_c)[0]
    print 'left hole eroded', lhole[0], lhole[-1]
    rhole = poly.erode(tool_radius/2.001, right_c);
    print 'right hole eroded', lhole[0], lhole[-1]
    rhole = rhole[0]
#    polyline = dxf.polyline(layer="OUTLINE")
#    polyline.add_vertices(lhole)
#    drawing.add(polyline)


    right_rough = poly.erode((tool_radius + 0.3)/2, right_c)[0]
    left_rough = poly.erode((tool_radius+0.3)/2, left_c)[0]
    print 'Calculated roughing passes'
    #lgroove = poly.erode(0.8, left_c)[0]
    #rgroove = poly.erode(0.8, right_c)[0]

    left_entry = poly.erode(5.0, left_c)[0][0];
    right_entry = poly.erode(5.0, right_c)[0][0];

    print 'calculated entry points'
    lhole = poly.reverse(lhole)
    rhole = poly.reverse(rhole)

    r = [
        "(Lens Holes)",
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(22000),
        cam.feedrate(2000),
        cam.rmh(right_entry + [-thickness - 1.0], 1.5, 0.5, 1.0),
        cam.contour(right_rough, True),
        cam.feedrate(1000),
        cam.contour(rhole, True),
        cam.feedrate(2000),
        cam.rmh(left_entry + [-thickness - 1.0], 1.5, 0.5, 1.0),
        cam.contour(left_rough, True),
        cam.feedrate(1000),
        cam.contour(lhole, True),

    ]
    return r

def bevel_entry_hole(temples, thickness):

    hole_radius = 8
    tool_radius = 3.175/2 # 1/2 inch
    helix_radius = hole_radius - tool_radius
    r = [
        cam.comment("Entry holes for temple bevelling operation"),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(22000),
        cam.feedrate(800),
        cam.dwell(5),

    ]

def index_holes(thickness):
# We put the index holes 1/2 between top and bottom, 160mm apart
    hole_radius = 4.85/2 # Measured from dowel pin
    tool_radius = 3.175/2
    helix_radius = hole_radius - tool_radius
    r = [
        cam.comment("Index holes for secondary operations"),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(22000),
        cam.feedrate(1000),
        cam.dwell(5),
        cam.rmh([0, 90, -thickness - 1.0], helix_radius, 0.5, 1),
        cam.rmh([0, -90, -thickness - 1.0], helix_radius, 0.5, 1),
        ]
    return r


def nose_pads(order, thickness):
    return []

def temple_dxf(o, filename, hinge, reverse):
    temple_contour = poly.mirror_y(o['ltemple_con'], False)

    # Location in y for the center of the hinge
    hinge_y = -(temple_contour[0][1]+temple_contour[-1][1])/2

    temple_contour = poly.mirror_y(temple_contour, True)
    # Only mirror the temple contour, not the hinge, because the hinge is modelled
    # for the RIGHT hinge and we have the LEFT temple
    hinge_contour = hinges.get_hinge(hinge)["temple_contour"]
    hinge_holes = hinges.get_hinge(hinge)["temple_holes"]

    # Temple should be at 6 degrees.  We tilt the hinge instead for manufacturing
    # ease.
#    hinge_contour = poly.rotate(hinge_contour, 6)
#    hinge_holes = poly.rotate(hinge_holes, 6)

    # 0 point of hinge is center of straight edge of hinge.  We want the
    # center of the hinge box.
    hinge_center_offset = (poly.top(hinge_contour)+poly.bottom(hinge_contour))/2
    y_offset = hinge_y-hinge_center_offset

    # Since we're tilted 6 degrees compared to design space we need to
    # compensate on the x axis
    x_offset = y_offset * math.tan(math.radians(6)) # This adjusts the hinge for the move to the center of the temple

    hinge_contour = poly.translate(hinge_contour, x_offset, y_offset)
    hinge_holes = poly.translate(hinge_holes, x_offset, y_offset)

    # Reverse for the right temple
    if reverse:
        temple_contour = poly.mirror_x(temple_contour, False)
        hinge_contour = poly.mirror_x(hinge_contour, False)
        hinge_holes = poly.mirror_x(hinge_holes, False)

    drawing = dxf.drawing(filename)
    drawing.add_layer("CUT", color=1)
    drawing.add_layer("POCKET", color=2)
    drawing.add_layer("TEXT", color=3)
    polyline = dxf.polyline(layer="CUT", color=1, thickness=0.1)
    polyline.add_vertices(temple_contour)
    drawing.add(polyline)
    polyline = dxf.polyline(layer="POCKET", color=2, thickness = 0.1)
    polyline.add_vertices(hinge_contour)
    drawing.add(polyline)
    for hole in hinge_holes:
        c = dxf.circle(0.5, hole, layer="CUT", color=1)
        drawing.add(c)
    if not reverse:
        # Left temple, add our company name
        # First find Y position 15mm from temple top corner
        xPos = temple_contour[0][0]-15
        yPos = -1
        yStartPos = -1
        for pt in reversed(temple_contour):
            if pt[0] < xPos and yPos < 0:
                yPos = pt[1]-3.5
            elif pt[0] < (xPos - 25):
                yStartPos = pt[1] - 3.5
                break
        angle = math.degrees(math.atan((yPos-yStartPos)/25))
        print 'text angle', angle
        text = dxf.text("GUILD eyewear", style="TIMES", color=3, rotation=angle, height=2.5, alignpoint=[xPos, yPos], halign=2, valign=2)
        drawing.add(text)


    drawing.save()



def arrange_temple_curves(left_temple_contour, hinge):
    right_temple_contour = poly.mirror_x(left_temple_contour, False)

    left_hinge = hinges.get_hinge(hinge)
    right_hinge = hinges.get_hinge(hinge, False)
    hinge_y = (left_temple_contour[0][1]+left_temple_contour[-1][1])/2

    left_holes = left_hinge['temple_holes']
    right_holes = right_hinge['temple_holes'] # opposite direction.  FIX someday
    right_hinge_contour = right_hinge['temple_contour']
    left_hinge_contour = left_hinge['temple_contour']
    print 'right hinge endpoints', right_hinge_contour[0], right_hinge_contour[-1], right_hinge_contour[-2]

    # y offset is mucked up because it's calculated from the center of the hinge boxing
    # square, not the center of the straight edge
    highpoint = poly.top(left_hinge_contour);
    lowpoint = poly.bottom(left_hinge_contour);
    y_offset = hinge_y-(lowpoint + (highpoint-lowpoint)/2.0) # Distance to move hinge to center it on temple

    # Since we're tilted 6 degrees compared to design space we need to
    # compensate on the x axis
    x_offset = y_offset * math.tan(math.radians(6)) # This adjusts the hinge for the move to the center of the temple

    x_offset = x_offset - 1 # Compensate for later bevelling of temple end

    left_hinge_contour = poly.translate(left_hinge_contour, x_offset, y_offset)
    left_holes = poly.translate(left_holes, x_offset, y_offset)
    right_hinge_contour = poly.translate(right_hinge_contour, -x_offset, y_offset)
    right_holes = poly.translate(right_holes, -x_offset, y_offset)

    left_temple_contour = poly.rotate_90(left_temple_contour)
    left_hinge_contour = poly.rotate_90(left_hinge_contour)
    left_holes = poly.rotate_90(left_holes)
    right_temple_contour = poly.rotate(right_temple_contour, 90)
    right_hinge_contour = poly.rotate(right_hinge_contour, 90)
    right_holes = poly.rotate(right_holes, 90)
    # Move them apart as much as possible
    height = poly.right(left_temple_contour) - poly.left(left_temple_contour)
    left_temple_contour = poly.translate(left_temple_contour, height, 0)
    left_holes = poly.translate(left_holes, height, 0)
    left_hinge_contour = poly.translate(left_hinge_contour, height, 0)
# Get the thing as horizontal as possible
    # When top-bottom distance is minimized, it's horizontal
#    height = poly.right(left_temple_contour) - poly.left(left_temple_contour)
#    opt_angle = 0
#
#    for angle in range(2, 41):
#        candidate_contour = poly.rotate(left_temple_contour, -angle)
#        candidate_height = poly.right(candidate_contour)-poly.left(candidate_contour)
#        if candidate_height < height:
#            height = candidate_height
#            opt_angle = angle
#        else:
#            break
##    print 'rotating temples by {0} degrees'.format(opt_angle)
#    left_temple_contour = poly.rotate(left_temple_contour, -opt_angle);
#    left_holes = poly.rotate(left_holes, -opt_angle)
#    left_hinge_contour = poly.rotate(left_hinge_contour, -opt_angle)

#    right_temple_contour = poly.rotate(right_temple_contour, opt_angle);
#    right_holes = poly.rotate(right_holes, opt_angle)
#    right_hinge_contour = poly.rotate(right_hinge_contour, opt_angle)

    temple_width = poly.top(left_temple_contour) - poly.bottom(left_temple_contour)
#    print 'temple width is {0}'.format(temple_width)
#
#    # Move left one left, right one right to offset them
    optimization_offset = (150-temple_width)/2.0
    left_temple_contour = poly.translate(left_temple_contour, 0, optimization_offset)
    left_holes = poly.translate(left_holes, 0, optimization_offset)
    left_hinge_contour = poly.translate(left_hinge_contour, 0, optimization_offset)

    right_temple_contour = poly.translate(right_temple_contour, 0, -optimization_offset)
    right_hinge_contour = poly.translate(right_hinge_contour, 0, -optimization_offset)
    right_holes = poly.translate(right_holes, 0, -optimization_offset)



    # Make sure they're not intersecting
    translation_step = 2.5
    intersection = poly.intersection(left_temple_contour, right_temple_contour)
    while(len(intersection) > 0):
        left_temple_contour = poly.translate(left_temple_contour, translation_step, 0)
        left_holes = poly.translate(left_holes, translation_step, 0)
        left_hinge_contour = poly.translate(left_hinge_contour, translation_step, 0)
        right_temple_contour = poly.translate(right_temple_contour, -translation_step, 0)
        right_holes = poly.translate(right_holes, -translation_step, 0)
        right_hinge_contour = poly.translate(right_hinge_contour, -translation_step, 0)
        intersection = poly.intersection(left_temple_contour, right_temple_contour)
    # Move them together until they touch (may will not have moved at all above because they are far)
    # then back off.
    translation_step = 1.0
    while(len(intersection) == 0):
        print "still intersecting", len(intersection)
        left_temple_contour = poly.translate(left_temple_contour, -translation_step, 0)
        left_holes = poly.translate(left_holes, -translation_step, 0)
        left_hinge_contour = poly.translate(left_hinge_contour, -translation_step, 0)
        right_temple_contour = poly.translate(right_temple_contour, translation_step, 0)
        right_holes = poly.translate(right_holes, translation_step, 0)
        right_hinge_contour = poly.translate(right_hinge_contour, translation_step, 0)
        intersection = poly.intersection(left_temple_contour, right_temple_contour)

    # We're just overlapping, so now back off
    translation_step = 3
    left_temple_contour = poly.translate(left_temple_contour, translation_step, 0)
    left_holes = poly.translate(left_holes, translation_step, 0)
    left_hinge_contour = poly.translate(left_hinge_contour, translation_step, 0)
    right_temple_contour = poly.translate(right_temple_contour, -translation_step, 0)
    right_holes = poly.translate(right_holes, -translation_step, 0)
    right_hinge_contour = poly.translate(right_hinge_contour, -translation_step, 0)


#    # sanity check that we fit on stock
    total_width =  poly.right(left_temple_contour) - poly.left(right_temple_contour)
    if total_width > 80:
        print 'Error! temples did not pack tight enough.', total_width
        raise 'Sizing error'

    # Now they're packed togegther OK but not centred on the stock
    #Midpoint of curves should be 0, so that's how much we'll shift it
    y_centering_shift = poly.bottom(right_temple_contour) + (poly.top(left_temple_contour) - poly.bottom(right_temple_contour))/2.0
    x_centering_shift = (poly.right(left_temple_contour) + poly.left(right_temple_contour))/2;



    left_temple_contour = poly.translate(left_temple_contour, -x_centering_shift, -y_centering_shift)
    right_temple_contour = poly.translate(right_temple_contour, -x_centering_shift, -y_centering_shift)
    left_hinge_contour = poly.translate(left_hinge_contour, -x_centering_shift, -y_centering_shift)
    right_hinge_contour = poly.translate(right_hinge_contour, -x_centering_shift, -y_centering_shift)
    left_holes = poly.translate(left_holes, -x_centering_shift, -y_centering_shift)
    right_holes = poly.translate(right_holes, -x_centering_shift, -y_centering_shift)

    # The temple contours aren't closed, close them here
    # NOTE: Taken out because we close the curves with a little arc to account for front curvature
    #left_temple_contour.append(left_temple_contour[0])
    #right_temple_contour.append(right_temple_contour[0])

    return {
        "pocket_depth": 1,
        "left_hinge_contour": left_hinge_contour,
        "right_hinge_contour": right_hinge_contour,
        "left_hinge_holes": left_holes,
        "right_hinge_holes": right_holes,
        "left_temple_contour": left_temple_contour,
        "right_temple_contour": right_temple_contour
            }



def create_dxf(filename, polys, close=False, close_with_arc=False, right_temple_text=None):
    drawing = dxf.drawing(filename)
    drawing.add_layer('OUTLINE', color=1)
    drawing.add_layer('TEXT', color=2)

    for p in polys:
        polyline = dxf.polyline(layer="OUTLINE", thickness=0.1)
        if close:
            p = p + [p[0]]  # Close the polygon to avoid a cusp
        elif close_with_arc and p[0] != p[-1]:
            # Temples get a little arc to account for frame curvature after bending.
            # Frame is bent at base 6 curvature, which is 88 mm radius sphere.  Plane
            # intersecting sphere at 30mm (about the offset of most temples) is 80mm radius.
            # First find the center of the arc

            center_point = (
                    p[-1][0] + (p[0][0]-p[-1][0])/2,
                    p[-1][1] + (p[0][1]-p[-1][1])/2)

            vect = (p[-1][0]-p[0][0], p[-1][1]-p[0][1])
            scale = (vect[0]**2 + vect[1]**2) ** 0.5
            unit_vect = (vect[0]/scale, vect[1]/scale)
            vect = (unit_vect[0]*88, unit_vect[1]*88)
            vect = (vect[1], -vect[0])
            center_point = (center_point[0]+vect[0], center_point[1]+vect[1])

            # We set it arbitrarily at 79 mm but the actual radius will be to the end points
            radius = ((center_point[0]-p[0][0])**2 + (center_point[1]-p[0][1])**2) **0.5
            print 'radius', radius
            angle1 = math.atan2(p[0][1]-center_point[1], p[0][0] - center_point[0])
            angle2 = math.atan2(p[-1][1]-center_point[1], p[-1][0] - center_point[0])
            print 'angle of p0', math.degrees(angle1)
            print 'angle of pn1', math.degrees(angle2)

            drawing.add(dxf.arc(radius, center_point, math.degrees(angle2), math.degrees(angle1), layer="OUTLINE", thickness=0.1))


        polyline.add_vertices(p)
        drawing.add(polyline)

    if right_temple_text:
        p = polys[1]
        insertion_point = (
                p[-1][0] + (p[0][0]-p[-1][0])/2,
                p[-1][1] + (p[0][1]-p[-1][1])/2)

        vect = (p[-1][0]-p[0][0], p[-1][1]-p[0][1])
        scale = (vect[0]**2 + vect[1]**2) ** 0.5
        unit_vect = (vect[0]/scale, vect[1]/scale)
        vect = (unit_vect[0]*15, unit_vect[1]*15)
        vect = (vect[1], -vect[0])
        insertion_point = (insertion_point[0]+vect[0], insertion_point[1]+vect[1]-1)

        bottom_vect = (p[100][0]-p[0][0], p[100][1]-p[0][1])
        text_angle = math.atan2(bottom_vect[1], bottom_vect[0])
        print 'text angle', text_angle
        txt = dxf.text(
                right_temple_text,
                insertion_point,
                rotation=math.degrees(text_angle),
                style="ARIAL",
                layer="TEXT",
                height=2.0
                )
        txt2 = dxf.text(
                "made with love by GUILD eyewear",
                (insertion_point[0]+0.5, insertion_point[1]-3),
                style="TIMES",
                rotation=math.degrees(text_angle),
                layer="TEXT",
                height=2.0
                )
        drawing.add(txt);
    #    drawing.add(txt2);
    drawing.save()

def check_frame_size(contour):
    # We're limited by our stock size and clamp clearances
    # We can be about 160mm wide max, and about 65mm high max.
    if abs(poly.top(contour)) > 86:
        return "Frame is too wide: %f mm" % (poly.top(contour) * 2)
    if poly.right(contour) - poly.left(contour) > 60:
        return "Frame is too tall: %f mm" % (poly.right(contour) - poly.left(contour))
    return None


def frame_offset(contour):
    """
    The origin of the glasses is in the middle of the pair (i.e. centered at Y) with the
    X origin at the line between the pupils. The fixture on the milling machine has its origin
    centered on the Y axis but with the X axis at the edge of the clamp.  We need to shift the
    machine origin toward the middle of the stock.
    """
    xoffset = poly.right(contour) + 12.5 # Offset by the furthest point, plus some extra for the tool
#Note: we don't actually need room for the tool since we're laser cutting the outside.
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
