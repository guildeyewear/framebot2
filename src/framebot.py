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
    if options.outname == None:
        print "extracting out name from ", options.infile
        outname = options.infile.split("/")[-1]
        outname = outname.split(".")[0]
        outname = datetime.date.today().isoformat() + "/" + outname
        print "out file is ", outname
    else:
        outname = options.outname
    # Create the output director
    #order_dir = "/Users/rodfrey/Dropbox/guild_orders/" + outname
    order_dir = "/Volumes/Untitled/" + outname
    if not os.path.exists(order_dir):
        os.makedirs(order_dir)


    # Create the milling program for the lens holes/groove, hinge pockets, etc.
    # and the dxf for the laser for the fronts
    mill_fronts(order_dir, o)
    # Rearrange the face curve a bit so the laser starts outside the curve.
    # Otherwise it leaves a little scar where it starts.
    laser_face = poly.new_start(o['face_con'], poly.leftmost_index(o['face_con']))
    laser_face.append(laser_face[0]) # Close the curve

    first_vector = [laser_face[1][0]-laser_face[0][0], laser_face[1][1]-laser_face[0][1]]
    unit = math.sqrt(first_vector[0]**2 + first_vector[1]**2)
    leadin_vector = [(4*first_vector[0])/unit,(4*first_vector[1])/unit]
    leadin_start = [laser_face[0][0] - leadin_vector[0], laser_face[0][1] - leadin_vector[1]]

    last_vector = [laser_face[-2][0]-laser_face[-1][0], laser_face[-2][1]-laser_face[-1][1]]
    unit = math.sqrt(last_vector[0]**2 + last_vector[1]**2)
    leadout_vector = [(4*last_vector[0])/unit,(4*last_vector[1])/unit]
    leadout_end = [laser_face[-1][0] + leadout_vector[0], laser_face[-1][1] + leadout_vector[1]]


#    rise = laser_face[0][0] - laser_face[1][0]
#    run = laser_face[0][1] - laser_face[1][1]
#    slope = rise/run
#    leadin = 4.0
#    leadin_y = math.sqrt((leadin*leadin)/(1+slope*slope))
#    leadin_x = slope*leadin_y
#    leadin_point =  [leadin_x+laser_face[0][0], leadin_y+laser_face[0][1]]

    laser_face.insert(0, leadin_start)
#    laser_face.append(leadout_end)

#    leadout = 4.0
#    rise = laser_face[-2][0] - laser_face[-1][0]
#    run = laser_face[-2][1] - laser_face[-1][1]
#    slope = rise/run
#    leadout_y = math.sqrt((leadin*leadin)/(1+slope*slope))
#    leadout_x = slope*leadin_y
#    leadout_point =  [laser_face[-1][0]-leadout_x, laser_face[-1][1]-leadout_y]
#    laser_face.append(leadout_point)
    create_dxf(order_dir + "/face_contour.dxf", [laser_face], close=False)


    # Arrange the temples to fit on the stock
    temples =  arrange_temple_curves(o['ltemple_con'], o['rtemple_con'], o['lhinge'], o['lhinge_y'], o['rhinge_y'])
    left_hinge = poly.rotate_90(temples['left_hinge_contour'])
    right_hinge = poly.rotate_90(temples['right_hinge_contour'])

    l_temple = poly.rotate_90(temples['left_temple_contour']);
    r_temple = poly.rotate_90(temples['right_temple_contour']);

    create_dxf(order_dir + "/temple_contour.dxf",
            [poly.rotate_90(temples['left_temple_contour']),
             poly.rotate_90(temples['right_temple_contour']),
#             left_hinge,
#             right_hinge,
            ],
#            close_with_arc=True,
#           right_temple_text='Custom made for Tim Sullivan',
            close=True)

    mill_temples(order_dir, temples)

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

#TODO: Replace with information in materials database
# Initial thickness of the forward-facing lamination.  0 if no lamination.
    total_thickness = 6
    thin_front = 0
    thin_back = 0

    #back_surface_thickness = 6


# The thickness of the highest point of the nosepad
    #nosepad_thickness = 4

# Final thickness of the main part of the frame.  Must be equal to or less than
# the back_surface_thickness.
    #final_back_thickness =6

    #thickness = final_front_thickness + final_back_thickness

    #front_surface_removal = final_front_thickness - front_surface_thickness
    #back_surface_removal = final_back_thickness - back_surface_thickness


# The machine has the stock clamp oriented 90 degrees to the way the
# software creates the contours.
    face_c = poly.rotate_90(order["face_con"])
    left_lens_c = poly.rotate_90(order["lhole_con"])
    right_lens_c = poly.rotate_90(order["rhole_con"])
    temple_height = abs(order["ltemple_con"][0][1] - order["ltemple_con"][-1][1])

    msg = check_frame_size(left_lens_c)
    if msg:
        print msg
        sys.exit(0)
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

#    groove_height = back_surface_removal + (thickness/2)
    groove_height = -(total_thickness - thin_front - 2)
    size_info = order['usersizes']
    #SAVI
#    size_info['nose_radius'] = 6
#    size_info['nose_height'] = 8
#    size_info['nose_splayangle'] = 25
#    size_info['nose_ridgeangle'] =36

    #ROD
#    size_info['nose_radius'] = 8
#    size_info['nose_height'] = 8
#    size_info['nose_splayangle'] = 38
#    size_info['nose_ridgeangle'] =36
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
#        cam.temporary_offset(offset),
# Note that X and Y parameters in the order are switched from our system
#TODO: replace the thickness offset with the thickness of the TEMPLE, not the fronts.
        index_holes([face_c], total_thickness),
        lens_holes(left_lens_c, right_lens_c, total_thickness),
        #lens_groove(left_lens_c, right_lens_c, back_surface_removal - (thickness/2.0)),
        lens_groove(left_lens_c, right_lens_c, groove_height),
        face_hinge_pockets(order["lhinge"], order["lhinge_y"], order["ltemple_x"], (-x_shift, -y_shift), thin_back),
        #nose_pads(order, thickness),
        nose_contour(
            float(size_info["nose_radius"]),
            float(size_info["nose_height"]),
            float(size_info["nose_splayangle"]),
            float(size_info["nose_ridgeangle"]),
            face_c, total_thickness,thin_back, -x_shift ),
        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        cam.change_tool("1/8in endmill"),
        #cam.rapid([face_c[0][0], face_c[0][1], -thin_back] ),
        #cam.contour(face_c, True),
        cam.end_program(),
    ]
    open(outdir + "/face_stage1.ngc", "w").write(to_string(program))

def nose_contour(nose_rad, nose_h, nose_sa, nose_ra, face_con, thickness, thin_back, centering_shift):
    """Creates the nose contour feature toolpath.  Angular arguments are in degrees."""#
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
    print 'nose height', nose_h
    print 'centering shift', centering_shift

    sa = math.radians(nose_sa)
    ra = math.radians(nose_ra)
    h = nose_h + centering_shift

    xfloor = poly.left(face_con) - 3.175  # bottom most point minus tool radius
    xfloor = max(xfloor, -27.0) # miminum safe distance without hitting clamp

    z_depth = 0.5 # Start a bit above the surface of the glasses
    nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, cutter_offset, z_depth)

    r = [
        "(Nose Contour)",
        cam.change_tool("1/4in ballmill"),
        cam.start_spindle(20000),
        cam.feedrate(2000),
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


def outline(face_con, thickness):
    """Creates the face outline feature toolpath.  face_con is the polygon representing the face outline."""

    rough =  poly.dilate(6.35/4.0 + 0.15, face_con)
    rough_ramp = ramped_outline(rough, 0.0, -thickness+0.5, 20.0)

    finish = poly.dilate(3.175/2, face_con)
    finish  = [p + [-thickness+0.5] for p in finish]

    tabs = poly.dilate(3.175/2+0.05, face_con)
    add_tabs(tabs, -thickness-1,-thickness+0.3, 12, 6)

    #return poly.ramp(poly.segment(polyline, poly.polyline_length(polyline) - lead_length, poly.polyline_length(polyline), True), top, bottom, False)
#    roughing_depth = -thickness + 0.5;
#    finish_depth = -thickness + 0.1;

#    finish_ramp = poly.ramp(poly.segment(finish, poly.polyline_length(finish) - 50.0, poly.polyline_length(finish), True), 0.0, finish_depth, False)

#    tab_lengths = poly.intercepts(finish, True, y=[-39.0, 39.0])
#    finish_and_tabs = poly.make_gaps(finish, tab_lengths, 12, finish_depth+1.5, finish_depth, True)

#    left_temple_rough = poly.dilate(3.175/2.0 + 0.1, l_temple)
#    right_temple_rough = poly.dilate(3.175/2.0 + 0.1, r_temple)

#    left_temple_ramp = ramped_outline(left_temple_rough, 0.0, -thickness+0.5, 10.0)
#    right_temple_ramp = ramped_outline(right_temple_rough, 0.0, -thickness+0.5, 10.0)
#
#    left_temple_finish = poly.dilate(3.175/2.0, l_temple)
#    right_temple_finish = poly.dilate(3.175/2.0, r_temple)
#    left_temple_finish = [p + [-thickness+0.5] for p in left_temple_finish]
#    right_temple_finish = [p + [-thickness+0.5] for p in right_temple_finish]

#    left_temple_tabs = poly.dilate(3.175/2.0+0.05, l_temple)
#    right_temple_tabs = poly.reverse(poly.dilate(3.175/2.0+0.05, r_temple))
#    add_tabs(left_temple_tabs, -thickness-1, -thickness+0.3, 10, 6)
#    add_tabs(right_temple_tabs, -thickness-1, -thickness+0.3, 10, 6)

#
    r = [
        "(Outline Rough)",
        cam.feedrate(1500),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(22000),
        cam.dwell(5),
        # Ramp in
        cam.rmp(rough_ramp[0]),
        cam.feedrate(1500),
        cam.contour(rough_ramp, False),
        cam.contour(rough, False),
        cam.feedrate(700),
        cam.contour(finish, True),
        cam.contour(tabs, True),
#        cam.contour(rough, True),

        "(Outline Finish and Tabs)",
        #cam.change_tool(1), # 4 flute 1/8" flat
        #cam.start_spindle(4500),

#        rmp(finishramp[0]),
        ["G1 Z-5"],
#        cam.contour(finish, False),
#        cam.contour(finish_and_tabs, True),
        cam.rapid([None, None, 20.0]),
    ]

    return r


def contour_face(body_removal, hinge_removal, nosepad_removal, temple_height, face_c, lens_c, x_pos):
    ''' Create the heightmap of the frame, surfacing the back and adding thickness for the
    hinge location and the nosepads.  '''
    if body_removal == hinge_removal == nosepad_removal == 0:
        return [] # Nothing to do

    cutter_radius = 6.35/2 # 3/4 inch cutter
    entry_point = [x_pos, 110, 0]

    facing_contour = poly.dilate(0.05, lens_c)

# Reshape the facing contour so the first point is near the hinge
    center_y = poly.bottom(facing_contour) + (poly.top(facing_contour) - poly.bottom(facing_contour))/2
    center_x = poly.right(facing_contour) + (poly.left(facing_contour) - poly.right(facing_contour))/2
    split_idx = -1
    for idx, pt in enumerate(facing_contour):
        if pt[1] > center_y and (idx+1) < len(facing_contour):
            if (pt[0] < x_pos and facing_contour[idx+1][0] > x_pos) or (pt[0] > x_pos and facing_contour[idx+1][0] < x_pos):
                split_idx = idx
                break
    if split_idx < 0:
        print 'Error contouring back of frame: could not locate entry point for surfacing cut'
        return []
    facing_contour = poly.new_start(facing_contour, split_idx)
# Ensure we're going clockwise, i.e. starting at the hinge and moving up over the frame
    if poly.is_ccw(facing_contour):
        facing_contour = poly.reverse(facing_contour)

# Calculate the Z values
# We'll need a few helper values.  nosepad_start is the inflection point of the nose bridge.
    nosepad_start = max([pt[0] for pt in face_c if pt[1] == 0]) + cutter_radius
    hinge_rampdown_start_x = x_pos + temple_height/2 + cutter_radius
    hinge_rampdown_start_y = facing_contour[0][1] - cutter_radius
    hinge_rampup_start_x = x_pos - temple_height/2 - cutter_radius
    hinge_rampup_start_y = facing_contour[0][1] - cutter_radius

    '''
    Arbitrary heuristic, adjusted for aesthetics.
    1. If we're past the center point of the lens hole, we're either on the body
    of the frame or over the raised hinge point.
    2. If we're before the center point we're either on the body or over the nosepiece.

    1a. If we're above the cutter-radius-adjusted top of the temple, we're ramping down
    1b. If we're below the cutter-radius-adjusted bottom of the temple, we're ramping up
    1c. Otherwise we're at body thickness

    2a. If we're above the top of the nose cutout, we're at body thickness
    2b. When we reach nose cutout, we do a s-curve over 3 mm to nosepad height
    2c. Continue for length of cutter diameter to get rear of cutter over highest point
    2d. Continue for 10mm
    2e. S-curve down over 10mm
    '''
    def add_hinge_heights(contour):
        heightmap = []
        over_hinge = True # Start over hinge

        items_to_skip = 0 # for fast-forwarding enumeration
        for idx, pt in enumerate(contour):
            if items_to_skip > 0:
                items_to_skip = items_to_skip - 1
                continue


            if pt[1] < center_y:
                heightmap = heightmap + [pt]
            # Going up and around: start ramping down when we're clear of X or Y
            elif pt[0] > x_pos:
                if pt[0] > hinge_rampdown_start_x or pt[1] < hinge_rampdown_start_y:
                    if(over_hinge): # starting transition
                        transition_length = poly.polyline_length(contour[:(idx+1)], False)
                        ramp_segment = poly.segment(contour, transition_length, transition_length+5, False)
                        ramp_segment = poly.ramp(ramp_segment, hinge_removal, body_removal, False)
                        heightmap = heightmap + ramp_segment[:-1]
                        items_to_skip = len(ramp_segment)
                        over_hinge = False
                    else: # past transition but still on hinge side of lens hole
                        heightmap = heightmap + [pt + [body_removal]]
                else: # We're on the top part but haven't reached the transition yet
                    heightmap = heightmap + [pt + [hinge_removal]]

            # Coming back up to the hinge: start ramping up if we encroach on both x and y
            elif pt[0] < x_pos and (pt[0] > hinge_rampup_start_x and pt[1] > hinge_rampdown_start_y):
                if(not over_hinge): # starting transition
                    transition_length = poly.polyline_length(contour[:(idx+1)], False)
                    ramp_segment = poly.segment(contour, transition_length, transition_length+5, False)
                    ramp_segment = poly.ramp(ramp_segment, body_removal, hinge_removal, False)
                    heightmap = heightmap + ramp_segment
                    items_to_skip = len(ramp_segment)
                    over_hinge = True
                else: # Over flat hinge area
                    heightmap = heightmap + [pt + [hinge_removal]]
            else: # We're over the body area but back on the hinge side
                heightmap = heightmap + [pt + [body_removal]]
        return heightmap

    def add_nosepad_heights(contour):
        heightmap = []
        over_nosepad = False
        past_nosepad = False
        nosepad_flat_idx = -1

        items_to_skip = 0 # for fast-forwarding the enumeration
        for idx, pt in enumerate(contour):
            if items_to_skip > 0:
                items_to_skip = items_to_skip-1
                continue
            if pt[1] >= center_y:
                heightmap = heightmap + [pt]
            elif not over_nosepad and not past_nosepad:
                if pt[0] < nosepad_start: # Transition
                    transition_length = poly.polyline_length(contour[:(idx+1)], False)
                    ramp_segment = poly.segment(contour, transition_length, transition_length+5, False)
                    ramp_segment = poly.ramp(ramp_segment, body_removal, nosepad_removal, False)
                    heightmap = heightmap + ramp_segment[:-1]
                    items_to_skip = len(ramp_segment)
                    nosepad_flat_idx = idx + items_to_skip # we'll need this to go down
                    over_nosepad = True
                else: # we're past the nosepad
                    heightmap = heightmap + [pt + [body_removal]]
            elif over_nosepad and not past_nosepad:
                if nosepad_flat_idx < 0:
                    print "ERROR! I think I'm on the nosepad but have not transitioned yet"
                    return []
                # We'll be cutting the far side with the back of the cutter, so need to move at
                # least the diameter to get any flat at all
                flat_length = poly.polyline_length(contour[nosepad_flat_idx:(idx+1)], False) - (cutter_radius*2)
                if flat_length < 5:
                    heightmap = heightmap + [pt + [nosepad_removal]]
                else: # ramp down
                    transition_length = poly.polyline_length(contour[:(idx+1)], False)
                    ramp_segment = poly.segment(contour, transition_length, transition_length+5, False)
                    ramp_segment = poly.ramp(ramp_segment, nosepad_removal, body_removal,  False)
                    heightmap = heightmap + ramp_segment[:-1]
                    items_to_skip = len(ramp_segment)
                    nosepad_flat_idx = idx + items_to_skip # we'll need this to go down
                    over_nosepad = False
                    past_nosepad = True
            else:
                heightmap = heightmap + [pt + [body_removal]]
        return heightmap


    facing_contour = add_hinge_heights(facing_contour)
    facing_contour = add_nosepad_heights(facing_contour)
    facing_contour = poly.reverse(facing_contour)
    right_facing = poly.mirror_y(facing_contour, True)

    passes = [1]
    heights = [p[2] for p in facing_contour]
    r = [
        cam.change_tool("1/4in ballmill"),
        cam.spindle_speed(22000),
        cam.feedrate(1000),
        cam.start_spindle(),
        cam.rmp(entry_point),
        cam.contour(facing_contour, True),
    ]

    for dilate in passes:
        dilated = poly.reverse(poly.dilate(dilate, facing_contour))
#        dilated = add_hinge_heights(dilated)
        dilated = add_nosepad_heights(dilated)
        r = r + [ cam.contour(dilated, True),]
    return r

def add_tabs(contour, base_height, tab_height, tab_width, count):
    for i, p in enumerate(contour):
        contour[i] = [p[0], p[1], base_height]
    length = poly.polyline_length(contour, True)
    interval = length/(count+1)

    def create_tab(tab_seg, height):
        total_dist = 0
        half = poly.polyline_length(tab_seg, False)/2
        for i in range(len(tab_seg)-2):
            p1 = tab_seg[i]
            p2 = tab_seg[i+1]
            d = poly.line_length([p1, p2])
            total_dist += d
            if total_dist < half:
                adj_height = 0.5*total_dist # 30 degree slope
                p2[2] = min(p2[2]+height, p2[2]+adj_height)
            else:
                adj_height = -0.5*(2*half-total_dist)
                p2[2] = min(p2[2] + height, p2[2] - adj_height)


    next_tab_start = interval
    current_dist = 0
    tabs = []
    for i in range(len(contour)-2):
       p1 = contour[i]
       p2 = contour[i+1]
       d = poly.line_length([p1, p2])
       if current_dist < next_tab_start and (current_dist+d) > next_tab_start:
           tabs.append([i])
       elif current_dist > next_tab_start and (current_dist+d) > (next_tab_start+tab_width):
           tabs[-1].append(i+1)
           next_tab_start = next_tab_start + interval
       current_dist += d
    for tab in tabs:
        tab_seg = contour[tab[0]:tab[1]]
        create_tab(tab_seg, tab_height-base_height)

def ramped_outline(polyline, top, bottom, lead_length):
    return poly.ramp(poly.segment(polyline, poly.polyline_length(polyline) - lead_length, poly.polyline_length(polyline), True), top, bottom, False)


def mill_temples(outdir, temples, thickness=4):
#TODO: Replace with information in materials database
    total_thickness = 4
    thin_front = 0
    thin_back = 0

    r_temple = temples['right_temple_contour']
    l_temple = temples['left_temple_contour']

#    def tabbed_outline(polyline, top, bottom, gap):
        # Try putting 6 equally spaced tabs
#        length = poly.polyline_length(polyline, True)
#        interval = length/7
#        tab_lengths = [interval*(i+1) for i in range(6)]
#        height = poly.right(polyline) - poly.left(polyline)
#        width = poly.top(polyline) - poly.bottom(polyline)
#        y_tabs = [poly.bottom(polyline) + width/3.0, poly.bottom(polyline) + 2.0*width/3.0]
#        x_tabs = [poly.left(polyline) + height/2.0]
#        tab_lengths = poly.intercepts(polyline, y=y_tabs, x=x_tabs, closed=True)
#        return poly.make_gaps(polyline, tab_lengths, gap, top, bottom, closed=True)



    #left_temple_rough = poly.dilate(3.175/2.0 + 0.1, l_temple)
    #right_temple_rough = poly.dilate(3.175/2.0 + 0.15, r_temple)

#    left_temple_rough = poly.dilate(3.175/2.0 + 0.2, l_temple)
#    right_temple_rough = poly.dilate(3.175/2.0 + 0.2, r_temple)
#
#    left_temple_ramp = ramped_outline(left_temple_rough, 0.0, -thickness+0.5, 10.0)
#    right_temple_ramp = ramped_outline(right_temple_rough, 0.0, -thickness+0.5, 10.0)
#
#    left_temple_finish = poly.dilate(3.175/2.0, l_temple)
#    right_temple_finish = poly.dilate(3.175/2.0, r_temple)
#
#    left_temple_finish = [p + [-thickness+0.5] for p in left_temple_finish]
#    right_temple_finish = [p + [-thickness+0.5] for p in right_temple_finish]
#
#    left_temple_tabs = poly.dilate(3.175/2.0, l_temple)
#    right_temple_tabs = poly.reverse(poly.dilate(3.175/2.0+0.05, r_temple))
#    add_tabs(left_temple_tabs, -thickness-1, -thickness+0.3, 10, 6)
#    add_tabs(right_temple_tabs, -thickness-1, -thickness+0.3, 10, 6)

#    left_temple_finish = tabbed_outline(poly.dilate(3.175/2.0, l_temple), -thickness+0.5, -thickness-0.1, 12)
#    right_temple_finish = tabbed_outline(poly.reverse(poly.dilate(3.175/2.0, r_temple)), -thickness+0.5, -thickness-0.1, 12)

    offset = frame_offset(l_temple)
    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.rapid([0,0]),
        cam.activate_pin("stock_clamp"),
        surface_front(thin_front),
        surface_back(thin_back),
        index_holes([l_temple, r_temple], total_thickness),
        cam.change_tool("1/16in endmill"),
        cam.rapid([0,0]),
#        cam.temporary_offset(offset), # Not needed if we center on stock and have G55 at center
        temple_hinge_pockets(temples, thin_back),

        #thin_temples([l_temple, r_temple], temple_length),

        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        #cam.rapid([l_temple[0][0],l_temple[0][1], -thin_back] ),
        #cam.contour(l_temple, True),
        #cam.rapid([r_temple[0][0],r_temple[0][1], -thin_back] ),
        #cam.contour(r_temple, True),
        cam.change_tool("1/8in endmill"),
        cam.end_program()
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

    # offset the path so our 1/2 mill cuts the whole thing
    # TODO: gauge width and make sure we're cutting the whole thing.
    flat_begin = left_side + temple_length -10 # 20 mm flat
    flat_end = flat_begin + 20
    front_slope = -2.0 / (halfway-flat_begin)

    def calc_thinning_z(pt):
        if pt[1] > flat_begin:
            return (abs(pt[1]-halfway) * front_slope)
        elif pt[1] > flat_end:
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
    if amount < 0.1:
        print 'Not surfacing front, returning'
        return None
    print "Surfacing front with amount", amount
    surface_amount = min(amount, 2)
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
    surface_amount = min(amount, 2)
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
    xposition = hinge_height;
    yposition = temple_position + 4;
    left_hinge = hinges.get_hinge(hinge_num)
    right_hinge = hinges.get_hinge(hinge_num, False)
    left_translate = [xposition, yposition]
    right_translate = [xposition, -yposition]
    #right_translate = [xposition, 0]
    # Adjust by pocket depth of hinge
    #pocket_depth = left_hinge['pocket_depth']+thin_back

    pocket_depth = 1 + thin_back
    drill_depth = -thin_back - 2.0

    left_contour = poly.rotate_90(left_hinge["face_contour"])
    right_contour = poly.rotate_90(right_hinge["face_contour"])
    left_holes = poly.rotate_90(left_hinge["face_holes"])
    right_holes = poly.rotate_90(right_hinge["face_holes"])
    left_contour = poly.translate(left_contour, left_translate[0], left_translate[1])
    right_contour = poly.translate(right_contour, right_translate[0], right_translate[1])
    left_holes = poly.translate(left_holes, left_translate[0], left_translate[1])
    right_holes = poly.translate(right_holes, right_translate[0], right_translate[1])

    # Now center everything on the stock
    left_contour = poly.translate(left_contour, centering_shift[0], centering_shift[1])
    right_contour = poly.translate(right_contour, centering_shift[0], centering_shift[1])
    left_holes = poly.translate(left_holes, centering_shift[0], centering_shift[1])
    right_holes = poly.translate(right_holes, centering_shift[0], centering_shift[1])

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
        cam.start_spindle(4500),
        cam.dwell(2),
        [cam.rmp(p + [drill_depth], retract=10.0) for p in right_holes],
        [cam.rmp(p + [drill_depth], retract=10.0) for p in left_holes],
        cam.rapid([None, None, 20.0]),
    ]
    return r

def temple_hinge_pockets(temples, thinned):
    # We're operating in a 90 degree rotated fixture
    #l_hinge = poly.rotate_90(temples["left_hinge_contour"])
    #r_hinge = poly.rotate_90(temples["right_hinge_contour"])

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
        while len(erode) > 0:
            if len(contours) > 2:
                break
            elif len(erode) >= 1:
                contours.append(erode[0])
            else:
                break
#                for path in erode:
#                    contours.append(path)
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
    rhole = rhole[0]
#    polyline = dxf.polyline(layer="OUTLINE")
#    polyline.add_vertices(lhole)
#    drawing.add(polyline)


    right_rough = poly.erode((tool_radius + 0.3)/2, right_c)[0]
    left_rough = poly.erode((tool_radius+0.3)/2, left_c)[0]

    #lgroove = poly.erode(0.8, left_c)[0]
    #rgroove = poly.erode(0.8, right_c)[0]

    left_entry = poly.erode(5.0, left_c)[0][0];
    right_entry = poly.erode(5.0, right_c)[0][0];

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

def index_holes(contours, thickness):
# We put the index holes 1/2 between top and bottom, 160mm apart
#    log(str(poly.right(face_c)))
#    log(str(poly.left(face_c)))
#    log(str(poly.right(face_c) - poly.left(face_c)))
#    log(str((poly.right(face_c) - poly.left(face_c))/2))
#    log(str(poly.right(face_c) - (poly.right(face_c) - poly.left(face_c))/2))

# EXPERIMENTAL:  Just mill the index holes on the 0 line.  It's the responsibility
# of the milling programs to make sure the contour is in the correct location
# with respect to the index holes.
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
# END experimental code

# START pre-experiment code
    rightmost = -1000000
    leftmost = 1000000
    for contour in contours:
        right = poly.right(contour)
        left = poly.left(contour)
        rightmost = max(right, rightmost)
        leftmost = min(left, leftmost)

#    x_offset = poly.right(face_c) - (poly.right(face_c) - poly.left(face_c))/2
    x_offset = rightmost - (rightmost - leftmost)/2

    r_hole = [x_offset, 90]
    l_hole = [x_offset, -90]

    r = [
        cam.comment("Index Holes for secondary operations"),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(20000),
        cam.feedrate(1000),
        cam.dwell(3),
        cam.rmh(r_hole + [-thickness - 1.0], helix_radius, 0.5, 1),
        cam.rmh(l_hole + [-thickness - 1.0], helix_radius, 0.5, 1),
    ]
    return r


def nose_pads(order, thickness):
    return []


def arrange_temple_curves(left_temple_contour, right_temple_contour, hinge, lhinge_y, rhinge_y):

    left_hinge = hinges.get_hinge(hinge)
    right_hinge = hinges.get_hinge(hinge, False)

    left_holes = left_hinge['temple_holes']
    right_holes = right_hinge['temple_holes'] # opposite direction.  FIX someday
    right_hinge_contour = right_hinge['temple_contour']
    left_hinge_contour = left_hinge['temple_contour']
    print 'right hinge endpoints', right_hinge_contour[0], right_hinge_contour[-1], right_hinge_contour[-2]

    # y offset is mucked up because it's calculated from the center of the hinge boxing
    # square, not the center of the straight edge
    highpoint = poly.top(left_hinge_contour);
    lowpoint = poly.bottom(left_hinge_contour);
    y_offset = lhinge_y-(lowpoint + (highpoint-lowpoint)/2.0) # Distance to move hinge to center it on temple

    # Since we're tilted 6 degrees compared to design space we need to
    # compensate on the x axis
    x_offset = y_offset * math.tan(math.radians(6)) # This adjusts the hinge for the move to the center of the temple

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
    # Get the thing as horizontal as possible
    # When top-bottom distance is minimized, it's horizontal
    height = poly.right(left_temple_contour) - poly.left(left_temple_contour)
    opt_angle = 0

    for angle in range(2, 41):
        candidate_contour = poly.rotate(left_temple_contour, -angle)
        candidate_height = poly.right(candidate_contour)-poly.left(candidate_contour)
        if candidate_height < height:
            height = candidate_height
            opt_angle = angle
        else:
            break
#    print 'rotating temples by {0} degrees'.format(opt_angle)
    left_temple_contour = poly.rotate(left_temple_contour, -opt_angle);
    left_holes = poly.rotate(left_holes, -opt_angle)
    left_hinge_contour = poly.rotate(left_hinge_contour, -opt_angle)

    right_temple_contour = poly.rotate(right_temple_contour, opt_angle);
    right_holes = poly.rotate(right_holes, opt_angle)
    right_hinge_contour = poly.rotate(right_hinge_contour, opt_angle)

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
            p = p + [p[0], p[1]]  # Close the polygon to avoid a cusp
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
