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

    (options, args) = parser.parse_args()

    if options.infile == None and options.url == None:
        log("Insufficient arguments.")
        parser.print_help()
        sys.exit(0)

    infile = open(options.infile) if options.infile else urllib.urlopen(options.url)
    o = json.loads(infile.read())

    # Create a dxf for cutting the outline on the laser cutter
    order_id = o.get("order_id") or "testorder"
    order_dir = "/Users/rodfrey/Dropbox/guild_orders/" + order_id
    if not os.path.exists(order_dir):
        os.makedirs(order_dir)

    create_dxf(order_dir + "/face_contour.dxf", [o['face_con']])

    # Create the milling program for the lens holes/groove, hinge pockets, etc.
    mill_fronts(order_dir, o)

    temples =  arrange_temple_curves(o['ltemple_con'], o['rtemple_con'], o['lhinge'], o['lhinge_y'], o['rhinge_y'])
    create_dxf(order_dir + "/temple_contour.dxf",
            [temples['left_temple_contour'],
            temples['right_temple_contour']])
    mill_temples(order_dir, temples);

    #create_dxf("/Users/rodfrey/Dropbox/temple_contour.dxf", [temples['left_temple_contour'], temples['right_temple_contour']])
    #mill_temples(outdir, temples, o['temple_length'])


def mill_fronts(outdir, order):
    """Creates the g-code for the first milling operation.  The
    first milling operation is done on an unregistered plastic blank,
    so includes creating registration holes for subsequent operations."""

#TODO: Replace with information in materials database
# Initial thickness of the forward-facing lamination.  0 if no lamination.
    front_surface_thickness = 0

# The final desired thickness of the fronto facing lamination.  Must
# be equal to or less than the front_surface_thickness.
    final_front_thickness = 0

# Initial thickness of the face side lamination. Use this thickness only if
# stock is solid
    back_surface_thickness = 5

# The final thickness of the left and right extremes of the frame where the
# hinge will go.  Must be equal to or less than back_surface_thickness.
    hinge_thickness = 4

# The thickness of the highest point of the nosepad
    nosepad_thickness = 4

# Final thickness of the main part of the frame.  Must be equal to or less than
# the back_surface_thickness.
    final_back_thickness =5

    thickness = final_front_thickness + final_back_thickness

    front_surface_removal = final_front_thickness - front_surface_thickness
    back_surface_removal = final_back_thickness - back_surface_thickness


# The machine has the stock clamp oriented 90 degrees to the way the
# software creates the contours.
    face_c = poly.rotate_90(order["face_con"])
    left_lens_c = poly.rotate_90(order["lhole_con"])
    right_lens_c = poly.rotate_90(order["rhole_con"])
    temple_height = abs(order["ltemple_con"][0][1] - order["ltemple_con"][-1][1])

    msg = check_frame_size(face_c)
    if msg:
        print msg
        sys.exit(0)

    offset = frame_offset(face_c)

    groove_height = back_surface_removal + (thickness/2)
    size_info = order['usersizes']

    #TEMPORARY
    size_info['nose_radius'] = 6
    size_info['nose_height'] = 0
    size_info['nose_splayangle'] = 20
    size_info['nose_ridgeangle'] = 23
    program = [
        cam.setup(),
        cam.select_fixture("blank_clamp"),
        cam.retract_spindle(),
        cam.activate_pin("stock_clamp"),
        surface_front(front_surface_removal),
#        surface_back(back_surface_removal),

        cam.change_tool("1/8in endmill"),
        cam.rapid([0, 0]),
        cam.temporary_offset(offset),
# Note that X and Y parameters in the order are switched from our system
#TODO: replace the thickness offset with the thickness of the TEMPLE, not the fronts.
#        index_holes([face_c], back_surface_thickness),
        lens_holes(left_lens_c, right_lens_c, back_surface_thickness),
        lens_groove(left_lens_c, right_lens_c, back_surface_removal - (thickness/2.0)),
#        contour_face(
#            back_surface_removal,
#            back_surface_thickness - hinge_thickness,
#            back_surface_thickness - nosepad_thickness,
#            temple_height,
#            face_c, left_lens_c, order['lhinge_y']),
        face_hinge_pockets(order["lhinge"], order["lhinge_y"], order["ltemple_x"]),
#        nose_pads(order, thickness),
        nose_contour(size_info["nose_radius"], size_info["nose_height"], size_info["nose_splayangle"], size_info["nose_ridgeangle"], face_c, back_surface_thickness),
#        cam.contour(face_c, True),
        outline(face_c),

        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
        cam.end_program(),
    ]
    open(outdir + "/face_stage1.ngc", "w").write(to_string(program))

def outline(face_con):
    """Creates the face outline feature toolpath.  face_con is the polygon representing the face outline."""

    #rough = poly.reverse(poly.dilate(3.175/2.0 + 0.25, face_con)) # Reversed to have climb milling
    rough =  poly.dilate(6.35/4.0 + 0.2, face_con)
    #finish = poly.reverse(poly.dilate(6.35/2.0, face_con)) # Reversed to have climb milling
    finish = poly.dilate(6.35/4.0, face_con)

    ramp = poly.ramp(poly.segment(rough, poly.polyline_length(rough, True) - 50.0, poly.polyline_length(rough), True), 0.0, -3.5, False)

    finishramp = poly.ramp(poly.segment(finish, poly.polyline_length(finish) - 50.0, poly.polyline_length(finish), True), 0.0, -4.5, False)
    tab_lengths = poly.intercepts(finish, True, y=[-39.0, 39.0])
#    tab_lengths = poly.intercepts(finish, y=[0])
    finish_and_tabs = poly.make_gaps(finish, tab_lengths, 10, -4.2, -5.2, True)

    r = [
        "(Outline Rough)",
        cam.feedrate(1500),
        cam.change_tool("1/8in endmill"),
        cam.start_spindle(22000),

        # Ramp in
        cam.rmp(ramp[0]),
        cam.contour(ramp, False),
        cam.contour(rough, True),

        "(Outline Finish and Tabs)",
        cam.feedrate(1500),
        #cam.change_tool(1), # 4 flute 1/8" flat
        #cam.start_spindle(4500),

#        rmp(finishramp[0]),
        ["G1 Z-5"],
#        cam.contour(finish, False),
        cam.contour(finish_and_tabs, True),
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






def mill_temples(outdir, temples):
#TODO: Replace with information in materials database
    front_surface_thickness = 0
    back_surface_thickness = 4
    final_front_thickness = 0
    final_back_thickness = 4
    thickness = final_front_thickness + final_back_thickness

    front_surface_removal = final_front_thickness - front_surface_thickness
    back_surface_removal = final_back_thickness - back_surface_thickness

    r_temple = temples['right_temple_contour']
    l_temple = temples['left_temple_contour']
    def ramped_outline(polyline, top, bottom, lead_length):
        return poly.ramp(poly.segment(polyline, poly.polyline_length(polyline) - lead_length, poly.polyline_length(polyline), True), top, bottom, False)

    def tabbed_outline(polyline, top, bottom, gap):
        height = poly.right(polyline) - poly.left(polyline)
        width = poly.top(polyline) - poly.bottom(polyline)
        y_tabs = [poly.bottom(polyline) + width/3.0, poly.bottom(polyline) + 2.0*width/3.0]
        x_tabs = [poly.left(polyline) + height/2.0]
        tab_lengths = poly.intercepts(polyline, y=y_tabs, x=x_tabs, closed=True)
        return poly.make_gaps(polyline, tab_lengths, gap, top, bottom, closed=True)



    left_temple_rough = poly.dilate(3.175/2.0 + 0.1, l_temple)
    right_temple_rough = poly.dilate(3.175/2.0 + 0.15, r_temple)

    left_temple_ramp = ramped_outline(left_temple_rough, 0.0, -3.5, 10.0)
    right_temple_ramp = ramped_outline(right_temple_rough, 0.0, -3.5, 10.0)

    left_temple_finish_1 = poly.dilate(3.175/2.0, l_temple)
    right_temple_finish_1 = poly.reverse(poly.dilate(3.175/2.0, r_temple))

    left_temple_finish = tabbed_outline(poly.dilate(3.175/2.0, l_temple), -3.5, -4.5, 6)
    right_temple_finish = tabbed_outline(poly.reverse(poly.dilate(3.175/2.0, r_temple)), -3.5, -4.5, 6)

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
        cam.change_tool("1/8in endmill"),
        #index_holes([l_temple, r_temple], thickness),
        #thin_temples([l_temple, r_temple], temple_length),

        cam.feedrate(1200),
        cam.start_spindle(22000),
        cam.rmp(left_temple_ramp[0]),
        cam.contour(left_temple_ramp, False),
        cam.contour(left_temple_rough, True),
        cam.contour(left_temple_finish_1, True),
        cam.feedrate(700),
        cam.contour(left_temple_finish, True),

        cam.feedrate(1000),
        cam.rmp(right_temple_ramp[0]),
        cam.contour(right_temple_ramp, False),
        cam.contour(right_temple_rough, True),
        cam.contour(right_temple_finish_1, True),
        cam.feedrate(700),
        cam.contour(right_temple_finish, True),

#        outline(r_temple),
#        outline(l_temple),
        cam.retract_spindle(),
        cam.deactivate_pin("stock_clamp"),
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




def face_hinge_pockets(hinge_num, hinge_height, temple_position):
    xposition = hinge_height;
    yposition = temple_position + 4;
    left_hinge = hinges.get_hinge(hinge_num)
    right_hinge = hinges.get_hinge(hinge_num, False)
    left_translate = [xposition, yposition]
    #left_translate = [xposition, 0]
    right_translate = [xposition, -yposition]
    #right_translate = [xposition, 0]
    # Adjust by pocket depth of hinge
    pocket_depth = left_hinge['pocket_depth']

    left_contour = poly.rotate_90(left_hinge["face_contour"])
    right_contour = poly.rotate_90(right_hinge["face_contour"])
    left_holes = poly.rotate_90(left_hinge["face_holes"])
    right_holes = poly.rotate_90(right_hinge["face_holes"])
    left_contour = poly.translate(left_contour, left_translate[0], left_translate[1])
    right_contour = poly.translate(right_contour, right_translate[0], right_translate[1])
    left_holes = poly.translate(left_holes, left_translate[0], left_translate[1])
    right_holes = poly.translate(right_holes, right_translate[0], right_translate[1])


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
        [cam.rmp(p + [-3.5], retract=10.0) for p in right_holes],
        [cam.rmp(p + [-3.5], retract=10.0) for p in left_holes],
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

    def pocket_contours(contour):
        contours = []
        erode = poly.erode(1.5875/2, contour)
        while len(erode) > 0:
            if len(erode) == 1:
                contours.append(erode[0])
            else:
                for path in erode:
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
#    print "Length of left hinge pocket contours", len(left_hinge_pocket_contours)
#    right_hinge_pocket_contours = [];
#    while len(r_hinge) > 0:
#            r_hinge = poly.erode(1.5875/2, r_hinge)
#            if len(rhinge_) == 1:
#                right_hinge_pocket
#            if len(r_hinge) > 0:
#                r_hinge = r_hinge[0]
#                right_hinge_pocket_contours.append(r_hinge)
#    print "Length of right hinge pocket contours", len(right_hinge_pocket_contours)
#    print right_hinge_pocket_contours
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

        #cam.move([None, None, 0]),
        #cam.contour(poly.rotate_90(temples['left_temple_contour']), True),
        #cam.contour(poly.rotate_90(temples['right_temple_contour']), True),
    ]
    return r



def lens_groove(left_c, right_c, height):
    """Generates the toolpath for the lens holes (holes, groove and tabs)."""
    print 'GROOVE HEIGHT', height
    if not poly.is_ccw(left_c):
        left_c = poly.reverse(left_c)
    if not poly.is_ccw(right_c):
        right_c = poly.reverse(right_c)
    print 'offsetting lens contours'
    lgroove = poly.erode(0.8, left_c)[0]
    rgroove = poly.erode(0.8, right_c)[0]

    print 'offsetting lens countours for entry'
    left_entry = poly.erode(7.0, left_c)[0][0];
    right_entry = poly.erode(7.0, right_c)[0][0];
    r = [
        "(Lens Grooves)",
        cam.change_tool("vgroove"),
        cam.start_spindle(20000),
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
    if not poly.is_ccw(left_c):
        left_c = poly.reverse(left_c)
    if not poly.is_ccw(right_c):
        right_c = poly.reverse(right_c)

    drawing = dxf.drawing('test.dxf')
    drawing.add_layer('OUTLINE', color=1)
    polyline = dxf.polyline(layer="OUTLINE")
    polyline.add_vertices(left_c)
    drawing.add(polyline)


    lhole = poly.erode(3.175/2.0, left_c)[0]
    rhole = poly.erode(3.175/2.0, right_c)[0]
    polyline = dxf.polyline(layer="OUTLINE")
    polyline.add_vertices(lhole)
    drawing.add(polyline)


    right_rough = poly.erode(0.1, rhole)[0]
    left_rough = poly.erode(0.1, lhole)[0]

    lgroove = poly.erode(0.8, left_c)[0]
    rgroove = poly.erode(0.8, right_c)[0]

    left_entry = poly.erode(5.0, left_c)[0][0];
    right_entry = poly.erode(5.0, right_c)[0][0];

    lhole = poly.reverse(lhole)
    rhole = poly.reverse(rhole)

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
    # Find the lower height of the bridge
    h = face_con[len(face_con)/2][0] + 1
    print 'face contour lower point', h
    print 'face contour upper point', face_con[0][0]
    # Offset so that the ridge slope ends somewhere around the front
    ridge_offset = thickness/math.tan(math.radians(90 - nose_ra))
    print 'ridge offset', ridge_offset, nose_ra
    h = h + ridge_offset + nose_h + 3.175
    print 'final height', h
    sa = math.radians(nose_sa * 2)
    ra = math.radians(nose_ra)
    xfloor = poly.left(face_con) - 3.175  # bottom most point minus tool radius
    xfloor = max(xfloor, -27.0) # miminum safe distance without hitting clamp
    nose_tool_radius = 3.175

    nextpoly = nose.nose_poly(nr, h, sa, ra, xfloor, nose_tool_radius, 0.0)

    r = [
        "(Nose Contour)",
        cam.change_tool("1/4in ballmill"),
        cam.start_spindle(20000),
        cam.feedrate(2000),
        cam.rmp(nextpoly[0] + [2.0])  # Start near our first contour
    ]

    direction = 1
    for i in range(-5, (thickness+2)*10 +5):
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

    left_holes = left_hinge['temple_holes']
    right_holes = right_hinge['temple_holes'] # opposite direction.  FIX someday
    right_hinge_contour = right_hinge['temple_contour']
    left_hinge_contour = left_hinge['temple_contour']

    # y offset is mucked up because it's calculated from the center of the hinge boxing
    # square, not the center of the straight edge
    highpoint = poly.top(left_hinge_contour);
    lowpoint = poly.bottom(left_hinge_contour);
    y_offset = lhinge_y-(lowpoint + (highpoint-lowpoint)/2.0)
# Since we're tilted 6 degrees compared to design space we need to
# compensate on the x axis
    x_offset = y_offset * math.tan(math.radians(6)) + 1.0
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

    left_temple_contour = poly.rotate(left_temple_contour, -opt_angle);
    left_holes = poly.rotate(left_holes, -opt_angle)
    left_hinge_contour = poly.rotate(left_hinge_contour, -opt_angle)

    right_temple_contour = poly.rotate(right_temple_contour, opt_angle);
    right_holes = poly.rotate(right_holes, opt_angle)
    right_hinge_contour = poly.rotate(right_hinge_contour, opt_angle)

    temple_width = poly.top(left_temple_contour) - poly.bottom(left_temple_contour)
    optimization_offset = (150-temple_width)/2.0
    left_temple_contour = poly.translate(left_temple_contour, 0, optimization_offset)
    left_holes = poly.translate(left_holes, 0, optimization_offset)
    left_hinge_contour = poly.translate(left_hinge_contour, 0, optimization_offset)
    right_temple_contour = poly.translate(right_temple_contour, 0, -optimization_offset)
    right_hinge_contour = poly.translate(right_hinge_contour, 0, -optimization_offset)
    right_holes = poly.translate(right_holes, 0, -optimization_offset)
    # Translate left one rightward as much as we can
    total_width = poly.left(right_temple_contour) - poly.right(left_temple_contour)
    print "Width before packing", total_width
    for i in range(1, 100):
        candidate = poly.translate(right_temple_contour, i, 0)
        intersection = poly.intersection(candidate, left_temple_contour)
        if len(intersection) > 0:
            right_temple_contour = poly.translate(right_temple_contour, i-15, 0)
            right_hinge_contour = poly.translate(right_hinge_contour, i-15, 0)
            right_holes = poly.translate(right_holes, i-15, 0)
            break

    # sanity check that we fit on stock
    total_width =  poly.right(left_temple_contour) - poly.left(right_temple_contour)
    print "Width after packing", total_width
    if total_width > 60:
        print 'Error! temples did not pack tight enough.'
        raise 'Sizing error'
    # Now they're packed togegther OK but not centred on the stock
    #Midpoint of curves should be 0, so that's how much we'll shift it
    centering_shift = poly.bottom(right_temple_contour) + (poly.top(left_temple_contour) - poly.bottom(right_temple_contour))/2.0

    left_temple_contour = poly.translate(left_temple_contour, 0, -centering_shift)
    right_temple_contour = poly.translate(right_temple_contour, 0, -centering_shift)
    left_hinge_contour = poly.translate(left_hinge_contour, 0, -centering_shift)
    right_hinge_contour = poly.translate(right_hinge_contour, 0, -centering_shift)
    left_holes = poly.translate(left_holes, 0, -centering_shift)
    right_holes = poly.translate(right_holes, 0, -centering_shift)
    centering_shift = (poly.top(left_temple_contour) - poly.bottom(right_temple_contour))/2.0
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
