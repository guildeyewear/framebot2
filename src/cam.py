"""Toolpath and gcode generation."""
TOOL_DEFINITIONS = {
    "1/16in endmill": 1,
    "1/8in endmill": 2,
    "1/4in endmill": 3,
    "dovetail": 4,
    "1/4in ballmill": 5,
    "vgroove" : 6,
    "1mm drill": 7,
    "3/8in endmill": 8,
    "engraver": 9,
    "empty10": 10,
    "tapered": 11,
};

FIXTURE_DEFINITIONS = {
    "default": "G54",
    "blank_clamp": "G55",
    "lens_clamp" : "G56"
};
PIN_DEFINITIONS = {
    "stock_clamp": "P03",
};

def xyz(point):
    # Format is (x, y) or (x, y, z).  (x, y) implies z is None
    x = point[0]
    y = point[1]
    z = None

    if len(point) == 3:
        z = point[2]
    r = ""
    if x is not None:
        r += " X%f" % x
    if y is not None:
        r += " Y%f" % y
    if z is not None:
        r += " Z%f" % z
    return r

def rapid(point):
    return ["G0" + xyz(point)]

def move(point):
    return ["G1" + xyz(point)]

def machine_rapid(point):
    return ["G0", "G53" + xyz(point)]

def machine_move(point):
    return ["G1", "G53" + xyz(point)]

def retract_spindle():
    return ["G0 G53 Z0"]

def feedrate(rate):
    return ["F%f" % rate]

def change_tool(n, retract=True, spindle_off=True):
    '''Change to tool n with an optional retract, included by default.'''

    if isinstance(n, basestring):
        n = TOOL_DEFINITIONS[n];

    return [
        #stop_spindle() if spindle_off else None,
        #machine_rapid((None, None, 0.0)) if retract else None,
        #dwell(5) if spindle_off else None,
        "T%d M6 G43 H%d" % (n, n)
        ]

def pause():
    return ["M0"]

def stop_spindle():
    return ["M5"]

def spindle_speed(rpm):
    return ["S%d" % rpm]

def start_spindle(rpm=None):
    if rpm != None:
        return spindle_speed(rpm) + start_spindle()
    else:
        return ["M3"]

def relative():
    return ["G91"]

def absolute():
    return ["G90"]

def xy_plane():
    return ["G17"]

def mm():
    return ["G21"]

def flip_stock():
    return [
        deactivate_pin("stock_clamp"),
        pause(),
        activate_pin("stock_clamp"),
    ]

def surface_along_y(start_x, start_y, end_x, end_y, tool_radius, depth):
    # Adjust start points for tool radius
    if start_y < end_y:
        tool_radius = abs(tool_radius)
    else:
        tool_radius = -(abs(tool_radius))
    #start_y = start_y - tool_radius
    #end_y = end_y + tool_radius

    direction = 1
    if start_x < end_x:
        tool_radius = abs(tool_radius)
    else:
        tool_radius = -(abs(tool_radius))
        direction = -1
    start_x = start_x + tool_radius
    end_x = end_x - tool_radius
    x_stepover = tool_radius * 1.5   # Sign of tool_radius depends on start_x compare just before

    start = [start_x, start_y]
    current_loc = [start[0], start[1]]

    path = [[current_loc[0], current_loc[1]],
            [current_loc[0], end_y]]
    current_loc[1] = end_y

    while (current_loc[0] + tool_radius)*direction  < end_x*direction:
        # Step over
        current_loc[0] = current_loc[0] + x_stepover
        path = path + [[current_loc[0], current_loc[1]]]

        # Cutting move
        if current_loc[1] == start_y:
            current_loc[1] = end_y
        else:
            current_loc[1] = start_y
        path = path + [[current_loc[0], current_loc[1]]]

    # Final stepover to finish
    current_loc[0] = end_x
#    path = path + [[current_loc[0], current_loc[1]]]
    if current_loc[1] == start_y:
        current_loc[1] = end_y
    else:
        current_loc[1] = start_y
#    path = path + [[current_loc[0], current_loc[1]]]
    return [
        rmp(start + [depth], retract=50),
        contour(path, False),
    ]

def cancel_offsets():
    return ["G92.1"]

def arc(point, i, j):
    return ["G3 X%f Y%f Z%f I%f J%F" % (point[0], point[1], point[2], i, j)]

def comment(s):
    return ["(%s)" % s]

def dwell(seconds):
    return ["G4 P%f" % seconds]

def rmp(p, start_height=10.0, retract=20.0):
    """Retract-move-plunge (RMP) to given point."""
    r = [
        rapid([None, None, retract]),
        rapid([p[0], p[1]]),
        rapid([None, None, start_height]),
        move([None, None, p[2]]),
    ]
    return r

def rmh(p, radius, start_height=1.0, pitch=0.1, retract=20.0):
    """Retract-move-helix (RMH) to given point."""
    r = [
        rapid([None, None, retract]),
        rapid([p[0], p[1]]),
        rapid([None, None, start_height]),
        helix(radius, pitch, [p[0], p[1], start_height], p[2]),
    ]
    return r


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

def select_fixture(fixture):
    code = None
    if isinstance(fixture, basestring):
        code = FIXTURE_DEFINITIONS[fixture]
    return [code]

def activate_pin(pinname):
    code = None
    if isinstance(pinname, basestring):
        code = "M64 %s" % PIN_DEFINITIONS[pinname]
    return [code]

def deactivate_pin(pinname):
    code = None
    if isinstance(pinname, basestring):
        code = "M65 %s" % PIN_DEFINITIONS[pinname]
    return [code]


def to_string(l):
    '''Converts a list of strings or a tree structure of strings to a single string for output.'''
    return "\n".join(flatten(l)) + "\n"  # Mach3 does not process the last line if it does not have a newline at the end

def setup():
    '''Standard preamble for most milling programs.'''
    return [
        xy_plane(),
        mm(),
        absolute(),
        cancel_offsets(),
        ["G64 Q0.02"],  # LinuxCNC trajectory planner parameter
    ]

def end_program():
    return [
        remove_temporary_offset(),
        retract_spindle(),
        select_fixture("default"),
        ["M30"], # End program
            ]

def helix(radius, max_pitch, start_point, end_z):
    # move from start in +x direction by radius of helix (at radius of helix)
    # do full turns at max_pitch until you canna do full turns no more
    # do a turn that has a partial z move to the end_height
    # move from position in -x direction by radius of helix (back to centre)

    max_pitch = -abs(max_pitch)  # We're always going down baby

    start_x = start_point[0]
    start_y = start_point[1]
    start_z = start_point[2]

    depth = start_z - end_z
    full_turns = abs(int(depth/float(max_pitch)))

    r = [
        move((start_x+radius, None)),
        [arc((start_x+radius, start_y, start_z+(max_pitch*n)), -radius, 0) for n in range(1, full_turns+1)],
        arc((start_x+radius, start_y, end_z), -radius, 0),
        move((start_x, None))
        ]

    return r

def contour(points, closed, cw=False):
    step = {False: 1, True: -1}[cw]

    r = [move(p) for p in points[::step]]

    if closed:
        r += move(points[0]) # have to move back to start

    return r

# Assumes cutter is in X
def polar_contour(points):
   return [
        ["G1 A%f F720" % points[0][0]],
        ["G1 A%f X%f" % (p[0], p[1]) for p in points],
           ]

def pocket(contours, bottom, retract=5.0):
    '''Poor man's pocket routine.  Runs a series of contours at depth to create a pocket.'''

    r = []

    for c in contours:
        r += rapid(c[0])  # Contour XY
        r += rapid([None, None, retract]) # Retract plane
        r += move([None, None, bottom]) # Plunge
        r += contour(c, True, cw=True) # Run contour
        r += rapid([None, None, retract]) # Retract plane

    return r

def work_offset(n):
    '''G59, Switches to work offset n.'''

    if n not in range(1, 255):
        raise Exception("Invalid work offset, must be from 1 to 254.")

    return ["G59 P%d" % n]

def temporary_offset(point):
    ''' This would best be done with G52, but that's not supported on Linuxcnc.
    Instead we use G92, which temporarily overwrites the current location of the
    named axis with the value passed in.  So if we're at the origin in G55, and
    pass temporary_offset([10, 10]), Linuxcnc will now think that the current
    location is 10,10, effectively shifting the origin by -10, -10.
    '''
    return ["G92 " + xyz(point)]

def remove_temporary_offset():
    return ["G92.1"]
