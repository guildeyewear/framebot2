Hinge Modelling Folder
======================

Repository of hinges that we can read and interpret.

Hinges are modelled using 2d CAD and saved as DXF.  create\_hinge.py is then run to create a json file that contains the polylines for the hinges. These json files are used by the website to render the hinge locations and by framebot to create the machine code for milling the hinge pockets.

The hinges are modelled using some location conventions.  These conventions are very important to ensure consistent rendering and modelling.

1. We model the **right** hinge.
2. The Y origin is at the middle of the **face** hinge, that is, the face hinge is split in half by the X axis.
3. The X axis for the **temple** hinge is calculated based on where the temple hinge is on the face hinge when assembled.
4. The Y axis (X=0) of the **face** hinge is the location of the exposed face of the temple hinge when the hinge is at 90 degrees, i.e. open. (This assumption might be revisited if we add hinges that are open at 180 degrees and closed at 90 degrees.)
5. The Y axis (X=0) of the **temple** hinge is the exposed face of the face hinge when the hinge is open, i.e. at 90 degrees. As a consequence of modelling the right hinge that means the entire temple hinge has its X coordinates < 0.
