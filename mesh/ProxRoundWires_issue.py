#DL240520
# two round wires model in python.
# 2d geometric model.
# to be used MFEM to simulate electromagnetic proximity effect.
# Control the mesh size with using points mesh size arguments.
#
#define the following physical
#contour, conductor.

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("proxroundwires_py")

dia=0.05;		#wire diameter.
dist=0.15;      #distance center to center`.`



gmsh.model.occ.addCircle(-dist/2,0,0,dia, 1)
gmsh.model.occ.addCurveLoop([1], 1)
gmsh.model.occ.addPlaneSurface([1], 1)

gmsh.model.occ.addCircle(dist/2,0,0,dia, 2)
gmsh.model.occ.addCurveLoop([2], 2)
gmsh.model.occ.addPlaneSurface([2], 2)

gmsh.model.occ.addCircle(0,0,0,2*dist, 3)
gmsh.model.occ.addCurveLoop([3], 3)
gmsh.model.occ.addPlaneSurface([3, -1, -2], 3)


gmsh.model.occ.addPoint(0,-2*dist,0, 4)
gmsh.model.occ.addPoint(0,2*dist,0, 5)
gmsh.model.occ.addLine(4, 5, 6)







#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(2, [3], 3, name="air")
gmsh.model.addPhysicalGroup(2, [2], 2, name="conductorright")
gmsh.model.addPhysicalGroup(2, [1], 1, name="conductorleft")
gmsh.model.addPhysicalGroup(1, [6], 4, name="centervertical")


# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

#gmsh.model.mesh.refine()

# glvis can read mesh version 2.2#
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("proxroundwires.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

