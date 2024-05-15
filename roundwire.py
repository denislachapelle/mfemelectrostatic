#DL240511
# round wire model in python.
# 2d geometric model.
#to be used MFEM to simulate electromagnetic skin effect.
#Control the mesh size with using points mesh size arguments.
#
#define the following physical
#contour, conductor.

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("roundwire_py")

dia=1.0;		#wire diameter.
tms=0.5;	    #target mesh size.


gmsh.model.occ.addCircle(0,0,0,dia, 1)
gmsh.model.occ.addCurveLoop([1], 1)

gmsh.model.occ.addPlaneSurface([1], 1)




#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(1, [1], 1, name="contour")
gmsh.model.addPhysicalGroup(2, [1], 2, name="conductor")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

gmsh.model.mesh.refine()

# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("roundwire.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

