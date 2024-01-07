#DL240106
#tst_ring model in python.
#2d geometric model used to test charge calculating software.
#to be used MFEM to simulate the electric field.
#Control the mesh size with using points mesh size arguments.
#
#define the following physical
#groundline, traceline, airsurface, dielectricsurface

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("tst_ring")

inner=99;		#
outter=100;		#


#
# Start by drawing the air points, line, loop, surface and physical.
#
gmsh.model.occ.addCircle(0, 0, 0, outter, 1)
gmsh.model.occ.addCurveLoop([1], 1)
gmsh.model.occ.addCircle(0, 0, 0, inner, 2)
gmsh.model.occ.addCurveLoop([2], 2)
gmsh.model.occ.addPlaneSurface([1, 2], 1)


#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(1, [2], 2, name="traceline")
gmsh.model.addPhysicalGroup(1, [1], 1, name="groundline")
gmsh.model.addPhysicalGroup(2, [1], 3, name="airsurface")
#gmsh.model.addPhysicalGroup(2, [2], 4, name="dielecsurface")
#gmsh.model.addPhysicalGroup(1, [13], 5, name="integcloseloop")


# We can then generate a 2D mesh...
gmsh.option.setNumber('Mesh.MeshSizeMax', 0.5)
gmsh.model.mesh.generate(2)

# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("tst_ring.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

