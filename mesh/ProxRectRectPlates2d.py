#DL241009
# two rect plates model in python.
# 2d geometric model.
# to be used MFEM proximity2d.cpp to simulate electromagnetic proximity effect.
#
#define the following physical
#contour, conductor.

import gmsh
import sys
import math

gmsh.initialize()

gmsh.model.add("proxrectrectplates2d_py")

skindept=0.0085         #skin dept in meter for copper at 60Hz.
width=20*skindept		#plates width.
thickness=4*skindept   #plate thickness.
dist=width*1.1



#left conductor
gmsh.model.occ.addRectangle(-dist/2-width/2, -thickness/2, 0, width, thickness, 1, thickness/4)

#right conductor
gmsh.model.occ.addRectangle(+dist/2-width/2, -thickness/2, 0, width, thickness, 2, thickness/4)


#domain limit
#domain limit
gmsh.model.occ.addCircle(0,0,0, 3*width, 100)
gmsh.model.occ.addCurveLoop([100], 100)
gmsh.model.occ.addPlaneSurface([100, -1, -2], 100)
#gmsh.model.occ.addRectangle(-3*(dist+width)/2, -5*thickness/2, 0, 3*(dist+width), 5*thickness, 3, 0)

#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(2, [2], 1, name="wire_1")
gmsh.model.addPhysicalGroup(2, [1], 2, name="wire_2")
gmsh.model.addPhysicalGroup(2, [100], 3, name="air")
gmsh.model.addPhysicalGroup(1, [100], 4, name="aircontour")



# We can then generate a 2D mesh...
gmsh.option.setNumber("Mesh.Algorithm", 6)

gmsh.model.mesh.generate(2)
gmsh.model.mesh.refine()
gmsh.model.mesh.refine()



# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("proxrectrectplates2d.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

