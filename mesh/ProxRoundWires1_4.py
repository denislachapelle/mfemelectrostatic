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

gmsh.model.add("proxroundwires_1_4_py")

ddia=0.3;       #domain diameter.
wdia=0.05;		#wire diameter.
dist=0.12;       #distance center to center`.`



gmsh.model.occ.addPoint(0,0,0,1)
gmsh.model.occ.addPoint(dist-wdia/2,0,0,2)
gmsh.model.occ.addPoint(dist+wdia/2,0,0,3)
gmsh.model.occ.addPoint(ddia,0,0,4)
gmsh.model.occ.addPoint(0,ddia,0,5)
gmsh.model.occ.addPoint(dist, 0, 0, 12)
gmsh.model.occ.addPoint(dist, wdia/2, 0, 13)

gmsh.model.occ.addLine(1, 2, 6)
gmsh.model.occ.addLine(2, 3, 7)
gmsh.model.occ.addLine(3, 4, 8)
gmsh.model.occ.addLine(5, 1, 9)

gmsh.model.occ.addCircleArc(3, 6, 7, 10)
gmsh.model.occ.addCircleArc(7, 6, 2, 14)
gmsh.model.occ.addCircleArc(5, 1, 4, 11)

gmsh.model.occ.addCurveLoop([7, 10, 14], 1)
gmsh.model.occ.addPlaneSurface([1], 1)

gmsh.model.occ.addCurveLoop([6, 14, 10, 8, 11, 9], 2)
gmsh.model.occ.addPlaneSurface([2], 2)

gmsh.model.occ.extrude([[2, 1], [2, 2]], 0, 0, -0.5)





#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(3, [1], 1, name="conductor")
gmsh.model.addPhysicalGroup(2, [1], 2, name="conductorsurfacefront")
gmsh.model.addPhysicalGroup(2, [6], 3, name="conductorsurfaceback")

gmsh.model.addPhysicalGroup(3, [2], 4, name="air")
gmsh.model.addPhysicalGroup(2, [2], 5, name="airsurfacefront")
gmsh.model.addPhysicalGroup(2, [11], 7, name="airsurfaceback")

gmsh.model.addPhysicalGroup(2, [10], 8, name="verticalsurface")
gmsh.model.addPhysicalGroup(2, [7, 3, 8], 9, name="horizontalsurface")
gmsh.model.addPhysicalGroup(2, [9], 10, name="circularsurface")


# We can then generate a 3D mesh...
gmsh.model.mesh.generate(3)

gmsh.model.mesh.refine()

# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("proxroundwires1_4.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

