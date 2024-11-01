#DL240520
# two round wires model in python.
# 3d geometric model.
# to be used MFEM proximity3d.cpp to simulate electromagnetic proximity effect.
#
#define the following physical
#contour, conductor.

import gmsh
import sys
import math

gmsh.initialize()

gmsh.model.add("proxrectroundwires3d_py")

wrad=0.05;		#wire radius.
dist=0.15;      #distance center to center`.`
lenght=0.5;

wside = 0.1;


#left conductor
gmsh.model.occ.addRectangle(-dist/2-wside/2, -wside/2, 0, wside, wside, 1, 0)
#gmsh.model.occ.addCurveLoop([1], 1)
#gmsh.model.occ.addPlaneSurface([1], 1)

#right conductor
gmsh.model.occ.addCircle(dist/2,0,0,wrad, 12)
gmsh.model.occ.addCurveLoop([12], 12)
gmsh.model.occ.addPlaneSurface([12], 12)

#domain limit
gmsh.model.occ.addCircle(0,0,0,2*dist, 13)
gmsh.model.occ.addCurveLoop([13], 13)
gmsh.model.occ.addPlaneSurface([13, -1, -12], 13)

#form the cable lenght
gmsh.model.occ.extrude([[2, 1], [2, 12], [2 ,13]], 0, 0, -lenght)
#gmsh.model.occ.revolve([[2, 1], [2, 12], [2 ,13]], 0, 0.3, 0, 0, 0.6, 0, math.pi/2 )

#gmsh.model.occ.addPoint(0, 0.3, 0, 0, 31)
#gmsh.model.occ.addPoint(0, 0.6, 0, 0, 32)
#gmsh.model.occ.addLine(31, 32, 33)






#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(3, [3], 1, name="air")
gmsh.model.addPhysicalGroup(3, [2], 2, name="conductorright")
gmsh.model.addPhysicalGroup(3, [1], 3, name="conductorleft")
gmsh.model.addPhysicalGroup(2, [13], 4, name="airsurffront")
gmsh.model.addPhysicalGroup(2, [22], 5, name="airsurfback")
gmsh.model.addPhysicalGroup(2, [1], 6, name="conductorleftsurffront")
gmsh.model.addPhysicalGroup(2, [18], 7, name="conductorleftsurfback")
gmsh.model.addPhysicalGroup(2, [12], 8, name="conductorrightsurffront")
gmsh.model.addPhysicalGroup(2, [20], 9, name="conductorrightsurfback")
gmsh.model.addPhysicalGroup(2, [21], 10, name="airsurfaround")






# We can then generate a 3D mesh...
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
#gmsh.option.setNumber("Mesh.CheckSurfaceNormalValidity", 1)
#gmsh.option.setNumber("Mesh.MinCircleNodes", 50)

gmsh.model.mesh.generate(3)
#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()



# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("proxrectroundwires3d.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

