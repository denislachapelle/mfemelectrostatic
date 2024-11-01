#DL240520
# two round wires model in python.
# 3d geometric model.
# to be used MFEM proximity3d.cpp to simulate electromagnetic proximity effect.
#
#define the following physical
#contour, conductor.

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("proxroundwires3d_py")

wrad=0.05;		#wire radius.
dist=0.15;      #distance center to center`.`
lenght=0.5;


#left conductor
gmsh.model.occ.addCircle(-dist/2,0,0,wrad, 1)
gmsh.model.occ.addCurveLoop([1], 1)
gmsh.model.occ.addPlaneSurface([1], 1)

#right conductor
gmsh.model.occ.addCircle(dist/2,0,0,wrad, 2)
gmsh.model.occ.addCurveLoop([2], 2)
gmsh.model.occ.addPlaneSurface([2], 2)

#domain limit
gmsh.model.occ.addCircle(0,0,0,2*dist, 3)
gmsh.model.occ.addCurveLoop([3], 3)
gmsh.model.occ.addPlaneSurface([3, -1, -2], 3)

#form the cable lenght
gmsh.model.occ.extrude([[2, 1], [2, 2], [2 ,3]], 0, 0, -lenght)





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
gmsh.model.addPhysicalGroup(2, [3], 4, name="airsurffront")
gmsh.model.addPhysicalGroup(2, [9], 5, name="airsurfback")
gmsh.model.addPhysicalGroup(2, [1], 6, name="conductorleftsurffront")
gmsh.model.addPhysicalGroup(2, [5], 7, name="conductorleftsurfback")
gmsh.model.addPhysicalGroup(2, [2], 8, name="conductorrightsurffront")
gmsh.model.addPhysicalGroup(2, [7], 9, name="conductorrightsurfback")
gmsh.model.addPhysicalGroup(2, [8], 10, name="airsurfaround")






# We can then generate a 3D mesh...
gmsh.model.mesh.generate(3)

#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()


# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("proxroundwires3d.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

