<<<<<<< HEAD
#DL231219
#microstrip model in python.
#microstrip 2d geometric model.
#to be used MFEM to simulate the electric field.
#Control the mesh size with using points mesh size arguments.
#
#define the following physical
#groundline, traceline, airsurface, dielectricsurface

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("microstrip_py")

th=0.1;		#trace height.
tw=10.0;		#trace width.
dh=5;		#dielectric heigh.
gw=20;		#ground width.
hspace=20;
rad=0.01;	#radius corner.
tms=0.5;	#target mesh size.


#
# Start by drawing the air points, line, loop, surface and physical.
#

Points=[    [-gw/2, 0, 0, tms, 1],     #air bottom-left
            [gw/2, 0, 0, tms,2],       #air bottom-right
            [gw/2, hspace, 0, 10*tms, 3], #air top-right
            [-gw/2, hspace, 0, 10*tms, 4], #air top-left
            [-tw/2, 0, 0, tms, 5],      #trace bottom-left
            [-tw/2, th, 0, tms, 6],     #trace top-left
            [tw/2, th, 0, tms, 7],      #trace top-right
            [tw/2, 0, 0, tms, 8]
            ]



for Point in Points:
   gmsh.model.occ.addPoint(Point[0], Point[1], Point[2], Point[3], Point[4])
   
   

Lines=[    [1, 5, 1],
            [5, 6, 2],
            [6, 7, 3],
            [7, 8, 4],
            [8, 2, 5],
            [2, 3, 6],
            [3, 4, 7],
            [4, 1, 8]
            ]

for Line in Lines:
   gmsh.model.occ.addLine(Line[0], Line[1], Line[2])
   
gmsh.model.occ.addCurveLoop([1, 2, 3, 4, 5, 6, 7, 8], 1)
gmsh.model.occ.addPlaneSurface([1], 1)
#
#draw the trace extra line, loop and 1d physical.
#
gmsh.model.occ.addLine(5, 8, 9)        #trace bottom line.
#
#draw the dielectric, extra two points, three lines, loop and surface and physical.
#
gmsh.model.occ.addPoint(-gw/2, -dh, 0, 5*tms, 10)        #dielectric bottom left
gmsh.model.occ.addPoint(gw/2, -dh, 0, 5*tms, 11)        #dielectric bottom right

gmsh.model.occ.addLine(10, 11, 10)       
gmsh.model.occ.addLine(11, 2, 11)        
gmsh.model.occ.addLine(1, 10, 12)        
gmsh.model.occ.addCurveLoop([10, 11, -5, -9, -1, 12], 3)
gmsh.model.occ.addPlaneSurface([3], 2)
#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(1, [9, -4, -3, -2], 2, name="traceline")
gmsh.model.addPhysicalGroup(1, [10, 11, 6, 7, 8, 12], 1, name="groundline")
gmsh.model.addPhysicalGroup(2, [1], 4, name="airsurface")
gmsh.model.addPhysicalGroup(2, [3], 3, name="dielectricsurface")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("microstrip_py.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

=======
#DL231219
#microstrip model in python.
#microstrip 2d geometric model.
#to be used MFEM to simulate the electric field.
#Control the mesh size with using points mesh size arguments.
#
#define the following physical
#groundline, traceline, airsurface, dielectricsurface

import gmsh
import sys

gmsh.initialize()

gmsh.model.add("microstrip_py")

th=1.0;		#trace height.
tw=2.0;		#trace width.
dh=5;		#dielectric heigh.
gw=20;		#ground width.
hspace=20;
tms=0.5;	#target mesh size.
tlgap=0.2;	#the gap around the trace for the int boundary.

#
# Start by drawing the air points, line, loop, surface and physical.
#

Points=[    [-gw/2, 0, 0, tms, 1],     #air bottom-left
            [gw/2, 0, 0, tms,2],       #air bottom-right
            [gw/2, hspace, 0, 10*tms, 3], #air top-right
            [-gw/2, hspace, 0, 10*tms, 4], #air top-left
            [-tw/2, 0, 0, tms, 5],      #trace bottom-left
            [-tw/2, th, 0, tms, 6],     #trace top-left
            [tw/2, th, 0, tms, 7],      #trace top-right
            [tw/2, 0, 0, tms, 8]	#trace bottom-right
            ]



for Point in Points:
   gmsh.model.occ.addPoint(Point[0], Point[1], Point[2], Point[3], Point[4])
   
   

Lines=[     [1, 5, 1], 	#air bottom-left
            [5, 6, 2],	#trace left
            [6, 7, 3],	#trace top
            [7, 8, 4],	#trace right
            [8, 2, 5],	#air horiz-right
            [2, 3, 6],	#air vert-right
            [3, 4, 7],	#air top
            [4, 1, 8]	#air vert-left
            ]

for Line in Lines:
   gmsh.model.occ.addLine(Line[0], Line[1], Line[2])
   
gmsh.model.occ.addCurveLoop([1, 2, 3, 4, 5, 6, 7, 8], 1)
gmsh.model.occ.addPlaneSurface([1], 1)
#
#draw the trace extra line, loop and 1d physical.
#
gmsh.model.occ.addLine(5, 8, 9)        #trace bottom line.
#
#draw the dielectric, extra two points, three lines, loop and surface and physical.
#
gmsh.model.occ.addPoint(-gw/2, -dh, 0, 5*tms, 10)        #dielectric bottom left
gmsh.model.occ.addPoint(gw/2, -dh, 0, 5*tms, 11)        #dielectric bottom right

gmsh.model.occ.addLine(10, 11, 10)	#diel bottom       
gmsh.model.occ.addLine(11, 2, 11)       #diel right
gmsh.model.occ.addLine(1, 10, 12)       #diel left
gmsh.model.occ.addCurveLoop([10, 11, -5, -9, -1, 12], 3)
gmsh.model.occ.addPlaneSurface([3], 2)


#gmsh.model.occ.addPoint(-5.0, 5, 0.0, tms, 12)
#gmsh.model.occ.addPoint(5.0, 5, 0.0, tms, 13)
#gmsh.model.occ.addLine(12, 13, 13)

#gmsh.model.occ.addCircle(0, th/2, 0, 2.5, 13)
#gmsh.model.occ.synchronize()
#gmsh.model.mesh.embed(1, [13], 2, 1)


#
#synchronize prior to add physical group.
#
gmsh.model.occ.synchronize()
#
#add physical groups.
#

gmsh.model.addPhysicalGroup(1, [9, -4, -3, -2], 2, name="traceline")
gmsh.model.addPhysicalGroup(1, [10, 11, 6, 7, 8, 12], 1, name="groundline")
gmsh.model.addPhysicalGroup(2, [1], 3, name="airsurface")
gmsh.model.addPhysicalGroup(2, [2], 4, name="dielecsurface")
#gmsh.model.addPhysicalGroup(1, [13], 5, name="integcloseloop")


# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# glvis can read mesh version 2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# ... and save it to disk
gmsh.write("microstrip_py.msh")

# start gmsh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#before leaving.
gmsh.finalize()

>>>>>>> 5a55bdaea108e5a3017ae82a5b065cdeda7010db
