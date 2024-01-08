Electrostatic.cpp is based on the MFEM library which is a finite element method simulation platform.
This software computes:
- the potential distribution between two boundaries including one or two physical area, such as air and dielectric.
- it then compute the D-field Gradient.
- then the close loop line integral of the gradient of the two boundaries, which give the total charge causing the field.

The model microstrip model are made using GMSH in python. microstrip.py.
