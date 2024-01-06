<<<<<<< HEAD
//DL230710
//microstrip 2d geometric model.
//to be used MFEM to simulate the electric field.

SetFactory("OpenCASCADE");

th=0.1;		//trace heigh.
tw=10;		//trace width.
dh=5;		//dielectric heigh.
gw=20;		//ground width.
hspace=20;
rad=0.01;	//radius corner.
tms=0.5;	//target mesh size.
div=1; //target mesh size divider used around trace.


// the square containig the entire domain.
//bottom left.
Point(41)={-gw/2+rad, 0, 0, tms/div};
Point(42)={-gw/2, 0+rad, 0, tms};
Point(43)={-gw/2+rad, 0+rad, 0, tms/div};
Circle(4)={42, 43, 41};
//Circle(4)={-gw/2+rad, rad, 0, rad, 2*Pi/2, 3*Pi/2};
//bottom right.
Point(31)={gw/2-rad, 0, 0, tms/div};
Point(32)={gw/2, 0+rad, 0, tms};
Point(33)={gw/2-rad, 0+rad, 0, tms/div};
Circle(3)={32, 33, 31};
//Circle(3)={gw/2-rad, rad, 0, rad, 3*Pi/2, 4*Pi/2};
//top right.
Point(22)={gw/2-rad, hspace, 0, tms/div};
Point(21)={gw/2, hspace-rad, 0, tms/div};
Point(23)={gw/2-rad, hspace-rad, 0, tms/div};
Circle(2)={21, 23, 22};
//Circle(2)={gw/2-rad, hspace-rad, 0, rad, 0*Pi/2, 1*Pi/2};
//top left.
Point(12)={-gw/2+rad, hspace, 0, tms/div};
Point(11)={-gw/2, hspace-rad, 0, tms/div};
Point(13)={-gw/2+rad, hspace-rad, 0, tms/div};
Circle(1)={11, 13, 12};
//Circle(1)={-gw/2+rad, hspace-rad, 0, rad, 1*Pi/2, 2*Pi/2};
// line bottom
Line(15)={41, 31};
//line top.
Line(7)={22, 12};


// the dielectric top.
Point(51)={-gw/2, dh, 0, tms/div};
Point(52)={gw/2, dh, 0, tms/div};

//line left.
Line(91)={51, 42};
Line(92)={51, 11};
//line right.
Line(111)={32, 52};
Line(112)={52, 21};

// trace bottom.
//left.
Point(91)={-tw/2+rad, dh, 0, tms/div};
Point(92)={-tw/2, dh+rad, 0, tms/div};
Point(93)={-tw/2+rad, dh+rad, 0, tms/div};
Circle(101)={92, 93, 91};
//Circle(101)={-tw/2+rad, dh+rad, 0, rad, 2*Pi/2, 3*Pi/2 };

Line(5)={91, 51};
//right.
Point(81)={tw/2-rad, dh, 0, tms/div};
Point(82)={tw/2, dh+rad, 0, tms/div};
Point(83)={tw/2-rad, dh+rad, 0, tms/div};
Circle(102)={81, 83, 82};
//Circle(102)={tw/2-rad, dh+rad, 0, rad, 3*Pi/2, 4*Pi/2 };
Line(93)={81,91};
Line(6)={52, 81};


// trace top.
//right.
Point(62)={tw/2-rad, dh+th, 0, tms/div};
Point(63)={tw/2, dh+th-rad, 0, tms/div};
Point(64)={tw/2-rad, dh+th-rad, 0, tms/div};
Circle(103)={62, 64, 63};
//Circle(103)={-tw/2+rad, dh+th-rad, 0, rad, 1*Pi/2, 2*Pi/2 };
//left.
Point(72)={-tw/2+rad, dh+th, 0, tms/div};
Point(73)={-tw/2, dh+th-rad, 0, tms/div};
Point(74)={-tw/2+rad, dh+th-rad, 0, tms/div};
Circle(104)={72, 74, 73};
//Circle(104)={tw/2-rad, dh+th-rad, 0, rad, 0*Pi/2, 1*Pi/2 };
Line(107)={62, 72};
//trace line left.
Line(108)={92, 73};
//trace line right.
Line(109)={63, 82};

// ground
Physical Line("groundline") = {-1, -92, 91, 4, 15, -3, 111, 112, 2, 7};

// trace
Physical Line("traceline") = {104, -108, 101, -93, 102, -109, -103, 107};

// the countour of the dielectric.
Curve Loop(1) = {15, 3, 111, 6, 93, 5, 91, 4};
Plane Surface(1) = {1};
//dielectric
Physical Surface("dielecsurface") = {1};

// contour of the trace.
//Curve Loop(2) = {-101, -93, 102, -109, -103, 107, 104, -108};
//Plane Surface(2) = {2};
//Physical Surface("tracesurface") = {2};

//contour of the air.
Curve Loop(4) = {5, 101, 108, -104, -107, 103, 109, -102, 6, 112, 2, -7, 1, -92};
Plane Surface(4) = {4};
Physical Surface("airsurface") = {4};




Mesh 2;
Mesh.MshFileVersion = 2.2;
Save "microstrip_rnd4.msh";



=======
//DL230710
//microstrip 2d geometric model.
//to be used MFEM to simulate the electric field.

SetFactory("OpenCASCADE");

th=0.1;		//trace heigh.
tw=1.0;		//trace width.
dh=5;		//dielectric heigh.
gw=20;		//ground width.
hspace=20;
rad=0.01;	//radius corner.
tms=0.5;	//target mesh size.
div=1; //target mesh size divider used around trace.


// the square containig the entire domain.
//bottom left.
Point(41)={-gw/2+rad, 0, 0, tms/div};
Point(42)={-gw/2, 0+rad, 0, tms};
Point(43)={-gw/2+rad, 0+rad, 0, tms/div};
Circle(4)={42, 43, 41};
//Circle(4)={-gw/2+rad, rad, 0, rad, 2*Pi/2, 3*Pi/2};
//bottom right.
Point(31)={gw/2-rad, 0, 0, tms/div};
Point(32)={gw/2, 0+rad, 0, tms};
Point(33)={gw/2-rad, 0+rad, 0, tms/div};
Circle(3)={32, 33, 31};
//Circle(3)={gw/2-rad, rad, 0, rad, 3*Pi/2, 4*Pi/2};
//top right.
Point(22)={gw/2-rad, hspace, 0, tms/div};
Point(21)={gw/2, hspace-rad, 0, tms/div};
Point(23)={gw/2-rad, hspace-rad, 0, tms/div};
Circle(2)={21, 23, 22};
//Circle(2)={gw/2-rad, hspace-rad, 0, rad, 0*Pi/2, 1*Pi/2};
//top left.
Point(12)={-gw/2+rad, hspace, 0, tms/div};
Point(11)={-gw/2, hspace-rad, 0, tms/div};
Point(13)={-gw/2+rad, hspace-rad, 0, tms/div};
Circle(1)={11, 13, 12};
//Circle(1)={-gw/2+rad, hspace-rad, 0, rad, 1*Pi/2, 2*Pi/2};
// line bottom
Line(15)={41, 31};
//line top.
Line(7)={22, 12};


// the dielectric top.
Point(51)={-gw/2, dh, 0, tms/div};
Point(52)={gw/2, dh, 0, tms/div};

//line left.
Line(91)={51, 42};
Line(92)={51, 11};
//line right.
Line(111)={32, 52};
Line(112)={52, 21};

// trace bottom.
//left.
Point(91)={-tw/2+rad, dh, 0, tms/div};
Point(92)={-tw/2, dh+rad, 0, tms/div};
Point(93)={-tw/2+rad, dh+rad, 0, tms/div};
Circle(101)={92, 93, 91};
//Circle(101)={-tw/2+rad, dh+rad, 0, rad, 2*Pi/2, 3*Pi/2 };

Line(5)={91, 51};
//right.
Point(81)={tw/2-rad, dh, 0, tms/div};
Point(82)={tw/2, dh+rad, 0, tms/div};
Point(83)={tw/2-rad, dh+rad, 0, tms/div};
Circle(102)={81, 83, 82};
//Circle(102)={tw/2-rad, dh+rad, 0, rad, 3*Pi/2, 4*Pi/2 };
Line(93)={81,91};
Line(6)={52, 81};


// trace top.
//right.
Point(62)={tw/2-rad, dh+th, 0, tms/div};
Point(63)={tw/2, dh+th-rad, 0, tms/div};
Point(64)={tw/2-rad, dh+th-rad, 0, tms/div};
Circle(103)={62, 64, 63};
//Circle(103)={-tw/2+rad, dh+th-rad, 0, rad, 1*Pi/2, 2*Pi/2 };
//left.
Point(72)={-tw/2+rad, dh+th, 0, tms/div};
Point(73)={-tw/2, dh+th-rad, 0, tms/div};
Point(74)={-tw/2+rad, dh+th-rad, 0, tms/div};
Circle(104)={72, 74, 73};
//Circle(104)={tw/2-rad, dh+th-rad, 0, rad, 0*Pi/2, 1*Pi/2 };
Line(107)={62, 72};
//trace line left.
Line(108)={92, 73};
//trace line right.
Line(109)={63, 82};

// ground
Physical Line("groundline") = {-1, -92, 91, 4, 15, -3, 111, 112, 2, 7};

// trace
Physical Line("traceline") = {104, -108, 101, -93, 102, -109, -103, 107};

// the countour of the dielectric.
Curve Loop(1) = {15, 3, 111, 6, 93, 5, 91, 4};
Plane Surface(1) = {1};
//dielectric
Physical Surface("dielecsurface") = {1};

// contour of the trace.
//Curve Loop(2) = {-101, -93, 102, -109, -103, 107, 104, -108};
//Plane Surface(2) = {2};
//Physical Surface("tracesurface") = {2};

//contour of the air.
Curve Loop(4) = {5, 101, 108, -104, -107, 103, 109, -102, 6, 112, 2, -7, 1, -92};
Plane Surface(4) = {4};
Physical Surface("airsurface") = {4};




Mesh 2;
Mesh.MshFileVersion = 2.2;
Save "microstrip_rnd4.msh";



>>>>>>> 5a55bdaea108e5a3017ae82a5b065cdeda7010db
