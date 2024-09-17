/*    Written by Denis Lachapelle July 2024.
                            MFEM
//
// Compile with: make proximity3d
//

The ultimate goal is to compute the proximity effect is a pair of wires in 3D.
1- compute the current density in DC (static) in the domain.
1a- compute the scalar voltage potential using ...
1b- project the gradient of the scalar potential on RT.
2- use the above current density to compute the vector magnetic potential A.
3- using the static current density and the current density caused by A compute the
   total current density.

-m <mesh file>, default "ProxRectRoundWires3d.msh".
-po <order potential mesh element order>, default 1.
-iro <integration order>, Default 1.
-rt <refine to>, default 1.
-rto <raviart thomas order>, default 1.

*/

/* These definitions fit with the mesh file ProxRoundWires3d.msh */
#define air 1 //domain
#define conductorright 2 //domain
#define conductorleft 3  //domain
#define airsurffront 4  //boundary
#define airsurfback 5
#define conductorleftsurffront 6
#define conductorleftsurfback 7
#define conductorrightsurffront 8
#define conductorrightsurfback 9
#define airsurfaround 10

#include <mfem.hpp>
#include <miniapps/common/fem_extras.hpp>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

//Compute the DC resistance of the wires.
double ComputeDCR(GridFunction &gf , int Boundary);

int main(int argc, char *argv[])
{
   const char *meshFile = "ProxRectRoundWires3d.msh"; //default mesh file.
   int potentialOrder = 2; // default order for potential elements.
   int eFieldOrder = potentialOrder; // raviart thomas element order.
   int refineTo = 1; // Default 1 cause no refinement.
   OptionsParser args(argc, argv);
   args.AddOption(&meshFile, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&potentialOrder, "-po", "--potOrder", "Finite element polynomial degree");
   args.AddOption(&refineTo, "-rerfto", "--refineTo", "Refine to _ elements");
   args.AddOption(&eFieldOrder, "-efo", "--eFOrder", "Raviat Thomas Element order");
   args.ParseCheck();

   //  Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device("cpu");
   device.Print();

   //***********************
   //***********************
   //  Read the mesh from the given mesh file. We can handle triangular,
   //  quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //  the same code.
   Mesh mesh(meshFile, 1, 1);

   //Refine the mesh, better to refine with the generator tool to improve section borders.
   int refineTimes = 0;
   while(mesh.GetNE()<refineTo) {
      mesh.UniformRefinement();
      refineTimes++;
   } 
   cout << "refineTimes =  " << refineTimes << endl; 
   
   int meshDim = mesh.Dimension();

   cout << "mesh.Dimension() = "<< mesh.Dimension() << endl;
   cout << "mesh.GetNE() = "<< mesh.GetNE() << endl;
   cout << "mesh.GetNBE() = "<< mesh.GetNBE() << endl;
   cout << "mesh.GetNEdges() = "<< mesh.GetNEdges() << endl;
   cout << "mesh.GetNFaces() = "<< mesh.GetNFaces() << endl;
   cout << "mesh.bdr_attributes.Max() = "<< mesh.bdr_attributes.Max() << endl;

   //***********************
   //***********************
   //  Define a finite element space on the mesh. Here we use 
   //  continuous Lagrange finite elements, since we seek a scalar.
    FiniteElementCollection *potFec = (FiniteElementCollection*)new H1_FECollection(potentialOrder, meshDim);
    FiniteElementSpace potFESpace(&mesh, potFec);
    int potFESpaceSize = potFESpace.GetTrueVSize();
    cout << "Number of finite element unknowns: " << potFESpaceSize << endl;
   
   // 4. Create "marker arrays" to define the portions of boundary associated
   //    with each type of boundary condition. These arrays have an entry
   //    corresponding to each boundary attribute.  Placing a '1' in entry i
   //    marks attribute i+1 as being active, '0' is inactive.
   //    in this case there are only dirichelet boundary.
      
   Array<int> potDirBdrMkr(mesh.bdr_attributes.Max());
   cout << "mesh.bdr_attributes.Max() = " << mesh.bdr_attributes.Max() << endl;
   assert(mesh.bdr_attributes.Max()==10);  // to fit the mesh file.
   potDirBdrMkr = 1;
   potDirBdrMkr[conductorright-1]=0;
   potDirBdrMkr[conductorleft-1]=0;
   potDirBdrMkr[air-1]=0;
   // both airsurffront and airsurfback should not be forced boundary value
   // to avoid very large gradient at the boundary.
   potDirBdrMkr[airsurffront-1]=0; 
   potDirBdrMkr[airsurfback-1]=0;
   

   Array<int> potEssTdofList(0);
   if (mesh.bdr_attributes.Size())
   {
      // For a continuous basis the linear system must be modified to enforce an
      // essential (Dirichlet) boundary condition. 
      potFESpace.GetEssentialTrueDofs(potDirBdrMkr, potEssTdofList);
   }

   //conduction air, copper, copper.
   double domainCoeffArray[mesh.bdr_attributes.Max()]={0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   Vector domainCoeffVector(domainCoeffArray, mesh.bdr_attributes.Max());
   PWConstCoefficient domainCoeff(domainCoeffVector);

   //  Define the solution vector u as a finite element grid function
   //  corresponding to fespace. Initialize u with initial guess of zero.
   GridFunction potGridFunction(&potFESpace);
   potGridFunction = 0.0;

   //  Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm potBilinearForm(&potFESpace);
   potBilinearForm.AddDomainIntegrator(new DiffusionIntegrator(domainCoeff));
   potBilinearForm.Assemble();

   //  Assemble the linear form for the right hand side vector.
   LinearForm potLinearForm(&potFESpace);

   // Set the Dirichlet values in the solution vector
   double BoundaryCoeffArray[mesh.bdr_attributes.Max()];
   BoundaryCoeffArray[air-1]=0.0;
   BoundaryCoeffArray[conductorright-1]=0.0;
   BoundaryCoeffArray[conductorleft-1]=0.0;
   BoundaryCoeffArray[airsurffront-1]=0.0;
   BoundaryCoeffArray[airsurfback-1]=0.0;
   BoundaryCoeffArray[conductorleftsurffront-1]=-1.0;  // -1.0 Volt.
   BoundaryCoeffArray[conductorleftsurfback-1]=1.0;    // +1.0 Volt.
   BoundaryCoeffArray[conductorrightsurffront-1]=1.0;
   BoundaryCoeffArray[conductorrightsurfback-1]=-1.0;
   BoundaryCoeffArray[airsurfaround-1]=0.0;

   Vector BoundaryCoeffVector(BoundaryCoeffArray, mesh.bdr_attributes.Max());
   PWConstCoefficient BoundaryCoeff(BoundaryCoeffVector);
   potGridFunction.ProjectBdrCoefficient(BoundaryCoeff, potDirBdrMkr);
   potLinearForm.Assemble();

   //  Construct the linear system.
   OperatorPtr A;
   Vector B, X;
   potBilinearForm.FormLinearSystem(potEssTdofList, potGridFunction, potLinearForm, A, X, B);
   
   //  Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system AX=B with PCG in the symmetric case, and GMRES in the
   //     non-symmetric one.
        GSSmoother M((SparseMatrix&)(*A));
         PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
 
//  Recover the grid function corresponding to U. This is the local finite
//     element solution.
   potBilinearForm.RecoverFEMSolution(X, potLinearForm, potGridFunction);


   // 15. Send the potential solution by socket to a GLVis server.
      string title_str = "H1";
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << potGridFunction
               << "window_title '" << title_str << " Solution'"
               << " keys 'mmc'" << flush;

// **********************
// **********************
   cout << "*****\ncompute gradient\n*****\n" << endl;
   eFieldOrder = potentialOrder-1;  //not sure if it is ok.
   // This section computes the gradient (field) using raviart thomas basis function.
   RT_FECollection eFieldFec(eFieldOrder, meshDim);
   FiniteElementSpace eFieldFESpace(&mesh, &eFieldFec);
   GridFunction eFieldGridFunction(&eFieldFESpace);  //trial space is RT.
   {
     LinearForm eFieldLinearForm(&eFieldFESpace);
     MixedBilinearForm eFieldMixedBilinearForm(&potFESpace, &eFieldFESpace);
     eFieldMixedBilinearForm.AddDomainIntegrator(new MixedVectorGradientIntegrator(domainCoeff));
     eFieldMixedBilinearForm.Assemble();
     eFieldMixedBilinearForm.Finalize();
     eFieldMixedBilinearForm.Mult(potGridFunction, eFieldLinearForm);

     BilinearForm eFieldBilinearForm(&eFieldFESpace);
     eFieldBilinearForm.AddDomainIntegrator(new VectorFEMassIntegrator);
     eFieldBilinearForm.Assemble();
     eFieldBilinearForm.Finalize();

     Array<int> eFieldEssTdofList;
     OperatorPtr A;
     Vector B, X;

     eFieldGridFunction = 0.0;
     eFieldBilinearForm.FormLinearSystem(eFieldEssTdofList, eFieldGridFunction, eFieldLinearForm, A, X, B);

     GSSmoother M((SparseMatrix&)(*A));
     PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
     eFieldBilinearForm.RecoverFEMSolution(X, eFieldLinearForm, eFieldGridFunction);
     eFieldGridFunction *= -1.0;
   }

  //Send the gradient solution by socket to a GLVis server.
   {
     string title_str = "E-Field";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << mesh << eFieldGridFunction
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   }

// from ex5.cpp  
// 14. Save data in the VisIt format
   VisItDataCollection visit_dc("two_wires", &mesh);
   visit_dc.RegisterField("potential", &potGridFunction);
   visit_dc.RegisterField("e-field", &eFieldGridFunction);
   visit_dc.Save();

   // 15. Save data in the ParaView format
   ParaViewDataCollection paraview_dc("two_wires_pw", &mesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(potentialOrder);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("e-field",&eFieldGridFunction);
   paraview_dc.RegisterField("potential",&potGridFunction);
   paraview_dc.Save();

   
   
// **********************
// **********************
   cout << "compute the currents in both wires at bothe extremities." << endl;
   double Rconductorrightsurffront = ComputeDCR(eFieldGridFunction, conductorrightsurffront);
   cout << "Rconductorrightsurffront = " << Rconductorrightsurffront << endl;
   double Rconductorrightsurfback = ComputeDCR(eFieldGridFunction, conductorrightsurfback);
   cout << "Rconductorrightsurfback = " << Rconductorrightsurfback << endl;
   double Rconductorleftsurffront = ComputeDCR(eFieldGridFunction, conductorleftsurffront);
   cout << "Rconductorleftsurffront = " << Rconductorleftsurffront << endl;
   double Rconductorleftsurfback = ComputeDCR(eFieldGridFunction, conductorleftsurfback);
   cout << "Rconductorleftsurfback = " << Rconductorleftsurfback << endl;
   // 16. Free the used memory.
   delete potFec;

// **********************
// **********************
   cout << "*****\compute the H-field\n*****\n" << endl;

   return 0;
}

/* Compute the DC resistance of the two wires.*/
double ComputeDCR(GridFunction &gf, int Boundary)
   {
      mfem::FiniteElementSpace *sp = gf.FESpace();
      mfem::Mesh *mesh = sp->GetMesh();
      mfem::Array<int> Current_bdr(mesh->bdr_attributes.Max());
      Current_bdr = 0;
      assert((mesh->bdr_attributes.Max()) >= Boundary);
      Current_bdr[Boundary-1]=1;
      mfem::ConstantCoefficient one(1.0);
      mfem::LinearForm Ifr_LF(sp);
      Ifr_LF.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(one), Current_bdr);      
      Ifr_LF.Assemble();
      double Current = (Ifr_LF)(gf);
      return 2.0/Current;
   }



/*
Function intgrad

// This code compute the line integral of the gradient
// of the potential computed with electrostatic.cpp.
// the initial goal was to compare to the result obtain with 
// the integration on raviart-thomas D-field to get the charge or capacitance.

compute the lenght of the line.
divide by the step to get the sample lenght.
in a while loop based on lenght.
check if there is one or more samples.
then find the middle of the sample.
get the grad and multiply by sample lenghtand and coeff and sum.
*/

double intgrad(double x0, double y0, double x1, double y1, double NbrStep, double delta, double *CoeffArray, GridFunction& u) {
   double TotalLenght=2.0*((x1-x0)+(y1-y0));
   double SampleLenght=TotalLenght/NbrStep;
   double pos=x0, RunningLenght=0.0, x, y;
   double CurrentSampleLenght;
   double grad=0.0;
   int j=0;
   int NbrPoints=2;
   Vector XYPoints(2*NbrPoints); //x, y.
   enum STATE {BOTTOM_X, RIGHT_Y, TOP_X, LEFT_Y};
   int state=BOTTOM_X;
   while(RunningLenght<TotalLenght) {
      switch(state) {

         case BOTTOM_X:  // from (x0, y0) to the right.
            y=y0;
            if(pos+SampleLenght<=x1) {
               x=pos+SampleLenght/2.0;
               pos+=SampleLenght;
               CurrentSampleLenght=SampleLenght;
               }
            else {
               x=pos+(x1-pos)/2;
               CurrentSampleLenght=x1-pos;
               pos=y0;
               state=RIGHT_Y;
               }
            XYPoints[0]=x;
            XYPoints[1]=y-delta/2; 
            XYPoints[2]=x;
            XYPoints[3]=y+delta/2;
            break;

         case RIGHT_Y:  // from (x1, y0) to the up.
            x=x1;
            if(pos+SampleLenght<=y1) {
               y=pos+SampleLenght/2.0;
               pos+=SampleLenght;
               CurrentSampleLenght=SampleLenght;
               }
            else {
               y=pos+(y1-pos)/2.0;
               CurrentSampleLenght=y1-pos;
               pos=x1;
               state=TOP_X;
               }
            XYPoints[2]=x-delta/2;
            XYPoints[3]=y;
            XYPoints[0]=x+delta/2;
            XYPoints[1]=y;
            break;
            
         case TOP_X:  //from (x1, y1) to the left.
            y=y1;
            if(pos-SampleLenght>=x0) {
               x=pos-SampleLenght/2.0;
               pos-=SampleLenght;
               CurrentSampleLenght=SampleLenght;
               }
            else {
               x=pos-(-x0+pos)/2.0;
               CurrentSampleLenght=(-x0+pos);
               pos=y1;
               state=LEFT_Y;
               }
            XYPoints[2]=x;
            XYPoints[3]=y-delta/2; 
            XYPoints[0]=x;
            XYPoints[1]=y+delta/2;
            break;

         case LEFT_Y:
            x=x0;
            if(pos-SampleLenght>=y0) {
               y=pos-SampleLenght/2.0;
               pos-=SampleLenght;
               CurrentSampleLenght=SampleLenght;
               }
            else {
               y=pos-(-y0+pos)/2.0;
               CurrentSampleLenght=(-y0+pos);
               pos=x0;
               RunningLenght=TotalLenght+1;
               }
            XYPoints[0]=x-delta/2;
            XYPoints[1]=y;
            XYPoints[2]=x+delta/2;
            XYPoints[3]=y;
            break;
         }
      RunningLenght+=CurrentSampleLenght;


   //Transfert the points in the matrix.
   DenseMatrix point_mat(XYPoints.GetData(), 2, 2);

   Array<int> elem_ids(NbrPoints); // el ement ids.
   Array<IntegrationPoint> ips(NbrPoints);  // the location within the element.
   //cout << "###, " << j++ << ", " << RunningLenght << " point_mat.Print(); = " << endl; 
   //point_mat.Print();
   Mesh * m = u.FESpace()->GetMesh();
   assert(m->FindPoints(point_mat, elem_ids, ips)==2); // find the element and the point in the element.
   double val[NbrPoints];

   int attr[2];
   attr[0]=m->GetAttribute(elem_ids[0]);
   attr[1]=m->GetAttribute(elem_ids[1]);
   

   // get the value of each point one by one.
   int i;
   for(i=0; i< NbrPoints; i++) {
      val[i] = u.GetValue(elem_ids[i], ips[i], 2);
      }

if(0) {
   cout << j++  << ", "  << RunningLenght << ", " <<  CurrentSampleLenght << ", " << XYPoints[0] << ", " << XYPoints[1]
        << ", " << val[0] << ", " << val[1]
        << ", " << val[1]-val[0] << endl;
        }

if(0) {
   cout << j++  << ", "  << XYPoints[1]
        << ", " << val[0] << ", " << attr[0] << ", " << CoeffArray[attr[0]] << endl;
        }

   // attr[0]-1, -1 is necessary to align dielectric coefficients;
   // but I do not fully understand why; i think it start to 1 instead of 0.
   grad+=CoeffArray[attr[0]-1]*CurrentSampleLenght*(val[1]-val[0])/delta; 
      }

   return grad;


   }
