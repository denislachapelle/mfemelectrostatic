/*    Written by Denis Lachapelle Feb 2024.
                            MFEM
//
// Compile with: make magnetostatic2d
//
// This code is about microstrip simulation.
// it compute the magnetic ...
// 
// The mesh is created with gmsh from the file microstrip.py.

-m <mesh file>, default "microstrip_py.msh".
-o <order potential mesh element order>, default 2.
-iro <integration order>, Default 1.
-rft <refine to>, default 1.
-dgi <>, default gradient integrator. default 1.
-ndo <nedelec order>, default potential order.

*/

#undef SECTION1

#include <mfem.hpp>
#include <miniapps/common/fem_extras.hpp>
#include <fstream>
#include <iostream>

#define TRACEATTRIB 5

using namespace std;
using namespace mfem;


double intgrad(double x0, double y0, double x1, double y1, double NbrStep, double delta, double *CoeffArray, GridFunction& u);


int main(int argc, char *argv[])
{
   double CurrentDensity;  //Current density in the trace causing 1 ampere.
   
   const char *mesh_file = "microstrip_py.msh"; //default mesh file.
   int order = 2; // default order for potential elements.
   int irorder = 1; // integration rule order, Default is (basis function order - 1)
   int nd_order = order; // nedelec element order.
   int refineTo = 1; // Default 1 cause no refinement.
   int dgi = 1; // 1 to use default gradient Integrator.
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.AddOption(&irorder, "-iro", "--irorder", "Gradient integration rule order");
   args.AddOption(&refineTo, "-rft", "--refineto", "Refine to _ elements");
   args.AddOption(&dgi, "-dgi", "--DefaultGradInt", "Default gradient integrator");
   args.AddOption(&nd_order, "-ndo", "--ndorder", "nedelec Element order");
   args.ParseCheck();

   Device device("cpu");
   device.Print();

   //  Read the mesh from the given mesh file, triangular elements.
   Mesh mesh(mesh_file, 1, 1);

   // Refine the mesh.
   cout << "mesh.GetNE() before refine = "<< mesh.GetNE() << endl;
   int RefineCount = 0;
   while(mesh.GetNE()<refineTo) {
      mesh.UniformRefinement();
      RefineCount++;
   } 
   cout << "Refine " << RefineCount << " times."<<endl; 
   
   int dim = mesh.Dimension();
   // in this program the mesh shall have dimension of 2.
   assert(dim == 2);

   cout << "mesh.Dimension() = "<< mesh.Dimension() << endl;
   cout << "mesh.GetNE() after refine = "<< mesh.GetNE() << endl;
   cout << "mesh.GetNBE() = "<< mesh.GetNBE() << endl;
   cout << "mesh.GetNEdges() = "<< mesh.GetNEdges() << endl;
   cout << "mesh.GetNFaces() = "<< mesh.GetNFaces() << endl;
   cout << "mesh.bdr_attributes.Max() = "<< mesh.bdr_attributes.Max() << endl;
   cout << "mesh.attributes.Max() = "<< mesh.attributes.Max() << endl;
   cout << "mesh.GetNV() = "<< mesh.GetNV() << endl;
   
   
   // the code is made for 5 attributes.
   assert(mesh.attributes.Max()==5);

   //compute the area of the physical group tracesurface,
   //and then compute the current density assuming current is 1A.
   {
      int i;
      double TraceArea=0.0;
         for(i=0; i<mesh.GetNE(); i++) {
            if(0) {cout << mesh.GetElementVolume(i) << endl;}
            if(mesh.GetAttribute(i)==TRACEATTRIB) TraceArea+=mesh.GetElementVolume(i);
         }
         CurrentDensity = 1.0/TraceArea; //Current is 1 ampere.
         cout << "TraceArea =" << TraceArea << endl;
         cout << "CurrentDensity =" << CurrentDensity << endl;
   }
   
   // setup a PWConstCoefficient to define the current density
   // in the trace. DL240129.
   // attribute 4 is the tracesurface.
   Vector CCJ(5); CCJ=0.0; CCJ[4]=CurrentDensity;
   PWConstCoefficient PWCJ(CCJ);

  // Define a finite element space on the mesh. Here we use 
  //    h1 finite elements. ex0.cpp.
    FiniteElementCollection *Afec = (FiniteElementCollection*)new H1_FECollection(order, dim);
    FiniteElementSpace Afespace(&mesh, Afec);
    int size = Afespace.GetTrueVSize();
    cout << "Number of finite element unknowns: " << size << endl;
   
   //  Create "marker arrays" to define the portions of boundary associated
   //    with each type of boundary condition. These arrays have an entry
   //    corresponding to each boundary attribute.  Placing a '1' in entry i
   //    marks attribute i+1 as being active, '0' is inactive.
   //    in this case there are only dirichelet boundary.
  
   assert(mesh.bdr_attributes.Max()==2);
   Array<int> dbc_bdr(mesh.bdr_attributes.Max());
   dbc_bdr = 0;
   dbc_bdr[0] = 1; //dirichelet for groundline physical line.
   dbc_bdr[1] = 0;

   Array<int> ess_tdof_list(0);
   if (mesh.bdr_attributes.Size())
   {
      // For a continuous basis the linear system must be modified to enforce an
      // essential (Dirichlet) boundary condition. 
      Afespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
   }


  // permeability pcb dielectric 1, copper 1 and air 1.
   double PermCoeffArray[]={0.0, 0.0, 1.0, 1.0, 1.0};
   Vector PermCoeffVector(PermCoeffArray, 5);
   PWConstCoefficient PermCoeff(PermCoeffVector);
   cout << "PermCoeffVector\n";
   PermCoeffVector.Print(cout, 5);


   // Define the solution vector u as a finite element grid function
   //    corresponding to fespace. Initialize u with initial guess of zero.
   GridFunction u(&Afespace);
   u = 0.0;

   // Set up the bilinear form a(.,.) on the finite element space
   //    and add the integrator.
   BilinearForm a(&Afespace);

   a.AddDomainIntegrator(new DiffusionIntegrator(PermCoeff));
   a.Assemble();

   // Assemble the linear form for the right hand side vector.
   LinearForm b(&Afespace);

   // Set the Dirichlet values in the solution vector
   // the trace is at 1V, the contour is 0V.
   double BoundaryCoeffArray[]={0.0, 0.0};
   Vector BoundaryCoeffVector(BoundaryCoeffArray, 2);
   PWConstCoefficient BoundaryCoeff(BoundaryCoeffVector);
   u.ProjectBdrCoefficient(BoundaryCoeff, dbc_bdr);


   b.AddDomainIntegrator(new DomainLFIntegrator(PWCJ));

   //set the newmann boundary conditions.
   assert(mesh.bdr_attributes.Max()==2);
   Array<int> nbc_marker(mesh.bdr_attributes.Max());
   nbc_marker = 0; nbc_marker[0] = 0; nbc_marker[1] = 0;

   double CurlCoeffArray[]={0.0, 0.0};
   Vector CurlCoeffVector(CurlCoeffArray, 2);
   VectorConstantCoefficient CurlCoeff(CurlCoeffVector);

   b.AddBoundaryIntegrator(new BoundaryTangentialLFIntegrator(CurlCoeff), nbc_marker);
   b.Assemble();
 
   // Construct the linear system.
   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);

   cout << "A->NumRows(): " << A->NumRows() << endl;
    cout << "A->NumCols(): " << A->NumCols() << endl;

   // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system AX=B with PCG in the symmetric case, and GMRES in the
   //     non-symmetric one.
        GSSmoother M((SparseMatrix&)(*A));
         PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
 cout << "step #12" << endl;
// 12. Recover the grid function corresponding to U. This is the local finite
//     element solution.
   a.RecoverFEMSolution(X, b, u);

//save the gridfuntion class for later use.
//this is not require, but was done by curiosity.
   u.Save("ESGridFile.txt");

//
// test get gradient.
//

   

// 14. Save the refined mesh and the solution. This output can be viewed
   //     later using GLVis: "glvis -m refined.mesh -g sol.gf".
   {
      ofstream mesh_ofs("refined.mesh");
      mesh_ofs.precision(8);
      mesh.Print(mesh_ofs);
      ofstream sol_ofs("sol.gf");
      sol_ofs.precision(8);
      u.Save(sol_ofs);
   }

   // 15. Send the magnetic vector potential solution by socket to a GLVis server.
      string title_str = "magnetic vector potential Z-axe";
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << u
               << "window_title '" << title_str << " Solution'"
               << " keys 'mmcvvh'" << flush;

/*compute the line integral along a close path; it should give the magnetomotrice force
  which are amp*turn. test along different paths and along a path enclosing 1/2 of the
  center rectangle. */
double AmpTurn;
AmpTurn=intgrad(-3.0, 9.0, 3.0, 12.0, 100, 0.01, PermCoeffArray, u);
cout << "expect 1, AmpTurn=" << AmpTurn <<endl;                 
AmpTurn=intgrad(-4.0, 9.0, 3.0, 12.0, 100, 0.01, PermCoeffArray, u);
cout << "expect 1, AmpTurn=" << AmpTurn <<endl;                 
AmpTurn=intgrad(-4.0, 9.0, 4.0, 12.0, 100, 0.01, PermCoeffArray, u);
cout << "expect 1, AmpTurn=" << AmpTurn <<endl;                 
AmpTurn=intgrad(-0.0, 9.0, 3.0, 12.0, 100, 0.01, PermCoeffArray, u);
cout << "expect 0.5, AmpTurn=" << AmpTurn <<endl;    


/*
Compute the integral over the entire domain.
Resulting in the Sum and The sum of square.
Get help from ex38 and ex27 line 639.
https://en.wikipedia.org/wiki/Gaussian_quadrature
*/
double Sum=0.0; //the running summation.
double SumSquare=0.0;
double Temp;
Element *El;   //a pointer to the current element.
IntegrationRule ir;  //integration rule for the current element.
IntegrationRules irs;
ElementTransformation *Trans;
for (int i = 0; i < mesh.GetNE(); i++)
    {
      Trans = mesh.GetElementTransformation(i);
      El = mesh.GetElement(i);
      ir = irs.Get(El->GetGeometryType(), order+1/2);
      Array<double> Weights(ir.GetNPoints());
      Weights = ir.GetWeights();  //quadrature points weights.
      if (0) cout << "ir.GetNPoints = " << ir.GetNPoints() << endl;
      for(int j=0; j<ir.GetNPoints(); j++)
      {
         Trans->SetIntPoint(&ir.IntPoint(j));
         // Trans->Weight() returns the determinant of the jaccobians.
         Temp=u.GetValue(i, ir.IntPoint(j))*Weights[j]*Trans->Weight();
         Sum+=Temp;
         SumSquare+=Temp*Temp;

         if(0) {
            cout << "i, j = " << i << " ," << j << endl;
            cout << "Trans.Weight() = " << Trans->Weight() << endl;         
            cout << "Point = " << ir.IntPoint(j).x << " ," << ir.IntPoint(j).y << endl;         
            cout << "u.GetValue = " << u.GetValue(i, ir.IntPoint(j)) << endl;
            cout << "Weights[j] = " << Weights[j] << endl;
         }
      }
    }
    cout << "SumSquare = " << SumSquare << endl;
    cout << "Sum = " << Sum << endl;
    


cout << "step compute gradient" << endl;

if(0){
   // This section computes the gradient (field).
   // use nedelec basis function.
   ND_FECollection nd_fec(nd_order, dim);
   FiniteElementSpace nd_fespace(&mesh, &nd_fec);
   GridFunction D(&nd_fespace);
   {
     LinearForm epsdT(&nd_fespace);
     MixedBilinearForm epsGrad(&Afespace, &nd_fespace);
     epsGrad.AddDomainIntegrator(new MixedVectorGradientIntegrator(PermCoeff));
     epsGrad.Assemble();
     epsGrad.Finalize();
     epsGrad.Mult(u, epsdT);

     BilinearForm nd_m(&nd_fespace);
     nd_m.AddDomainIntegrator(new VectorFEMassIntegrator);
     nd_m.Assemble();
     nd_m.Finalize();

     Array<int> ess_tdof_rt_list;
     OperatorPtr A;
     Vector B, X;

     D = 0.0;
     nd_m.FormLinearSystem(ess_tdof_rt_list, D, epsdT, A, X, B);

     GSSmoother M((SparseMatrix&)(*A));
     PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
     nd_m.RecoverFEMSolution(X, epsdT, D);
     D *= -1.0;
   }
}

  /* From example 24, gradient can be computed using GradientInterpolator
     as suggested in https://github.com/mfem/mfem/issues/4118.
  */
   ND_FECollection nd_fec(nd_order, dim);
   FiniteElementSpace nd_fespace(&mesh, &nd_fec);
   GridFunction D(&nd_fespace);
   DiscreteLinearOperator dlo(&Afespace, &nd_fespace);
   dlo.AddDomainInterpolator(new GradientInterpolator());
   dlo.Assemble();
   dlo.Mult(u, D);


   D.Save("GradOfPotGridFile.txt");

  //Send the gradient solution by socket to a GLVis server.
   {
     string title_str = "Gradient of magentic potential";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << mesh << D
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   }
   
   /* the magnetis potential computed above is rotated 90 debree compared to the 
   magentic field. Copying the dofs from nedelec to raviart thomas fe space.
   Note the order of RT FE are nd_order-1.
   Get help from https://github.com/mfem/mfem/issues/4118
   */
   Vector Data( nd_fespace.GetVSize());
   D.GetTrueDofs(Data);
   RT_FECollection rt_fec(nd_order-1, dim);
   FiniteElementSpace rt_fespace(&mesh, &rt_fec);
   //Create a GridFunction for rt fe space with the dofa of the nedelec.
   GridFunction ROT90(&rt_fespace, Data);

   //Send the msgnetic solution by socket to a GLVis server.
   {
     string title_str = "suoposed magnetic field";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << mesh << ROT90
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   }

   if(0)   {
      int Size = Data.Size();
      int vdim = D.VectorDim();
      assert(vdim==2);
      assert(Size%vdim==0);
      double temp;
      if(0) for(int i=0; i<Size/2; i++) {
         temp=Data[i];
         Data[i]=Data[i+Size/2];
         Data[i+Size/2]=temp;
      }
      if(1) for(int i=0; i<Size; i+=2) {
         temp=Data[i];
         Data[i]=-Data[i+1];
         Data[i+1]=temp;
      }
   //  D.SetFromTrueDofs(Data);
   }

// from ex5.cpp  
// 14. Save data in the VisIt format
 //  VisItDataCollection visit_dc("microstrip", &mesh);
//   visit_dc.RegisterField("potential", &u);
//   visit_dc.RegisterField("gradient", &dT);
//   visit_dc.Save();

#ifdef NOTDEF
 cout << "step compute integral" << endl;
   //Compute the integral of the gradient on the boundary,
   //this will be the charge.
   // With the help of https://github.com/mfem/mfem/issues/993
   // and from from volta_solver.cpp and 
   // https://github.com/mfem/mfem/issues/3753

   LinearForm *rt_surf_int_;
   {
      rt_surf_int_ = new LinearForm(&fespace_rt);
      cout << "integration on Raviart Thomas"<<endl;
      Array<int> bdr_marker(2); bdr_marker[0]=1; bdr_marker[1]=0;

     if(dgi==1) {
        rt_surf_int_->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator,
					bdr_marker);
     }
     else {
      //get access to the integrator to change the order.
        LinearFormIntegrator *linint = new VectorFEBoundaryFluxLFIntegrator;
        const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, irorder);
        linint->SetIntRule(&ir);
        rt_surf_int_->AddBoundaryIntegrator(linint, bdr_marker);
        cout<<"integration rule order = "<<irorder<<endl;
     }    
    
     rt_surf_int_->Assemble();
     double charge_D1 = (*rt_surf_int_)(D);
     cout << endl<<"charge_D1 = "<<charge_D1<<endl;
     delete rt_surf_int_;
   }

   {
     rt_surf_int_ = new LinearForm(&fespace_rt);
     cout << "integration on Raviart Thomas"<<endl;
          Array<int> bdr_marker(2); bdr_marker[0]=0; bdr_marker[1]=1;
     if(dgi==1) {
        rt_surf_int_->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator,
					bdr_marker);
     }
     else {
      //get access to the integrator to change the order.
        LinearFormIntegrator *linint = new VectorFEBoundaryFluxLFIntegrator;
        const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, irorder);
        linint->SetIntRule(&ir);
        rt_surf_int_->AddBoundaryIntegrator(linint, bdr_marker);
        cout<<"integration rule order = "<<irorder<<endl;
     }    
     rt_surf_int_->Assemble();
     double charge_D2 = (*rt_surf_int_)(D);
     cout << endl<<"charge_D2 = "<<charge_D2<<endl;
     delete rt_surf_int_;
   }

   if(0) {
      ifstream infile;
      infile.open ("ESGridFile.txt");
      GridFunction g(&mesh, infile);
      double D3=intgrad(-9.0, 1.0, 9.0, 19.0, 760, 0.01, CoeffArray, g);
      cout << "charge_D3 = " << D3 << endl;
      }
#endif
   // 16. Free the used memory.
   delete Afec;
   mesh.Save("microstrip1_mfem1.mesh");
   return 0;
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
