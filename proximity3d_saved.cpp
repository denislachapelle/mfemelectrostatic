//                                proximity effect
//                                inspired from MFEM Example 22
//                                prob 2, case 1.
//
// Compile with: make proximity3d
//
// Sample runs:  proximity3d -m ProxRoundWires3d.msh
//
// Description:  This example code demonstrates the use of MFEM to define and
//               solve simple complex-valued linear systems. It implements three
//               variants of a damped harmonic oscillator:
//
//               1) A scalar H1 field
//                  -Div(a Grad u) - omega^2 b u + i omega c u = 0
//
//               2) A vector H(Curl) field
//                  Curl(a Curl u) - omega^2 b u + i omega c u = 0
//
//               3) A vector H(Div) field
//                  -Grad(a Div u) - omega^2 b u + i omega c u = 0
//
//               In each case the field is driven by a forced oscillation, with
//               angular frequency omega, imposed at the boundary or a portion
//               of the boundary.
//
//               In electromagnetics, the coefficients are typically named the
//               permeability, mu = 1/a, permittivity, epsilon = b, and
//               conductivity, sigma = c. The user can specify these constants
//               using either set of names.
//
//               The example also demonstrates how to display a time-varying
//               solution as a sequence of fields sent to a single GLVis socket.
//
//               We recommend viewing examples 1, 3 and 4 before viewing this
//               example.

#include "mfem.hpp"
// #include "element.hpp"

#include <fstream>
#include <iostream>

#define airsurfaround 10
#define airsurffront 4
#define airsurfback 5
#define conductorrightsurffront 8
#define conductorleftsurffront 6
#define air 1
#define conductorright 2
#define conductorleft 3
#define conductorrightsurfback 9
#define conductorleftsurfback 7



using namespace std;
using namespace mfem;

static double mu_ = 1.257e-6;
static double epsilon_ = 8.854E-12;
static double sigma_ = 1.0/16.78e-9;
static double omega_ = 2.0*M_PI*60;

static double mu0_ = 1.25663706212E-6;
static double epsilon0_ = 8.8541878188E-12;

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "ProxRoundWires3d.msh";
   int ref_levels = 0;
   int order = 1;
   double freq = -1.0;
   double a_coef = 0.0;
   bool visualization = 1;
   bool herm_conv = true;
   bool pa = false;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&a_coef, "-a", "--stiffness-coef",
                  "Stiffness coefficient (spring constant or 1/mu).");
   args.AddOption(&epsilon_, "-b", "--mass-coef",
                  "Mass coefficient (or epsilon).");
   args.AddOption(&sigma_, "-c", "--damping-coef",
                  "Damping coefficient (or sigma).");
   args.AddOption(&mu_, "-mu", "--permeability",
                  "Permeability of free space (or 1/(spring constant)).");
   args.AddOption(&epsilon_, "-eps", "--permittivity",
                  "Permittivity of free space (or mass constant).");
   args.AddOption(&sigma_, "-sigma", "--conductivity",
                  "Conductivity (or damping constant).");
   args.AddOption(&freq, "-f", "--frequency",
                  "Frequency (in Hz).");
   args.AddOption(&herm_conv, "-herm", "--hermitian", "-no-herm",
                  "--no-hermitian", "Use convention for Hermitian operators.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   if ( freq > 0.0 )
   {
      omega_ = 2.0 * M_PI * freq;
   }

   ComplexOperator::Convention conv =
      herm_conv ? ComplexOperator::HERMITIAN : ComplexOperator::BLOCK_SYMMETRIC;

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes
   //    with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 4. Refine the mesh to increase resolution. In this example we do
   //    'ref_levels' of uniform refinement where the user specifies
   //    the number of levels with the '-r' option.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }

   cout << "mesh->Dimension() = "<< mesh->Dimension() << endl;
   cout << "mesh->GetNE() = "<< mesh->GetNE() << endl;
   cout << "mesh->GetNBE() = "<< mesh->GetNBE() << endl;
   cout << "mesh->GetNEdges() = "<< mesh->GetNEdges() << endl;
   cout << "mesh->GetNFaces() = "<< mesh->GetNFaces() << endl;
   cout << "mesh->bdr_attributes.Max() = "<< mesh->bdr_attributes.Max() << endl;
   
   
   FiniteElementCollection *fec = NULL;
   fec = new ND_FECollection(order /*for RT -1*/, dim); // switch to nedelec.
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl;


   // 8. Define the solution vector u as a complex finite element grid function
   //    corresponding to fespace. Initialize u with initial guess of 1+0i or
   //    the exact solution if it is known.
   ComplexGridFunction u(fespace);
   u=0.0;
 

   //compute the boundary area of the physical group conductorleftsurfback an
   //conductorrightback and then compute the current density assuming current is 1A.
   double CondLeftArea=0.0, CondRightArea=0.0;
   for (int it = 0; it < fespace->GetNBE(); it++)
   {
      if(mesh->GetBdrAttribute(it)==conductorleftsurfback) 
      {
         Vector normal(dim);
         ElementTransformation *Trans = fespace->GetBdrElementTransformation(it);
         Trans->SetIntPoint(&Geometries.GetCenter(Trans->GetGeometryType()));
         CalcOrtho(Trans->Jacobian(), normal);
         if(0) Trans->Jacobian().Print(cout);
         CondLeftArea += normal.Norml2();
         if (0) cout << "CondLeftArea = " << CondLeftArea << endl;
      }
   if(mesh->GetBdrAttribute(it)==conductorrightsurfback) 
      {
         Vector normal(dim);
         ElementTransformation *Trans = fespace->GetBdrElementTransformation(it);
         Trans->SetIntPoint(&Geometries.GetCenter(Trans->GetGeometryType()));
         CalcOrtho(Trans->Jacobian(), normal);
         if(0) Trans->Jacobian().Print(cout);
         CondRightArea += normal.Norml2();
         if (0) cout << "CondRightArea = " << CondRightArea << endl;
      }
   }
   CondLeftArea /= 2.0; CondRightArea /= 2.0; // the result is twice the area.
   // Compute current density in each conductor.
   double CondLeftJ=1.0/CondLeftArea;
   double CondRightJ=-1.0/CondRightArea;
   cout << "CondLeftArea = " << CondLeftArea << endl;
   cout << "CondRightArea = " << CondRightArea << endl;

   // setup a VectorRestrictedCoefficient to define the current density
   // in the trace.
  /*
   Array<int> AttrArray(5);
   //i put 1 marker all over for debug.
   AttrArray=0; AttrArray[4]=1;  //attribute 4 is the tracesurface.
   // i put 1.0 in x y and z just to debug.
   double ArrayJCoeff[]={0.0, 0.0, 1.0};  //current density vector in 3D.
   Vector VectJCoeff((double *)ArrayJCoeff, 3);
   VectorConstantCoefficient VCCJ(VectJCoeff);
   VectorRestrictedCoefficient VRCJ(VCCJ, AttrArray);
*/

   Vector zeroVec(dim); zeroVec = 0.0;
   VectorConstantCoefficient zeroVecCoef(zeroVec);
   Vector oneVec(dim); oneVec = 1.0;
   VectorConstantCoefficient oneVecCoef(oneVec);
   
   
// 6a. Determine the list of conductors boundaries.
// clb: conductor left back, f: front, r: right.
   Array<int> clb_brd_list, clb_bdr;
   Array<int> crb_brd_list, crb_bdr;
   Array<int> clf_brd_list, clf_bdr;
   Array<int> crf_brd_list, crf_bdr;
   if (mesh->bdr_attributes.Size())
   {
      clb_bdr.SetSize(mesh->bdr_attributes.Max());
      crb_bdr.SetSize(mesh->bdr_attributes.Max());
      clf_bdr.SetSize(mesh->bdr_attributes.Max());
      crf_bdr.SetSize(mesh->bdr_attributes.Max());
      assert(clb_bdr.Size() == 10);
      clb_bdr = 0;
      clb_bdr[conductorleftsurfback-1] = 1;
      crb_bdr = 0;
      crb_bdr[conductorrightsurfback-1] = 1;
      clf_bdr = 0;
      clf_bdr[conductorleftsurffront-1] = 1;
      crf_bdr = 0;
      crf_bdr[conductorrightsurffront-1] = 1;
   }

   

   // 6b. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined based on the type
   //    of mesh and the problem type.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   cout << "mesh->bdr_attributes.Size() = " << mesh->bdr_attributes.Size() << endl;
   cout << "mesh->bdr_attributes.Max() = " << mesh->bdr_attributes.Max() << endl;
   
   if (mesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      assert(ess_bdr.Size() == 10);
      ess_bdr = 0;
      ess_bdr[airsurffront-1] = 1;
      ess_bdr[airsurfback-1] = 1;
      ess_bdr[airsurfaround-1] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   //dl240529 
   //shall be zero
   u.ProjectBdrCoefficientTangent(zeroVecCoef, zeroVecCoef, ess_bdr);

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   ComplexLinearForm b(fespace, conv);
   b.Vector::operator=(0.0);
   
   
   
   Vector vp(3); vp[0]=1.0; vp[1]=1.0;vp[2]=1.0; VectorConstantCoefficient Coefp(vp);
   Vector vn(3); vn[0]=1.0; vn[1]=1.0;vn[2]=-1.0; VectorConstantCoefficient Coefn(vn);
   
   b.AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(Coefp),
                           new VectorFEBoundaryTangentLFIntegrator(zeroVecCoef),
                           clb_bdr);
   b.AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(Coefn),
                           new VectorFEBoundaryTangentLFIntegrator(zeroVecCoef),
                           clf_bdr);
   b.AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(Coefn),
                           new VectorFEBoundaryTangentLFIntegrator(zeroVecCoef),
                           crb_bdr);
   b.AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(Coefp),
                           new VectorFEBoundaryTangentLFIntegrator(zeroVecCoef),
                           crf_bdr);
   
  
 
   // 9. Set up the sesquilinear form a(.,.) on the finite element space
   //    corresponding to the damped harmonic oscillator operator of the
   //    appropriate type:
   //
   //    0) A scalar H1 field
   //       -Div(a Grad) - omega^2 b + i omega c
   //
   //    1) A vector H(Curl) field
   //       Curl(a Curl) - omega^2 b + i omega c
   //
   //    2) A vector H(Div) field
   //       -Grad(a Div) - omega^2 b + i omega c
   //

   ConstantCoefficient a_Coef(1.0);
   //b and c coefficient are dependant on the element attributes.
   // 0: copper, 1: copper, 2: air.
   int NbrAttr = mesh->attributes.Size();
   assert(NbrAttr == 3);

   Vector b_Vector(NbrAttr);
   b_Vector[2] = omega_ * omega_ * mu_ * epsilon_;  //copper.
   b_Vector[1] = omega_ * omega_ * mu_ * epsilon_;  //copper.
   b_Vector[0] = omega_ * omega_ * mu0_ * epsilon0_;  //air.
   PWConstCoefficient b_Coef(b_Vector);

   Vector c_Vector(NbrAttr);
   c_Vector[2] = -omega_ * mu_ * sigma_;   //copper.
   c_Vector[1] = -omega_ * mu_ * sigma_;   //copper.
   c_Vector[0] = -omega_ * mu_ * epsilon0_;      //air.
   PWConstCoefficient c_Coef(c_Vector);

   SesquilinearForm *a = new SesquilinearForm(fespace, conv);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
         a->AddDomainIntegrator(new CurlCurlIntegrator(a_Coef),
                                NULL);
         a->AddDomainIntegrator(new VectorFEMassIntegrator(b_Coef),
                                new VectorFEMassIntegrator(c_Coef));

   // 9a. Set up the bilinear form for the preconditioner corresponding to the
   //     appropriate operator
   //
   //      0) A scalar H1 field
   //         -Div(a Grad) - omega^2 b + omega c
   //
   //      1) A vector H(Curl) field
   //         Curl(a Curl) + omega^2 b + omega c
   //
   //      2) A vector H(Div) field
   //         -Grad(a Div) - omega^2 b + omega c
   //
   BilinearForm *pcOp = new BilinearForm(fespace);
   if (pa) { pcOp->SetAssemblyLevel(AssemblyLevel::PARTIAL); }

         pcOp->AddDomainIntegrator(new CurlCurlIntegrator(a_Coef));
         pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(b_Coef));
         pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(c_Coef));

   // 10. Assemble the form and the corresponding linear system, applying any
   //     necessary transformations such as: assembly, eliminating boundary
   //     conditions, conforming constraints for non-conforming AMR, etc.
   a->Assemble();
   pcOp->Assemble();

   OperatorHandle A;
   Vector B, U;

   a->FormLinearSystem(ess_tdof_list, u, b, A, U, B);

   cout << "Size of linear system: " << A->Width() << endl << endl;

   // 11. Define and apply a GMRES solver for AU=B with a block diagonal
   //     preconditioner based on the appropriate sparse smoother.
   {
      Array<int> blockOffsets;
      blockOffsets.SetSize(3);
      blockOffsets[0] = 0;
      blockOffsets[1] = A->Height() / 2;
      blockOffsets[2] = A->Height() / 2;
      blockOffsets.PartialSum();

      BlockDiagonalPreconditioner BDP(blockOffsets);

      Operator * pc_r = NULL;
      Operator * pc_i = NULL;

      if (pa)
      {
         pc_r = new OperatorJacobiSmoother(*pcOp, ess_tdof_list);
      }
      else
      {
         OperatorHandle PCOp;
         pcOp->SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
         pcOp->FormSystemMatrix(ess_tdof_list, PCOp);
        pc_r = new DSmoother(*PCOp.As<SparseMatrix>());
      }
      double s = 1.0;
      pc_i = new ScaledOperator(pc_r,
                                (conv == ComplexOperator::HERMITIAN) ?
                                s:-s);

      BDP.SetDiagonalBlock(0, pc_r);
      BDP.SetDiagonalBlock(1, pc_i);
      BDP.owns_blocks = 1;

      GMRESSolver gmres;
      gmres.SetPreconditioner(BDP);
      gmres.SetOperator(*A.Ptr());
      gmres.SetRelTol(1e-12);
      gmres.SetMaxIter(1000);
      gmres.SetPrintLevel(1);
      gmres.Mult(B, U);
   }

   // 12. Recover the solution as a finite element grid function and compute the
   //     errors if the exact solution is known.
   a->RecoverFEMSolution(U, b, u);


   // 13. Save the refined mesh and the solution. This output can be viewed
   //     later using GLVis: "glvis -m mesh -g sol".
   {
      ofstream mesh_ofs("refined.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      ofstream sol_r_ofs("sol_r.gf");
      ofstream sol_i_ofs("sol_i.gf");
      sol_r_ofs.precision(8);
      sol_i_ofs.precision(8);
      u.real().Save(sol_r_ofs);
      u.imag().Save(sol_i_ofs);
   }

   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock_r(vishost, visport);
      socketstream sol_sock_i(vishost, visport);
      sol_sock_r.precision(8);
      sol_sock_i.precision(8);
      sol_sock_r << "solution\n" << *mesh << u.real()
                 << "window_title 'Solution: Real Part'" << flush;
      sol_sock_i << "solution\n" << *mesh << u.imag()
                 << "window_title 'Solution: Imaginary Part'" << flush;
   }
  

   // 15. Free the used memory.
   delete a;
   
   delete pcOp;
   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}

