//                                skineffect2d.cpp, inspired from MFEM Example 22
//
// Compile with: make ex22
//
// Sample runs:  ex22 -m ../data/inline-segment.mesh -o 3
//               ex22 -m ../data/inline-tri.mesh -o 3
//               ex22 -m ../data/inline-quad.mesh -o 3
//               ex22 -m ../data/inline-quad.mesh -o 3 -p 1
//               ex22 -m ../data/inline-quad.mesh -o 3 -p 1 -pa
//               ex22 -m ../data/inline-quad.mesh -o 3 -p 2
//               ex22 -m ../data/inline-tet.mesh -o 2
//               ex22 -m ../data/inline-hex.mesh -o 2
//               ex22 -m ../data/inline-hex.mesh -o 2 -p 1
//               ex22 -m ../data/inline-hex.mesh -o 2 -p 2
//               ex22 -m ../data/inline-hex.mesh -o 2 -p 2 -pa
//               ex22 -m ../data/inline-wedge.mesh -o 1
//               ex22 -m ../data/inline-pyramid.mesh -o 1
//               ex22 -m ../data/star.mesh -r 1 -o 2 -sigma 10.0
//

// Device sample runs:
//               ex22 -m ../data/inline-quad.mesh -o 3 -p 1 -pa -d cuda
//               ex22 -m ../data/inline-hex.mesh -o 2 -p 2 -pa -d cuda
//               ex22 -m ../data/star.mesh -r 1 -o 2 -sigma 10.0 -pa -d cuda
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
#include <fstream>
#include <iostream>


using namespace std;
using namespace mfem;

typedef double real_t;

static real_t mu_ = 1.0;
static real_t epsilon_ = 1.0;
static real_t sigma_ = 20.0;
static real_t omega_ = 1.0;

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "star.mesh";
   int ref_levels = 0;
   int order = 1;
   real_t freq = -1.0;
   real_t a_coef = 0.0;
   bool visualization = 1;
   bool herm_conv = true;
   const char *device_config = "cpu";

   OptionsParser args(argc, argv);
   args.AddOption(&omega_, "-w", "--omega",
                  "Omega used.");
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
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   if ( a_coef != 0.0 )
   {
      mu_ = 1.0 / a_coef;
   }
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

   cout << "mesh->GetNE() before = " << mesh->GetNE() << endl;

   // 4. Refine the mesh to increase resolution. In this example we do
   //    'ref_levels' of uniform refinement where the user specifies
   //    the number of levels with the '-r' option.
   for (int l = 0; l < ref_levels; l++)
   {
      mesh->UniformRefinement();
   }
   cout << "mesh->GetNE() after = " << mesh->GetNE() << endl;

getchar();
   

   FiniteElementCollection *fec = NULL;
   
   fec = new ND_FECollection(order, dim);

   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined based on the type
   //    of mesh and the problem type.

cout << "mesh->bdr_attributes.Size() =" << mesh->bdr_attributes.Size() << endl;

   Array<int> ess_tdof_list;
   Array<int> ess_bdr;
   if (mesh->bdr_attributes.Size())
   {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   ComplexLinearForm b(fespace, conv);
   b.Vector::operator=(0.0);

   // 8. Define the solution vector u as a complex finite element grid function
   //    corresponding to fespace. Initialize u with initial guess of 1+0i or
   //    the exact solution if it is known.
   ComplexGridFunction u(fespace);
   

   ConstantCoefficient zeroCoef(0.0);
   ConstantCoefficient oneCoef(1.0);

   Vector zeroVec(dim); zeroVec = 0.0;
   Vector  oneVec(dim);  oneVec = 0.0; oneVec[0] = 1.0;
   VectorConstantCoefficient zeroVecCoef(zeroVec);
   VectorConstantCoefficient oneVecCoef(oneVec);

   
u.ProjectBdrCoefficientTangent(oneVecCoef, zeroVecCoef, ess_bdr);
  
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
   ConstantCoefficient stiffnessCoef(1.0/mu_);
   ConstantCoefficient massCoef(-omega_ * omega_ * epsilon_);
   ConstantCoefficient lossCoef(omega_ * sigma_);
   ConstantCoefficient negMassCoef(omega_ * omega_ * epsilon_);

   SesquilinearForm *a = new SesquilinearForm(fespace, conv);
   a->AddDomainIntegrator(new CurlCurlIntegrator(stiffnessCoef),
                                NULL);
   a->AddDomainIntegrator(new VectorFEMassIntegrator(massCoef),
                                new VectorFEMassIntegrator(lossCoef));

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
   pcOp->AddDomainIntegrator(new CurlCurlIntegrator(stiffnessCoef));
   pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(negMassCoef));
   pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(lossCoef));
      
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

      OperatorHandle PCOp;
      pcOp->SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
      pcOp->FormSystemMatrix(ess_tdof_list, PCOp);
      pc_r = new GSSmoother(*PCOp.As<SparseMatrix>());
               
      

      pc_i = new ScaledOperator(pc_r,
                  (conv == ComplexOperator::HERMITIAN) ? -1.0:1.0);

      BDP.SetDiagonalBlock(0, pc_r);
      BDP.SetDiagonalBlock(1, pc_i);
      BDP.owns_blocks = 1;

      GMRESSolver gmres;
      gmres.SetPreconditioner(BDP); //works fine even if not called.
      gmres.SetOperator(*A.Ptr());
      gmres.SetRelTol(1e-12);
      gmres.SetMaxIter(1000);
      gmres.SetPrintLevel(1);
      gmres.Mult(B, U);
   }

   // 12. Recover the solution as a finite element grid function.
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
                  << "window_title 'Solution: Real Part'" << " keys 'c'" << flush;
      sol_sock_i << "solution\n" << *mesh << u.imag()
                 << "window_title 'Solution: Imaginary Part'" << " keys 'c'" << flush;
   }
   
   

   // 15. Free the used memory.
   delete a;
   delete pcOp;
   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}
