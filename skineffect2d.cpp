//                                skineffect2d
//                                based on MFEM Example 22 prob 1 (case 0)
//
// Compile with: make skineffect2d
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

static double mu_ = 1.257e-6;
static double epsilon_ = 8.854E-12;
static double sigma_ = 1.0/16.78e-9;
static double omega_ = 2.0*M_PI*60;

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "roundwire.msh";
   int ref_levels = 0;
   int order = 1;
   int prob = 0;
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

   FiniteElementCollection *fec = NULL;
   fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined based on the type
   //    of mesh and the problem type.
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


   u.ProjectBdrCoefficient(oneCoef, zeroCoef, ess_bdr);
         
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
   ConstantCoefficient a_(-1.0);
   ConstantCoefficient b_(-omega_ * omega_ * epsilon_ * mu_);
   ConstantCoefficient c_(omega_ * mu_ * sigma_);
   ConstantCoefficient negc_(omega_ * omega_ * epsilon_* mu_);

   SesquilinearForm *a = new SesquilinearForm(fespace, conv);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   
         a->AddDomainIntegrator(new DiffusionIntegrator(a_),
                                NULL);
         a->AddDomainIntegrator(new MassIntegrator(b_),
                                new MassIntegrator(c_));
         

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

   
         pcOp->AddDomainIntegrator(new DiffusionIntegrator(a_));
         pcOp->AddDomainIntegrator(new MassIntegrator(b_));
         pcOp->AddDomainIntegrator(new MassIntegrator(c_));
         

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
      double s = (prob != 1) ? 1.0 : -1.0;
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

bool check_for_inline_mesh(const char * mesh_file)
{
   string file(mesh_file);
   size_t p0 = file.find_last_of("/");
   string s0 = file.substr((p0==string::npos)?0:(p0+1),7);
   return s0 == "inline-";
}

complex<double> u0_exact(const Vector &x)
{
   int dim = x.Size();
   complex<double> i(0.0, 1.0);
   complex<double> alpha = (epsilon_ * omega_ - i * sigma_);
   complex<double> kappa = std::sqrt(mu_ * omega_* alpha);
   return std::exp(-i * kappa * x[dim - 1]);
}

double u0_real_exact(const Vector &x)
{
   return u0_exact(x).real();
}

double u0_imag_exact(const Vector &x)
{
   return u0_exact(x).imag();
}

void u1_real_exact(const Vector &x, Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[0] = u0_real_exact(x);
}

void u1_imag_exact(const Vector &x, Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[0] = u0_imag_exact(x);
}

void u2_real_exact(const Vector &x, Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[dim-1] = u0_real_exact(x);
}

void u2_imag_exact(const Vector &x, Vector &v)
{
   int dim = x.Size();
   v.SetSize(dim); v = 0.0; v[dim-1] = u0_imag_exact(x);
}
