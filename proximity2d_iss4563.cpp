//                                proximity2d
//                                inspired from skineffect2d
//                                based on MFEM Example 22 prob 1 (case 0)
//
// Compile with: make proximity2d, need MFEM version 4.7 and GLVIS-4.3.
//
// Sample runs:  ./proximity2d -m ../mesh/ProxRectRoundWires2d.msh
//
/*
Description:  

-1/u Div(grad(Az)) - e w^2 Az + i w s Az= J0z
Jinduced = -i w s Az
Jtot = J0 + Jinduced, note J0 is real, Jind is imaginary (90 degree)

a (-div(grad((Az))) - b Az + i c Az= J0z



u: permeability.
e: permitivity.
s: conductivity.
i: sqrt(-1)
*/

#include "mfem.hpp"
#include <fstream>
#include <iostream>

#define wire_1 1
#define wire_2 2
#define air 3
#define aircontour 4

using namespace std;
using namespace mfem;

double IntScalar3(FiniteElementSpace &fes, Coefficient &coeff, int Attr)
{
   QuadratureSpace qs(fes.GetMesh(), 2*fes.GetMaxElementOrder());
   Array<int> attrs;
   attrs=0;
   if (fes.GetMesh()->attributes.Size())
   {
      attrs.SetSize(fes.GetMesh()->attributes.Max());
      attrs[Attr-1] = 1;
   }
   RestrictedCoefficient restr_coeff(coeff, attrs);
   return qs.Integrate(restr_coeff);
}

double IntScalar4(FiniteElementSpace &fes, Coefficient &coeff, int Attr)
{
   QuadratureSpace qs(fes.GetMesh(), 2*fes.GetMaxElementOrder());
   Array<int> attrs(2); attrs=0;
   assert(fes.GetMesh()->attributes.Max() >= Attr);
   attrs[0] = Attr;
   RestrictedCoefficient restr_coeff(coeff, attrs);
   return qs.Integrate(restr_coeff);
}

/// Function calculating the integral of a scalar Coefficient, Coeff, over a domain
/// identified by an Attribute, Attr, on the FiniteElementSpace fes.
double IntegrateScalar(FiniteElementSpace &fes, Coefficient &coeff, int Attr)
{
   double integral_value = 0.0;
   for (int i = 0; i < fes.GetMesh()->GetNE(); i++)  // Loop over elements
   {
      if(fes.GetAttribute(i) == Attr)
      {   
         ElementTransformation *trans = fes.GetMesh()->GetElementTransformation(i);
         const FiniteElement &fe = *(fes.GetFE(i));
         // Use a quadrature rule that matches the element's order
         const IntegrationRule &ir = IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder());
         for (int j = 0; j < ir.GetNPoints(); j++)  // Loop over quadrature points
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            trans->SetIntPoint(&ip);
            // Evaluate scalar function at the quadrature point in physical coordinates
            double scalar_value = coeff.Eval(*trans, ip);
            // Accumulate the integral (scalar value * weight * Jacobian determinant)
            integral_value += scalar_value * ip.weight * trans->Weight();
         }
      }
   }
   return integral_value;
}

   void Glvis(Mesh *m, GridFunction *gf, string title, int precision = 8, string keys=" keys 'mmc'")
      {
         char vishost[] = "localhost";
         int  visport   = 19916;
         socketstream sol_sock(vishost, visport);
         sol_sock.precision(precision);
         sol_sock << "solution\n" << *m << *gf
                 << "window_title '" + title +"'"
                 << keys << flush;   
      }


static double mu_ = 1.257E-6;
static double epsilon_ = 8.854E-12;
static double sigma_ = 58E6;
static double omega_ = 2.0*M_PI*60;

int main(int argc, char *argv[])
{
   tic();
   // 1. Parse command-line options.
   const char *mesh_file = "mesh/proxrectroundwires2d.msh";
   int ref_levels = 0;
   int order = 1;
   int prob = 0;
   double freq = -1.0;
   double a_coef = 0.0;
   bool visualization = 1;
   bool herm_conv = true;
   bool pa = false;
   bool noDispl = false;
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
   args.AddOption(&noDispl, "-nd", "--no-displ", "-dincl", "--displ-include",
                  "Do not include displacement current in the computation.");
                  
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   if ( freq >= 0.0 )
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
   Mesh *mesh = new Mesh(mesh_file, 1, 1); // 
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
   assert(mesh->bdr_attributes.Max()==4);
   ess_bdr.SetSize(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[aircontour-1]=1;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   

   // current density in conductor left and conductor right.
   PWConstCoefficient *CurrentDensityCoeff = new PWConstCoefficient;
   {
      double CoeffArray[]={1.0, -1.0, 0.0, 0.0};  // wire_1 1A/m, wire_2 -1A/m.
      Vector CoeffVector(CoeffArray, 4);
      CurrentDensityCoeff->UpdateConstants(CoeffVector);
   }
   
   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   ComplexLinearForm b(fespace, conv);
   b.Vector::operator=(0.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(*CurrentDensityCoeff), nullptr);
   b.Assemble();

   

   // 8. Define the solution vector u as a complex finite element grid function
   //    corresponding to fespace. Initialize u with initial guess of 1+0i or
   //    the exact solution if it is known.
   ComplexGridFunction u(fespace);
   //u=0.0;
      
   // 9. Set up the sesquilinear form a(.,.) 
   //
   // Wrote this way to show that diffusion integrator is negative, -div(grad)).
   ConstantCoefficient a_(-(-1.0/mu_));
   ConstantCoefficient b_(-omega_ * omega_ * epsilon_);

   PWConstCoefficient *c_ = new PWConstCoefficient;
   {
      double CoeffArray[]={omega_*sigma_, omega_*sigma_, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 4);
      c_->UpdateConstants(CoeffVector);
   }
   
   PWConstCoefficient *negc_ = new PWConstCoefficient;
   {
      double CoeffArray[]={-omega_*sigma_, -omega_*sigma_, 0.0, 0.0};
      Vector CoeffVector(CoeffArray, 4);
      negc_->UpdateConstants(CoeffVector);
   }

   SesquilinearForm *a = new SesquilinearForm(fespace, conv);
   if (pa) { a->SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   
         a->AddDomainIntegrator(new DiffusionIntegrator(a_),
                                NULL);
         a->AddDomainIntegrator(noDispl? nullptr : new MassIntegrator(b_),
                                new MassIntegrator(*c_));

   // 9a. Set up the bilinear form for the preconditioner 
   //
   BilinearForm *pcOp = new BilinearForm(fespace);   
   pcOp->AddDomainIntegrator(new DiffusionIntegrator(a_));
   if (!noDispl) pcOp->AddDomainIntegrator(new MassIntegrator(b_));
   pcOp->AddDomainIntegrator(new MassIntegrator(*c_));     

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
      pc_r = new DSmoother(*PCOp.As<SparseMatrix>());
               
      
      
      pc_i = new ScaledOperator(pc_r,
                                (conv == ComplexOperator::HERMITIAN) ?
                                1:-1);

      BDP.SetDiagonalBlock(0, pc_r);
      BDP.SetDiagonalBlock(1, pc_i);
      BDP.owns_blocks = 1;

      GMRESSolver gmres;
      gmres.SetPreconditioner(BDP);
      gmres.SetOperator(*A.Ptr());
      gmres.SetRelTol(1e-12);
      gmres.SetMaxIter(2000);
      gmres.SetPrintLevel(1);
      gmres.Mult(B, U);
   }

   // 12. Recover the solution as a finite element grid function and compute the
   //     errors if the exact solution is known.
   a->RecoverFEMSolution(U, b, u);



// see https://github.com/mfem/mfem/issues/4543

   ConstantCoefficient Zero(0.0);
   GridFunctionCoefficient ACoeffReal(&(u.real())), ACoeffImag(&(u.imag()));
   ProductCoefficient JindCoeffReal(ACoeffImag, *c_);
   ProductCoefficient JindCoeffImag(ACoeffReal, *negc_);
   PowerCoefficient JindRealSquare(JindCoeffReal, 2.0), JindImagSquare(JindCoeffImag, 2.0);
   SumCoefficient JindSquare(JindRealSquare, JindImagSquare);
   PowerCoefficient JindCoeff(JindSquare, 0.5);

   SumCoefficient JtotCoeffReal(*CurrentDensityCoeff, JindCoeffReal);
   SumCoefficient JtotCoeffImag(Zero, JindCoeffImag);
   PowerCoefficient JtotRealSquare(JtotCoeffReal, 2.0), JtotImagSquare(JtotCoeffImag, 2.0);
   SumCoefficient JtotSquare(JtotRealSquare, JtotImagSquare);
   PowerCoefficient JtotCoeff(JtotSquare, 0.5);
   
   
   FiniteElementCollection *JIndFec = new DG_FECollection(order, dim);
   FiniteElementSpace *JIndFESpace = new FiniteElementSpace(mesh, JIndFec);

   GridFunction J0r(JIndFESpace);
   J0r.ProjectCoefficient(*CurrentDensityCoeff);

   GridFunction JIndGridFunctionReal(JIndFESpace), JIndGridFunctionImag(JIndFESpace);
   JIndGridFunctionReal.ProjectCoefficient(JindCoeffReal);
   JIndGridFunctionImag.ProjectCoefficient(JindCoeffImag);
   GridFunction JindGridFunction(JIndFESpace);
   JindGridFunction.ProjectCoefficient(JindCoeff);


  GridFunction Jtr(JIndFESpace), Jti(JIndFESpace);
  Jtr.ProjectCoefficient(JtotCoeffReal);
  Jti.ProjectCoefficient(JindCoeffImag);
 
   GridFunction Jtot(JIndFESpace);
   Jtot.ProjectCoefficient(JtotCoeff); 

   

   Glvis(mesh, &J0r, "J0r");
   
   Glvis(mesh, &u.real(), "A-field: Real Part" );
   Glvis(mesh, &u.imag(), "A-field: Imag Part" );

   Glvis(mesh, &JIndGridFunctionReal, "JIndReal" );
   Glvis(mesh, &JIndGridFunctionImag, "JIndImag" );
   Glvis(mesh, &JindGridFunction, "JIndtot" );   

   Glvis(mesh, &Jtr, "Jtot: Real Part" );
   Glvis(mesh, &Jti, "Jtot: Imaginary Part" );
   Glvis(mesh, &Jtot, "Jtot" );


cout << "\n ****************************** \n";
cout << "\n IntegrateScalar J0r 2 = " << IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_1) << endl;
cout << "\n IntegrateScalar J0r 3 = " << IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_2) << endl;
cout << "\n IntegrateScalar J0r 1 = " << IntegrateScalar(*fespace, *CurrentDensityCoeff, air) << endl;

cout << "\n IntegrateScalar Jtot 2 = " << IntegrateScalar(*fespace, JtotCoeff, wire_1) << endl;
cout << "\n IntegrateScalar Jtot 3 = " << IntegrateScalar(*fespace, JtotCoeff, wire_2) << endl;
cout << "\n IntegrateScalar Jtot 1 = " << IntegrateScalar(*fespace, JtotCoeff, air) << endl;

cout << "\n Ratio wire_1 = " << IntegrateScalar(*fespace, JtotCoeff, wire_1)/IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_1) << endl;
cout << "\n Ratio wire_2 = " << IntegrateScalar(*fespace, JtotCoeff, wire_2)/IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_2) << endl;

cout << "\n Ratio wire_1 real= " << IntegrateScalar(*fespace, JtotCoeffReal, wire_1)/IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_1) << endl;
cout << "\n Ratio wire_1 imag= " << IntegrateScalar(*fespace, JtotCoeffImag, wire_1)/IntegrateScalar(*fespace, *CurrentDensityCoeff, wire_1) << endl;

cout << "\n ****************************** \n";
cout << "\n IntScalar3 J0r 1 = " << IntScalar3(*fespace, *CurrentDensityCoeff, wire_1) << endl;
cout << "\n IntScalar3 J0r 2 = " << IntScalar3(*fespace, *CurrentDensityCoeff, wire_2) << endl;
cout << "\n IntScalar3 J0r air = " << IntScalar3(*fespace, *CurrentDensityCoeff, air) << endl;




   // 15. Free the used memory.
   delete a;
   delete pcOp;
   delete fespace;
   delete fec;
   delete mesh;

   cout << "\n time elapsed = " << toc() << endl;

   return 0;
}
