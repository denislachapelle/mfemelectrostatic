/*    Written by Denis Lachapelle December 2023.
                            MFEM
//
// Compile with: make intgrad
//
// This code is about microstrip simulation.
// 
// DL231229.
//
*/

#include <mfem.hpp>
#include <miniapps/common/fem_extras.hpp>
#include <fstream>
#include <iostream>


using namespace std;
using namespace mfem;