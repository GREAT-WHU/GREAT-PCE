
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include "gutils/gmatrixconv.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gmodels/geop.h"

using namespace std;

namespace gnut {   

//Constructor
t_geop::t_geop()
{
}

//Destructor
t_geop::~t_geop()
{
}

// Nutation Matrix 
// --------------------------------------
Matrix t_geop::nutMatrix(double mjd)
{
//  boost::mutex::scoped_lock lock(_mutex);      

  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
     
  return I;
 
}

// Precession Matrix
// ------------------------------------------
Matrix t_geop::precMatrix (double mjd_1)
{
//  boost::mutex::scoped_lock lock(_mutex);
   
  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
     
  return I;   
}

// Normalize angle into interval 0 - 2pi
// -----------------------------------
double t_geop::_normangle(double x)
{
   double norm;
   norm = fmod(x, 2*G_PI);
   if (norm < 0 ) norm += 2*G_PI;
   
   return norm;
}

} // namespace
