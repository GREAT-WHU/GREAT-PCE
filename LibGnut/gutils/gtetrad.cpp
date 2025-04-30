
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <cmath>
#include <iomanip>

#include "gutils/gtetrad.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// null constructor
// ----------
t_gtetrad::t_gtetrad()
{
  _crd[0] = _crd[1] = _crd[2] = _crd[3] = 0.0;
}


// constructor
// ----------
t_gtetrad::t_gtetrad(double x, double y, double z, double t)
{
  _crd[0] = x;
  _crd[1] = y;
  _crd[2] = z;
  _crd[3] = t;
}


// constructor
// ----------
t_gtetrad::t_gtetrad(double crd[])
{
  _crd[0] = crd[0];
  _crd[1] = crd[1];
  _crd[2] = crd[2];
  _crd[3] = crd[3];
}


// constructor
// ----------
t_gtetrad::t_gtetrad(const ColumnVector& crd)
{
  _crd[0] = crd(1);
  _crd[1] = crd(2);
  _crd[2] = crd(3);
  _crd[3] = crd(4);
}


// destructor
// ----------
t_gtetrad::~t_gtetrad(){}


// get a reference of element
// ----------
double& t_gtetrad::operator[](const size_t idx)
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad);
  if( idx > 3 ){
    cerr << "Not valid tetrad index [used 0]\n";
    return _crd[0];
  }
  return _crd[idx];
}


// get a value of element
// ----------
double t_gtetrad::operator[](const size_t idx) const
{
//  boost::mutex::scoped_lock lock(
//  _mutex_tetrad);
  if( idx < 4 ) return _crd[idx];
   
  return 0.0;
}

// operator +
// ------------------------
t_gtetrad t_gtetrad::operator+(const t_gtetrad& other) const
{
  t_gtetrad tmp(*this);
  tmp[0] += other[0];
  tmp[1] += other[1];
  tmp[2] += other[2];   
  tmp[3] += other[3];   

  return tmp;   
}

   
// operator for sorting
// ----------
bool t_gtetrad::operator<(const t_gtetrad& tr) const
{
  return ( ( _crd[0]  < tr.crd(0) ) ||
           ( _crd[0] == tr.crd(0) && _crd[1]  < tr.crd(1) ) ||
           ( _crd[0] == tr.crd(0) && _crd[1] == tr.crd(1) && _crd[2]  < tr.crd(2) ) ||
           ( _crd[0] == tr.crd(0) && _crd[1] == tr.crd(1) && _crd[2]  < tr.crd(2) && _crd[3]  < tr.crd(3) ) );
}


// get single element
// ----------
double t_gtetrad::crd(int idx) const
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad);
  if( idx >= 0 && idx < 4 ) return _crd[static_cast<unsigned int>(idx)];
   
  return 0.0;
}


// set single element
// ------------------------
void t_gtetrad::set(int idx, double newValue)
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad);
  if( idx >= 0 && idx < 4 ) _crd[static_cast<unsigned int>(idx)] = newValue;
     
}


// copy operator
// ----------
t_gtetrad& t_gtetrad::operator=(const t_gtetrad& other)
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad);
  if( this != &other ){
    _crd[0] = other.crd(0);
    _crd[1] = other.crd(1);
    _crd[2] = other.crd(2);
    _crd[3] = other.crd(3);
  }
  return *this;
}

// equal operator
// ----------
bool t_gtetrad::operator==(const t_gtetrad& tr) const
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad);
  return ( _crd[0] == tr.crd(0) &&
           _crd[1] == tr.crd(1) &&
           _crd[2] == tr.crd(2) &&
           _crd[3] == tr.crd(3)     );
}


// get array
// ----------
double* t_gtetrad::crd_array()
{
//  boost::mutex::scoped_lock lock(_mutex_tetrad); 
  return _crd;
}


// get ColumnVector[3]
// ----------
ColumnVector t_gtetrad::crd_cvect()
{
 ColumnVector tmp(3);
 tmp(1) = _crd[0];
 tmp(2) = _crd[1];
 tmp(3) = _crd[2];
 tmp(4) = _crd[3];
 return tmp;
}


// get tetrad
// ----------
t_gtetrad& t_gtetrad::crd_tetrad()
{
 return *this;
}

// set array by ColumnVector
// ------------------------------
void t_gtetrad::set(const ColumnVector& crd)
{
  _crd[0] = crd(1);
  _crd[1] = crd(2);
  _crd[2] = crd(3);   
  _crd[3] = crd(4);   
}

// set array by array
// ------------------------------
void t_gtetrad::set(double crd[4])
{
  _crd[0] = crd[0];
  _crd[1] = crd[1];
  _crd[2] = crd[2];
  _crd[3] = crd[3];
}

// get unit ColumnVector
// -----------------------------
ColumnVector t_gtetrad::unitary()
{
  ColumnVector tmp(4);
  tmp = this->crd_cvect();
  double s = tmp.norm_Frobenius();
  tmp /= s;
   
  return tmp;
   
}

// overloading << operator
// -----------------------------
ostream& operator<<(ostream& os, const t_gtetrad& x)
{
   os << fixed << setprecision(3) 
      << dbl2str(x[0]) + " "
       + dbl2str(x[1]) + " " 
       + dbl2str(x[2]) + " " 
       + dbl2str(x[3]);
   return os;
}

// test if all elements are zero
// -----------------------------
bool t_gtetrad::zero()
{
   if ( double_eq(_crd[0], 0.0) &&
        double_eq(_crd[1], 0.0) &&
        double_eq(_crd[2], 0.0) &&
        double_eq(_crd[3], 0.0) ) return true;
   else return false;
}

} // namespace

