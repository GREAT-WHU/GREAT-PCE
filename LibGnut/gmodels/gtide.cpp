
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gutils/gsysconv.h"
#include "gmodels/gtide.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gtide::t_gtide(t_glog* l)
  : _gotl(0),
    _log(l)
{}


// destructor
// ----------
t_gtide::~t_gtide(){}


// solid earth tides
// ----------
t_gtriple t_gtide::tide_searth(const t_gtime& epoch, t_gtriple& crd)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _mutex.lock();
  t_gtriple dxyz(0.0,0.0,0.0);
  _mutex.unlock(); return dxyz;
}


// pole tides
// ----------
t_gtriple t_gtide::tide_pole()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _mutex.lock();
  t_gtriple dxyz(0.0,0.0,0.0);
  _mutex.unlock(); return dxyz;
}


// ocean tide loading
// ----------
t_gtriple t_gtide::load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _mutex.lock();
  t_gtriple dxyz(0.0,0.0,0.0);
  _mutex.unlock(); return dxyz;
}


// atmospheric tide loading
// ----------
t_gtriple t_gtide::load_atmosph()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _mutex.lock();
  t_gtriple dxyz(0.0,0.0,0.0);  
  _mutex.unlock();
  return dxyz;
}

// set OTL pointer
// -----------------
void t_gtide::setOTL(t_gallotl* gallotl)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _mutex.lock();
  _gotl = gallotl;  
  _mutex.unlock();
   return;

}

} // namespace
