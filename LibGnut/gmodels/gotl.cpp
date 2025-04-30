
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include "gmodels/gotl.h"

namespace gnut {   

//Constructor
t_gotl::t_gotl()
{
  id_type(t_gdata::OTL);   
}

//Destructor
t_gotl::~t_gotl()
{
}

string t_gotl::site()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   string tmp = _site;
   _gmutex.unlock(); return tmp;
}

Matrix t_gotl::data()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   Matrix tmp = _data;
   _gmutex.unlock(); return tmp;
}

double t_gotl::lat()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   double tmp = _lat;
   _gmutex.unlock(); return tmp;
}

double t_gotl::lon()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   double tmp = _lon;
   _gmutex.unlock(); return tmp;
}

// Set data
// -------------
void t_gotl::setdata(const string& site, const double& lon, const double& lat, const Matrix& data)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   _site = site;
   _lon  = lon;
   _lat  = lat;
   _data = data;
   _gmutex.unlock();
   return;
}

} // namespace
