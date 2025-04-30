
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include "gall/gallotl.h"
#include "gutils/gconst.h"

#include <math.h>
namespace gnut {  

//Constructor
t_gallotl::t_gallotl()
{
//   cout << "vytvarim t_gallotl" << endl;
   id_type(  t_gdata::ALLOTL );
}

//Destructor
t_gallotl::~t_gallotl()
{
   _mapotl.clear();
}

int t_gallotl::data(Matrix& otldata, const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  if( _mapotl.find(site) == _mapotl.end() || _mapotl.size() == 0 ){
    string site_short = site.substr(0,4);
    if( _mapotl.find(site_short) == _mapotl.end() || _mapotl.size() == 0 ){    
      _gmutex.unlock(); return -1;
    }else otldata = _mapotl[site_short].data();
  }else otldata = _mapotl[site].data();

  _gmutex.unlock(); return 1;
}

int t_gallotl::data(Matrix& otldata, double lon,double lat)
{
#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif


	const double dlon_eps = 1E4 / 6378235.0 * R2D, dlat_eps = dlon_eps / cos(lat*D2R) ;

	if (lon > 180) lon -= 360.0;

	for (auto otl_iter = _mapotl.begin(); otl_iter != _mapotl.end(); otl_iter++)
	{
		double dlon, dlat;
		dlon = fabs(otl_iter->second.lon() - lon);
		dlat = fabs(otl_iter->second.lat() - lat);
		if (dlon < dlon_eps && dlat < dlat_eps) {
			otldata = otl_iter->second.data();
			return 1;
		}
	}

	return -1;
}

double t_gallotl::lat(const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double tmp = 0.0;
  if( _mapotl.find(site) != _mapotl.end() ) tmp = _mapotl[site].lat();
   
  _gmutex.unlock(); return tmp;
}

double t_gallotl::lon(const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  double tmp = 0.0;
  if( _mapotl.find(site) != _mapotl.end() ) tmp = _mapotl[site].lon();

  _gmutex.unlock(); return tmp;
}


// Add element to map container
// -----------------------------
void t_gallotl::add(t_gotl& otl)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  string site = otl.site();
  _mapotl[site] = otl;
   
  _gmutex.unlock(); return;
}

// Print all data
// --------------------
void t_gallotl::print()
{
   map<string, t_gotl>::iterator it;
   for (it = _mapotl.begin(); it != _mapotl.end(); it++)
   {
      t_gotl otl = it->second;
      cout << "Site: " << otl.site() << " Lon: " << otl.lon() << " Lat: " << otl.lat() << endl;
      cout << "Data: \n" << otl.data() << endl;
   }
   
}

} // namespace
