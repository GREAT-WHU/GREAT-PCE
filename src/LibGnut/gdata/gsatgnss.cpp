/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gdata/gsatgnss.h"
#include "gmodels/gephplan.h"

using namespace std;

namespace gnut {


	t_gsatgnss::t_gsatgnss()
	{
	}

	// -----------
t_gsatgnss::t_gsatgnss(t_gobsgnss* obs)
 :  _satcrd(0.0,0.0,0.0),
    _clk(0.0),
	_dclk(0.0),
    _ele(0.0),
    _azi(0.0),
    _rho(0.0),
    _eclipse(false)
{   
  gobs = obs;
  id_type(t_gdata::SATDATA);
  id_group(t_gdata::GRP_OBSERV);
  _lastEcl = FIRST_TIME;
}


// -----------
t_gsatgnss::~t_gsatgnss()
{

}


// -----------
void t_gsatgnss::addcrd( const t_gtriple& crd )
{
  gtrace("t_gsatgnss::addcrd");
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _satcrd = crd;
  _gmutex.unlock();  
  return;
}

//----------
int t_gsatgnss::addprd( t_gallnav* gnav, bool msk_health )
{
  gtrace("t_gsatgnss::addprd");
   
  if( ! gobs ) return -1;

  double P3 = gobs->P3();
  double L3 = gobs->L3();

  if( P3 == 0 || L3 == 0 || gnav == 0 )  return -1;

  string satname( gobs->sat() );
  if( satname.substr(0,1) != "G" &&
      satname.substr(0,1) != "R" )       return -1;

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  double xyz[3]  = {0.0, 0.0, 0.0};
  double vel[3]  = {0.0, 0.0, 0.0};
  double var[3]  = {0.0, 0.0, 0.0};
  double clk =    0.0;
  double dclk =   0.0;
  double clkrms = 0.0;   
   
  t_gtime epoT(t_gtime::GPS);
  double satclk = 0.0;
  double satclk2 = 1.0;

  while( fabs(satclk-satclk2) > 1.e-3/CLIGHT )
  {
    satclk2 = satclk;
    epoT = gobs->epoch() - P3/CLIGHT - satclk;
    int irc = gnav->clk( satname, epoT, &clk, &clkrms, &dclk, msk_health );

    if( irc < 0 ){

      if( _log && _log->verb() >= 2 )
       _log->comment(2, "gppp"," satelite " + satname 
		       + gobs->epoch().str_ymdhms(" clocks not calculated for epoch "));
//     cout << satname + gobs->epoch().str_ymdhms(" clocks not calculated for epoch ") << endl;
       _gmutex.unlock(); return -1;
    }
    satclk = clk;
  }

  int irc = gnav->pos( satname, epoT, xyz, var, vel, msk_health );

  if( irc < 0 ) {
    if( _log && _log->verb() >= 2 )
       _log->comment(2, "gppp"," satelite " + satname
			 + gobs->epoch().str_ymdhms(" coordinates not calculated for epoch "));
//	 cout << satname + gobs->epoch().str_ymdhms(" coordinates not calculated for epoch ") << endl;
      _gmutex.unlock(); return -1;
  }  

  t_gtriple txyz(xyz);   


  // relativistic correction
  satclk -= 2.0 * ( txyz[0]*vel[0] + txyz[1]*vel[1] + txyz[2]*vel[2] ) /CLIGHT /CLIGHT;

  // filling gsatgnss
  _satcrd = txyz;
  _clk    = satclk*CLIGHT;
  _dclk = dclk*CLIGHT;
  addecl(txyz, epoT);


#ifdef BMUTEX
  lock.unlock();
#endif
  _gmutex.unlock(); return 1;
}


// -----------
void t_gsatgnss::addclk( double clk )
{
  gtrace("t_gsatgnss::addclk");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _clk = clk;
  _gmutex.unlock();
   return;
}


// -----------
void t_gsatgnss::addele( double ele )
{
  gtrace("t_gsatgnss::addele");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _ele = ele;
  _gmutex.unlock();
   return;
}


// -----------
void t_gsatgnss::addazi( double azi )
{
  gtrace("t_gsatgnss::addazi");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _azi = azi;
  _gmutex.unlock();
   return;
}


// -----------
void t_gsatgnss::addrho( double rho )
{
  gtrace("t_gsatgnss::addrho");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _rho = rho;
  _gmutex.unlock();
  return;
}


// -----------
t_gtriple t_gsatgnss::satcrd()
{
  gtrace("t_gsatgnss::satcrd");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  t_gtriple crd(_satcrd[0], _satcrd[1], _satcrd[2] );
  _gmutex.unlock();  

  return crd;
}


// -----------
double t_gsatgnss::clk()
{
  gtrace("t_gsatgnss::clk");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _clk;
  _gmutex.unlock();  
  return tmp;
}

double t_gsatgnss::dclk()
{
	gtrace("t_gsatgnss::dclk");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	double tmp = _dclk;
	_gmutex.unlock();
	return tmp;
}


// -----------
double t_gsatgnss::ele()
{
  gtrace("t_gsatgnss::ele");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _ele;
  _gmutex.unlock(); 
  return tmp;
}

// -----------
double t_gsatgnss::ele_deg()
{
  gtrace("t_gsatgnss::ele_deg");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _ele*180.0/G_PI;
  _gmutex.unlock();
  return tmp;
}


// -----------
double t_gsatgnss::azi()
{
  gtrace("t_gsatgnss::azi");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _azi;
  _gmutex.unlock();
  return tmp;
}


// -----------
double t_gsatgnss::rho()
{
  gtrace("t_gsatgnss::rho");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _rho;
  _gmutex.unlock();
  return tmp;
}


// valid
// ----------
bool t_gsatgnss::valid()
{
  gtrace("t_gsatgnss::valid");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  bool tmp = _valid();
  _gmutex.unlock();
  return tmp;
}


// clean data
// ----------
void t_gsatgnss::clear()
{ 
  gtrace("t_gsatgnss::clear");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _clear();
  _gmutex.unlock();
   return;
}


// clean internal function
// ----------
void t_gsatgnss::_clear()
{
  gtrace("t_gsatgnss::_clear");   
   
  _ele = 0.0;
  _azi = 0.0;
  _rho = 0.0;
  _clk = 0.0;
  _dclk = 0.0;
}


// clean internal function
// ----------
bool t_gsatgnss::_valid() const
{   
  gtrace("t_gsatgnss::_valid");   
   
  // single validity identification for gsatgnss!
  if( _rho == 0.0 ) return false;
//if( gobs && ! gobs->_valid() ) return false;

  return true; 
}

// get eclipsing
// ---------------------
bool t_gsatgnss::ecl()
{
  gtrace("t_gsatgnss::ecl");   
   
  return _eclipse;
}

// determi ne wheather eclipsed or not
// -------------------------------------
void t_gsatgnss::addecl(t_gtriple& sat, t_gtime& epoch)
{
  gtrace("t_gsatgnss::addecl");   
   
  t_gephplan plan;
  double mjd = epoch.dmjd();
   
  t_gtriple sun = plan.sunPos(mjd);
  
  // Unit vector Earth - Sun 
  ColumnVector eSun = sun.unitary();
     
  // Unit vector Earth - Satellite
  ColumnVector eSat = sat.unitary();
   
  double cosin = dotproduct(eSun, eSat);
  if( cosin < 0 )
  {
    double r = sat.crd_cvect().norm_Frobenius();
    double val = r*sqrt(1-cosin*cosin);
    if( val < Aell )
    {
      _eclipse = true;
      _lastEcl = epoch;
       return;
    }else{
       _eclipse = false;
    }
  }

  double tdiff = epoch.diff(_lastEcl);
  if (tdiff <= 1800) _eclipse = true;
}

} // namespace
