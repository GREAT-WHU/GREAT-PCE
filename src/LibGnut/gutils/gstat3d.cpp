
/**
* @verbatim
    History
    2012-05-11  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gstat3d.cpp
* @brief       Purpose: statistical function (3D)
* @author      JD
* @version     1.0.0
* @date        2012-05-11
*
*/

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gutils/gstat3d.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// Constructor
// ----------
t_gstat3d::t_gstat3d(vector<t_gtriple>& data, double cint)
 : t_gstat(cint)
{
  gtrace("t_gstat3d::constructor");

  _add_data(data);
  _statistics();
}


// Destructor
// ----------
t_gstat3d::~t_gstat3d()
{
  gtrace("t_gstat3d::destructor");
}


// Calculate statistics
// ----------
int t_gstat3d::calc_stat(double sig)
{
  gtrace("t_gstat3d::calc_stat");
   
  _statistics(sig); // set validity status !
    
  return 0;
}

// =================
// INTERNAL METHODS
// =================

// reset stat
// ----------
void t_gstat3d::_clear()
{
  t_gstat::_clear();
  
  _data3d.clear();

  t_gtriple null3d(0.0, 0.0, 0.0);
  _medi3d = null3d;
  _mean3d = null3d;
  _sdev3d = null3d;
  _rms3d  = null3d;
  _var3d  = null3d;
}


// add data
// ----------
void t_gstat3d::_add_data(vector<t_gtriple>& data)
{
  gtrace("t_gstat3d::_add_data");
   
  // 1D initialization
  vector<double> data1d;
  for( size_t i = 0; i < data.size(); i++ ){
     data1d.push_back( data[i].norm() );
  }
  t_gstat::_add_data(data1d);

  // 3D initialization
  _data3d = data;
}


// statistics
// ----------
int t_gstat3d::_statistics(double sig)
{
  gtrace("t_gstat3d::_statistics");
  
//if( t_gstat::_statistics(sig) > 0 ) return -1;  // calculated along with 3D

  if( _data3d.size() == 0 ) return -1;

  do{ _calc_median(); _mean3d = _medi3d;
      _calc_sdev();
  }while( _chk_resid3d( sig ) > 0 );

  if( _data3d.size() == 0 ) return -1;

  _calc_mean();
  _calc_sdev();
  _calc_rms();

  _valid = true;

  return 0;
}

// check_residuals
// ----------
int t_gstat3d::_chk_resid3d(double sig)
{
  gtrace("t_gstat3d::_chk_resid3d");
  
  if( _data3d.size() == 0 ) return -1;
    
  double threshold = 0.0;
  double max_res = 0;
  double max_idx = 0;

  if( !double_eq(sig, 0.0) ){ threshold = _cint *   sig; }
  else                      { threshold = _cint * _sdev; }
              
  for( size_t i = 0; i < _data3d.size(); ++i ){
    double res = (_data3d[i] - _mean3d).norm();
    if( fabs( res ) > max_res ){
      max_res = fabs( res );
      max_idx = i;
    }
  }

  if( max_res > threshold ){
    _outl3d.push_back( _data3d[max_idx] );
    _data3d.erase(     _data3d.begin() + (long)max_idx );
    return 1; // residual found
  }
                                                           
  return 0;
}
   

// median
// ----------
int t_gstat3d::_calc_median()
{
  gtrace("t_gstat3d::_calc_median");
   
  // 1D
  t_gstat::_calc_median();
     
  // 3D -- is not a real 3D median, just a collection of three 1D medians..
  // .. consider re-implementation ..
  t_gstat stat1d;
  vector<double> a, b, c;
  for( size_t i = 0; i < _data3d.size(); ++i ){
    a.push_back(_data3d[i][0]);
    b.push_back(_data3d[i][1]);
    c.push_back(_data3d[i][2]);
  }
  stat1d.add_data(a); if( stat1d.valid() ) _medi3d[0] = stat1d.get_median();
  stat1d.add_data(b); if( stat1d.valid() ) _medi3d[1] = stat1d.get_median();
  stat1d.add_data(c); if( stat1d.valid() ) _medi3d[2] = stat1d.get_median();
   
  return 0;
}


// simple mean
// ----------
int t_gstat3d::_calc_mean()
{
  gtrace("t_gstat3d::_calc_mean");

  if( t_gstat::_calc_mean() > 0 || _data3d.size() == 0 ) return -1;

  for( size_t i = 0; i < _data3d.size(); ++i ){ _mean3d += _data3d[i]; }
  _mean3d /= _data3d.size();

  return 0;
}


// simple rms/var
// ----------
int t_gstat3d::_calc_rms()
{
  gtrace("t_gstat3d::_calc_rms");

  if( t_gstat::_calc_rms() > 0 || _data3d.size() == 0 ) return -1;

  t_gtriple sum3d(0.0, 0.0, 0.0);
  for( size_t i = 0; i < _data3d.size(); ++i )
  {
    sum3d[0] += pow(_data3d[i][0], 2);
    sum3d[1] += pow(_data3d[i][1], 2);
    sum3d[2] += pow(_data3d[i][2], 2);
  }

   sum3d /= _data3d.size();
  _rms3d /= _data3d.size();

  return 0;
}


// simple rms/var
// ----------
int t_gstat3d::_calc_sdev()
{
  gtrace("t_gstat3d::_calc_sdev");

  if( t_gstat::_calc_sdev() > 0 || _data3d.size() == 0 ) return -1;

  _var3d = t_gtriple(0.0, 0.0, 0.0);
  for( size_t i = 0; i < _data3d.size(); ++i )
  {
     _var3d[0] += pow(_data3d[i][0] - _mean3d[0], 2);
     _var3d[1] += pow(_data3d[i][1] - _mean3d[1], 2);
     _var3d[2] += pow(_data3d[i][2] - _mean3d[2], 2);
   }

  _var3d   /= _data3d.size();
  _sdev3d[0] = sqrt( _var3d[0] );
  _sdev3d[1] = sqrt( _var3d[1] );
  _sdev3d[2] = sqrt( _var3d[2] );

  return 0;
}

} // namespace
