
/**
* @verbatim
    History
    2012-05-11  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gstat.cpp
* @brief       Purpose: statistical function (1D)
* @author      JD
* @version     1.0.0
* @date        2012-05-11
*
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "gutils/gstat.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// Constructor
// ----------
t_gstat::t_gstat(double cint)
: _valid(false),
  _cint(cint)
{
  gtrace("t_gstat::constructor");
  _clear();  
}


// Constructor
// ----------
t_gstat::t_gstat(vector<double>& data, double cint)
: _cint(cint)
{  
  gtrace("t_gstat::constructor");

  _add_data(data);
  _statistics();    // set validity status!
}


// Destructor
// ----------
t_gstat::~t_gstat()
{
  gtrace("t_gstat::destructor");
}


// Add new data (and reset stats)
// ----------
void t_gstat:: add_data(vector<double>& data)
{
  gtrace("t_gstat::add_data");

  _clear();
  _add_data(data);
  _statistics();
  _valid = false;
}


// Calculate statistics
// ----------
int t_gstat::calc_stat(double sig)
{
  gtrace("t_gstat::calc_stat");

  _statistics(sig); // set validity status !

  return 0;
}


// calculate quartiles (upper/lower)
// ---------
int t_gstat::calc_quartiles(double& low, double& upp)
{     
  gtrace("t_stat::calc_quartiles");
   
  vector<double> vec = _data;
  sort(vec.begin(), vec.end());
   
  int n = vec.size(); if(n==0) return -1;
   
  if(n%4 == 0){
    int idx_1 = n/4;
    int idx_2 = n/4 + 1;
    low = ( vec[idx_1] + vec[idx_2] )/2;
  }else if(n%4 != 0){
    double num = n/4.0;
    int idx = (int)ceil(num);
    low = vec[idx-1];
  }else{
    cout << "warning: Problem with lower quartile calculation!" << endl;
  }
   
  if(n%4 == 0){
    int idx_1 = 3*n/4;
    int idx_2 = 3*n/4 + 1;
    upp = ( vec[idx_1] + vec[idx_2] )/2;
  }else if(n%4 != 0){
    double num = n/4.0;
    int idx = (int)ceil(3*num);
    upp = vec[idx-1];
  }else{
    cout << "warning: Problem with upper quartile calculation!" << endl;
  }

  return 0;
}
   

// interquartile limits
// --------------------
int t_gstat::calc_iqrlimits(double& low, double& upp)
{
  gtrace("t_gstat::calc_iqrlimits");

  double min, max;
  calc_quartiles(min, max);
   
  double condition = _cint*(max-min);
  low = min - condition;
  upp = max + condition;
   
  return 0;
}   

// =================
// INTERNAL METHODS
// =================

// clear statistics
// ----------
void t_gstat::_clear()
{
  _valid = false;

  _data.clear();
  _outl.clear();

  _sdev = _rms  =  _var =  0.0;
  _medi = _mean =  0.0;
  _min  = _max  =  0.0;
//_mad  =  0.0;
}


// add data
// ----------
void t_gstat::_add_data(vector<double>& data)
{
  gtrace("t_gstat::_add_data");

  _data   = data;
}


// calculate statistics
// ----------
int t_gstat::_statistics(double sig)
{
  gtrace("t_gstat::statistics");

  _valid = false;
   
  if( _data.size() == 0 ) return -1;

  do{ _calc_median(); _mean = _medi;
      _calc_sdev();
  }while( _chk_resid( sig ) > 0 );

  if( _data.size() == 0 ) return -1;

  _calc_mean();
  _calc_sdev();
  _calc_rms();    // could be merged with above
  _calc_minmax(); // use both _data + _outl

  _valid = true;
  
  return 0;
}

   
// check_residuals
// ----------
int t_gstat::_chk_resid(double sig)
{
  gtrace("t_gstat::_chk_resid");
   
  if( _data.size() == 0 ) return -1;

  double threshold = 0.0;
  double max_res = 0;
  double max_idx = 0;

  if( !double_eq(sig, 0.0) ){ threshold = _cint *   sig; }
  else                      { threshold = _cint * _sdev; }

  for( size_t i = 0; i < _data.size(); ++i ){
    if( fabs( _data[i] - _mean ) > max_res ){
      max_res = fabs( _data[i] - _mean );
      max_idx = i;
    }
  }

  if( max_res > threshold ){
    _outl.push_back( _data[max_idx] );
    _data.erase(     _data.begin() + (long)max_idx );
    return 1; // residual found
  }

  return 0;
}


// simple mean
// ----------
int t_gstat::_calc_mean()
{
  gtrace("t_gstat::_calc_mean");
   
  if(_data.size() == 0){ return -1; }

  for( size_t i = 0; i < _data.size(); ++i ){ _mean  += _data[i]; }
  _mean /= _data.size();

  return 0;
}


// simple rms/var
// ----------
int t_gstat::_calc_rms()
{
  gtrace("t_gstat::_calc_rms");
   
  if(_data.size() == 0){ return -1; }

  double sum = 0.0;
  for( size_t i = 0; i < _data.size(); ++i ){ sum += pow(_data[i], 2); }
  sum /= _data.size();
  _rms = sqrt(sum);

  return 0;
}


// simple sdev
// ----------
int t_gstat::_calc_sdev()
{
  gtrace("t_gstat::_calc_sdev");

  if(_data.size() == 0){ return -1; }

  _var = 0.0;
  for( size_t i = 0; i < _data.size(); ++i ){ 
    _var += pow(_data[i] - _mean, 2);
  }
   _var /= _data.size();
  _sdev  = sqrt(_var);

  return 0;
}


// median
// ----------
int t_gstat::_calc_median()
{
  gtrace("t_gstat::_calc_median");
  _medi = 0.0;

  vector<double> vec = _data;   
  sort( vec.begin(), vec.end() );
  
  int n = vec.size(); if( n == 0 ) return -1;

#ifdef DEBUG   
  for(int i = 0; i<n; i++) cout << i << "  " << vec[i] << endl;
#endif
   
  if( n%2 == 0 ){
    int i = n/2;
    int j = i+1;
    _medi = ( vec[i-1] + vec[j-1] )/2;

#ifdef DEBUG
     cout << i << "  " << j << "  " << vec[i-1] << "  " << vec[j-1] << endl;
#endif

  }else{ 
    int i = n/2 + 1; 
    _medi = vec[i-1];
  }
   
  return 0;
}


// get min & max val
// -----------------
int t_gstat::_calc_minmax()
{
  gtrace("t_gstat::_calc_minmax");

  _min = 0.0;
  _max = 0.0;
   
  if( _data.size() == 0 ) return -1;

  vector<double> vec = _data;
  vec.insert( vec.end(), _outl.begin(), _outl.end() ); // add outliers
  sort( vec.begin(), vec.end() );
  
  _min = *(vec. begin());
  _max = *(vec.rbegin());

  return 0;
}


// histogram
// ----------
t_gstat::t_hist t_gstat::histogram(vector<double>& data, set<double>& bound)
{
  gtrace("t_gstat::histogram");
  double mult = 10.0;  // to specify left/right inf boundary
  
  t_hist hist;
  size_t i = 1;
  for( auto itBOUND = bound.begin(); itBOUND != bound.end(); itBOUND++ )
  {
    set<double>::iterator itBOUND_PRE;
    if( itBOUND != bound.begin() ){
      itBOUND_PRE = itBOUND; itBOUND_PRE--;
    }

    for( auto itDATA = data.begin(); itDATA != data.end(); itDATA++ )
    {
      if(itBOUND == bound.begin())
      {
        t_gpair p(*itBOUND*mult, *itBOUND);          
        if( *itDATA <= *itBOUND ){
               hist[p]++;
        }else{ hist[p] += 0; }              // allocate boundary
      }else{
        t_gpair p(*itBOUND_PRE, *itBOUND);
        if( *itDATA <= *itBOUND && *itDATA > *itBOUND_PRE ){
               hist[p]++;
        }else{ hist[p] += 0; }              // allocate boundary

        if( i == bound.size() ){
          t_gpair p(*itBOUND, *itBOUND*mult);
          if( *itDATA >= *itBOUND ){          // right outliers
     	         hist[p]++;
          }else{ hist[p] += 0; }              // allocate boundary
        }
      }
    }
    i++;
  }
  return hist;
}

/*
// weights
// ----------
double t_gstat::_p(double v)
{
  return 1/sqrt(1+v*v/2);
}
*/
   
} // namespace
