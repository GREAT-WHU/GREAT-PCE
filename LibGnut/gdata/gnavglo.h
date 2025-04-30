/**
*
* @verbatim
    History
    2011-02-14  JD: created
    2018-08-13  JD: updated

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gnavglo.h
* @brief      GLONASS navigation
*.
* @author     JD
* @version    1.0.0
* @date       2011-02-14
*
*/
#ifndef GNAVGLO_H
#define GNAVGLO_H 
 

#include <vector>

#include "newmat/newmat.h"
#include "gdata/gnav.h"
#include "gutils/gtime.h"
#include "gutils/gconst.h"
#include "gutils/gcommon.h"
#include "gutils/gtriple.h"

#define MAX_RINEXN_REC_GLO 15

#define MIN_GLO_RADIUS 25300       // km
#define MAX_GLO_RADIUS 25700       // km
#define MAX_GLO_RADDIF 50          // km
#define MAX_GLO_CLKDIF 50          // ns
#define RAD_GLO_FACTOR 0.001       // m->km
#define CLK_GLO_FACTOR 1000000000  // sec->ns

using namespace std;

namespace gnut {


	class LibGnut_LIBRARY_EXPORT  t_gnavglo : public t_gnav {

 public:
  t_gnavglo();
  virtual ~t_gnavglo();

  // pointers to support NULL if not requested!
  virtual int nav( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true );
  virtual int pos( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true );
  virtual int clk( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true );

  int channel() const;
  int chk(set<string>& msg);
  int ura( double acc ) const;
  int data2nav( string  sat, const t_gtime& ep, const t_gnavdata& data );
  int nav2data( t_gnavdata& data );
   
  int iod() const { return _iodc; }
  int rec() const { return MAX_RINEXN_REC_GLO; }
   
  t_timdbl param( const NAVDATA& n );
  int      param( const NAVDATA& n, double val );

  string line() const;
  string linefmt() const;

  int freq_num() const { return _freq_num; }

 protected:
   int  _iod() const;
   bool _healthy() const;

 private:
   ColumnVector _deriv(const ColumnVector& xx, const t_gtriple& acc);
   ColumnVector _RungeKutta(double step, int nsteps, const ColumnVector& yy, const t_gtriple& acc);
   
   double _maxEphAge;    // max age of ephemerises [s]

   int      _pos( const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_health = true );
   int      _clk( const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL, bool chk_health = true );
   
   int      _iodc;    // issue of clocks

   double   _x;       // position X [km]
   double   _x_d;     // velocity X [km/s]
   double   _x_dd;    // acceleration X [km/s^2]
   double   _y;       // position Y [km]
   double   _y_d;     // velocity Y [km/s]
   double   _y_dd;    // acceleration Y [km/s^2]
   double   _z;       // position Z [km]
   double   _z_d;     // velocity Z [km/s]
   double   _z_dd;    // acceleration Z [km/s^2]
   double   _E;       // Age of oper. information [days]
   int      _freq_num;// frequency number (-7 ... 13)
   double   _health;  // health 0 = OK
   t_gtime  _toc;     // Epoch of clocks [s]
   double   _gamma;   // SV relative frequency bias []
   double   _tau;     // SV clock bias [s]
   double   _tki;     // message frame time [0 ... 86400 s]
   int      _min_step;// mininal step length for Runge Kutta
   
};

} // namespace

#endif
