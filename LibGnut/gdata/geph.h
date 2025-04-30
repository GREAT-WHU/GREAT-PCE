/**
*
* @verbatim
    History
    2011-01-10  JD: created
  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        geph.h
* @brief       implements ephemerides and navigation classes
*.
* @author      JD
* @version     1.0.0
* @date        2011-01-10
*
*/

#ifndef GEPH_H
#define GEPH_H 
 

#include <memory>

#include "gio/gio.h"
#include "gdata/gdata.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"

using namespace std;
    
namespace gnut {
    /** @brief navigation data. */
    enum NAVDATA {
            NAV_UNDEF,
            NAV_A, NAV_E, NAV_M, NAV_I, NAV_IDOT, NAV_OMEGA, NAV_OMG, NAV_OMGDOT, NAV_DN,
            NAV_CRC, NAV_CIC, NAV_CUC, NAV_CRS, NAV_CIS, NAV_CUS, NAV_F0, NAV_F1, NAV_F2,
            NAV_X, NAV_XD, NAV_XDD, NAV_Y, NAV_YD, NAV_YDD, NAV_Z, NAV_ZD, NAV_ZDD,
            NAV_IOD, NAV_HEALTH,
            NAV_TGD0, NAV_TGD1, NAV_TGD2, NAV_TGD3
          };
    /** @brief time + double. */
    typedef pair<t_gtime, double>  t_timdbl;
   
    /** @brief class for navigation data. */
    class LibGnut_LIBRARY_EXPORT t_geph : public t_gdata {

     public:
      /** @brief default constructor. */
      explicit t_geph();
      virtual ~t_geph();
   
      virtual int clk( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true ){ return -1; }  // [s]
      virtual int pos( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ){ return -1; }  // [m]
      virtual int nav( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ){ 
                  return this->pos(t, xyz, var, vel, chk_health); }  // [m]
   
      virtual bool healthy() const { return true; }
      virtual string health_str() const { return ""; }
      virtual int chk() const { return -1; }

      virtual GNAVTYPE gnavtype(bool full = true) const { return NAV; }
      virtual int src(bool full = true) const { return  0; }
   
      virtual string linefmt() const { return ""; }
      virtual string line() const { return ""; }
      virtual void print() const { cout << linefmt(); }

      virtual t_timdbl param( const NAVDATA& n );
      virtual int      param( const NAVDATA& n, double val );
      virtual bool     param_cyclic( const NAVDATA& n );

      virtual void gio(shared_ptr<t_gio> p){ _gio_ptr = p; }
      /** @brief get the value of gio. */
      shared_ptr<t_gio> gio(){ return _gio_ptr; }
      /** @brief clear data. */
      void clear();
      /** @brief valid. */
      bool valid();
      /** @brief override valid. */
      void valid(bool validity) {_validity = validity;}   
      /** @brief get the name of GNSS system. */
      GSYS    gsys() const;                                          // GNSS system
      /** @brief get the name of satellite. */
      string  gsat() const;                                          // satellite number
      /** @brief get the value of _sat. */
      string  sat()      const{ return _sat; }                       // satellite number  
      /** @brief get the value of validity interval. */
      double  interval() const{ return _interval; }                  // get validity interval
      /** @brief get the value of reference epoch. */
      t_gtime epoch()    const{ return _epoch; }                     // reference epoch
      /** @brief get the begin time of validity. */
      t_gtime begin()    const{ return _epoch - _interval/2; }       // beg of validity
      /** @brief get the end time of validity. */
      t_gtime end()      const{ return _epoch + _interval/2; }       // end of validity
      /** @brief chktot. */
      virtual bool chktot(const t_gtime& t) { return true; }
   
     protected:
         /** @brief clear. */
      virtual void _clear();
      /** @brief valid. */
      virtual bool _valid() const;

      string            _sat;         // satellite number
      t_gtime           _epoch;       // reference epoch
      double            _interval;    // validity interval
      bool              _validity;    // validity
   
      shared_ptr<t_gio> _gio_ptr;     // gio pointer

     private:

    };

} // namespace

#endif
