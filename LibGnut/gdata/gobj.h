/**
*
* @verbatim
    History
    2011-01-10  JD: created
    2019-04-05  HJ_Zheng: add partial function for orbit

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gobj.h
* @brief      Purpose: implements gobject class
*.
* @author     JD
* @version    1.0.0
* @date       2011-01-10
*
*/
#ifndef GOBJ_H
#define GOBJ_H


#include <stdio.h>
#include <string>
#include <memory>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gexport/ExportLibGnut.h"
#include "gdata/gdata.h"
#include "gdata/grnxhdr.h"
#include "gall/gallpcv.h"
#include "gmodels/gpcv.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {
    /** @brief class for t_gobj. */
    class LibGnut_LIBRARY_EXPORT t_gobj : public t_gdata {

     public:
         /** @brief default constructor. */
       t_gobj();
       virtual ~t_gobj();

       typedef map<t_gtime,shared_ptr<t_gpcv> >          t_mappcv;
       typedef map<t_gtime,t_gtriple>                   t_mapecc_xyz;
       typedef map<t_gtime,t_gtriple>                   t_mapecc_neu;   
       typedef map<t_gtime,pair<t_gtriple, t_gtriple> >  t_mapcrd;
       //   typedef map<t_gtime,t_gtriple>              t_mapstd;
       typedef map<t_gtime,string>                      t_mapant;

    //   typedef map<t_gtime,t_gbias*>    t_mapbias;

       virtual void id( string str );                    // set/get id (uniq internal id)
       virtual string id() const;
   
       virtual void name( string str );                  // set/get name (name)
       virtual string name() const;
   
       virtual void domes( string str );                 // set/get name (domes)
       virtual string domes() const;
   
       virtual void desc( string str );                  // set/get description (full name)
       virtual string desc() const;
   
       virtual void eccxyz(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       virtual t_gtriple eccxyz(const t_gtime& t) const;        // set/get ecc offsets (>=t) w.r.t. center of mass/reference point
       virtual void eccxyz_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const;
   
       virtual void eccneu(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       virtual t_gtriple eccneu(const t_gtime& t) const;        // set/get ecc offsets (>=t) w.r.t. center of mass/reference point
       virtual void eccneu_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const;
   
       virtual void pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       virtual shared_ptr<t_gpcv> pcv(const t_gtime& t) const;  // set/get external pointer to pcv element (>=t)
   
       virtual void ant(string ant, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       virtual string ant(const t_gtime& t) const;              // set/get antenna
       virtual void ant_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const;
   
       virtual void channel(int chk){};
       virtual int  channel() const {return 255;};


       virtual void crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end = LAST_TIME,bool overwrite = false);
       virtual void crd(const t_gtriple& crd, const t_gtriple& std);
       virtual t_gtriple crd(const t_gtime& t) const;      // get crd of Marker point
       virtual t_gtriple std(const t_gtime& t) const;      // get std of crd
       virtual t_gtriple crd_arp(const t_gtime& t) const;  // get crd of ARP (ARP = MARKER + ECC)
       virtual void crd_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const;
       virtual bool get_recent_crd(const t_gtime& t,const double& ref_std,t_gtriple& crd,t_gtriple& std);
       virtual bool get_adjacent_crd(const t_gtime& t, const double& ref_std, t_gtriple& crd, t_gtriple& std);
   
       virtual vector<t_gtime> pcv_id() const;
       virtual vector<t_gtime> ant_id() const;
       virtual vector<t_gtime> crd_id() const;

       virtual void compare(shared_ptr<t_gobj> gobj, const t_gtime& tt, string source);

       virtual bool isrec() = 0;
       virtual bool istrn() = 0;
   
               bool overwrite();
               void overwrite(bool overwrite);
   
       virtual void sync_pcv(t_gallpcv* pcvs);

       virtual bool    operator<(const t_gobj& t) const;
       virtual bool    operator==(const t_gobj& t) const;
   
     protected:  
       string           _id;                      // object id (internal)
       string           _name;                    // object name 
       string           _domes;                   // object domes
       string           _desc;                    // object description (full name)
       t_mapcrd         _mapcrd;                  // object position
       t_mappcv         _mappcv;                  // map of pco+pcv
       t_mapecc_xyz     _mapeccxyz;               // map of xyz eccentricities (to center of mass or reference point)
       t_mapecc_neu     _mapeccneu;               // map of neu eccentricities (to center of mass or reference point)   
       t_mapant         _mapant;                  // map of antennas + dome
   
       bool             _overwrite;
   
       // source for public (mutexed) interfaces
       shared_ptr<t_gpcv>  _pcv(const t_gtime& t) const;   
       void                _pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       void                _ant(string ant, const t_gtime& beg, const t_gtime& end = LAST_TIME );
       string              _ant(const t_gtime& t) const;   
       vector<t_gtime>     _ant_id() const;   
       void                _crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end = LAST_TIME, bool overwrite = false); // overwrite for solution must over write; by ZHJ
       t_gtriple           _crd(const t_gtime& t) const;
       t_gtriple           _std(const t_gtime& t) const;   
       void                _eccxyz(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end = LAST_TIME );   
       t_gtriple           _eccxyz(const t_gtime& t) const;
       void                _eccneu(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end = LAST_TIME );   
       t_gtriple           _eccneu(const t_gtime& t) const;
       shared_ptr<t_gpcv>  _pcvnull;

     private:
    };

} // namespace

#endif
