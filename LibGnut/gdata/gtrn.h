/**
*
* @verbatim
    History
    2013-06-20  PV: created
    2018-08-05  JD: updated

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gtrn.h
* @brief      about trn
*.
* @author     PV
* @version    1.0.0
* @date       2013-06-20
*
*/

#ifndef GTRN_H
#define GTRN_H


#include <stdio.h>
#include <string>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gdata/gobj.h"
#include "gdata/gdata.h"
#include "gdata/grxnhdr.h"

using namespace std;

namespace gnut {
    /** @brief class for gtrn. */
    class t_gtrn : public t_gobj {

     public:
       /** @brief default constructor. */
       t_gtrn();
       virtual ~t_gtrn();

       typedef pair<string,t_rxnhdr>            t_header_pair;
       typedef vector<t_header_pair>            t_header;
      
       virtual void     header(const t_rxnhdr& hdr, string path);
       virtual t_rxnhdr header(string path) const;
       virtual t_header headers() const;

       virtual bool isrec() {return true;}
       virtual bool istrn() {return false;}
   
       virtual void channel(int chk);
       virtual int  channel() const;
   
   
     protected:
       t_rxnhdr        _header(const t_gtime& epo);

       t_header        _headers;               // map of rinex header information
       int _channel;   // temporarily only one value. Must be enhance via _mapchk 
     private:
    };

} // namespace

#endif
