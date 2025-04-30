/**
*
* @verbatim
    History
    2011-04-20  JD: created
*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     rinexn.h
* @brief    Purpose: brdm RINEX encoder/decoder
*
* @author   JD
* @version  1.0.0
* @date     2011-04-20
*
*/

#ifndef RINEXN_H
#define RINEXN_H
 

#define ID_ALLSAT "ALL_SATELLITES"

#include <vector> 
#include <sstream>

#include "gcoders/gcoder.h"
#include "gdata/grxnhdr.h"
#include "gset/gsetinp.h"
#include "gdata/gnav.h"
using namespace std;

namespace gnut {
    /**
    *@brief Class for t_rinexn derive from t_gcoder
    */
class LibGnut_LIBRARY_EXPORT t_rinexn : public t_gcoder {

 public:
   t_rinexn( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  virtual ~t_rinexn(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
   
  virtual  int encode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

  void gnsssys(char s ){ _gnsssys = s; }
  char gnsssys(){ return _gnsssys; }

  void rxnhdr_all(bool b){ _rxnhdr_all = b; }

 protected:
   virtual int _fill_head();                                                // fill header information
   
   bool _filter_gnav(shared_ptr<t_gnav> geph, const string& prn);
   t_rxnhdr _consolidate_header();

 private:
   char            _gnsssys;
   bool            _rxnhdr_all; // use all RINEX headers in encoder
   t_gtime         _check_dt;   // to validate the message
   t_rxnhdr        _rxnhdr;     // RINEX header
   vector<string>  _comment;    // RINEX comments
   
};

} // namespace

#endif
