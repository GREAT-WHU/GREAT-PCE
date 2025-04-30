/**
*
* @verbatim
    History
      2011-11-04  JD: created
      2019-11-29  ZHJ: add encoder
*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     rinexc.h
* @brief    Purpose: Clock RINEX encoder/decoder
*
* @author   JD
* @version  1.0.0
* @date     2011-11-04
*
*/

#ifndef RINEXC_H
#define RINEXC_H


#include <vector> 

#include "gcoders/gcoder.h"
#include "gutils/gtime.h"
#include "gall/gallobj.h"
#include "gall/gallprec.h"

#define RINEXC_BUFFER_LEN 81

using namespace std;

namespace gnut {
    /**
    *@brief Class for t_rinexc derive from t_gcoder
    */
class LibGnut_LIBRARY_EXPORT t_rinexc : public t_gcoder {

 public:
   t_rinexc( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  virtual ~t_rinexc(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

  /**
  * @brief encode header of  resfile
  * @param[in]  buff        buffer of the data
  * @param[in]  sz          buffer size of the data
  * @param[in]  errmsg      error message of the data decoding
  * @return  size of header encoding
  */
  virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

  /**
  * @brief encode data of  resfile
  * @param[in]  buff        buffer of the data
  * @param[in]  sz          buffer size of the data
  * @param[in]  errmsg      error message of the data decoding
  * @return  size of data body encoding
  */
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;


  void gnsssys(char s ){ _gnsssys = s; }
  char gnsssys(){ return _gnsssys; }
 
 protected:
  struct t_meta {
    t_gtriple  xyz;
    string     name,domes,desc,ant,rec;
    t_gtime    begCLK, endCLK;
  };
     
  t_gtime                         _file_beg;              // file first epoch
  t_gtime                         _file_end;              // file last epoch
  t_gtime                         _file_run;              // file created

  t_gallobj*                      _allobj;
  map<string,shared_ptr<t_gobj>>  _mapsat;
  map<string,shared_ptr<t_gobj>>  _maprec;
  map<string,shared_ptr<t_gobj>>::const_iterator itOBJ;

  set<string>                     _sites;
  set<string>                     _satellites;

 private:
   char  _gnsssys;

   t_gallprec* _get_prec_data();

};

} // namespace

#endif
