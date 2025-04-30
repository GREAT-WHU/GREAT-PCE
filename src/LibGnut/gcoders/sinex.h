/**
*
* @verbatim
    History
    2015-11-24  JD: created
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     sinex.h
* @brief    Purpose: SINEX
*    Todo:    header information completed
            encoders implementation
* @author   JD
* @version  1.0.0
* @date     2015-11-24
*
*/

#ifndef SINEX_H
#define SINEX_H

#include "gcoders/gcoder.h"
#include "gall/gallobj.h"
#include "gall/gallprod.h"
#include "gutils/gtime.h"
#include "gmodels/ggpt.h"

using namespace std;

namespace gnut {

enum t_snx_type { SINEX_GNS, TROSNX_GNS, TROSNX_NWM, TROSNX_OTH };
    /**
    *@brief Class for t_sinex derive from t_gcoder
    */
class LibGnut_LIBRARY_EXPORT t_sinex : public t_gcoder {

 public:
   t_sinex( t_gsetbase* s, string version, int sz = DEFAULT_BUFFER_SIZE, string id = "sinex" );
 virtual ~t_sinex(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
   
  virtual  int encode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

  virtual void technique(char c);
  virtual char technique();

 protected:
  struct t_meta {
    t_gtriple  ell,xyz,ecc,apr,std,rms,var,cov,idx;
    t_gtriple  gps_neu_L1,gps_neu_L2;
    vector<string> par;
    string     id,name,domes,desc,ant,rec,snx_code;
    t_gtime    begPOS, endPOS;
  };

  virtual t_gtime  _site_beg(string site);
  virtual t_gtime  _site_end(string site);
   
  virtual  shared_ptr<t_grec> _get_grec(string site);
  virtual void  _add_data(string id, t_gdata* pt_data);
  virtual  int  _initialize_data();
  virtual  int  _decode_vers();
  virtual  int  _decode_data();
  virtual  int  _decode_comm();
  virtual  int  _decode_block();

  virtual  int  _fill_site_INI(string site, t_meta& meta);
  virtual void  _fill_site_IDE(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_SOL(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_XYZ(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_ECC(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_ANT(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_REC(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_PCO(string site, t_meta& meta, ostringstream& os);
   
  virtual void  _fill_site_EST(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_APR(string site, t_meta& meta, ostringstream& os);
  virtual void  _fill_site_COV(string site, t_meta& meta, ostringstream& os);

  virtual void  _fill_head_INI();
  virtual void  _fill_head_IDE(ostringstream& os);
  virtual void  _fill_head_SOL(ostringstream& os);
  virtual void  _fill_head_XYZ(ostringstream& os);
  virtual void  _fill_head_ECC(ostringstream& os);
  virtual void  _fill_head_ANT(ostringstream& os);
  virtual void  _fill_head_REC(ostringstream& os);
  virtual void  _fill_head_PCO(ostringstream& os);

  virtual void  _fill_data_STT();
  virtual void  _fill_data_EST(ostringstream& os);
  virtual void  _fill_data_APR(ostringstream& os);
  virtual void  _fill_data_COV(ostringstream& os);

  virtual void  _complete_obj(shared_ptr<t_grec> obj, const t_gtime& epo);

  t_snx_type                      _snx_type;              // TRO-SINEX TYPE  
  char                            _technique;             // TECHNIQUE TYPE
  int                             _parindex;              // PARAMETER INDEX
  int                             _tmpsize;               // working amount of bytes processed
  int                             _consume;               // working total amount of bytes (return)
  bool                            _complete;              // working flag for completed epoch decoding
  bool                            _estimation;            // estimated parameters
  string                          _code_label;            // site code label (CODE:4 vs STATION__:9)
  string                          _list_gnss;             // list of GNSSs (GREC)
  string                          _line;                  // working line read from
  string                          _site;                  // cache
  string                          _block;                 // working block name
  string                          _pco_mod;               // PCO model
  string                          _ac;                    // analyses center abbr
  vector<string>                  _comment;               // vector of comments

  t_gpt                           _ggpt;
  t_gallprod*                     _pt_prod;

  t_gtime                         _file_beg;              // file first epoch
  t_gtime                         _file_end;              // file last epoch
  t_gtime                         _file_run;              // file created
   
  t_gallobj*                      _allobj;
  map<string,shared_ptr<t_gobj>>  _mapobj;
  map<string,shared_ptr<t_gobj>>::const_iterator itOBJ;
   
  map<string,set<t_gtime>>        _epo_pos;               // site/product epochs (POS)
  map<string,pair<int,int>>       _mapidx;
  set<string>                     _set_types;
  set<string>                     _sites;
  set<string>::const_iterator      itSET;

  // Decoding estimated coordinates
  map<string, set<t_gtime>>                       _sol_epoch;   // For decoding SOLUTION/EPOCH
  bool                                            _end_epoch;   // The block completed
  map<string, map<t_gtime, map<string, t_gpair>>> _sol_estim;   // For decoding SOLUTION/ESTIMATE
  bool                                            _end_estim;   // The block completed
  void _set_rec_crd();
};

} // namespace

#endif
