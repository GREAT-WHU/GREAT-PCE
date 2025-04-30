/**
 * @file         glsqproc.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GPROCMODEL_H
#define GPROCMODEL_H

#include "gexport/ExportLibGREAT.h"
#include "gdata/gdata.h"
#include "gall/gallproc.h"
#include "gall/gallprod.h"
#include "gall/gallobj.h"
#include "gall/gallobs.h"
#include "gproc/glsq.h"
#include "gmodels/gmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gcombmodel.h"
#include "gproc/glsqmatrix.h"
#include "gdata/gsatdata.h"
#include "gutils/gcycleslip.h"
#include "gset/gsetamb.h"
#include "gproc/gqualitycontrol.h"

namespace great {
typedef vector<tuple<string, string, string> > t_lsq_equ_info;

/**
* @brief virtual class for process
*/
class LibGREAT_LIBRARY_EXPORT t_glsqproc {
 public:
  t_glsqproc();
  t_glsqproc(t_gsetbase *set, t_gallproc *data, t_glog *log);
  virtual ~t_glsqproc();

  void add_coder(const vector<t_gcoder *> &coder);
  /** @brief ProcessBatch */
  virtual bool ProcessBatch(t_gallproc *data, const t_gtime &beg, const t_gtime &end);
  /** @brief Process One Epoch Data */
  bool processOneEpoch(const t_gtime &crt_epoch, std::vector<t_gsatdata> &crt_obs) {
    return _processOneEpoch(crt_epoch, crt_obs);
  };
  /** @brief generate product */
  virtual bool GenerateProduct();
  /** @brief get sat list */
  set<string> &sat_list() { return _sat_list; }
  /** @brief get rec list */
  set<string> &rec_list() { return _rec_list; }
  /** @brief get sys list */
  set<string> &sys_list() { return _sys_list; }
  /** @brief get obs interval */
  double get_obs_intv() const { return _obs_intv; };
  set<string> get_rec_list() { return _rec_list; };
  t_glog *get_glog() { return _glog; };
  /** @brief get begin time */
  t_gtime get_beg_time() { return _beg_time; };
  /** @brief get end time */
  t_gtime get_end_time() { return _end_time; };
  /** @brief get current time */
  t_gtime get_crt_time() { return _crt_time; };
  /** @brief get all obs data */
  t_gallobs *get_gall_obs() { return _gall_obs; };
  /** @brief get lsq eqution */
  t_glsq *get_lsq() { return _lsq; };
  t_gallrecover *get_gall_rcv() { return _gall_rcv; };
  /** @brief get current lsq eqution */
  t_glsqEquationMatrix get_crt_equ() { return _crt_equ; };
  /** @brief set current time */
  void set_crt_time(const t_gtime& time) { _crt_time = time; };
  /** @brief set end time */
  void set_end_time(const t_gtime& time) { _end_time = time; };
  /** @brief set begin time */
  void set_beg_time(const t_gtime& time) { _beg_time = time; };
  /** @brief set all obs data */
  void set_gall_obs(t_gallobs *obs) { _gall_obs = obs; };

 protected:

  t_glog *_glog = nullptr;
  t_gsetbase *_gset = nullptr;

  set<string> _rec_list;                  ///< all sites of proc
  set<string> _sat_list;                  ///< all sats of proc
  set<string> _sys_list;                  ///< sat systems of proc
  map<string, t_gtriple> _rec_crds;
  map<string, t_gtriple> _rec_stds;

  int _frequency = 0;
  double _obs_intv = 0.0;      /// obs interval
  double _pwc_intv = 0.0;      /// pwc interval
  double _grd_pwc_intv = 0.0;  /// grd pwc interval
  double _isb_pwc_intv = 0.0;  /// isb pwc interval
  t_gtime _beg_time;           /// begin time
  t_gtime _end_time;           /// end time
  t_gtime _crt_time;           /// crt time

  set<pair<FREQ_SEQ, FREQ_SEQ>>                _freq_pair;
  map< GSYS, map<FREQ_SEQ, GOBSBAND> >         _band_index;
  map< GSYS, map<GOBSBAND, FREQ_SEQ> >         _freq_index;
  map< GSYS, int>                              _freq_number;
  map<string, int> _glofrq_num;

  /** @brief set init mode */
  LSQMODE _lsq_mode = LSQMODE::LSQMODE_UNDEF;
  OBSCOMBIN _obs_mode = OBSCOMBIN::DEF_OBSCOMBIN;
  CONSTRPAR _crd_est = CONSTRPAR::CONSTRPAR_UNDEF;
  SLIPMODEL _slip_model = SLIPMODEL::SLIPMODEL_UNDEF;
  ZTDMODEL _ztd_model = ZTDMODEL::DEF_ZTDMODEL;
  SYSBIASMODEL _sysbias_model = SYSBIASMODEL::AUTO_WHIT;
  IFCB_MODEL _ifcb_model = IFCB_MODEL::IFCB_MODEL_UNDEF;
  IONO_ORDER _iono_order = IONO_ORDER::NONE;
  IFB_MODEL _ifb_model = IFB_MODEL::NONE;
  bool _lite_turboedit = false;

  std::string _crt_rec;                    ///< current site of proc
  std::string _crt_sat;                    ///< current sat  of proc

  t_gallobs *_gall_obs = nullptr;                ///< all observ info
  t_gallobj *_gall_obj = nullptr;                ///< all obj info include pcv and intial crd
  t_gallrecover *_gall_rcv = nullptr;
  t_gallbias *_gall_bias = nullptr;

  t_gallproc *_gall_proc = nullptr;
  t_gpoleut1 *_gall_erp = nullptr;
  t_gallnav *_gall_nav = nullptr;
  t_gallprec *_gall_prec = nullptr;

  t_gpoleut1 *_gdata_erp = nullptr;

  t_glsq *_lsq = nullptr;        ///< least square estimator
  shared_ptr<t_gcycleslip> _slip12 = nullptr;
  shared_ptr<t_gbasemodel> _base_model = nullptr;
  shared_ptr<t_gbiasmodel> _bias_model = nullptr;
  t_gprecisemodel *_bias_pcemodel = nullptr;

  t_glsqEquationMatrix _crt_equ;
  map<string, map<t_gtime, t_gtriple> > _kin_xyz;


  map<int, string> _slr_list;
  map<string, vector<string>> _sat_sta;

  double _recclk_threshold = 0.0;

  int _num_threads = 1;
  bool _matrix_remove = false;
  bool _cmb_equ_multi_thread = false;

  int _cmb_equ_msec{};
  int _remove_par_msec{};

 protected:
  double _maxres_norm = 0.0;
  t_glog *_clk_log = nullptr;
  t_giof *_satclkfile = nullptr;
  string _biasdir;
  string _biasname;
  map<GSYS, double> _maxres_C;
  map<GSYS, double> _maxres_L;
  map<string, t_gallrecover *> _siteres;
  vector<t_gcoder *> _gcoder;
  bool _write_equ = true;

  int _obs_crt_num = 0;    

  double _ambcon_weight = 1e7;
  shared_ptr<t_gqualitycontrol> _quality;

 protected:
  /** @brief init receiver trop pars */
  virtual bool _init_rec_trop_pars(t_glsq *lsq);
  /** @brief init receiver crd pars */
  virtual bool _init_rec_crd_pars(t_glsq *lsq, bool isPPP = true);
  /** @brief init receiver clk pars */
  virtual bool _init_rec_clk_pars(t_glsq *lsq);
  /** @brief init satellite trop pars */
  virtual bool _init_sat_clk_pars(t_glsq *lsq);
  /** @brief init lsq product data */
  virtual bool _initLsqProdData(t_gallprod *data);
  /** @brief init lsq process pars */
  virtual bool _initLsqProcPars(t_glsq *lsq);
  /** @brief init lsq process data */
  virtual bool _initLsqProcData(t_gallproc *data);
  /** @brief select current epoch obs */
  bool _select_obs(const t_gtime& epoch, vector<t_gsatdata>& all_obs);
  /** @brief select receiver current epoch obs */
  bool _select_rec_obs(const string &rec,
                       const t_gtime &epoch,
                       vector<t_gsatdata> &all_obs,
                       vector<t_gsatdata> &rec_obs);
  /** @brief init satellite clock */
  bool _init_sat_clk_IF(vector<t_gsatdata> &crt_obs, ColumnVector &l, t_glsq *lsq);
  bool _init_sat_clk();
  /** @brief init xml settings */
  bool _init_xml_settings(t_gsetbase *set);
  /** @brief init lsq settings */
  bool _init_lsq_settings(t_gsetbase *set, t_gallproc *data, t_glog *log);
  /** @brief extract clock file */
  bool _extract_clkfile(t_gallrecover &recover_data, t_gallprec::clk_type type);
  /** @brief extract residual file */
  bool _extract_resfile(t_gallrecover &recover_data, string resfile = "");

  bool _myLog = false;
  bool _create_log(t_glog *log);
  static bool _delete_log(t_glog *log);
  /** @brief process one epoch data */
  virtual bool _processOneEpoch(const t_gtime &crt_epoch, std::vector<t_gsatdata> &crt_obs);
  virtual bool _processOneEpoch_thread_safe(const t_gtime &crt_epoch, std::vector<t_gsatdata> &crt_obs);
  /** @brief process one receiver data */
  virtual bool _processOneRec(const t_gtime &crt_epoch, const std::string &crt_rec, std::vector<t_gsatdata> &crt_obs);
  virtual bool _processOneRec_thread_safe(const t_gtime &crt_epoch,
                                          const std::string &crt_rec,
                                          std::vector<t_gsatdata> &crt_obs,
                                          t_glsqEquationMatrix &equ_result);
  /** @brief process one receiver and one satellite data */
  virtual bool _processOneRecOneSat(const t_gtime &crt_epoch,
                                    const std::string &rec,
                                    const std::string &sat,
                                    t_gsatdata &crt_obs);
};

}
#endif