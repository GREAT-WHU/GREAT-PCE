/**
*
* @verbatim
    History

    @endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

*
* @file        gsetproc.h
* @brief       implements process setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/
#ifndef GSETPROC_H
#define GSETPROC_H

#define XMLKEY_PROC "process" ///< The defination of process node in XML
#define XMLKEY_ROM "read_ofile_mode"

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gutils/gobs.h"
#include "gset/gsetbase.h"
#include "gexport/ExportLibGnut.h"

using namespace std;
using namespace pugi;

namespace gnut {

enum class CONSTRPAR { EST, FIX, KIN, SIMU_KIN, CONSTRPAR_UNDEF };
enum HMTCONSTR { NNT, NNR, NNTNNR, HMTCONSTR_UNDEF };
enum class GRDMPFUNC { DEF_GRDMPFUNC, TILTING, CHEN_HERRING, BAR_SEVER };
enum class ZTDMPFUNC { DEF_ZTDMPFUNC, COSZ, GMF, NO_MF };
enum class IONMPFUNC { DEF_IONMPFUNC, ICOSZ, QFAC, NONE };
enum class OBSWEIGHT { DEF_OBSWEIGHT, EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH, PARTELE, SNR };
enum class TROPMODEL { DEF_TROPMODEL, SAASTAMOINEN, DAVIS, HOPFIELD, MOPS, GPTW, GPT2W, GAL27, GALTROPO27, EXTERN };
enum class ZTDMODEL { DEF_ZTDMODEL, PWC, STO };
enum class RESIDTYPE { DEF_RESIDTYPE, RES_ORIG, RES_NORM, RES_ALL };
enum class OBSCOMBIN {
  DEF_OBSCOMBIN,
  IONO_FREE,
  RAW_SINGLE,
  RAW_DOUBLE,
  RAW_ALL,
  MW_COMBIN,
  EWL_COMBIN,
  WL_COMBIN,
  RAW_MIX,
  IF_P1
};
enum class PHASEBIAS { DEF_POSTPRD, RTCM_PRD, CNES_PRD, ESTI_PRD };
enum class ATTITUDES { DEF_YAWMODEL, YAW_NOMI, YAW_RTCM };
enum class CBIASCHAR { DEF_CBIASCHAR, ORIG, CHAR2, CHAR3 };
enum class SYSBIASMODEL { AUTO_CON, AUTO_RWK, AUTO_WHIT, AUTO_PWC, ISB_CON, ISB_RWK , ISB_WHIT, ISB_PWC };
enum class GRDMODEL { DEF_GRDMODEL, GRD_PWC, GRD_STO };
enum class BASEPOS { RINEXH, SPP, CFILE };

/** @brief  enum for lsq model */
enum class LSQMODE {
  EPO,    ///< EPO model:solve every epoch
  LSQ,    ///< LSQ model:Overall solution
  LSQMODE_UNDEF
};

enum class SLIPMODEL {
  DEF_DETECT_MODEL, ///<
  TURBO_EDIT,
  SLIPMODEL_UNDEF
};

enum class IONMODEL {
  VION,
  SION,
  DEF_ION,
  IONMODEL_UNDEF
};

enum class IFCB_MODEL {
  EST,
  COR,
  DEF,
  IFCB_MODEL_UNDEF
};

enum class IFB_MODEL {
  EST_REC_IFB,
  EST_SAT_IFB,
  EST_IFB,
  COR_REC_IFB,
  COR_SAT_IFB,
  COR_IFB,
  NONE,
  IFB_MODEL_UNDEF
};

enum class DCB_MODEL {
  CODE,
  CAS,
  NONE,
  DCB_MODEL_UNDEF
};

enum class FREQ_MODEL {
  DOUBLE,
  TRIPLE,
  QUADRUPLE, //
  QUINTUPLE,
  FREQ_MODEL_UNDEF
};

enum class IONO_ORDER {
  FIRST_ORDER,
  SECOND_ORDER,
  THIRD_ORDER,
  NONE,
  IONO_ORDER_UNDEF
};


// The Define for LC PC
class LibGnut_LIBRARY_EXPORT t_gobscombtype {

 public:
  t_gobscombtype();
  t_gobscombtype(const t_gobscombtype &other);
  explicit t_gobscombtype(const string &obscombtype);
  t_gobscombtype(const t_gobs &obstype, OBSCOMBIN combtype);
  t_gobscombtype(const t_gobs &obstype, GOBSBAND b1, FREQ_SEQ freq_1, OBSCOMBIN combtype);
  t_gobscombtype(const t_gobs &obstype, GOBSBAND b1, GOBSBAND b2, FREQ_SEQ freq_1, FREQ_SEQ freq_2, OBSCOMBIN combtype);
  t_gobscombtype(GOBSTYPE t, GOBSBAND b, OBSCOMBIN obscomb);

  string convert2str() const;

  bool operator==(const t_gobscombtype &g) const;
  bool operator<(const t_gobscombtype &g) const;
  bool is_freq(const FREQ_SEQ& freq_1, const FREQ_SEQ& freq_2) const {
	  return (_obs_freq_1 == freq_1 && _obs_freq_2 == freq_2);
  }
  bool is_freq12() const { return (_obs_freq_1 == FREQ_1 && _obs_freq_2 == FREQ_2); }
  bool is_freq_raw1() const { return (_obs_freq_1 == FREQ_1 && _obs_freq_2 == FREQ_X); }

  GOBSBAND getBand_1() const { return _obs_band_1; };
  GOBSBAND getBand_2() const { return _obs_band_2; };
  FREQ_SEQ getFreq_1() const { return _obs_freq_1; };
  FREQ_SEQ getFreq_2() const { return _obs_freq_2; };

  bool is_phase() const;
  bool is_code() const;

protected:
  GOBSTYPE  _obs_type   = GOBSTYPE::TYPE;
  GOBSBAND  _obs_band_1 = BAND;
  GOBSBAND  _obs_band_2 = BAND;
  FREQ_SEQ  _obs_freq_1 = FREQ_X;
  FREQ_SEQ  _obs_freq_2 = FREQ_X;
  OBSCOMBIN _obs_combine= OBSCOMBIN::DEF_OBSCOMBIN;
};

/// The class for setting of process modelue in XML file
class LibGnut_LIBRARY_EXPORT t_gsetproc : public virtual t_gsetbase {
 public:
  /// constructor
  t_gsetproc();
  /// destructor
  ~t_gsetproc() override;

  /// settings check
  void check() override;
  /// settings help
  void help() override;

  /**
   * @brief         get the decimals for sampling interval (for high-rate) in process
   * @return bool : decimals for sampling interval (for high-rate) in process
   */
  bool tropo();
  bool iono();
  bool tropo_slant();
  bool tropo_grad();
  bool phase();
  bool doppler();
  bool pos_kin();
  bool auto_band();

  int frequency();
  int num_threads();
  bool matrix_remove();
  bool cmb_equ_multi_thread();
  /**@brief initial sigma */
  double sig_init_ztd();
  double sig_init_vion();
  double sig_init_grd();
  double sig_init_crd();
  double sig_init_vel();
  double sig_init_amb();
  double sig_init_glo();
  double sig_init_gal();
  double sig_init_bds();
  double sig_init_qzs();

  double sig_ref_clk();
  double minimum_elev();
  double max_res_norm();
  int minsat();

  string ref_clk();
  string ref_rec();
  LSQMODE lsq_mode();
  SLIPMODEL slip_model();
  IONMODEL ion_model();
  int lsq_buffer_size();
  /**@brief set process */
  TROPMODEL tropo_model();
  ZTDMODEL ztd_model(double *dt = nullptr);
  ZTDMPFUNC tropo_mf();
  IONMPFUNC iono_mf();
  GRDMPFUNC grad_mf();
  CONSTRPAR crd_est();
  OBSWEIGHT weighting();
  RESIDTYPE residuals();
  OBSCOMBIN obs_combin();
  ATTITUDES attitudes();
  CBIASCHAR cbiaschar();

  SYSBIASMODEL sysbias_model(double& dt);
  SYSBIASMODEL ifbbias_model();
  GRDMODEL grd_model(double *dt = nullptr);
  /**@brief str convert */
  GRDMPFUNC str2grdmpfunc(const string &mf);
  ZTDMPFUNC str2ztdmpfunc(const string &mf);
  IONMPFUNC str2ionmpfunc(const string &mf);
  OBSWEIGHT str2obsweight(const string &wg);
  TROPMODEL str2tropmodel(const string &tm);
  static ZTDMODEL str2ztdmodel(const string &ztd);
  RESIDTYPE str2residtype(const string &rs);
  OBSCOMBIN str2obscombin(const string &oc);
  ATTITUDES str2attitudes(const string &at);
  CBIASCHAR str2cbiaschar(const string &cb);
  SYSBIASMODEL str2sysbiasmodel(const string &sys);
  static GRDMODEL str2grdmodel(const string &grd);
  /**@brief str convert */
  static string grdmpfunc2str(GRDMPFUNC MF);
  static string ztdmpfunc2str(ZTDMPFUNC MF);
  static string ionmpfunc2str(IONMPFUNC MF);
  static string obsweight2str(OBSWEIGHT WG);
  static string tropmodel2str(TROPMODEL TM);
  static string ztdmodel2str(ZTDMODEL ZTD);
  static string residtype2str(RESIDTYPE RS);
  static string obscombin2str(OBSCOMBIN OC);
  static string attitude2str(ATTITUDES AT);
  static string cbiaschar2str(CBIASCHAR CB);

  bool simulation();

  IFCB_MODEL ifcb_model();
  IFB_MODEL ifb_model();
  bool trimcor();
  bool write_equ();

  IONO_ORDER iono_order();

  BASEPOS basepos();
  static BASEPOS str2basepos(const string &str);
  static string basepos2str(BASEPOS BP);
  // smooth range
  string range_smooth_mode(int *smt_windows = nullptr);
  // correct bds code bias
  bool bds_code_bias_correction();

  int clk_init_epo();
  bool ambupd();

 protected:

  bool _phase;                                ///< use phase data [true=1, false=0]
  bool _tropo;                                ///< estimate troposphere [true=1, false=0]
  bool _iono;                                 ///< estimate ionosphere [true=1, false=0]
  bool _tropo_grad;                           ///< estimate troposphere gradinet [true=1, false=0]
  bool _tropo_slant;                          ///< estimate tropo slant delays
  TROPMODEL _tropo_model;                       ///< tropospheric model [SAASTAMOINEN, DAVIS, HOPFIELD, ...]
  ZTDMODEL _ztd_model;                          ///< ztdmodle [PWC STO]
  ZTDMPFUNC _tropo_mf;                          ///< tropo mapping function [COSZ, NIELL, GMF, VMF1, ... ]
  IONMPFUNC _iono_mf;                           ///< iono mapping function [COSZ, QFAC, NONE, ... ]
  GRDMPFUNC _grad_mf;                           ///< grad mapping function [tilt, CHH]
  OBSWEIGHT _obs_weight;                        ///< observations weighting
  RESIDTYPE _res_type;                          ///< types of residuals
  OBSCOMBIN _obs_combin;                        ///< observation combination
  ATTITUDES _attitudes;                          ///< satellite attitude model
  CBIASCHAR _cbiaschar;                          ///< forcing code bias signals

  double _sig_init_ztd;                         ///< accuracy of initial zenith path delays [m]
  double _sig_init_vion;                        ///< accuracy of initial vertical iono path delays [m]
  double _sig_init_grd;                         ///< accuracy of initial tropo gradients [m]
  double _sig_init_crd;                         ///< accuracy of initial coordinates [m]
  double _sig_init_vel;                         ///< accuracy of initial velocities [m/s]
  double _sig_init_amb;                         ///< accuracy of initial ambiguity [m]
  double _sig_init_glo;                         ///< accuracy of initial GLONASS system time difference
  double _sig_init_gal;                         ///< accuracy of initial Galileo system time difference
  double _sig_init_bds;                         ///< accuracy of initial BeiDou system time difference
  double _sig_init_qzs;                         ///< accuracy of initial QZSS system time difference
  double _minimum_elev;                         ///< elevation angle cut-off [degree]
  double _max_res_norm;                         ///< normalized residuals threshold
  string _crd_est;                              ///< FIX or estimate CRD
  bool _pos_kin;                              ///< static/kinematic receiver (true == kinematic)
  bool _auto_band;                            ///< automatic band selection (mainly for Anubis purpose)
  int _dynamics;
  int _frequency;
  bool _simulation;                            ///< judge whether simulation data

  int _minsat;                                 ///<  minimum satellite number
  BASEPOS _basepos;
  bool _sd_sat;                              ///< single differented between sat and sat_ref
 private:
};

} // namespace

#endif
