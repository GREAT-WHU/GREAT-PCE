/*
*
* @verbatim
     (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file        gsetproc.cpp
* @brief       implements process setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#include <iostream>
#include <algorithm>

#include "gset/gsetproc.h"
#include "gproc/gsppflt.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace pugi;

#define XMLKEY_PROC_IFCB "ifcb_model"
#define XMLKEY_PROC_ION_ORDER  "iono_order"
#define XMLKEY_PROC_IFB  "ifb_model"
#define XMLKEY_PROC_READ_OFILE "read_ofile_mode"
#define XMLKEY_PROC_TRIMCOR "trimcor"
#define XMLKEY_PROC_EQU "write_equ"
using namespace great;

namespace gnut {
// Constructor
// ----------
t_gsetproc::t_gsetproc()
    : t_gsetbase() {
  _set.insert(XMLKEY_PROC);
  _phase = true;
  _tropo = true;
  _iono = true;
  _tropo_grad = false;
  _tropo_slant = false;
  _tropo_model = TROPMODEL::SAASTAMOINEN;
  _ztd_model = ZTDMODEL::PWC;
  _tropo_mf = ZTDMPFUNC::GMF;
  _iono_mf = IONMPFUNC::ICOSZ;
  _grad_mf = GRDMPFUNC::TILTING;
  _obs_weight = OBSWEIGHT::SINEL2;
  _res_type = RESIDTYPE::RES_NORM;
  _obs_combin = OBSCOMBIN::IONO_FREE;
  _attitudes = ATTITUDES::DEF_YAWMODEL;
  _cbiaschar = CBIASCHAR::ORIG;

  _sig_init_ztd = 0.1;
  _sig_init_vion = 10;
  _sig_init_grd = 0.0005;
  _sig_init_amb = 1000.0;
  _sig_init_crd = 100.0;
  _sig_init_vel = 10.0;
  _sig_init_glo = 1000.0;
  _sig_init_gal = 1000.0;
  _sig_init_bds = 1000.0;
  _sig_init_qzs = 1000.0;
  _minimum_elev = 10;
  _max_res_norm = 10;
  _crd_est = "EST";
  _pos_kin = false;
  _auto_band = false;
  _frequency = 2;
  _sd_sat = false;
  _simulation = false;
  _dynamics = 0;
  _basepos = BASEPOS::SPP;
  _minsat = static_cast<size_t>(SPP_MINSAT);

}

// Destructor
// ----------
t_gsetproc::~t_gsetproc() = default;


// Return value
// ----------
bool t_gsetproc::tropo() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::iono() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::tropo_slant() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_slant").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::tropo_grad() {
  _gmutex.lock();
  // temporary for legacy mode with keyword "gradient"
  bool tmp = false;
  bool tmp1 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_grad").as_bool();
  bool tmp2 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("gradient").as_bool();   // legacy mode
  if (tmp1 || tmp2) tmp = true;
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::phase() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("phase").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::doppler() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("doppler").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::pos_kin() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("pos_kin").as_bool();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
bool t_gsetproc::auto_band() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("auto_band").as_bool();
  _gmutex.unlock();
  return tmp;
}

BASEPOS t_gsetproc::basepos() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("basepos").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  BASEPOS type = str2basepos(tmp);
  _gmutex.unlock();
  return type;
}

BASEPOS t_gsetproc::str2basepos(const string &str) {
  if (str == "RINEXH") return BASEPOS::RINEXH;
  else if (str == "CFILE") return BASEPOS::CFILE;
  else return BASEPOS::SPP;
}

// Return value
// ----------
double t_gsetproc::sig_init_ztd() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_ztd").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_vion() {
  _gmutex.lock();

  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_vion").as_double();

  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_grd() {
  _gmutex.lock();

  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_grd").as_double();

  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_crd() {
  _gmutex.lock();

  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_crd").as_double();

  _gmutex.unlock();
  return tmp;
}

double t_gsetproc::sig_init_vel() {
  _gmutex.lock();

  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_vel").as_double();

  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_amb() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_amb").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_glo() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_glo").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_gal() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_gal").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_bds() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_bds").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_qzs() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_qzs").as_double();
  _gmutex.unlock();
  return tmp;
}

double t_gsetproc::sig_ref_clk() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_ref_clk").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::minimum_elev() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("minimum_elev").as_double();
  _gmutex.unlock();
  return tmp;
}

// Return value
// ----------
double t_gsetproc::max_res_norm() {
  _gmutex.lock();
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("max_res_norm").as_double();
  _gmutex.unlock();
  return tmp;
}

int t_gsetproc::minsat() {
  _gmutex.lock();
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("min_sat").as_int();
  _gmutex.unlock();
  return tmp;
}

string t_gsetproc::ref_clk() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ref_clk").value();
  _gmutex.unlock();
  return tmp;
}

string t_gsetproc::ref_rec() {
  return string();
}

LSQMODE t_gsetproc::lsq_mode() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("lsq_mode").value();
  _gmutex.unlock();
  if (tmp == "LSQ") return LSQMODE::LSQ;
  else return LSQMODE::EPO;
}

int t_gsetproc::lsq_buffer_size() {
  _gmutex.lock();
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("lsq_buffer_size").as_int(10);
  _gmutex.unlock();
  return tmp;
}

SLIPMODEL t_gsetproc::slip_model() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("slip_model").value();
  _gmutex.unlock();
  if (tmp == "turboedit") return SLIPMODEL::TURBO_EDIT;
  else return SLIPMODEL::DEF_DETECT_MODEL;
}

IONMODEL t_gsetproc::ion_model() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ion_model").value();
  _gmutex.unlock();
  if (tmp == "SION") return IONMODEL::SION;
  else if (tmp == "VION") return IONMODEL::VION;
  else return IONMODEL::DEF_ION;
}

// Return value
// ----------
TROPMODEL t_gsetproc::tropo_model() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_model").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  TROPMODEL TM = str2tropmodel(tmp);
  if (TM == TROPMODEL::DEF_TROPMODEL) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = tropmodel2str(_tropo_model);
    _default_attr(node, "tropo_model", tmp);
    TM = _tropo_model;
  }
  _gmutex.unlock();
  return TM;
}

ZTDMODEL t_gsetproc::ztd_model(double *dt) {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ztd_model").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  *dt = 0.0;
  ZTDMODEL TM = str2ztdmodel(tmp);

  if (TM == ZTDMODEL::PWC) {
    int loc = tmp.find(':');
    *dt = str2dbl(tmp.substr(loc + 1));
  } else if (TM == ZTDMODEL::DEF_ZTDMODEL) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = ztdmodel2str(_ztd_model);
    _default_attr(node, "tropo_model", tmp);

    TM = _ztd_model;
    *dt = 60;
  }

  _gmutex.unlock();
  return TM;
}

// Return value
// ----------
ATTITUDES t_gsetproc::attitudes() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("attitudes").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  _gmutex.unlock();
  return str2attitudes(tmp);
}

// Return value
// ----------
CBIASCHAR t_gsetproc::cbiaschar() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("cbiaschar").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  _gmutex.unlock();
  return str2cbiaschar(tmp);
}

// Return value
// ----------
SYSBIASMODEL t_gsetproc::sysbias_model(double& dt) {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sysbias_model").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  _gmutex.unlock();
  SYSBIASMODEL model = str2sysbiasmodel(trim(tmp));
  if (model == SYSBIASMODEL::AUTO_PWC || model == SYSBIASMODEL::ISB_PWC) {
      int loc = tmp.find(':');
      if (loc > 0) dt = str2dbl(tmp.substr(loc + 1));
      else dt = 24;
  }
  return model;
}

// Return value
// ----------
SYSBIASMODEL t_gsetproc::ifbbias_model() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ifbbias_model").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  _gmutex.unlock();
  return str2sysbiasmodel(trim(tmp));

}

// Return value
// ----------
GRDMPFUNC t_gsetproc::grad_mf() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("grad_mf").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  GRDMPFUNC MF = str2grdmpfunc(tmp);
  if (MF == GRDMPFUNC::DEF_GRDMPFUNC) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = grdmpfunc2str(_grad_mf);
    _default_attr(node, "grad_mf", tmp);

    MF = _grad_mf;
  }

  _gmutex.unlock();
  return MF;
}

// Return value
// ----------
OBSWEIGHT t_gsetproc::weighting() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_weight").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  OBSWEIGHT WG = str2obsweight(tmp);
  if (WG == OBSWEIGHT::DEF_OBSWEIGHT) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = obsweight2str(_obs_weight);
    _default_attr(node, "obs_weight", tmp);
    WG = _obs_weight;
  }

  _gmutex.unlock();
  return WG;
}

// Return value
// ----------
RESIDTYPE t_gsetproc::residuals() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("residuals").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  RESIDTYPE RS = str2residtype(tmp);
  if (RS == RESIDTYPE::DEF_RESIDTYPE) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = residtype2str(_res_type);
    _default_attr(node, "residuals", tmp);
    RS = _res_type;
  }

  _gmutex.unlock();
  return RS;
}

// Return value
// ----------
OBSCOMBIN t_gsetproc::obs_combin() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_combination").value();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  OBSCOMBIN OC = str2obscombin(tmp);
  if (OC == OBSCOMBIN::DEF_OBSCOMBIN) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = obscombin2str(_obs_combin);
    _default_attr(node, "obs_combination", tmp);
    OC = _obs_combin;
  }

  _gmutex.unlock();
  return OC;
}

// Return value
// ----------
ZTDMPFUNC t_gsetproc::tropo_mf() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_mf").value();
  _gmutex.unlock();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  ZTDMPFUNC MF = str2ztdmpfunc(tmp);
  if (MF == ZTDMPFUNC::DEF_ZTDMPFUNC) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = ztdmpfunc2str(_tropo_mf);
    _default_attr(node, "tropo_mf", tmp);
    MF = _tropo_mf;
  }

  return MF;
}

// Return value
// ----------
IONMPFUNC t_gsetproc::iono_mf() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono_mf").value();
  _gmutex.unlock();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  IONMPFUNC MF = str2ionmpfunc(tmp);
  if (MF == IONMPFUNC::DEF_IONMPFUNC) {
    xml_node parent = _doc.child(XMLKEY_ROOT);
    xml_node node = _default_node(parent, XMLKEY_PROC);
    tmp = ionmpfunc2str(_iono_mf);
    _default_attr(node, "iono_mf", tmp);
    MF = _iono_mf;
  }

  return MF;
}

// Return value
// ----------
CONSTRPAR t_gsetproc::crd_est() {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("crd_constr").value();
  _gmutex.unlock();

  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  if (tmp.empty()) tmp = _crd_est; // default

  if (tmp == "EST") return CONSTRPAR::EST;
  else if (tmp == "KIN") return CONSTRPAR::KIN;
  else if (tmp == "SIMU_KIN") return CONSTRPAR::SIMU_KIN;
  else return CONSTRPAR::FIX;
}

// settings check
// ----------
void t_gsetproc::check() {
  _gmutex.lock();

  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_PROC);

  // check existence of attributes
  _default_attr(node, "phase", _phase);
  _default_attr(node, "tropo", _tropo);
  _default_attr(node, "iono", _iono);
  _default_attr(node, "tropo_grad", _tropo_grad);
  _default_attr(node, "tropo_model", tropmodel2str(_tropo_model));
  _default_attr(node, "tropo_slant", _tropo_slant);
  _default_attr(node, "sig_init_crd", _sig_init_crd);
  _default_attr(node, "sig_init_vel", _sig_init_vel);
  _default_attr(node, "sig_init_ztd", _sig_init_ztd);
  _default_attr(node, "sig_init_grd", _sig_init_grd);
  _default_attr(node, "sig_init_vion", _sig_init_vion);
  _default_attr(node, "sig_init_amb", _sig_init_amb);
  _default_attr(node, "sig_init_glo", _sig_init_glo);
  _default_attr(node, "sig_init_gal", _sig_init_gal);
  _default_attr(node, "sig_init_bds", _sig_init_bds);
  _default_attr(node, "sig_init_qzs", _sig_init_qzs);
  _default_attr(node, "minimum_elev", _minimum_elev);
  _default_attr(node, "max_res_norm", _max_res_norm);
  _default_attr(node, "crd_constr", _crd_est);
  _default_attr(node, "pos_kin", _pos_kin);
  _default_attr(node, "sd_sat", _sd_sat);
  _default_attr(node, "min_sat", _minsat);
  _default_attr(node, "auto_band", _auto_band);
  _default_attr(node, "simulation", _simulation);
  _default_attr(node, "basepos", basepos2str(_basepos));

  _gmutex.unlock();
}

// help body
// ----------
void t_gsetproc::help() {
  _gmutex.lock();

  cerr << " <process \n"
       << "   phase=\"" << _phase << "\" \n"
       << "   tropo=\"" << _tropo << "\" \n"
       << "   iono=\"" << _iono << "\" \n"
       << "   tropo_grad=\"" << _tropo_grad << "\" \n"
       << "   tropo_model=\"" << tropmodel2str(_tropo_model) << "\" \n"
       << "   tropo_slant=\"" << _tropo_slant << "\" \n"
       << "   tropo_mf=\"" << ztdmpfunc2str(_tropo_mf) << "\" \n"
       << "   iono_mf=\"" << ionmpfunc2str(_iono_mf) << "\" \n"
       << "   grad_mf=\"" << grdmpfunc2str(_grad_mf) << "\" \n"
       << "   obs_weight=\"" << obsweight2str(_obs_weight) << "\" \n"
       << "   residuals=\"" << residtype2str(_res_type) << "\" \n"
       << "   obs_combination=\"" << obscombin2str(_obs_combin) << "\" \n"
       << "   sig_init_crd=\"" << _sig_init_crd << "\" \n"
       << "   sig_init_ztd=\"" << _sig_init_ztd << "\" \n"
       << "   sig_init_vion=\"" << _sig_init_vion << "\" \n"
       << "   minimum_elev=\"" << _minimum_elev << "\" \n"
       << "   max_res_norm=\"" << _max_res_norm << "\" \n"
       << "   basepos=\"" << basepos2str(_basepos) << "\" \n"
       << " />\n";

  cerr << "\t<!-- process description:\n"
       << "\t phase [0,1]   .. use carrier phase data\n"
       << "\t tropo [0,1]   .. estimate troposphere\n"
       << "\t tropo_grad    .. tropospheric horizontal gradient models\n"
       << "\t tropo_model   .. tropospheric model (SAAS, GPT, ...)\n"
       << "\t tropo_slant   .. tropospheric slant delays produced\n"
       << "\t tropo_mf      .. tropospheric mapping function (COSZ, GMF, ...)\n"
       << "\t grad_mf       .. tropo gradient mapping function (TILTING,CHEN_HERRING, BAR_SEVER ...)\n"
       << "\t obs_weight    .. observation elevation dependant weighting (EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH)\n"
       << "\t residuals     .. type of residuals (RES_ORIG, RES_NORM, RES_ALL)\n"
       << "\t sig_init_crd  .. accuracy of initial coordinates [m]\n"
       << "\t sig_init_ztd  .. accuracy of initial zenith path delay [m]\n"
       << "\t sig_init_vion .. accuracy of initial vertical iono path delay [m]\n"
       << "\t minimum_elev  .. elevation angle cut-off [degree]\n"
       << "\t max_res_norm  .. maximal normalized residuals\n"
       << "\t basepos  .. base site coordinate\n"
       << "\t -->\n\n";

  _gmutex.unlock();
}

// convert str to ATTITUDES enum
// ----------
ATTITUDES t_gsetproc::str2attitudes(const string &ati) {
  if (ati == "NOMINAL") {
    return ATTITUDES::YAW_NOMI;
  } else if (ati == "RTCM") {
    return ATTITUDES::YAW_RTCM;
  } else {
    stringstream ostr;
    ostr << "Unsupported attitude model (" << attitude2str(_attitudes) << ")! Used default value ( YAW_MODEL )";
    _add_log("gsetproc", ostr.str());
    return ATTITUDES::DEF_YAWMODEL;
  }
}

// convert str to GRDMPFUNC enum
// ----------
GRDMPFUNC t_gsetproc::str2grdmpfunc(const string &mf) {
  if (mf == "TILTING") {
    return GRDMPFUNC::TILTING;
  } else if (mf == "CHEN_HERRING") {
    return GRDMPFUNC::CHEN_HERRING;
  } else if (mf == "BAR_SEVER") {
    return GRDMPFUNC::BAR_SEVER;
  } else {
    stringstream ostr;
    ostr << "Unsupported GRD mapping function (" << grdmpfunc2str(_grad_mf) << ")! Used default value ("
         << grdmpfunc2str(_grad_mf) << ")";
    _add_log("gsetproc", ostr.str());
    return GRDMPFUNC::DEF_GRDMPFUNC;
  }
}

// convert str to CBIASCHAR enum
// ----------
CBIASCHAR t_gsetproc::str2cbiaschar(const string &cb) {
  if (cb == "2CHAR") {
    return CBIASCHAR::CHAR2;
  } else if (cb == "3CHAR") {
    return CBIASCHAR::CHAR3;
  } else if (cb == "ORIG") {
    return CBIASCHAR::ORIG;
  } else {
    stringstream ostr;
    ostr << "Unsupported forcing code bias signals! Used default value (" << cbiaschar2str(_cbiaschar) << ")";
    _add_log("gsetproc", ostr.str());
    return CBIASCHAR::DEF_CBIASCHAR;
  }
}

// convert str to SYSBIASMODEL enum
// ----------
SYSBIASMODEL t_gsetproc::str2sysbiasmodel(const string &sys) {
  if (sys == "AUTO+CON") {
    return SYSBIASMODEL::AUTO_CON;
  } else if (sys == "ISB+CON") {
    return SYSBIASMODEL::ISB_CON;
  } else if (sys == "AUTO+RWK") {
    return SYSBIASMODEL::AUTO_RWK;
  } else if (sys == "ISB+RWK") {
    return SYSBIASMODEL::ISB_RWK;
  } else if (sys == "AUTO+WHIT") {
      return SYSBIASMODEL::AUTO_WHIT;
  } else if (sys == "ISB+WHIT") {
      return SYSBIASMODEL::ISB_WHIT;
  } else if (sys.substr(0, 8) == "AUTO+PWC") {
      return SYSBIASMODEL::AUTO_PWC;
  } else if (sys.substr(0, 7) == "ISB+PWC") {
      return SYSBIASMODEL::ISB_PWC;
  } else {
    stringstream ostr;
    ostr << "Unsupported forcing code bias signals! Used default value (AUTO+WHIT)";
    _add_log("gsetproc", ostr.str());
    return SYSBIASMODEL::AUTO_WHIT;
  }
}

// convert str to ZTDMPFUNC enum
// ----------
ZTDMPFUNC t_gsetproc::str2ztdmpfunc(const string &mf) {
  if (mf == "COSZ") {
    return ZTDMPFUNC::COSZ;
  } else if (mf == "GMF") {
    return ZTDMPFUNC::GMF;
  } else if (mf == "NO_MF") {
    return ZTDMPFUNC::NO_MF;
  } else {
    stringstream ostr;
    ostr << "Unsupported ZTD mapping function (" << mf << ")! Used default value (" << ztdmpfunc2str(_tropo_mf) << ")";
    _add_log("gsetproc", ostr.str());
    return ZTDMPFUNC::DEF_ZTDMPFUNC;
  }
}

// convert str to IONMPFUNC enum
// ----------
IONMPFUNC t_gsetproc::str2ionmpfunc(const string &mf) {
  if (mf == "COSZ") {
    return IONMPFUNC::ICOSZ;
  } else if (mf == "QFAC") {
    return IONMPFUNC::QFAC;
  } else if (mf == "NONE") {
    return IONMPFUNC::NONE;
  } else {
    stringstream ostr;
    ostr << "Unsupported ION mapping function (" << mf << ")! Used default value (" << ionmpfunc2str(_iono_mf) << ")";
    _add_log("gsetproc", ostr.str());
    return IONMPFUNC::DEF_IONMPFUNC;
  }
}

// convert str to OBSWEIGHT enum
// ----------
OBSWEIGHT t_gsetproc::str2obsweight(const string &wg) {
  if (wg == "EQUAL") {
    return OBSWEIGHT::EQUAL;
  } else if (wg == "SINEL" || wg == "SINEL1") {
    return OBSWEIGHT::SINEL;
  } else if (wg == "SINEL2") {
    return OBSWEIGHT::SINEL2;
  } else if (wg == "SINEL4") {
    return OBSWEIGHT::SINEL4;
  } else if (wg == "CODPHA") {
    return OBSWEIGHT::CODPHA;
  } else if (wg == "MLTPTH") {
    return OBSWEIGHT::MLTPTH;
  } else if (wg == "PARTELE") {
    return OBSWEIGHT::PARTELE;
  } else if (wg == "SNR") {
    return OBSWEIGHT::SNR;
  } else {
    stringstream ostr;
    ostr << "Unsupported observation weighting model (" << wg << ")! Used default value (" << obsweight2str(_obs_weight)
         << ")";
    _add_log("gsetproc", ostr.str());
    return OBSWEIGHT::DEF_OBSWEIGHT;
  }
}

// convert str to TROPMODEL enum
// ----------
TROPMODEL t_gsetproc::str2tropmodel(const string &tm) {
  if (tm == "SAASTAMOINEN") {
    return TROPMODEL::SAASTAMOINEN;
  } else if (tm == "DAVIS") {
    return TROPMODEL::DAVIS;
  } else if (tm == "HOPFIELD") {
    return TROPMODEL::HOPFIELD;
  } else if (tm == "MOPS") {
    return TROPMODEL::MOPS;
  } else if (tm == "GPTW") {
    return TROPMODEL::GPTW;
  } else if (tm == "GPT2W") {
    return TROPMODEL::GPT2W;
  } else if (tm == "GAL27") {
    return TROPMODEL::GAL27;
  } else if (tm == "GALTROPO27") {
    return TROPMODEL::GALTROPO27;
  } else if (tm == "EXTERN") {
    return TROPMODEL::EXTERN;
  } else {
    stringstream ostr;
    ostr << "Unsupported tropospheric model (" << tm << ")! Used default value (" << tropmodel2str(_tropo_model) << ")";
    _add_log("gsetproc", ostr.str());
    return TROPMODEL::DEF_TROPMODEL;
  }
}

ZTDMODEL t_gsetproc::str2ztdmodel(const string &ztd) {
  if (ztd.substr(0, 3) == "PWC") return ZTDMODEL::PWC;
  else if (ztd.substr(0, 3) == "STO") return ZTDMODEL::STO;
  else return ZTDMODEL::DEF_ZTDMODEL;
}

// convert str to RESIDTYPE enum
// ----------
RESIDTYPE t_gsetproc::str2residtype(const string &rs) {
  if (rs == "RES_ORIG") {
    return RESIDTYPE::RES_ORIG;
  } else if (rs == "RES_NORM") {
    return RESIDTYPE::RES_NORM;
  } else if (rs == "RES_ALL") {
    return RESIDTYPE::RES_ALL;
  } else {
    stringstream ostr;
    ostr << "Unsupported type of residuals (" << rs << ")! Used default value (" << residtype2str(_res_type) << ")";
    _add_log("gsetproc", ostr.str());
    return RESIDTYPE::DEF_RESIDTYPE;
  }
}

// convert str to OBSCOMB enum
// ----------
OBSCOMBIN t_gsetproc::str2obscombin(const string &oc) {
  if (oc == "IONO_FREE") {
    return OBSCOMBIN::IONO_FREE;
  } else if (oc == "RAW_SINGLE") {
    return OBSCOMBIN::RAW_SINGLE;
  } else if (oc == "RAW_DOUBLE") {
    return OBSCOMBIN::RAW_DOUBLE;
  } else if (oc == "RAW_ALL") {
    return OBSCOMBIN::RAW_ALL;
  } else if (oc == "MW_COMBIN") {
    return OBSCOMBIN::MW_COMBIN;
  } else if (oc == "WL_COMBIN") {
    return OBSCOMBIN::WL_COMBIN;
  } else if (oc == "RAW_MIX") {
    return OBSCOMBIN::RAW_MIX;
  } else if (oc == "IF_P1") {
    return OBSCOMBIN::IF_P1;
  } else {
    stringstream ostr;
    ostr << "Unsupported observations combination (" << oc << ")! Used default value (" << obscombin2str(_obs_combin)
         << ")";
    _add_log("gsetproc", ostr.str());
    return OBSCOMBIN::DEF_OBSCOMBIN;
  }
}

// convert GRDMPFUNC enum to str
// ---------
string t_gsetproc::grdmpfunc2str(GRDMPFUNC MF) {
  switch (MF) {
    case GRDMPFUNC::TILTING: return "TILTING";
    case GRDMPFUNC::CHEN_HERRING: return "CHEN_HERRING";
    case GRDMPFUNC::BAR_SEVER: return "BAR_SEVER";
    case GRDMPFUNC::DEF_GRDMPFUNC: return "NOT DEFINED";
    default: return "";
  }
}

// convert ZTDMPFUNC enum to str
// ---------
string t_gsetproc::ztdmpfunc2str(ZTDMPFUNC MF) {
  switch (MF) {
    case ZTDMPFUNC::COSZ: return "COSZ";
    case ZTDMPFUNC::GMF: return "GMF";
    case ZTDMPFUNC::NO_MF: return "NO_MF";
    case ZTDMPFUNC::DEF_ZTDMPFUNC: return "NOT DEFINED";
    default: return "";
  }
}

// convert IONMPFUNC enum to str
// ---------
string t_gsetproc::ionmpfunc2str(IONMPFUNC MF) {
  switch (MF) {
    case IONMPFUNC::ICOSZ: return "ICOSZ";
    case IONMPFUNC::QFAC: return "QFAC";
    case IONMPFUNC::NONE: return "NONE";
    case IONMPFUNC::DEF_IONMPFUNC: return "NOT DEFINED";
    default: return "";
  }
}

// convert OBSWEIGHT enum to str
// ---------
string t_gsetproc::obsweight2str(OBSWEIGHT WG) {
  switch (WG) {
    case OBSWEIGHT::EQUAL: return "EQUAL";
    case OBSWEIGHT::SINEL: return "SINEL";
    case OBSWEIGHT::SINEL2: return "SINEL2";
    case OBSWEIGHT::SINEL4: return "SINEL4";
    case OBSWEIGHT::CODPHA: return "CODPHA";
    case OBSWEIGHT::MLTPTH: return "MLTPTH";
    case OBSWEIGHT::PARTELE: return "PARTELE";
    case OBSWEIGHT::SNR: return "SNR";
    case OBSWEIGHT::DEF_OBSWEIGHT: return "NOT DEFINED";
    default: return "";
  }
}

// convert TROPMODEL enum to str
// ---------
string t_gsetproc::tropmodel2str(TROPMODEL TM) {
  switch (TM) {
    case TROPMODEL::SAASTAMOINEN: return "SAASTAMOINEN";
    case TROPMODEL::DAVIS: return "DAVIS";
    case TROPMODEL::HOPFIELD: return "HOPFIELD";
    case TROPMODEL::MOPS: return "MOPS";
    case TROPMODEL::GPTW: return "GPTW";
    case TROPMODEL::GPT2W: return "GPT2W";
    case TROPMODEL::GAL27: return "GAL27";
    case TROPMODEL::GALTROPO27: return "GALTROPO27";
    case TROPMODEL::EXTERN: return "EXTERN";
    case TROPMODEL::DEF_TROPMODEL: return "NOT DEFINED";
    default: return "";
  }
}

// convert ZTD MODEL enum to str
// ---------
string t_gsetproc::ztdmodel2str(ZTDMODEL ZTD) {
  switch (ZTD) {
    case ZTDMODEL::DEF_ZTDMODEL: return "NOT DEFINED";
    case ZTDMODEL::PWC: return "PWC";
    case ZTDMODEL::STO: return "STO";
    default: return "";
  }
}

// convert RESIDTYPE enum to str
// ---------
string t_gsetproc::residtype2str(RESIDTYPE RS) {
  switch (RS) {
    case RESIDTYPE::RES_ORIG: return "RES_ORIG";
    case RESIDTYPE::RES_NORM: return "RES_NORM";
    case RESIDTYPE::RES_ALL: return "RES_ALL";
    case RESIDTYPE::DEF_RESIDTYPE: return "NOT DEFINED";
    default: return "";
  }
}

string t_gsetproc::attitude2str(ATTITUDES AT) {
  switch (AT) {
    case ATTITUDES::YAW_NOMI: return "YAW_NOMI";
    case ATTITUDES::YAW_RTCM: return "YAW_RTCM";
    case ATTITUDES::DEF_YAWMODEL: return "DEF_YAWMODEL";
    default: return "";
  }
}

// convert OBSCOMB enum to str
// ---------
string t_gsetproc::obscombin2str(OBSCOMBIN OC) {
  switch (OC) {
    case OBSCOMBIN::IONO_FREE: return "IONO_FREE";
    case OBSCOMBIN::RAW_SINGLE: return "RAW_SINGLE";
    case OBSCOMBIN::RAW_DOUBLE: return "RAW_DOUBLE";
    case OBSCOMBIN::RAW_ALL: return "RAW_ALL";
    case OBSCOMBIN::DEF_OBSCOMBIN: return "NOT DEFINED";
    case OBSCOMBIN::MW_COMBIN: return "MW_COMBIN";
    case OBSCOMBIN::RAW_MIX: return "RAW_MIX";
    case OBSCOMBIN::IF_P1: return "IF_P1";
    case OBSCOMBIN::WL_COMBIN: return "WL_COMBIN";
    case OBSCOMBIN::EWL_COMBIN: return "EWL_COMBIN";
    default: return "";
  }
}

// convert CBIASCHAR enum to str
// ---------
string t_gsetproc::cbiaschar2str(CBIASCHAR CB) {
  switch (CB) {
    case CBIASCHAR::CHAR2: return "2CHAR";
    case CBIASCHAR::CHAR3: return "3CHAR";
    case CBIASCHAR::ORIG: return "ORIG";
    case CBIASCHAR::DEF_CBIASCHAR: return "NOT DEFINED";
    default: return "";
  }
}

string t_gsetproc::basepos2str(BASEPOS BP) {
  switch (BP) {
    case BASEPOS::RINEXH: return "RINEXH";
    case BASEPOS::CFILE: return "CFILE";
    case BASEPOS::SPP: return "SPP";
    default: return "";
  }
}

t_gobscombtype::t_gobscombtype() :
    _obs_type(GOBSTYPE::TYPE),
    _obs_band_1(GOBSBAND::BAND),
    _obs_combine(OBSCOMBIN::DEF_OBSCOMBIN) {
}

t_gobscombtype::t_gobscombtype(const t_gobscombtype &other) :
    _obs_type(other._obs_type),
    _obs_band_1(other._obs_band_1),
    _obs_band_2(other._obs_band_2),
    _obs_freq_1(other._obs_freq_1),
    _obs_freq_2(other._obs_freq_2),
    _obs_combine(other._obs_combine) {
}

t_gobscombtype::t_gobscombtype(const string &obscombtype) :
    _obs_type(GOBSTYPE::TYPE),
    _obs_band_1(GOBSBAND::BAND),
    _obs_combine(OBSCOMBIN::DEF_OBSCOMBIN) {
  if (obscombtype.length() < 2) {
    return;
  }
  string str_type = obscombtype.substr(0, 1);
  if (str_type == "P") str_type = "C";
  _obs_type = str2gobstype(str_type);

  string str_band = obscombtype.substr(1);
  if (str_band.substr(0, 1) == "C") {
    _obs_combine = OBSCOMBIN::IONO_FREE;
    _obs_freq_1 = str2gnssfreq(obscombtype.substr(2, 1));
    _obs_freq_2 = str2gnssfreq(obscombtype.substr(3, 1));
  } else {
    _obs_freq_1 = str2gnssfreq(obscombtype.substr(1, 1));
  }
}

t_gobscombtype::t_gobscombtype(const t_gobs &obstype, OBSCOMBIN combtype) :
    _obs_type(obstype.type()),
    _obs_band_1(obstype.band()),
    _obs_combine(combtype) {
}

t_gobscombtype::t_gobscombtype(const t_gobs& obstype, GOBSBAND b1, FREQ_SEQ freq_1, OBSCOMBIN combtype) :
    _obs_type(obstype.type()),
    _obs_band_1(b1),
    _obs_freq_1(freq_1),
    _obs_combine(combtype)
{
}

t_gobscombtype::t_gobscombtype(const t_gobs &obstype,
                               GOBSBAND b1,
                               GOBSBAND b2,
                               FREQ_SEQ freq_1,
                               FREQ_SEQ freq_2,
                               OBSCOMBIN combtype) :
    _obs_type(obstype.type()),
    _obs_band_1(b1),
    _obs_band_2(b2),
    _obs_freq_1(freq_1),
    _obs_freq_2(freq_2),
    _obs_combine(combtype) {
}

t_gobscombtype::t_gobscombtype(GOBSTYPE t, GOBSBAND b, OBSCOMBIN obscomb) :
    _obs_type(t),
    _obs_band_1(b),
    _obs_combine(obscomb) {
}

string t_gobscombtype::convert2str() const {
  string ans(gobstype2str(_obs_type == TYPE_C ? TYPE_P : _obs_type));
  if (_obs_combine == OBSCOMBIN::IONO_FREE) {
    ans += "C";
    if (_obs_freq_1 != FREQ_X) ans += gfreqseq2str(_obs_freq_1);
    if (_obs_freq_2 != FREQ_X) ans += gfreqseq2str(_obs_freq_2);
  } else {
    if (_obs_freq_1 != FREQ_X) ans += gfreqseq2str(_obs_freq_1);
  }
  return ans;
}

bool t_gobscombtype::operator==(const t_gobscombtype &g) const {
  return this->convert2str() == g.convert2str();
}

bool t_gobscombtype::operator<(const t_gobscombtype &g) const {
  return this->convert2str() < g.convert2str();
}

bool t_gobscombtype::is_phase() const {
  return _obs_type == TYPE_L;
}

bool t_gobscombtype::is_code() const {
  return (_obs_type == TYPE_C || _obs_type == TYPE_P);
}

bool t_gsetproc::simulation() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("simulation").as_bool();
  _gmutex.unlock();
  return tmp;
}

IFCB_MODEL t_gsetproc::ifcb_model() {
  //get rkf model node
  string ifcb = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_IFCB);

  if (ifcb.empty()) return IFCB_MODEL::COR;

  // delete spaces
  ifcb.erase(remove(ifcb.begin(), ifcb.end(), ' '), ifcb.end());
  transform(ifcb.begin(), ifcb.end(), ifcb.begin(), ::toupper);

  if (ifcb == "EST") return IFCB_MODEL::EST;
  if (ifcb == "COR") return IFCB_MODEL::COR;

  return IFCB_MODEL::DEF;
}

IONO_ORDER t_gsetproc::iono_order() {
  string order = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_ION_ORDER);

  if (order.empty()) return IONO_ORDER::NONE;

  order.erase(remove(order.begin(), order.end(), ' '), order.end());
  transform(order.begin(), order.end(), order.begin(), ::toupper);

  if (order == "FIRST") return IONO_ORDER::FIRST_ORDER;
  if (order == "SECOND") return IONO_ORDER::SECOND_ORDER;
  if (order == "THIRD") return IONO_ORDER::THIRD_ORDER;
  if (order == "NONE") return IONO_ORDER::NONE;

  throw logic_error("ERROR type of iono order, you should use FIRST/SECOND/NONE");
}

string t_gsetproc::range_smooth_mode(int *smt_windows) {
  *smt_windows = 20;
  xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child("smooth_range");

  string tmp_mode = tmp_set.attribute("mode").value();
  if (tmp_mode.empty()) return "NONE";

  tmp_mode.erase(remove(tmp_mode.begin(), tmp_mode.end(), ' '), tmp_mode.end());
  transform(tmp_mode.begin(), tmp_mode.end(), tmp_mode.begin(), ::toupper);

  if (tmp_mode != "DOPPLER" && tmp_mode != "PHASE") return "NONE";

  string tmp_window = tmp_set.attribute("window").value();
  if (tmp_window.empty()) return tmp_mode;

  // delete spaces
  tmp_window.erase(remove(tmp_window.begin(), tmp_window.end(), ' '), tmp_window.end());
  int win_int = str2int(tmp_window);
  *smt_windows = win_int <= 0 ? 20 : win_int;

  return tmp_mode;
}

bool t_gsetproc::bds_code_bias_correction() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("bds_code_bias_corr").as_bool();
  _gmutex.unlock();
  return tmp;
}

int t_gsetproc::clk_init_epo()
{
    _gmutex.lock();
    int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("clk_init_epo").as_int();
    _gmutex.unlock();
    return tmp;
}

bool t_gsetproc::ambupd()
{
    _gmutex.lock();
    bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("ambupd").as_bool();
    _gmutex.unlock();
    return tmp;
}

IFB_MODEL t_gsetproc::ifb_model() {
  string ifb = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_IFB);

  if (ifb.empty()) return IFB_MODEL::NONE;

  ifb.erase(remove(ifb.begin(), ifb.end(), ' '), ifb.end());
  transform(ifb.begin(), ifb.end(), ifb.begin(), ::toupper);

  if (ifb == "EST_REC_IFB") return IFB_MODEL::EST_REC_IFB;
  if (ifb == "EST_SAT_IFB") return IFB_MODEL::EST_SAT_IFB;
  if (ifb == "EST_IFB") return IFB_MODEL::EST_IFB;
  if (ifb == "COR_REC_IFB") return IFB_MODEL::COR_REC_IFB;
  if (ifb == "COR_SAT_IFB") return IFB_MODEL::COR_SAT_IFB;
  if (ifb == "COR_IFB") return IFB_MODEL::COR_IFB;
  if (ifb == "NONE") return IFB_MODEL::NONE;

  throw logic_error(
      "ERROR type of ifb model, you should use EST_REC_IFB/EST_SAT_IFB/EST_IFB/COR_REC_IFB/COR_SAT_IFB/COR_IFB/NONE");
}

bool t_gsetproc::trimcor() {
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_TRIMCOR);

  if (tmp.empty()) return true;

  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  if (tmp == "TRUE") return true;
  if (tmp == "FALSE") return false;

  return false;
}

bool t_gsetproc::write_equ() {
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).child_value(XMLKEY_PROC_EQU);

  if (tmp.empty()) return true;

  tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  if (tmp == "TRUE") return true;
  if (tmp == "FALSE") return false;

  return true;
}

int t_gsetproc::frequency() {
  _gmutex.lock();
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("frequency").as_int(2);
  _gmutex.unlock();
  return tmp;
}

int t_gsetproc::num_threads() {
  _gmutex.lock();
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("num_threads").as_int(1);
  _gmutex.unlock();

#ifdef USE_OPENMP
  if (tmp == 0) {
      tmp = 1;
  }
  else if (tmp==-1) {
      tmp = omp_get_num_procs();
  }
#else
  // if no use openmp number of threads is 1
  tmp = 1;
#endif
  return tmp;
}

bool t_gsetproc::matrix_remove() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("matrix_remove").as_bool();
  _gmutex.unlock();
  return tmp;
}

bool t_gsetproc::cmb_equ_multi_thread() {
  _gmutex.lock();
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("cmb_equ_multi_thread").as_bool(true);
  _gmutex.unlock();

#ifdef USE_OPENMP
  return tmp;
#else
  return false;
#endif
}

GRDMODEL t_gsetproc::grd_model(double *dt) {
  _gmutex.lock();
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("grd_model").value();
  _gmutex.unlock();
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  *dt = 0.0;
  GRDMODEL TM = str2grdmodel(tmp);
  if (TM == GRDMODEL::GRD_PWC) {
    int loc = tmp.find(':');
    *dt = str2dbl(tmp.substr(loc + 1));
  } else if (TM == GRDMODEL::DEF_GRDMODEL) {
    TM = GRDMODEL::GRD_PWC;
    *dt = 120;
  }

  return TM;
}

GRDMODEL t_gsetproc::str2grdmodel(const string &grd) {
  string tmp = trim(grd);
  if (tmp.substr(0, 3) == "PWC") return GRDMODEL::GRD_PWC;
  else if (tmp.substr(0, 3) == "STO") return GRDMODEL::GRD_STO;
  else return GRDMODEL::DEF_GRDMODEL;
}

} // namespace
