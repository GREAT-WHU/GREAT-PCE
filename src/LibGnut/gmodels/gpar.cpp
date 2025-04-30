
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <cmath>

#include "newmat/newmat.h"

#include "gmodels/gpar.h"
#include "gmodels/ggmf.h"

using namespace std;

namespace gnut {

const string t_gpar::PAR_STR_SEP = "_";

const vector<string> t_gpar::PAR_STR =
    {
        "CRD_X", "CRD_Y", "CRD_Z",                                          // coordinates
        "TRP", "GRD_N", "GRD_E", "SION", "VION",                            // atmospheric parameters
        "CLK", "CLK_SAT",                                               // clocks
        "CLK_G", "CLK_E", "CLK_C", "CLK_R",
        "CLK_J",
        "CLK_ICB", "CLUSTERB",
        "AMB_IF",
        "AMB_WL", "AMB_L1", "AMB_L2","AMB_L3","AMB_L4","AMB_L5",
        "GLO_ISB", "GLO_IFB", "GLO_IFCB", "GLO_IFPB", "GAL_ISB", "BDS_ISB", "BD2_ISB", "QZS_ISB",
        "IFB_GPS", "IFB_BDS", "IFB_QZS", "IFB_GAL", "IFB_GAL_2", "IFB_GAL_3",
        "IFB_BDS_2", "IFB_BDS_3",
        "GPS_REC_IFB_C3",
        "VEL_X", "VEL_Y", "VEL_Z",
        "ATT_X", "ATT_Y", "ATT_Z",
        "eb_X", "eb_Y", "eb_Z",
        "db_X", "db_Y", "db_Z",
        "CAM_ATT_X", "CAM_ATT_Y", "CAM_ATT_Z",
        "CAM_CRD_X", "CAM_CRD_Y", "CAM_CRD_Z",
        "EXTRINSIC_ATT_X", "EXTRINSIC_ATT_Y", "EXTRINSIC_ATT_Z",
        "EXTRINSIC_CRD_X", "EXTRINSIC_CRD_Y", "EXTRINSIC_CRD_Z",
        "GLO_IFB","GLO_IFCB", "GLO_IFPB",
        "VTEC", "DCB_REC", "DCB_SAT", 
        "P1P2G_REC", "P1P2E_REC", "P1P2R_REC", "P1P2C_REC",
        "VEL_X", "VEL_Y", "VEL_Z",
        "ACC_X", "ACC_Y", "ACC_Z",
        "CLK_RAT",
        "PXSAT", "PYSAT", "PZSAT", "VXSAT", "VYSAT", "VZSAT",
        "SION_3"
    };

t_gparhead::t_gparhead(const par_type& type, const string& site, const string& sat, const int& channel):
    type(type),
    site(site),
    sat(sat),
    str_type(ptype2str(type)),
    channel(channel)
{
}
t_gparhead::t_gparhead(const t_gparhead &Other) :
    type(Other.type),
    sat(Other.sat),
    site(Other.site),
    str_type(Other.str_type),
    channel(Other.channel)
{
}
t_gparhead::~t_gparhead() {
}
bool t_gparhead::operator==(const t_gparhead &Other) const {
  if (this->type == Other.type &&
      this->sat == Other.sat &&
      this->site == Other.site &&
      this->channel == Other.channel) {
    return true;
  } else {
    return false;
  }
}
bool t_gparhead::operator<(const t_gparhead &Other) const {
  if (this->type < Other.type ||
      this->type == Other.type && this->site < Other.site ||
      this->type == Other.type && this->site == Other.site && this->sat < Other.sat ||
      this->type == Other.type && this->site == Other.site && this->sat == Other.sat && this->channel < Other.channel) {
    return true;
  } else {
    return false;
  }
}
bool t_gparhead::operator<=(const t_gparhead &Other) const {
  if (*this == Other || *this < Other) {
    return true;
  } else {
    return false;
  }
}
bool t_gparhead::operator>(const t_gparhead &Other) const {
  if (*this <= Other) {
    return false;
  } else {
    return true;
  }
}
bool t_gparhead::operator>=(const t_gparhead &Other) const {
  if (*this < Other) {
    return false;
  } else {
    return true;
  }
}

size_t t_gparhead_hash::operator()(const t_gparhead& a) const
{
    return std::hash<string>()(a.site + a.sat + a.str_type + int2str(a.channel));
}

t_gtimearc::t_gtimearc(const t_gtime &beg, const t_gtime &end) :
    begin(beg),
    end(end) {
}
t_gtimearc::~t_gtimearc() = default;
bool t_gtimearc::operator!=(const t_gtimearc &Other) const {
  return !(*this == Other);
}
bool t_gtimearc::operator==(const t_gtimearc &Other) const {
  return (this->begin == Other.begin && this->end == Other.end);
}
bool t_gtimearc::operator<(const t_gtimearc &Other) const {
  return (this->begin < Other.begin || this->begin == Other.begin && this->end < Other.end);
}
bool t_gtimearc::operator<=(const t_gtimearc &Other) const {
  return (*this == Other || *this < Other);
}
bool t_gtimearc::operator>(const t_gtimearc &Other) const {
  return !(*this <= Other);
}
bool t_gtimearc::operator>=(const t_gtimearc &Other) const {
  return !(*this < Other);
}
bool t_gtimearc::inside(const t_gtimearc &Other) const {
  return (this->begin >= Other.begin && this->end <= Other.end);
}

// constructor
// --------------------------------------------------------
t_gpar::t_gpar(const string &site, par_type t, unsigned i, const string &p, bool remove) {
  beg = FIRST_TIME;
  end = LAST_TIME;

  this->site = site;
  lremove = remove;
  parType = t;
  index = i;
  prn = p;
  amb_ini = false;
  //indPWC = 0;
  nPWC = 1;
  value(0.0);
  apriori(0.0);
}

t_gpar::t_gpar() {
  beg = FIRST_TIME;
  end = LAST_TIME;
}

// t_gpar destructor
// ---------------------------------------------------
t_gpar::~t_gpar() {
}

// setting mapping function for ZTD
// ---------------------------------------------------
void t_gpar::setMF(ZTDMPFUNC MF) {
  _mf_ztd = MF;
}

// setting mapping function for GRD
// ---------------------------------------------------
void t_gpar::setMF(GRDMPFUNC MF) {
  _mf_grd = MF;
}

// Partial derivatives 
// ----------------------------------------------------
double t_gpar::partial(t_gsatdata &satData, t_gtime &epoch, t_gtriple ground, t_gobs &gobs) {
  double mfw, dmfw, mfh, dmfh;
  mfw = dmfw = mfh = dmfh = 0.0;
  switch (parType) {
  case par_type::CRD_X:
      if (satData.site() == this->site) {
        return (value() - satData.satcrd().crd(0)) / satData.rho();
      } else return 0.0;
    case par_type::CRD_Y:
      if (satData.site() == this->site) {
        return (value() - satData.satcrd().crd(1)) / satData.rho();
      } else return 0.0;
    case par_type::CRD_Z:
      if (satData.site() == this->site) {
        return (value() - satData.satcrd().crd(2)) / satData.rho();
      } else return 0.0;
    case par_type::CLK: {
      if (satData.site() == this->site) return 1.0;
      else return 0.0;
    }
    case par_type::CLK_SAT: {
      if (satData.sat() == this->prn) return -1.0;
      else return 0.0;
    }
    case par_type::TRP   :
      if (satData.site() == this->site) {
        _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);

        return mfw;
      } else return 0.0;
    case par_type::SION  : {
      if (satData.site() == this->site) {
        double f1 = satData.frequency(t_gsys::band_priority(satData.gsys(), FREQ_1));
        double fk = satData.frequency(gobs.band());
        double alfa = 0.0;
        if (gobs.is_phase() && prn == satData.sat()) { alfa = -(f1 * f1) / (fk * fk); }
        if (gobs.is_code() && prn == satData.sat()) { alfa = (f1 * f1) / (fk * fk); }
        return alfa;
      } else return 0.0;
    }
    case par_type::VION  : {
      if (satData.site() == this->site) {
        double f1 = satData.frequency(t_gsys::band_priority(satData.gsys(), FREQ_1));
        double fk = satData.frequency(gobs.band());
        double mf = 1.0 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - satData.ele()), 2));
        double alfa = 0.0;
        if (gobs.is_phase() && prn == satData.sat()) { alfa = -(f1 * f1) / (fk * fk); }
        if (gobs.is_code() && prn == satData.sat()) { alfa = (f1 * f1) / (fk * fk); }
        return alfa * mf;
      } else return 0.0;
    }
    case par_type::P1P2G_REC : {
      if (satData.site() == this->site) {
        double f1 = G01_F;
        double fk = satData.frequency(gobs.band());
        double alfa = (f1 * f1) / (fk * fk);
        double beta = (G02_F * G02_F) / (G01_F * G01_F - G02_F * G02_F);
        FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());
        if (satData.gsys() == GPS && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
        else return 0.0;
      } else return 0.0;
    }
    case par_type::P1P2E_REC : {
      if (satData.site() == this->site) {
        double f1 = E01_F;
        double fk = satData.frequency(gobs.band());
        double alfa = (f1 * f1) / (fk * fk);
        double beta = (E05_F * E05_F) / (E01_F * E01_F - E05_F * E05_F);
        FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());
        if (satData.gsys() == GAL && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
        else return 0.0;
      } else return 0.0;
    }
    case par_type::GRD_N :
      if (satData.site() == this->site) {
        _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
        if (_mf_grd == GRDMPFUNC::CHEN_HERRING) {
          double sinel = sin(satData.ele());
          double tanel = tan(satData.ele());
          double cosaz = cos(satData.azi());
          return (1.0 / (sinel * tanel + 0.0032)) * cosaz;
        } else if (_mf_grd == GRDMPFUNC::TILTING) {
          double cosaz = cos(satData.azi());
          return dmfw * cosaz;
        } else if (_mf_grd == GRDMPFUNC::BAR_SEVER) {
          double tanel = tan(satData.ele());
          double cosaz = cos(satData.azi());
          return mfw * (1.0 / tanel) * cosaz;
        } else cerr << "Grad N mapping function is not set up correctly!!!" << endl;
      } else return 0.0;
    case par_type::GRD_E :
      if (satData.site() == this->site) {
        _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
        if (_mf_grd == GRDMPFUNC::CHEN_HERRING) {
          double sinel = sin(satData.ele());
          double tanel = tan(satData.ele());
          double sinaz = sin(satData.azi());
          return (1.0 / (sinel * tanel + 0.0032)) * sinaz;
        } else if (_mf_grd == GRDMPFUNC::TILTING) {
          double sinaz = sin(satData.azi());
          return dmfw * sinaz;
        } else if (_mf_grd == GRDMPFUNC::BAR_SEVER) {
          double tanel = tan(satData.ele());
          double sinaz = sin(satData.azi());
          return mfw * (1 / tanel) * sinaz;
        } else cerr << "Grad E mapping function is not set up correctly!!!" << endl;
      } else return 0.0;

    case par_type::AMB_IF  :
      if (satData.site() == this->site && gobs.is_phase() && prn == satData.sat()) return 1.0;
      else return 0.0;
    case par_type::AMB_L1  :
      if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_1
          && prn == satData.sat())
        return 1.0;
      else return 0.0;
    case par_type::AMB_L2  :
      if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_2
          && prn == satData.sat())
        return 1.0;
      else return 0.0;
    case par_type::AMB_L3:
        if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3
            && prn == satData.sat())
            return 1.0;
        else return 0.0;
    case par_type::AMB_L4:
        if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4
            && prn == satData.sat())
            return 1.0;
        else return 0.0;
    case par_type::AMB_L5:
        if (satData.site() == this->site && gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5
            && prn == satData.sat())
            return 1.0;
        else return 0.0;
    case par_type::GLO_ISB :
      if (satData.site() == this->site && satData.gsys() == GLO) return 1.0;
      else return 0.0;
    case par_type::GLO_IFB:
      if (satData.site() == this->site && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
      else return 0.0;
    case par_type::GLO_IFCB:
        if (satData.site() == this->site && !gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
        else return 0.0;
    case par_type::GLO_IFPB:
        if (satData.site() == this->site && gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
        else return 0.0;
    case par_type::GAL_ISB :
      if (satData.site() == this->site && satData.gsys() == GAL) return 1.0;
      else return 0.0;
    case par_type::BDS_ISB :
      if (satData.site() == this->site && satData.gsys() == BDS) return 1.0;
      else return 0.0;
    case par_type::BD2_ISB:
        if (satData.site() == this->site && satData.gsys() == BDS && satData.sat() < "C17") return 1.0;
        else return 0.0;
    case par_type::QZS_ISB :
      if (satData.site() == this->site && satData.gsys() == QZS) return 1.0;
      else return 0.0;
  }

  return 0.0;
}

double t_gpar::partial_doppler(t_gsatdata& satData, t_gtriple& groundXYZ, t_gtriple& groundVEL) {
    t_gtriple satcrd = satData.satcrd();
    t_gtriple satvel = satData.satvel();
    ColumnVector cSat = satcrd.crd_cvect();
    ColumnVector vSat = satvel.crd_cvect();
    ColumnVector cRec = groundXYZ.crd_cvect();
    ColumnVector vRec = groundVEL.crd_cvect();
    ColumnVector conf_crd = dotproduct(cSat - cRec, vSat - vRec) * (cSat - cRec) / pow(satData.rho(), 3);
    conf_crd -= (vSat - vRec) / satData.rho();
    ColumnVector e = (cSat - cRec) / satData.rho();
    satData._conf_crd.set(conf_crd);
    satData._e.set(e);
    switch (parType) {
    case par_type::CRD_X: return conf_crd(1);
    case par_type::CRD_Y: return conf_crd(2);
    case par_type::CRD_Z: return conf_crd(3);
    case par_type::VEL_X: return -e(1);
    case par_type::VEL_Y: return -e(2);
    case par_type::VEL_Z: return -e(3);
    case par_type::CLK_RAT: return 1.0;
    }
    return 0.0;
}

double t_gpar::partial_orb(t_gsatdata &satdata) {
  Matrix funct = satdata.orbfunct();

  int index_begin = satdata.satindex();

  // get unit_vector in TRS
  ColumnVector unit_rec2sat(3);

  unit_rec2sat = (satdata.satcrd() - satdata.reccrd()).crd_cvect();
  // change TRS to CRS
  unit_rec2sat = (satdata.rotmat()) * unit_rec2sat;
  unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

  funct = funct * unit_rec2sat;
  funct = funct * 1E3;


  return funct(index - index_begin + 1, 1);

}


// Operators for t_gpar
// -------------------------------------------------
bool t_gpar::operator==(const t_gpar &par) const {
  if (parType == par.parType &&
      beg == par.beg &&
      end == par.end) {
    return true;
  } else {
    return false;
  }

}

// Doplnil gabo aby par mohla byt pouzita ako klic v mape
bool t_gpar::operator<(const t_gpar &par) const {
  if ((parType < par.parType) ||
      (parType == par.parType && prn < par.prn) ||
      (parType == par.parType && prn == par.prn && site < par.site) ||
      (parType == par.parType && prn == par.prn && site == par.site && beg < par.beg) ||
      (parType == par.parType && prn == par.prn && site == par.site && beg == par.beg && end < par.end)
      ) {
    return true;
  }
  return false;
}

bool t_gpar::operator>(const t_gpar &par) const {
  if ((parType > par.parType) ||
      (parType == par.parType && prn < par.prn) ||
      (parType == par.parType && prn == par.prn && site < par.site) ||
      (parType == par.parType && prn == par.prn && site == par.site && beg < par.beg) ||
      (parType == par.parType && prn == par.prn && site == par.site && beg == par.beg && end < par.end)
      ) {
    return true;
  }
  return false;
}

t_gpar t_gpar::operator-(const t_gpar &p) const {
  t_gpar par = (*this);
  par.value(value() - p.value());
  return par;
}

t_gpar t_gpar::operator+(const t_gpar &p) const {
  t_gpar par = (*this);
  par.value(value() + p.value());
  return par;
}

// Setting begin and end time for validity object
// -------------------------------------------------
void t_gpar::setTime(const t_gtime &t1, const t_gtime &t2) {
  this->beg = t1;
  this->end = t2;
}

// get data type
// ----------
string t_gpar::str_type() const {
  string type;
  switch (parType) {
    case par_type::CRD_X    : type = "CRD_X";
      break;
    case par_type::CRD_Y    : type = "CRD_Y";
      break;
    case par_type::CRD_Z    : type = "CRD_Z";
      break;
    case par_type::TRP      : type = "TRP";
      break;
    case par_type::VEL_X: type = "VEL_X";
        break;
    case par_type::VEL_Y: type = "VEL_Y";
        break;
    case par_type::VEL_Z: type = "VEL_Z";
        break;
    case par_type::CLK_RAT: type = "CLK_RAT";
      break;
    case par_type::SION     : type = "SION_" + prn;
      break;
    case par_type::VION     : type = "VION_" + prn;
      break;
    case par_type::CLK      : type = "CLK";
      break;
    case par_type::CLK_G: type = "CLK_G";
        break;
    case par_type::CLK_C: type = "CLK_C";
        break;
    case par_type::CLK_E: type = "CLK_E";
        break;
    case par_type::CLK_R: type = "CLK_R";
        break;
    case par_type::CLK_J: type = "CLK_J";
        break;
    case par_type::CLK_SAT  : type = "CLK_SAT_" + prn;
      break;
      break;
    case par_type::AMB_IF   : type = "AMB_IF_" + prn;
      break;
      break;
    case par_type::AMB_WL   : type = "AMB_WL_" + prn;
      break;
    case par_type::AMB_L1   : type = "AMB_L1_" + prn;
      break;
    case par_type::AMB_L2   : type = "AMB_L2_" + prn;
      break;
    case par_type::AMB_L3: type = "AMB_L3_" + prn;
        break;
    case par_type::AMB_L4: type = "AMB_L4_" + prn;
        break;
    case par_type::AMB_L5: type = "AMB_L5_" + prn;
        break;
    case par_type::CLK_ICB  : type = "CLK_ICB";
      break;
    case par_type::CLUSTERB : type = "CLUSTERB";
      break;
    case par_type::GRD_N    : type = "GRD_N";
      break;
    case par_type::GRD_E    : type = "GRD_E";
      break;
    case par_type::GLO_ISB  : type = "GLO_ISB";
      break;
    case par_type::GLO_IFB : type = "GLO_IFB_" + int2str(channel, 2);
      break;
    case par_type::GLO_IFCB: type = "GLO_IFCB";
        break;
    case par_type::GLO_IFPB: type = "GLO_IFPB";
        break;
    case par_type::GAL_ISB  : type = "GAL_ISB";
      break;
    case par_type::BDS_ISB  : type = "BDS_ISB";
      break;
    case par_type::BD2_ISB : type = "BD2_ISB";
      break;
    case par_type::QZS_ISB  : type = "QZS_ISB";
      break;
    case par_type::IFB_GPS  : type = "IFB_GPS";
      break;
    case par_type::IFB_BDS: type = "IFB_BDS";
      break;
    case par_type::IFB_GAL: type = "IFB_GAL";
      break;
    case par_type::IFB_QZS: type = "IFB_QZS";
      break;
    case par_type::GPS_REC_IFB_C3: type = "GPS_REC_IFB_C3";
        break;
    case par_type::P1P2G_REC : type = "P1P2G_REC";
      break;
    case par_type::P1P2E_REC : type = "P1P2E_REC";
      break;
    case par_type::PXSAT: type = "PXSAT_" + prn;
        break;
    case par_type::PYSAT: type = "PYSAT_" + prn;
        break;
    case par_type::PZSAT: type = "PZSAT_" + prn;
        break;
    default        : type = "UNDEF";
  }

  return type;
}

t_gparhead t_gpar::get_head() const {
  return t_gparhead(this->parType, this->site, this->prn, this->channel);
}

t_gtimearc t_gpar::get_timearc() const {
  return t_gtimearc(this->beg, this->end);
}

// get ZTD mf according to settings   
void t_gpar::_getmf(t_gsatdata &satData,
                    const t_gtriple &crd,
                    const t_gtime &epoch,
                    double &mfw,
                    double &mfh,
                    double &dmfw,
                    double &dmfh) {
  if (parType != par_type::TRP && parType != par_type::GRD_N && parType != par_type::GRD_E) return;

  double ele = satData.ele();

  if (_mf_ztd == ZTDMPFUNC::COSZ) {
    mfw = mfh = 1.0 / sin(ele);
  } else if (_mf_ztd == ZTDMPFUNC::GMF) {
    t_gmf mf;
    mf.gmf(epoch.mjd(), crd[0], crd[1], crd[2], G_PI / 2.0 - ele,
           mfh, mfw, dmfh, dmfw);
  } else cerr << "ZTD mapping function is not set up correctly!!!" << endl;
}

string gpar2str(const t_gpar &par) {
  string partype = par.str_type();
  partype = par.site + t_gpar::PAR_STR_SEP + partype;
  if (par.prn == "") {
    partype += t_gpar::PAR_STR_SEP;
  }
  return partype;
}

string LibGnut_LIBRARY_EXPORT ptype2str(const par_type &parType) {
  // get data type
  // ----------
  string type;
  switch (parType) {
    case par_type::CRD_X: type = "CRD_X";
      break;
    case par_type::CRD_Y: type = "CRD_Y";
      break;
    case par_type::CRD_Z: type = "CRD_Z";
      break;
    case par_type::TRP: type = "TRP";
      break;
    case par_type::SION: type = "SION_";
      break;
    case par_type::VION: type = "VION_";
      break;
    case par_type::CLK: type = "CLK";
      break;
    case par_type::CLK_G: type = "CLK_G";
        break;
    case par_type::CLK_E: type = "CLK_E";
        break;
    case par_type::CLK_C: type = "CLK_C";
        break;
    case par_type::CLK_R: type = "CLK_R";
        break;
    case par_type::CLK_J: type = "CLK_J";
        break;
    case par_type::CLK_SAT: type = "CLK_SAT_";
      break;
    case par_type::AMB_IF: type = "AMB_IF_";
      break;
    case par_type::AMB_WL: type = "AMB_WL_";
      break;
    case par_type::AMB_L1: type = "AMB_L1_";
      break;
    case par_type::AMB_L2: type = "AMB_L2_";
      break;
    case par_type::AMB_L3: type = "AMB_L3_";
        break;
    case par_type::AMB_L4: type = "AMB_L4_";
        break;
    case par_type::AMB_L5: type = "AMB_L5_";
        break;
    case par_type::CLK_ICB: type = "CLK_ICB";
      break;
    case par_type::CLUSTERB: type = "CLUSTERB";
      break;
    case par_type::GRD_N: type = "GRD_N";
      break;
    case par_type::GRD_E: type = "GRD_E";
      break;
    case par_type::GLO_ISB: type = "GLO_ISB";
      break;
    case par_type::GLO_IFB: type = "GLO_IFB";
      break;
    case par_type::GLO_IFCB: type = "GLO_IFCB";
        break;
    case par_type::GLO_IFPB: type = "GLO_IFPB";
        break;
    case par_type::GAL_ISB: type = "GAL_ISB";
      break;
    case par_type::BDS_ISB: type = "BDS_ISB";
      break;
    case par_type::BD2_ISB: type = "BD2_ISB";
        break;
    case par_type::QZS_ISB: type = "QZS_ISB";
      break;
    case par_type::IFB_GPS: type = "IFB_GPS";
      break;
    case par_type::IFB_BDS: type = "IFB_BDS";
      break;
    case par_type::IFB_GAL: type = "IFB_GAL";
      break;
    case par_type::IFB_QZS: type = "IFB_QZS";
      break;
    case par_type::GPS_REC_IFB_C3: type = "IFB_REC_F3";
        break;
    case par_type::PXSAT: type = "PXSAT";
        break;
    case par_type::PYSAT: type = "PYSAT";
        break;
    case par_type::PZSAT: type = "PZSAT";
        break;
    case par_type::P1P2G_REC: type = "P1P2G_REC";
      break;
    case par_type::P1P2E_REC: type = "P1P2E_REC";
      break;
    default: type = "UNDEF";
  }
  return type;
}

t_gpar str2gpar(const string &str_par) {
  int sep_first, sep_last = 0;
  string site, sat, str_type;
  sep_first = str_par.find(t_gpar::PAR_STR_SEP);
  site = str_par.substr(0, sep_first);

  sep_last = str_par.rfind(t_gpar::PAR_STR_SEP);
  sat = str_par.substr(sep_last + 1);
  int channel = DEF_CHANNEL;
  if (sat.size() == 2 && (sat.substr(0, 1) == "-" || sat.substr(0, 1) == "0")) {
      channel = str2int(sat);
      sat = "";
  }

  str_type = str_par.substr(sep_first + 1, sep_last - sep_first - 1);

  par_type par_tp(par_type::NO_DEF);

  for (int i = 0; i < t_gpar::PAR_STR.size(); i++) {
    if (t_gpar::PAR_STR[i] == str_type) {
        par_tp = par_type(i);
      break;
    }
  }

  t_gpar par(site, par_tp, 0, sat);
  par.channel = channel;
  return par;
}

} // namespace
