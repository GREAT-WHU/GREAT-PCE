
/**
*
* @verbatim
    History
    2011-04-18  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gpar.h
* @brief      Purpose: parametres class
*.
* @author     PV
* @version    1.0.0
* @date       2011-04-18
*
*/
#ifndef GPAR_H
#define GPAR_H

#include <string>

#include "gexport/ExportLibGnut.h"
#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gset/gsetproc.h"

using namespace std;

namespace gnut {

    enum class LibGnut_LIBRARY_EXPORT par_type {
        CRD_X,
        CRD_Y,
        CRD_Z,  // coordinates
        TRP,
        GRD_N,
        GRD_E,
        SION,
        VION,  // atmospheric parameters
        CLK,
        CLK_SAT,// clocks
        CLK_G,
        CLK_E,
        CLK_C,
        CLK_R,
        CLK_J,
        CLK_ICB,
        CLUSTERB,
        IFCB_F3,
        IFCB_F4,
        IFCB_F5,
        AMB_IF,
        AMB_WL,
        AMB_L1,
        AMB_L2,
        AMB_L3,
        AMB_L4,
        AMB_L5, // ambiguities for indiv. freq. (number indicates freq not band)
        GLO_ISB,
        GLO_IFB,
        GLO_IFCB,
        GLO_IFPB,
        GAL_ISB,
        BDS_ISB,
        BD2_ISB,
        QZS_ISB,
        GPS_REC_IFB_C3,
        IFB_GPS,
        IFB_BDS,
        IFB_QZS,
        IFB_GAL,
        PXSAT,
        PYSAT,
        PZSAT,
        VEL_X,
        VEL_Y,
        VEL_Z, // velocity
        VTEC,
        DCB_REC,
        DCB_SAT,
        P1P2G_REC,
        P1P2E_REC,
        P1P2R_REC,
        P1P2C_REC,
        CLK_RAT,
        SION_3,
        NO_DEF
    };

class LibGnut_LIBRARY_EXPORT t_gparhead;
/** @brief class for t_gtimearc. */
class LibGnut_LIBRARY_EXPORT t_gtimearc {
 public:

  t_gtimearc(const t_gtime &beg, const t_gtime &end);
  ~t_gtimearc();
  /** @brief override operator. */
  bool operator!=(const t_gtimearc &Other) const;
  bool operator==(const t_gtimearc &Other) const;
  bool operator<(const t_gtimearc &Other) const;
  bool operator<=(const t_gtimearc &Other) const;
  bool operator>(const t_gtimearc &Other) const;
  bool operator>=(const t_gtimearc &Other) const;
  /** @brief inside. */
  bool inside(const t_gtimearc &Other) const;

  t_gtime begin;
  t_gtime end;
};
/** @brief class for t_gpar. */
class LibGnut_LIBRARY_EXPORT t_gpar {
 public:

  const static string PAR_STR_SEP;
  const static vector<string> PAR_STR;

  static bool is_amb(const par_type& tp) {
      return (tp == par_type::AMB_IF || tp == par_type::AMB_L1 || tp == par_type::AMB_L2);
  }

  static bool is_crd(const par_type& tp) {
      return (tp == par_type::CRD_X || tp == par_type::CRD_Y || tp == par_type::CRD_Z);
  }

  static bool  is_sysbias(const par_type& tp) {
      return (tp == par_type::GAL_ISB || tp == par_type::BD2_ISB || tp == par_type::BDS_ISB || tp == par_type::QZS_ISB ||
          tp == par_type::GLO_ISB || tp == par_type::GLO_IFB);
  }

  static string sysbias2cgsys(const par_type& tp, const int& chn = DEF_CHANNEL) {
      if (tp == par_type::BDS_ISB) return "C";
      else if (tp == par_type::BD2_ISB) return "C2";
      else if (tp == par_type::GAL_ISB) return "E";
      else if (tp == par_type::QZS_ISB) return "J";
      else if (tp == par_type::GLO_ISB) return "R";
      else if (tp == par_type::GLO_IFB) return "R" + int2str(chn, 2);
  }

  t_gpar(const string &site, par_type t, unsigned i, const string &p, bool remove = true);
  t_gpar();
  virtual ~t_gpar();
  /** @brief override operator. */
  bool operator==(const t_gpar &) const;
  bool operator<(const t_gpar &) const;
  bool operator>(const t_gpar &) const;
  t_gpar operator-(const t_gpar &) const;
  t_gpar operator+(const t_gpar &) const;

  void setTime(const t_gtime &, const t_gtime &);

  double partial(t_gsatdata &, t_gtime &, t_gtriple, t_gobs &gobs);
  double partial_doppler(t_gsatdata& satData,
      t_gtriple& groundXYZ,
      t_gtriple& groundVEL);
  double partial_orb(t_gsatdata &);

  void value(double val) { _value = val; }
  double value() const { return _value; }

  void apriori(double apr) { _apriori = apr; }
  double apriori() const { return _apriori; }

  par_type parType;///< par type
  int index;///< index
  string prn;///< satellite name
  string site;///< site name
  t_gtime beg;///< begin time
  t_gtime end;///< end time
  t_gtime stime;///< s time
  double aprval;///< apr value
  double pred;///< pred
  double zhd;///< for ztd par
  double zwd; ///< for zwd par
  int channel = DEF_CHANNEL; ///< channel
  bool amb_ini = false;///<amb to be initialized
  int nPWC;     // number of parameters for this kind of PWC parameter
  int piecewise_time;// time
  bool lremove;
  
  string str_type() const;
  t_gparhead get_head() const;
  t_gtimearc get_timearc() const;

  void setMF(ZTDMPFUNC MF);
  void setMF(GRDMPFUNC MF);

 protected:
  double _value;///< value
  ZTDMPFUNC _mf_ztd;///< mapping function for ZTD
  GRDMPFUNC _mf_grd;///< mapping function for GRD
  double _apriori;///< apriori
  void _getmf(t_gsatdata &satData,
              const t_gtriple &crd,
              const t_gtime &epoch,
              double &mfw,
              double &mfh,
              double &dmfw,
              double &dmfh);
};
/** @brief class for t_gprhead. */
class LibGnut_LIBRARY_EXPORT t_gparhead {
 public:
  t_gparhead(const par_type& type, const string& site, const string& sat, const int& channel = DEF_CHANNEL);
  t_gparhead(const t_gparhead &Other);
  ~t_gparhead();
  /** @brief override operator. */
  bool operator==(const t_gparhead &Other) const;
  bool operator<(const t_gparhead &Other) const;
  bool operator<=(const t_gparhead &Other) const;
  bool operator>(const t_gparhead &Other) const;
  bool operator>=(const t_gparhead &Other) const;

  par_type type;
  string site;
  string sat;
  string str_type;
  int channel = DEF_CHANNEL;
};

class t_gparhead_hash
{
 public:
    size_t operator()(const t_gparhead& a) const;
};
/** @brief convert par to string. */
string LibGnut_LIBRARY_EXPORT gpar2str(const t_gpar &par);
/** @brief convert partype to string. */
string LibGnut_LIBRARY_EXPORT ptype2str(const par_type &par);
/** @brief convert string to par. */
t_gpar LibGnut_LIBRARY_EXPORT str2gpar(const string &str_par);

}

#endif
