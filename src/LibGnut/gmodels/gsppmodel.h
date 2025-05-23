
/**
*
* @verbatim
    History
    2014-11-27  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gsppmodel.h
* @brief      Purpose: various SPP models
*.
* @author     PV
* @version    1.0.0
* @date       2014-11-27
*
*/

#ifndef GSPPMODEL_H
#define GSPPMODEL_H

#include <string>
#include <map>
#include <cmath>

#include "gmodels/gmodel.h"

#include "newmat/newmat.h"
#include "newmat/newmatio.h"

#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gutils/gsysconv.h"
#include "gutils/gnss.h"
#include "gmodels/ggmf.h"
#include "gutils/gnss.h"

using namespace std;

namespace gnut
{
	/** @brief class for t_gsppmodel derive from t_godel. */
	class LibGnut_LIBRARY_EXPORT t_gsppmodel : public t_gmodel
	{
	public:
		/** @brief constructor 1. */
		t_gsppmodel(string site, t_glog* glog, t_gsetbase* settings);
		t_gsppmodel();
		virtual ~t_gsppmodel();
		/** @brief Outliers detection. */
		virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&);
		virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, vector<t_gsatdata>::iterator &);
		virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&);
		virtual int outlierDetect_chi(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&);
		/** @brief model for computed range value. */
		virtual double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata&, t_gobs& gobs, bool com = false);
		/** @brief model computed D range value (phase/code). */
		virtual double cmpObsD(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs);
		/** @brief get trop Delay. */
		double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata);
		/** @brief get iono Delay. */
		virtual double ionoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata, t_gobs& gobs);
		/** @brief is Correction. */
		double isbCorrection(t_gallpar& param, string& sat, string&rec, t_gobs&gobs);
		/** @brief is reset observation */
		virtual void reset_observ(OBSCOMBIN observ) override;
		/** @brief is reset tropmf */
		virtual void reset_tropmf(ZTDMPFUNC mf) { _tropo_mf = mf; }
		/** @brief is set rec */
		virtual void setrec(shared_ptr<t_gobj> rec);
		/** @brief get outlier satellite */
		vector<string> get_outlier_sat() { return _outlier_sat; }

	protected:
		/** @brief Find maximal residual */
		double               _maxres(const ColumnVector& v, GSYS gs, bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA);  // OLD
		double               _maxres(bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA, RESIDTYPE res_type, GSYS gs = GNS);
		/** @brief check maximal residual */
		bool                 _check_outl(bool phase, double& maxres, vector<t_gsatdata>::iterator& itData, vector<t_gsatdata>& data); // OLD
		bool                 _check_outl(bool phase, double& maxresNORM, vector<t_gsatdata>::iterator& itDataNORM,
		double& maxresORIG, vector<t_gsatdata>::iterator& itDataORIG,
		vector<t_gsatdata>::iterator& itDataErase, vector<t_gsatdata>& data);
		/** @brief logging outlier */
		void                 _logOutl(bool phase, string prn, int data_size, double maxres, double ele, t_gtime epo, RESIDTYPE resid_type);
		vector<ColumnVector> _devide_res(const ColumnVector& v_orig);

		shared_ptr<fstream> debug_output;
		map<GSYS, double>   _maxres_C;                  ///< code maximal residual
		map<GSYS, double>   _maxres_L;                  ///< phase maximal residual
		double              _maxres_norm;               ///< normal maximal residual
		shared_ptr<t_gobj>  _grec;                      ///< grec
		t_gtriple           _antennal_height;           ///< antennal height
		TROPMODEL           _trpModStr;                 ///< trop mod(str)
		ZTDMPFUNC           _tropo_mf;                  ///< trop mf
		RESIDTYPE           _resid_type;                ///< residual type
		OBSCOMBIN           _observ;                    ///< observation
		CBIASCHAR           _cbiaschar;                 ///< bias(char)
		map<GSYS, map<FREQ_SEQ, GOBSBAND>> _band_index; ///< band
		map<GSYS, map<GOBSBAND, FREQ_SEQ>> _freq_index; ///< freq
		vector<string> _outlier_sat;                    ///< outlier satellite

		string strEst;
		bool emptyBase = false;
	};

} // namespace

#endif //  GSPPMODEL_H
