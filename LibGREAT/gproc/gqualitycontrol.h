/**
 * @file         gqualitycontrol.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        main about quality control
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef QUALITYCONTROL_H
#define QUALITYCONTROL_H

#include "gexport/ExportLibGREAT.h"
#include "gall/gallnav.h"
#include "gall/gallobs.h"

#include "gdata/gsatdata.h"

#include "gset/gsetgen.h"
#include "gset/gsetrec.h"
#include "gset/gsetgnss.h"
#include "gset/gsetproc.h"
#include "gset/gsetturboedit.h"

using namespace gnut;

namespace great
{
	enum SMOOTH_MODEL
	{
		SMT_DOPPLER,
		SMT_PHASE,
		SMT_NONE
	};

	class LibGREAT_LIBRARY_EXPORT t_gsmooth
	{
	public:
		/**
		 * @brief Construct a new t gsmooth object
		 * @param[in]  settings  setbase control
		 */
		t_gsmooth() {};
		t_gsmooth(t_gsetbase* settings);
		/**
		 * @brief Destroy the t gsmooth object
		 */
		virtual ~t_gsmooth() {};

		void smooth_range_obs(vector<t_gsatdata>& obsdata, const t_gtime& now);

	private:

		void _doppler_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now);
		void _phase_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now);

		SMOOTH_MODEL                    _smoothModel;
		int                             _smoothWindow;
		double                          _smoothFactor;
		double                          _sampling;
		// site/sat/range_obs/pre_smt_value
		map<string, map<string, map<GOBS, double> > >                   _pre_smt_range;
		map<string, map<string, t_gtime > >                             _smt_beg_time;
		map<string, map<string, map<GOBSBAND, pair<GOBS, double> > > >  _pre_orig_val;

	};

	/**
	 * @brief class for BeiDou satellite-induced code pseudorange variations correct
	 */
	class LibGREAT_LIBRARY_EXPORT t_gbds_codebias_cor
	{
	public:
		t_gbds_codebias_cor(t_gsetbase* settings);
		virtual ~t_gbds_codebias_cor() {};

		void apply_IGSO_MEO(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata);
		void apply_IGSO_MEO_obs(t_gtriple& rec_crd, t_gallnav* gnav, const shared_ptr<t_gobsgnss>& obsdata);
		void apply_IGSO_MEO_ele(const double& elev, t_gsatdata& obsdata);

	private:

		t_gsetbase*                          _set;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		bool                                 _correct_bds_code_bias;
		// Wanninger & Beer : BeiDou satellite-induced code pseudorange variations: diagnosis and therapy [unit:m]
		map<GOBSBAND, map<string, map<int, double> > > _IGSO_MEO_Corr;
		/**
		 * @brief Approximate location of the receiver
		 */
		bool _recAprCoordinate(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata);

	};

	/**
	 * @brief Class for outliers process
	 */
	class LibGREAT_LIBRARY_EXPORT t_goutliers_process
	{
	public:
		t_goutliers_process(t_gsetbase* settings, t_glog* glog);

		virtual ~t_goutliers_process();

		void setLog(t_glog* glog) { _log = glog; };

		void excludeBadObs(vector<t_gsatdata>& obsdata); // Low SNR < 10dB

		void flagRangeOutliers(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs, double sampling);

	private:
		t_glog*            _log;
		t_gsetbase*        _set;                    ///< setbase control
		t_giof*            _debug_outliers=nullptr; ///< debug outliers

	};
	/**
	 * @brief Class for qualitycontrol
	 */
	class LibGREAT_LIBRARY_EXPORT t_gqualitycontrol
	{
	public:
		t_gqualitycontrol(t_gsetbase* settings, t_gallnav* gnav);
		virtual ~t_gqualitycontrol();

		int processOneEpoch(const t_gtime& now, const string& rec, t_gtriple& rec_crd, vector<t_gsatdata>& obsdata);
		void bds_code_bias(t_gtriple& rec_crd, const shared_ptr<t_gobsgnss>& obsdata);
		void setNav(t_gallnav* gnav) { _gnav = gnav; };

	protected:

		t_glog*            _log;  ///< logbase control
		t_gsetbase*        _set;  ///< setbase control
		t_gallnav*         _gnav; ///< navigation data

		t_gbds_codebias_cor   _bds_codebias_cor; ///< bds codebias correction
		t_gsmooth             _smooth_range;     ///< smooth range
		t_goutliers_process   _outliers_proc;    ///< outliers process


	};


} // namespace

#endif
