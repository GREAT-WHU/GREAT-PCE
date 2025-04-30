/**
 * @file         gqualitycontrol.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        main about quality control
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include <iomanip>
#include <memory>

#include "gproc/gqualitycontrol.h"
#include "gmodels/gbancroft.h"

namespace great {

	t_gsmooth::t_gsmooth(t_gsetbase* settings)
	{
		string smoothModel = dynamic_cast<t_gsetproc*>(settings)->range_smooth_mode(&_smoothWindow);
		if      (smoothModel == "DOPPLER") _smoothModel = SMOOTH_MODEL::SMT_DOPPLER;
		else if (smoothModel == "PHASE")   _smoothModel = SMOOTH_MODEL::SMT_PHASE;
		else                               _smoothModel = SMOOTH_MODEL::SMT_NONE;

		_smoothFactor = 1.0 / _smoothWindow;
		_sampling = dynamic_cast<t_gsetgen*>(settings)->sampling();

	}

	void t_gsmooth::smooth_range_obs(vector<t_gsatdata>& obsdata, const t_gtime& now)
	{
		if (_smoothModel == SMOOTH_MODEL::SMT_NONE) return;

		// doppler
		if (_smoothModel == SMOOTH_MODEL::SMT_DOPPLER)
		{
			_doppler_smt_range(obsdata, now);
		}
		// phase
		else
		{
			_phase_smt_range(obsdata, now);
		}
	}

	void t_gsmooth::_doppler_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now)
	{		
		for (auto &oneObs : obsdata)
		{
			if (oneObs.getrangestate("smooth_range")) continue;
			string site = oneObs.site();
			GSYS   gsys = oneObs.gsys();
			string gsat = oneObs.sat();

			if (_smt_beg_time.find(site) == _smt_beg_time.end() || _smt_beg_time[site].find(gsat) == _smt_beg_time[site].end())
			{
				_smt_beg_time[site][gsat] = now;
			}
			else if (now.diff(_smt_beg_time[site][gsat]) > _smoothWindow * _sampling)
			{
				_smt_beg_time[site][gsat] = now;
			}

			vector<GOBS> vec_obs = oneObs.obs();
			for (GOBS obs_type : vec_obs)
			{
				// skip not code obs
				if (!t_gobs(obs_type).is_code())
				{
					continue;
				}

				t_gobs gobs(obs_type);
				GOBSBAND band = gobs.band();

				GOBS GOBSP = obs_type;
				string strGOBS = gobs2str(GOBSP);
				strGOBS.replace(0, 1, "D");
				GOBS GOBSD = str2gobs(strGOBS);

				// get range/doppler obs.
				double obsP = oneObs.obs_C(GOBSP);  // m
				double obsD = oneObs.obs_D(GOBSD);  // m/s

				if (_smt_beg_time[site][gsat] == now || GOBSD != _pre_orig_val[site][gsat][band].first
					|| oneObs.getoutliers(GOBSD) > 0)
				{
					_pre_orig_val[site][gsat][band] = make_pair(GOBSD, obsD);
					_pre_smt_range[site][gsat][GOBSP] = obsP;
					_smt_beg_time[site][gsat] == now;
				}
				else
				{
					double smt = _smoothFactor * obsP + (1.0 - _smoothFactor) * (_pre_smt_range[site][gsat][GOBSP] - obsD);
					_pre_orig_val[site][gsat][band].first = GOBSD;
					_pre_orig_val[site][gsat][band].second = obsD;
					_pre_smt_range[site][gsat][GOBSP] = smt;
					oneObs.resetobs(GOBSP, smt);
					oneObs.setrangestate("smooth_range", true);
				}
			}
		}

	}
	
	void t_gsmooth::_phase_smt_range(vector<t_gsatdata>& obsdata, const t_gtime& now)
	{
		for (auto &oneObs : obsdata)
		{
			if (oneObs.getrangestate("smooth_range")) continue;
			string site = oneObs.site();
			GSYS   gsys = oneObs.gsys();
			string gsat = oneObs.sat();

			if (_smt_beg_time.find(site) == _smt_beg_time.end() || _smt_beg_time[site].find(gsat) == _smt_beg_time[site].end())
			{
				_smt_beg_time[site][gsat] = now;
			}
			else if (now.diff(_smt_beg_time[site][gsat]) > _smoothWindow * _sampling)
			{
				_smt_beg_time[site][gsat] = now;
			}

			vector<GOBS> vec_obs = oneObs.obs();
			for (GOBS obs_type : vec_obs)
			{
				// skip not code obs
				if (!t_gobs(obs_type).is_code())
				{
					continue;
				}

				t_gobs gobs(obs_type);
				GOBSBAND band = gobs.band();

				GOBS GOBSP = obs_type;
				GOBS GOBSL = oneObs.select_phase(band, true);

				double obsP = oneObs.obs_C(GOBSP);
				double obsL = oneObs.obs_L(GOBSL);

				if (_smt_beg_time[site][gsat] == now || GOBSL != _pre_orig_val[site][gsat][band].first 
					|| oneObs.getlli(GOBSL) > 0)
				{
					_pre_orig_val[site][gsat][band] = make_pair(GOBSL, obsL);
					_pre_smt_range[site][gsat][GOBSP] = obsP;
					_smt_beg_time[site][gsat] == now;
				}
				else
				{
					double smt = _smoothFactor * obsP + (1.0 - _smoothFactor) * 
						(_pre_smt_range[site][gsat][GOBSP] + obsL - _pre_orig_val[site][gsat][band].second);
					_pre_orig_val[site][gsat][band].first = GOBSL;
					_pre_orig_val[site][gsat][band].second = obsL;
					_pre_smt_range[site][gsat][GOBSP] = smt;
					oneObs.resetobs(GOBSP, smt);
					oneObs.setrangestate("smooth_range", true);
				}
			}
		}
	}

	t_gbds_codebias_cor::t_gbds_codebias_cor(t_gsetbase* settings)
	{
		_correct_bds_code_bias = dynamic_cast<t_gsetproc*>(settings)->bds_code_bias_correction();
		set<string> sys_list = dynamic_cast<t_gsetgen*>(settings)->sys();
		if (sys_list.find("BDS") == sys_list.end())
		{
			_correct_bds_code_bias = false;
		}

		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(settings)->band_index(gnut::BDS);
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(settings)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(settings)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(settings)->band_index(gnut::GLO);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(settings)->band_index(gnut::QZS);
		_set = settings;

		// set Wanninger & Beer 
		_IGSO_MEO_Corr[BAND_2]["IGSO"] = { {0,-0.55}, {1,-0.40}, {2,-0.34}, {3,-0.23}, {4,-0.15}, {5,-0.04}, {6,0.09}, {7,0.19}, {8,0.27}, {9,0.35} }; // B1
		_IGSO_MEO_Corr[BAND_7]["IGSO"] = { {0,-0.71}, {1,-0.36}, {2,-0.33}, {3,-0.19}, {4,-0.14}, {5,-0.03}, {6,0.08}, {7,0.17}, {8,0.24}, {9,0.33} }; // B2
		_IGSO_MEO_Corr[BAND_6]["IGSO"] = { {0,-0.27}, {1,-0.23}, {2,-0.21}, {3,-0.15}, {4,-0.11}, {5,-0.04}, {6,0.05}, {7,0.14}, {8,0.19}, {9,0.32} }; // B3

		_IGSO_MEO_Corr[BAND_2]["MEO"]  = { {0,-0.47}, {1,-0.38}, {2,-0.32}, {3,-0.23}, {4,-0.11}, {5,0.06}, {6,0.34}, {7,0.69}, {8,0.97}, {9,1.05} }; //B1
		_IGSO_MEO_Corr[BAND_7]["MEO"]  = { {0,-0.40}, {1,-0.31}, {2,-0.26}, {3,-0.18}, {4,-0.06}, {5,0.09}, {6,0.28}, {7,0.48}, {8,0.64}, {9,0.69} }; //B2
		_IGSO_MEO_Corr[BAND_6]["MEO"]  = { {0,-0.22}, {1,-0.15}, {2,-0.13}, {3,-0.10}, {4,-0.04}, {5,0.05}, {6,0.14}, {7,0.27}, {8,0.36}, {9,0.47} }; //B3
	}

	void t_gbds_codebias_cor::apply_IGSO_MEO(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata)
	{
		if (!this->_correct_bds_code_bias || obsdata.size() == 0 || !gnav) return;

		if (rec_crd.zero())
		{
			if (!this->_recAprCoordinate(rec, rec_crd, gnav, obsdata)) return;
		}

		// Calculate elevation and Correct Obs
		map<GOBSBAND, double> Band_cor;
		string sat_type;

		vector<t_gsatdata>::iterator data;
		for (data = obsdata.begin(); data != obsdata.end();)
		{
			if (data->getrangestate("bds_code_bias"))
			{ 
				data++; 
				continue;
			}
			if (data->site() != rec || data->gsys() != BDS)
			{
				data++;
				continue;
			} 
			string sat = data->sat();
			if (sat <= "C05" || sat > "C16") // BDS2 - IGSO/MEO
			{
				data++;
				continue;
			}
			if (sat == "C11" || sat == "C12" || sat == "C14") sat_type = "MEO";
			else                                              sat_type = "IGSO";


			if (data->satcrd().zero())
			{
				if (data->addprd(gnav) < 0) 
				{ 
					data = obsdata.erase(data); 
					continue; 
				}
			}
			t_gtriple sat_crd = data->satcrd();
			t_gtriple rec_sat_vector = sat_crd - rec_crd;
			double distance = rec_sat_vector.norm();
			double elev = (rec_crd[0] * rec_sat_vector[0] + rec_crd[1] * rec_sat_vector[1] + rec_crd[2] * rec_sat_vector[2]) / rec_crd.norm() / distance;
			elev = 90.0 - acos(elev) * 180.0 / G_PI;

			// get correction
			double elev0 = elev / 10.0;
			int elev0_int = floor(elev0);
			if (elev0_int < 0)
			{
				Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][0];
				Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][0];
				Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][0];
			}
			else if (elev0_int >= 9)
			{
				Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][9];
				Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][9];
				Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][9];
			}
			else
			{
				Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
					+ _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int + 1] * (elev0 - elev0_int);
				Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
					+ _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int + 1] * (elev0 - elev0_int);
				Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
					+ _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int + 1] * (elev0 - elev0_int);
			}

			vector<GOBS> obs_vec = data->obs();
			for (auto obs_type : obs_vec)
			{
				// skip not code obs
				if (!t_gobs(obs_type).is_code())
				{
					continue;
				}

				GOBSBAND b = t_gobs(obs_type).band();
				double obs_P = data->obs_C(obs_type);

				obs_P += Band_cor[b];
				data->resetobs(obs_type, obs_P);
				data->setrangestate("bds_code_bias", true);
			}

			data++;
		}
	}

	void t_gbds_codebias_cor::apply_IGSO_MEO_obs(t_gtriple& rec_crd, t_gallnav* gnav, const shared_ptr<t_gobsgnss>& obsdata)
	{
		if (!_correct_bds_code_bias || !obsdata || !gnav || 
			rec_crd.zero() || obsdata->getrangestate("bds_code_bias")) return;
		const string& sat = obsdata->sat();
		if (sat.substr(0, 1) != "C" || sat > "C17" || sat < "C06") return;
		
		const t_gtime& epoch = obsdata->epoch();
		double xyz[3] = { 0.0, 0.0, 0.0 };
		if (gnav->pos(sat, epoch, xyz, nullptr, nullptr, false) < 0) return;
		t_gtriple sat_crd = t_gtriple(xyz);

		t_gtriple rec_sat_vector = sat_crd - rec_crd;
		double distance = rec_sat_vector.norm();
		double elev = (rec_crd[0] * rec_sat_vector[0] + rec_crd[1] * rec_sat_vector[1] + rec_crd[2] * rec_sat_vector[2]) / rec_crd.norm() / distance;
		elev = 90.0 - acos(elev) * 180.0 / G_PI;
		double elev0 = elev / 10.0;
		int elev0_int = floor(elev0);

		// Calculate elevation and Correct Obs
		string sat_type;
		if (sat == "C11" || sat == "C12" || sat == "C14") sat_type = "MEO";
		else sat_type = "IGSO";
		map<GOBSBAND, double> Band_cor;
		if (elev0_int < 0) {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][0];
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][0];
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][0];
		}
		else if (elev0_int >= 9) {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][9];
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][9];
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][9];
		}
		else {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int + 1] * (elev0 - elev0_int);
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int + 1] * (elev0 - elev0_int);
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int + 1] * (elev0 - elev0_int);
		}

		vector<GOBS> obs_vec = obsdata->obs();
		for (auto obs_type : obs_vec)
		{
			if (!t_gobs(obs_type).is_code()) continue;

			GOBSBAND b = t_gobs(obs_type).band();
			double obs_P = obsdata->obs_C(obs_type);

			obs_P += Band_cor[b];
			obsdata->resetobs(obs_type, obs_P);
			obsdata->setrangestate("bds_code_bias", true);
		}
	}

	void t_gbds_codebias_cor::apply_IGSO_MEO_ele(const double& elev, t_gsatdata& obsdata)
	{
		if (!_correct_bds_code_bias || obsdata.obs_empty() || obsdata.getrangestate("bds_code_bias")) return;
		const string& sat = obsdata.sat();
		if (sat.substr(0, 1) != "C" || sat > "C17" || sat < "C06") return;
		double elev0 = elev / 10.0;
		int elev0_int = floor(elev0);

		// Calculate elevation and Correct Obs
		string sat_type;
		if (sat == "C11" || sat == "C12" || sat == "C14") sat_type = "MEO";
		else sat_type = "IGSO";
		map<GOBSBAND, double> Band_cor;
		if (elev0_int < 0) {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][0];
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][0];
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][0];
		}
		else if (elev0_int >= 9) {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][9];
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][9];
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][9];
		}
		else {
			Band_cor[BAND_2] = _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_2][sat_type][elev0_int + 1] * (elev0 - elev0_int);
			Band_cor[BAND_7] = _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_7][sat_type][elev0_int + 1] * (elev0 - elev0_int);
			Band_cor[BAND_6] = _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int] * (1.0 - elev0 + elev0_int)
				+ _IGSO_MEO_Corr[BAND_6][sat_type][elev0_int + 1] * (elev0 - elev0_int);
		}

		vector<GOBS> obs_vec = obsdata.obs();
		for (auto obs_type : obs_vec)
		{
			if (!t_gobs(obs_type).is_code()) continue;

			GOBSBAND b = t_gobs(obs_type).band();
			double obs_P = obsdata.obs_C(obs_type);

			obs_P += Band_cor[b];
			obsdata.resetobs(obs_type, obs_P);
			obsdata.setrangestate("bds_code_bias", true);
		}
	}

	bool t_gbds_codebias_cor::_recAprCoordinate(const string& rec, t_gtriple& rec_crd, t_gallnav* gnav, vector<t_gsatdata>& obsdata)
	{
		Matrix BB;

		BB.ReSize(obsdata.size(), 4);
		BB = 0.0;
		int iobs = 0;

		vector<t_gsatdata>::iterator iter;
		for (iter = obsdata.begin(); iter != obsdata.end();)
		{
			if (iter->site() != rec)
			{
				iter++;
				continue;
			}

			GSYS gs = iter->gsys();

			GOBSBAND b1 = _band_index[gs][FREQ_1];
			GOBSBAND b2 = _band_index[gs][FREQ_2];

			GOBS l1 = iter->select_phase(b1);
			GOBS l2 = iter->select_phase(b2);
			GOBS p1 = iter->select_range(b1);
			GOBS p2 = iter->select_range(b2);

			if (p1 == X || l1 == X || p2 == X || l2 == X)
			{
				iter++;
				continue;
			}

			double P3 = iter->P3(p1, p2);
			double L3 = iter->L3(l1, l2);

			if (double_eq(L3, 0.0) || double_eq(P3, 0.0))
			{
				iter++;
				continue;
			}


			if (iter->addprd(gnav) < 0) 
			{ 
				iter = obsdata.erase(iter);
				continue; 
			}

			iobs++;
			BB(iobs, 1) = iter->satcrd().crd(0);
			BB(iobs, 2) = iter->satcrd().crd(1);
			BB(iobs, 3) = iter->satcrd().crd(2);
			BB(iobs, 4) = P3 + iter->clk();

			iter++;
		}

		if (iobs < 4) return false;

		BB = BB.Rows(1, iobs);    // delete zero rows

		ColumnVector vBanc;

		gbancroft(BB, vBanc);

		rec_crd[0] = vBanc(1);
		rec_crd[1] = vBanc(2);
		rec_crd[2] = vBanc(3);

		return true;
	}



	t_goutliers_process::t_goutliers_process(t_gsetbase* settings, t_glog* glog):
		_log(glog)
	{
		if (1) {
			if (_debug_outliers)
			{
				if (_debug_outliers->is_open()) _debug_outliers->close();
				delete _debug_outliers;
				_debug_outliers = nullptr;
			}
		}
		else {
			_debug_outliers = new t_giof;
			_debug_outliers->tsys(t_gtime::GPS);
			_debug_outliers->mask("debug_outliers.log");
			_debug_outliers->append(false);
		}
	}

	t_goutliers_process::~t_goutliers_process()
	{
		if (_debug_outliers)
		{
			if (_debug_outliers->is_open()) _debug_outliers->close();
			delete _debug_outliers;
			_debug_outliers = nullptr;
		}
	}

	void t_goutliers_process::flagRangeOutliers(shared_ptr<t_gobsgnss> ObsPre, shared_ptr<t_gobsgnss> Obs, double sampling)
	{
		ostringstream os;
		os.str("");
		t_gtime crt_time = Obs->epoch();
		string  crt_site = Obs->site();

		vector<GOBS> obs_vec = Obs->obs();

		// Gao Y et al. Modeling and estimation of C1-P1 bias in GPS receivers[J]. Journal of Geodesy, 2001.
		GSYS gsys = Obs->gsys();
		GOBSBAND b1 = GNSS_BAND_PRIORITY.at(gsys)[1];
		GOBS ref_obsP = Obs->select_range(b1);
		double ref_P = Obs->obs_C(ref_obsP);  // m

		for (auto obs_type : obs_vec)
		{
			// skip not code obs
			if (!t_gobs(obs_type).is_code()) continue;
			// skip ref obs(add case: ref_obsP == X)
			if (obs_type == ref_obsP || ref_obsP == X) continue;

			t_gobs gobs(obs_type);
			GOBSBAND bX = gobs.band();

			double obs_P = Obs->obs_C(obs_type);
			double diff_P = abs(obs_P - ref_P);

			double value = (b1 == bX) ? 10 : 30;  // meters

			os << fixed << setw(10) << " PP[m] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
				"  " << gobs2str(obs_type) << "  " << gobs2str(ref_obsP) << setw(15) << setprecision(4) << diff_P << endl;

			if (diff_P > value)
			{
				Obs->addoutliers(obs_type, 1);
				Obs->addoutliers(ref_obsP, 1);

				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(ref_obsP)
						<< "   " << gobs2str(obs_type) << "  value  " << setw(16) << setprecision(3) << diff_P << "  >  " << setw(16) << setprecision(3) << value;
					_log->comment(1, crt_site, msg.str());
				}
			}
		}

		if (crt_time.diff(ObsPre->epoch()) > sampling) return;

		for (auto obs_doppler : obs_vec)
		{
			// skip not doppler obs
			if (!t_gobs(obs_doppler).is_doppler()) continue;

			double crt_doppler = Obs->getobs(obs_doppler);      // cycle / s
			double pre_doppler = ObsPre->getobs(obs_doppler);   // cycle / s

			// check Obstype Gap
			if (double_eq(pre_doppler, NULL_GOBS) || double_eq(crt_doppler, NULL_GOBS)) continue;

			// cycle/s
			double delta_doppler = crt_doppler - pre_doppler;

			// Codeless
			if (obs_doppler == D1N)
			{
				if (abs(delta_doppler) > 10)  Obs->addoutliers(obs_doppler, 1);
				continue;
			}

			string strGOBS = gobs2str(obs_doppler);
			strGOBS.replace(0, 1, "C");
			GOBS obs_range = str2gobs(strGOBS);

			if (Obs->getoutliers(obs_range) > 0)  continue;

			GOBSBAND b = t_gobs(obs_doppler).band();
			crt_doppler *= Obs->wavelength(b);     // m/s
			pre_doppler *= ObsPre->wavelength(b);  // m/s

			double crt_range = Obs->obs_C(obs_range);    // m
			double pre_range = ObsPre->obs_C(obs_range); // m

			if (double_eq(crt_range, NULL_GOBS) || double_eq(pre_range, NULL_GOBS)) continue;

			double delta_PD = crt_range - pre_range + 0.5 * (crt_doppler + pre_doppler) * sampling;  // m

			os << fixed << setw(10) << " D[c/s] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
			    "  " << gobs2str(obs_doppler) << setw(15) << setprecision(4) << delta_doppler << endl;

			os << fixed << setw(10) << " PD[m] " << Obs->sat() << "  " << setw(6) << ObsPre->epoch().sod() << setw(6) << crt_time.sod() <<
				"  " << gobs2str(obs_range) << "  " << gobs2str(obs_doppler) << setw(15) << setprecision(4) << delta_PD << endl;

			if (abs(delta_doppler) > 10 && abs(delta_PD) > 2)     // cycle/s, m
			{
				Obs->addoutliers(obs_doppler, 1);       // doppler
				Obs->addoutliers(obs_range, 1);         // range
				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Doppler and Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_doppler)
						<< "   " << gobs2str(obs_range) << "  delta_doppler[cycle/s] =  " << setw(16) << setprecision(3) << delta_doppler
						<< "  delta_PD[m] =  " << setw(16) << setprecision(3) << delta_PD;
					_log->comment(1, crt_site, msg.str());
				}
			}
			else if (abs(delta_PD) > 5)  // m
			{
				Obs->addoutliers(obs_range, 1);         // range
				if (_log)
				{
					ostringstream msg;
					msg << fixed << setw(10) << " Outliers in Range : " << crt_time.str_hms() << "   " << Obs->sat() << "   " << gobs2str(obs_range)
						<< "  delta_PD[m] =  " << setw(16) << setprecision(3) << delta_PD;
					_log->comment(1, crt_site, msg.str());
				}
			}
		}

		// Print flt results
		if (_debug_outliers) {
			_debug_outliers->write(os.str().c_str(), os.str().size());
			_debug_outliers->flush();
		}
	}

	void t_goutliers_process::excludeBadObs(vector<t_gsatdata>& obsdata)
	{
		vector<t_gsatdata>::iterator it;
		bool valid_SNR, valid_Code, valid_Phase;
		bool lack_Phase, single_freq;
	
		for (it = obsdata.begin(); it != obsdata.end();)
		{
			 valid_SNR = true;
			 valid_Code = true;
			 valid_Phase = true;
			 lack_Phase = true;
			 single_freq = false;
			GSYS gsys = it->gsys();

			vector<GOBS> obs_vec = it->obs();
			for (auto obs_type : obs_vec)
			{   
				GOBSTYPE type = t_gobs(obs_type).type();
				
				// check_SNR
				if ( type == TYPE_S)
				{
					double obs_snr = it->getobs(obs_type);
					if (obs_snr <= 10.0) valid_SNR = false;
					else valid_SNR = true;
				}
				// check Range
				else if (type == TYPE_C || type == TYPE_P)
				{
					double obs_code = it->getobs(obs_type);
					if (obs_code <= 1.9e7 ) valid_Code = false;
					else valid_Code = true;
				}
				else if (type == TYPE_L)
				{
					double obs_phase = it->obs_L(obs_type);
					if (obs_phase <= 1.9e7) valid_Phase = false;
					else if (lack_Phase)    lack_Phase = false;
					else  valid_Phase = true; 
				}		
			}

			if (!valid_SNR || !valid_Phase || !valid_Code || lack_Phase)
			{
				it = obsdata.erase(it); 
				continue;
			}

			it++;
		}
	}

    // Constructor
    // ----------
	t_gqualitycontrol::t_gqualitycontrol(t_gsetbase* settings, t_gallnav* gnav):
		_gnav(gnav),
		_bds_codebias_cor(settings),
		_smooth_range(settings),
		_outliers_proc(settings, nullptr)
    {
    }
    
    // Destructor
    // ----------
	t_gqualitycontrol::~t_gqualitycontrol()
    { 
    }

	int t_gqualitycontrol::processOneEpoch(const t_gtime& now, const string& rec, t_gtriple& rec_crd, vector<t_gsatdata>& obsdata)
	{
		this->_bds_codebias_cor.apply_IGSO_MEO(rec, rec_crd, _gnav, obsdata);

		this->_smooth_range.smooth_range_obs(obsdata, now);

		return 1;
	}

	void t_gqualitycontrol::bds_code_bias(t_gtriple& rec_crd, const shared_ptr<t_gobsgnss>& obsdata)
	{
		_bds_codebias_cor.apply_IGSO_MEO_obs(rec_crd, _gnav, obsdata);
	}

} // namespace
