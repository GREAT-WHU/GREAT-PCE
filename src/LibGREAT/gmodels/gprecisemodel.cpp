/**
 * @file         gprecisemodel.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        mainly about how to cacultae B P l single
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/gprecisemodel.h"

#include "gutils/ginfolog.h"
#include "gutils/gstring.h"
#include <algorithm>
#include <assert.h>

namespace great
{
	t_gprecisemodel::t_gprecisemodel(t_gallproc* data, t_glog* log, t_gsetbase* setting) :
		t_gpppmodel("", log, setting),
		_minElev(7),
		_gdata_erp(0),
		_rec_obj_flag(make_pair("", t_gtime()))
	{
		// Get Settings
		_tide        = shared_ptr<t_gtide>(new t_gtideIERS(log));
		_minElev     = dynamic_cast<t_gsetproc*>(setting)->minimum_elev();
		_weight      = dynamic_cast<t_gsetproc*>(setting)->weighting();

		_mf_ztd      = dynamic_cast<t_gsetproc*>(setting)->tropo_mf();
		_mf_grd      = dynamic_cast<t_gsetproc*>(setting)->grad_mf();
		_crd_est = dynamic_cast<t_gsetproc*>(setting)->crd_est();
		_ion_model = dynamic_cast<t_gsetproc*>(setting)->ion_model();

		log->logInfo("t_gprecisemodel", "min Elev",     format("%16.8f", _minElev));
		log->logInfo("t_gprecisemodel", "mf ztd", dynamic_cast<t_gsetproc*>(setting)->ztdmpfunc2str(_mf_ztd));
		log->logInfo("t_gprecisemodel", "mf grd", dynamic_cast<t_gsetproc*>(setting)->grdmpfunc2str(_mf_grd));

		_sigCodeGPS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GPS);
		_sigCodeGLO = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GLO);
		_sigCodeGAL = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(GAL);
		_sigCodeBDS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(BDS);
		_sigCodeQZS = dynamic_cast<t_gsetgnss*>(setting)->sigma_C(QZS);
		_sigPhaseGPS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GPS);
		_sigPhaseGLO = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GLO);
		_sigPhaseGAL = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(GAL);
		_sigPhaseBDS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(BDS);
		_sigPhaseQZS = dynamic_cast<t_gsetgnss*>(setting)->sigma_L(QZS);
		
		log->logInfo("t_gprecisemodel", "sigCodeGPS ", format("%16.8f", _sigCodeGPS));
		log->logInfo("t_gprecisemodel", "sigCodeGLO ", format("%16.8f", _sigCodeGLO));
		log->logInfo("t_gprecisemodel", "sigCodeGAL ", format("%16.8f", _sigCodeGAL));
		log->logInfo("t_gprecisemodel", "sigCodeBDS ", format("%16.8f", _sigCodeBDS));
		log->logInfo("t_gprecisemodel", "sigCodeQZS ", format("%16.8f", _sigCodeQZS));
		log->logInfo("t_gprecisemodel", "sigPhaseGPS", format("%16.8f", _sigPhaseGPS));
		log->logInfo("t_gprecisemodel", "sigPhaseGLO", format("%16.8f", _sigPhaseGLO));
		log->logInfo("t_gprecisemodel", "sigPhaseGAL", format("%16.8f", _sigPhaseGAL));
		log->logInfo("t_gprecisemodel", "sigPhaseBDS", format("%16.8f", _sigPhaseBDS));
		log->logInfo("t_gprecisemodel", "sigPhaseQZS", format("%16.8f", _sigPhaseQZS));

		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(setting)->freq_index(gnut::QZS);

		// Get Data
		setOTL(dynamic_cast<t_gallotl*>((*data)[t_gdata::ALLOTL]));
		setOBJ(dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]));
		_gall_nav    = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gdata_erp   = dynamic_cast<t_gpoleut1*>((*data)[t_gdata::ALLPOLEUT1]);
		_gdata_navde = dynamic_cast<t_gnavde*>((*data)[t_gdata::ALLDE]);

		base_size = dynamic_cast<t_gsetgen*>(_settings)->list_base().size();
	}

	t_gprecisemodel::~t_gprecisemodel()
	{

	}

	void t_gprecisemodel::setOTL(t_gallotl* otl)
	{
		_tide->setOTL(otl);
	}
	double t_gprecisemodel::windUp(t_gsatdata& satdata, const ColumnVector& rRec)
	{
		gtrace("t_gprecisemodel::windUp");

		t_gtime epoch = satdata.epoch();
		string  prn   = satdata.sat();
		ColumnVector rSat = _trs_sat_crd.crd_cvect();

		double Mjd = epoch.mjd() + epoch.sod() / 86400.0;

		// Compute the correction for new time
		// -----------------------------------
		if (_phase_windup[_site][prn].find(epoch) == _phase_windup[_site][prn].end()) {

			// the last epoch
			double dphi0 = 0;
			if (_phase_windup[_site][prn].size() != 0) {
				dphi0 = _phase_windup[_site][prn].rbegin()->second;
				_phase_windup[_site][prn].clear();
			}
			ColumnVector rho = rRec - rSat;
			rho /= rho.norm_Frobenius();

			ColumnVector i(3);
			ColumnVector j(3);
			ColumnVector k(3);

			// attitude model
			string antype = "";
			if (_gallobj != 0) {
				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(satdata.sat());
				shared_ptr<t_gpcv>  sat_pcv;
				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(satdata.epoch());
				if (sat_pcv != 0)  antype = sat_pcv->anten();
			}
			if (_attitudes == ATTITUDES::YAW_NOMI)		   attitude(satdata, "", i, j, k);
			else if (_attitudes == ATTITUDES::YAW_RTCM) attitude(satdata, satdata.yaw(), i, j, k);
			else							   attitude(satdata, antype, i, j, k);

			if (antype.find("BLOCK IIR") != string::npos)
			{
				i = (-1.0)*i;
				j = (-1.0)*j;
			}
			double rlength = DotProduct(rho, i);
			ColumnVector dipSat = i - rlength * rho - crossproduct(rho, j);

			// Receiver unit Vectors rx, ry
			// ----------------------------
			ColumnVector rx(3);
			ColumnVector ry(3);

			double recEll[3]; xyz2ell(rRec.data(), recEll, false);
			double neu[3];

			neu[0] = 1.0;
			neu[1] = 0.0;
			neu[2] = 0.0;
			neu2xyz(recEll, neu, rx.data());

			neu[0] = 0.0;
			neu[1] = -1.0;
			neu[2] = 0.0;
			neu2xyz(recEll, neu, ry.data());

			// Effective Dipole of the Receiver Antenna
			// ----------------------------------------
			rlength = DotProduct(rho, rx);
			ColumnVector dipRec = rx - rlength * rho + crossproduct(rho, ry);

			// Resulting Effect
			// ----------------
			double alpha = DotProduct(dipSat, dipRec) / (dipSat.norm_Frobenius() * dipRec.norm_Frobenius());

			if (alpha > 1.0) alpha = 1.0;
			if (alpha < -1.0) alpha = -1.0;

			double dphi = acos(alpha) / 2.0 / G_PI;  // in cycles

			if (DotProduct(rho, crossproduct(dipSat, dipRec)) < 0.0) dphi = -dphi;

			_phase_windup[_site][prn][epoch] = floor(dphi0 - dphi + 0.5) + dphi;
		}

		satdata.addwind(_phase_windup[_site][prn][epoch]);
		return _phase_windup[_site][prn][epoch];
	}

	void t_gprecisemodel::set_multi_debug_output(string filename)
	{
		debug_output = make_shared<fstream>(filename,fstream::out);
	}

	bool t_gprecisemodel::_prepare_obs(const t_gtime& crt_epo, t_gallpar& pars)
	{
		_factorP = 1;
		if (!this->_gall_nav || !this->_gallobj)
		{
			write_log_info(_log, 1, "NOTE", "no navgation data or atx data for epoch " + crt_epo.str_ymdhms());
			return false;
		}

		// compute reciver time 
		_rec_epo = crt_epo - _crt_rec_clk;
		_crt_obs.addrecTime(_rec_epo);

		// get Rec crd
		bool apply_obj_valid;

		apply_obj_valid = _apply_rec(crt_epo, _rec_epo, pars);

		if (!apply_obj_valid)
		{
			write_log_info(_log, 1, "NOTE", "can not apply site in " + crt_epo.str_ymdhms());
			return false;
		}

		// get sat crd
		bool apply_sat_valid = _apply_sat(_rec_epo, _sat_epo);
		if (!apply_sat_valid)
		{
			write_log_info(_log, 1, "NOTE", "can not apply sat in " + crt_epo.str_ymdhms());
			return false;
		}
		_crt_obs.addsatTime(_sat_epo);

		// get f01 PCO
		this->_crs_rec_pco = this->_crs_rec_crd;
		this->_crs_sat_pco = this->_crs_sat_crd;

		shared_ptr<t_gobj>  sat_obj = this->_gallobj->obj(_crt_sat);
		shared_ptr<t_gobj>  rec_obj = this->_gallobj->obj(_crt_rec);

		shared_ptr<t_gpcv>  sat_pcv = (sat_obj != 0) ? sat_obj->pcv(crt_epo) : nullptr;
		shared_ptr<t_gpcv>  rec_pcv = (rec_obj != 0) ? rec_obj->pcv(crt_epo) : nullptr;

		if (!_isCalPCO)
		{
			sat_pcv = nullptr;
			rec_pcv = nullptr;
		}

		if (sat_pcv)
		{
			// Satellite phase center offset
			t_gtriple pco(0, 0, 0);
			if (sat_pcv->pcoS(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _rec_epo, sat_obj, true);
				this->_crs_sat_pco += t_gtriple(_rot_matrix*pco.crd_cvect());
			}

		}

		if (rec_pcv)
		{
			// Receiver phase center offset
			t_gtriple pco(0.0, 0.0, 0.0);
			if (rec_pcv->pcoR(_crt_obs, pco, _observ, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
			{
				Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _rec_epo, rec_obj, true);
				this->_crs_rec_pco += t_gtriple(_rot_matrix*pco.crd_cvect());
			}
		}
		
		// compute satclk[s]
		double sat_clk = _crt_sat_clk;

		// compute reldelay[m]
		double reldelay = 0.0;
		bool rel_valid = apply_reldelay(this->_crs_rec_pco, this->_crs_rec_vel, this->_crs_sat_pco, this->_crs_sat_vel, reldelay);
		if (!rel_valid)
		{
			return false;
		}
		_crt_obs.addclk(sat_clk * CLIGHT - reldelay); // include reldelay
		_crt_obs.addreldelay(reldelay);

		// set sat parameter PXSAT in lsq
		int sat_index = pars.getParam("", par_type::PXSAT, _crt_sat);
		if (sat_index >= 0) 
		{
			_sat_index = pars.getParIndex(sat_index);

			// get unit_vector in TRS
			ColumnVector unit_rec2sat(3);

			// change TRS to CRS
			unit_rec2sat = _crs_sat_pco.crd_cvect() - _crs_rec_pco.crd_cvect();
			unit_rec2sat = unit_rec2sat / unit_rec2sat.NormFrobenius();

			// order: P X0 Y0 Z0 VX0..VZ0 solarPar
			_sat_partial *= unit_rec2sat * 1E3;
		}

		// addrho
		double tmp = (_crs_sat_crd - _crs_rec_crd).norm();
		_crt_obs.addrho(tmp);

		// add drate
		_crt_obs.adddrate((DotProduct((_crs_sat_vel - _crs_rec_vel).crd_cvect(),
			(_crs_sat_pco - _crs_rec_pco).crd_cvect())) / (CLIGHT * tmp));

		// add azim && elev
		t_gtriple xyz_rho = _crs_sat_pco - _crs_rec_pco;
		t_gtriple ell_r, neu_s;

		Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, _rec_epo, _crt_obj, true);
		neu_s = t_gtriple(_rot_matrix.t()*xyz_rho.crd_cvect());

		double NE2 = neu_s[0] * neu_s[0] + neu_s[1] * neu_s[1];
		double ele = acos(sqrt(NE2) / _crt_obs.rho());

		if (sqrt(NE2) / _crt_obs.rho() > 1.0)
		{
			_crt_obs.addele(0.0);
		}
		else
		{
			_crt_obs.addele(ele);
		}

		double offnadir = dotproduct(xyz_rho.crd_cvect(), _crs_sat_pco.crd_cvect())/xyz_rho.norm()/_crs_sat_pco.norm();
		offnadir = acos(offnadir);
		_crt_obs.addnadir(offnadir);
		
		double azi = atan2(neu_s[1], neu_s[0]);
		if (azi < 0)
		{
			azi += 2 * G_PI;
		}
		_crt_obs.addazi(azi);

		Matrix _rot_matrix_scf2crs = _RotMatrix_Ant(_crt_obs, _sat_epo, sat_obj, true);
		t_gtriple xyz_s2r = t_gtriple(_rot_matrix_scf2crs.t()*((-1)*xyz_rho.crd_cvect())); // from sat. to rec. in SCF XYZ
		double azi_sat = atan2(xyz_s2r[0], xyz_s2r[1]);
		if (azi_sat < 0) azi_sat += 2 * G_PI;
		_crt_obs.addazi_sat(azi_sat);

		///add for another elev and azi used in calculating weight matric
		t_gtriple xyz_rh = _trs_sat_crd - _trs_rec_crd;
		t_gtriple ell_(0, 0, 0), neu_sa(0, 0, 0), xRec(0, 0, 0), xyz_s(0, 0, 0);
		xyz2ell(_trs_rec_crd, ell_, false);
		xyz2neu(ell_, xyz_rh, neu_sa);
		double rho0 = sqrt(pow(_trs_rec_crd[0] - xyz_s[0], 2) + pow(_trs_rec_crd[1] - xyz_s[1], 2) + pow(_trs_rec_crd[2] - xyz_s[2], 2));
		double dPhi = OMEGA * rho0 / CLIGHT;
		xRec[0] = _trs_rec_crd[0] * cos(dPhi) - _trs_rec_crd[1] * sin(dPhi);
		xRec[1] = _trs_rec_crd[1] * cos(dPhi) + _trs_rec_crd[0] * sin(dPhi);
		xRec[2] = _trs_rec_crd[2];
		double NE2_ = neu_sa[0] * neu_sa[0] + neu_sa[1] * neu_sa[1];
		double ele_ = acos(sqrt(NE2_) / _crt_obs.rho());
		_crt_obs.addele_leo(ele_);

		// check elevation cut-off
		if (_crt_obj->id_type() == t_gdata::REC && _crt_obs.ele_deg() < _minElev)
		{
			if (_log)
			{
				_log->comment(1, "t_gprecisemodel::prepareObs", "Prepare fail! the elevation is too small");
			}
			return false;
		}

		return true;
	}

	bool t_gprecisemodel::_update_obs_info(const t_gtime & epoch, t_gsatdata & obsdata, t_gallpar& pars)
	{
		if (!_gall_nav || !_gallobj)
		{
			write_log_info(_log, 1, "NOTE", "no pars");
			return false;
		}

		// Some frequently used variables are also defined here.
		_crt_epo = epoch;
		_crt_obs = obsdata;
		_crt_obs.tb12(obsdata.tb12());
		_crt_sat = obsdata.sat();
		_crt_rec = obsdata.site();
		_crt_sys = obsdata.gsys();
		_crt_obj = _gallobj->obj(_crt_rec);
		_crt_obj_type = _crt_obj->id_type();
		_site = _crt_rec;
		_grec = _crt_obj;
		_trs_sat_crd = obsdata.satcrd();

		bool rec_clk_valid = _update_obj_clk("rec"+_crt_rec, epoch, pars, _crt_rec_clk, _obj_clk);
		bool sat_clk_valid = _update_obj_clk("sat"+_crt_sat, epoch, pars, _crt_sat_clk, _obj_clk);
		if (sat_clk_valid) _crt_obs.addnavclk(_crt_sat_clk);

		if (!rec_clk_valid || !sat_clk_valid)
		{
			write_log_info(_log, 1, "NOTE", _crt_rec + " " + _crt_sat + " no rec or sat clk for epoch " + epoch.str_ymdhms());
			return false;
		}
		return true;
	}

	bool t_gprecisemodel::_update_obs_info(t_gsatdata & obsdata)
	{
		// smgcz
		obsdata = _crt_obs;
		return true;
	}

	void t_gprecisemodel::update_obj_clk(const string & obj, const t_gtime & epo, double clk)
	{
		_obj_clk[obj].first = epo;
		_obj_clk[obj].second = clk;
		_rec_clk[obj] = clk;
	}

	bool t_gprecisemodel::_update_obj_clk(const string& obj, const t_gtime& crt_epoch,
		t_gallpar& pars, double& clk, map<string, pair<t_gtime, double>>& obj_clk)
	{
		if (!_gall_nav || obj.substr(3).empty())
		{
			write_log_info(_log, 1, "NOTE", "no navigation files or obj is NONE");
			return false;
		}

		// clk[s] from clk files
		double clk_rms = 0.0;
		double dclk = 0.0;
		string type = obj.substr(0, 3);
		string name = obj.substr(3);

		// get clk from clk files
		int idx_clk = -1;
		if (type == "rec")
		{
			idx_clk = pars.getParam(name, par_type::CLK, "");
		}
		else if (type == "sat")
		{
			idx_clk = pars.getParam("", par_type::CLK_SAT, name);
		}

		// update clk
		if (obj_clk[name].first != crt_epoch)
		{
			int clk_valid = _gall_nav->clk(name, crt_epoch, &clk, &clk_rms, &dclk);
			if (clk_valid < 0) {
				if (type == "sat" && idx_clk >= 0) return false;
				clk = 0.0;
			}
			if (type == "sat" && (idx_clk < 0) &&
				(base_size == 0))
			{
				write_log_info(_log, 1, "NOTE", "no sat" + name + "  clk correction for " + crt_epoch.str_ymdhms());
				return false;
			}

			obj_clk[name].first = crt_epoch;
			obj_clk[name].second = clk;
		}
		else
		{
			clk = obj_clk[name].second;
		}
	

		if (double_eq(clk, 0.0) && idx_clk < 0)
		{
			write_log_info(_log, 1, "NOTE", "no rec clk for " + crt_epoch.str_ymdhms());
			return false;
		}
		//-- no clk file but clk par exist
		if (double_eq(clk, 0.0) && idx_clk >= 0)
		{
			clk = pars[idx_clk].value() / CLIGHT;
		}
		//-- update clk par
		if (idx_clk >= 0)
		{
			pars[idx_clk].value(clk * CLIGHT);
		}

		return true;
	}

	bool t_gprecisemodel::_omc_obs_ALL(const t_gtime& crt_epo, t_gallpar& pars, t_gobs& gobs, double& omc)
	{
		double obs = 0.0;
		if (gobs.is_code())
		{
			obs = _crt_obs.obs_C(gobs);
		}
		else if (gobs.is_phase())
		{
			obs = _crt_obs.obs_L(gobs);
		}
		else
		{
			throw logic_error("only support code and phase observation.");
		}

		if (double_eq(obs, 0.0))
		{
			_log->logDebug("t_gprecisemodel", "_omc_obs_ALL", "obs type is empty, the obs is zero!");
			return false;
		}

		t_gtime epo = crt_epo;
		double ModelObs = cmpObs(epo, pars, _crt_obs, gobs);
		if (ModelObs < 0)
		{
			_log->logDebug("t_gprecisemodel", "_omc_obs_ALL", "Compute ModelObs error!");
			return false;
		}

		double Amb = 0.0;
		omc = (obs - ModelObs - Amb);

		return true;
	}

	bool t_gprecisemodel::_wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt)
	{
		GOBSTYPE type = gobs1.type();
		if (type != TYPE_C &&
			type != TYPE_P &&
			type != TYPE_L)
		{
			_log->logError("t_gprecisemodel", "_wgt_obs", " type should be TYPE_C/TYPE_P/TYPE_L");
			return false;
		}

		double factor = 0.0;
		double sigRange = 0.0;
		double sigPhase = 0.0;

		factor = factorP;
		if (_crt_sys == GSYS::BDS && t_gsys::bds_geo(_crt_sat)) factor = 5.0;          

		// get the sys sigRange
		if (type == TYPE_C || type == TYPE_P)
		{
			if (_crt_obj_type == t_gdata::REC)
			{
				switch (_crt_sys)
				{

				case GPS: sigRange = _sigCodeGPS; break;
				case GLO: sigRange = _sigCodeGLO; break;
				case GAL: sigRange = _sigCodeGAL; break;
				case BDS: sigRange = _sigCodeBDS; break;
				case QZS: sigRange = _sigCodeQZS; break;
				default: sigRange = 0.0; return false;
				}
			}
			else
			{
				sigRange = 0.0;
				return false;
			}
		}
		else if (type == TYPE_L)
		{
			if (_crt_obj_type == t_gdata::REC)
			{
				switch (_crt_sys)
				{
				case GPS: sigPhase = _sigPhaseGPS; break;
				case GLO: sigPhase = _sigPhaseGLO; break;
				case GAL: sigPhase = _sigPhaseGAL; break;
				case BDS: sigPhase = _sigPhaseBDS; break;
				case QZS: sigPhase = _sigPhaseQZS; break;
				default: sigPhase = 0.0; return false;
				}
			}
			else
			{
				sigPhase = 0.0;
				return false;
			}
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "no obs type");
			return	false;
		}

		if (_crt_obj_type == t_gdata::REC)
		{
			double leo_ele = _crt_obs.ele_leo();
			double leo_ele_deg = _crt_obs.ele_leo_deg();
			double sin_leo_ele = sin(leo_ele);

			// get the _weight factor
			switch (_weight)
			{
			case OBSWEIGHT::DEF_OBSWEIGHT:
				cerr << "gspplsq: WeightObs (default) should not happened!\n";
				break;
			case OBSWEIGHT::EQUAL:
				factor *= 1;
				break;
			case OBSWEIGHT::SINEL:
				factor *= 1.0 / 2.0 / sin(leo_ele);
				break;
			case OBSWEIGHT::SINEL2:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 2);
				break;
			case OBSWEIGHT::SINEL4:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 4);
				break;
			case OBSWEIGHT::PARTELE:
				factor = (leo_ele_deg <= 30.0) ? factor * (1.0 / 2.0 / sin_leo_ele) : factor;
				break;
			default:
				cerr << "gspplsq: we can't deal with this WeightObs method!";
				return false;
			}
		}

		// get the combine obs factor
		auto b1 = gobs1.band();
		if (type == TYPE_C || type == TYPE_P)  
		{
			sigRange = sigRange /** 3*/;
		}
		else if (type == TYPE_L)
		{
			sigPhase = sigPhase * _crt_obs.wavelength(b1) /** 3*/;
		}
		else
		{
			sigRange = 0.0;
			sigPhase = 0.0;
		}

		if (sigPhase == 0.0 && sigRange == 0.0) return false;

		switch (type)
		{
		case TYPE_L:
			if (sigPhase == 0.0) return false;
			wgt = 1.0 / pow(factor*sigPhase, 2);
			break;
		case TYPE_C:
		case TYPE_P:
			if (sigRange == 0.0) return false;
			wgt = 1.0 / pow(factor * sigRange, 2);
			if (_crt_obs.getoutliers(gobs1.gobs()) >= 1) wgt *= 0.2;
			break;
		default:
			return false;
		}
		return true;
	}

	bool t_gprecisemodel::_wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt, double snr_value)
	{
		GOBSTYPE type = gobs1.type();
		if (type != TYPE_C &&
			type != TYPE_P &&
			type != TYPE_L)
		{
			_log->logError("t_gprecisemodel", "_wgt_obs", " type should be TYPE_C/TYPE_P/TYPE_L");
			return false;
		}

		double factor = 0.0;
		double sigRange = 0.0;
		double sigPhase = 0.0;

		factor = factorP;
		if (_crt_sys == GSYS::BDS) factor = 2.0;
		if (_crt_sys == GSYS::BDS && t_gsys::bds_geo(_crt_sat)) factor = 5.0;

		// get the sys sigRange
		if (type == TYPE_C || type == TYPE_P)
		{
			if (_crt_obj_type == t_gdata::REC)
			{
				switch (_crt_sys)
				{

				case GPS: sigRange = _sigCodeGPS; break;
				case GLO: sigRange = _sigCodeGLO; break;
				case GAL: sigRange = _sigCodeGAL; break;
				case BDS: sigRange = _sigCodeBDS; break;
				case QZS: sigRange = _sigCodeQZS; break;
				default: sigRange = 0.0; return false;
				}
			}
			else
			{
				sigRange = 0.0;
				return false;
			}
		}
		else if (type == TYPE_L)
		{
			if (_crt_obj_type == t_gdata::REC)
			{
				switch (_crt_sys)
				{
				case GPS: sigPhase = _sigPhaseGPS; break;
				case GLO: sigPhase = _sigPhaseGLO; break;
				case GAL: sigPhase = _sigPhaseGAL; break;
				case BDS: sigPhase = _sigPhaseBDS; break;
				case QZS: sigPhase = _sigPhaseQZS; break;
				default: sigPhase = 0.0; return false;
				}
			}
			else
			{
				sigPhase = 0.0;
				return false;
			}
		}
		else
		{
			write_log_info(_log, 1, "NOTE", "no obs type");
			return	false;
		}

		if (_crt_obj_type == t_gdata::REC)
		{
			double leo_ele = _crt_obs.ele_leo();
			double leo_ele_deg = _crt_obs.ele_leo_deg();
			double sin_leo_ele = sin(leo_ele);
			// get the _weight factor
			switch (_weight)
			{
			case OBSWEIGHT::DEF_OBSWEIGHT:
				cerr << "gspplsq: WeightObs (default) should not happened!\n";
				break;
			case OBSWEIGHT::EQUAL:
				factor *= 1;
				break;
			case OBSWEIGHT::SINEL:
				factor *= 1.0 / 2.0 / sin(leo_ele);
				break;
			case OBSWEIGHT::SINEL2:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 2);
				break;
			case OBSWEIGHT::SINEL4:
				factor *= 1.0 / 2.0 / pow(sin_leo_ele, 4);
				break;
			case OBSWEIGHT::PARTELE:
				factor = (leo_ele_deg <= 30.0) ? factor * (1.0 / 2.0 / sin_leo_ele) : factor;
				break;
			case OBSWEIGHT::SNR:
				factor *= sqrt(134.02 * pow(10, -(snr_value / 17.91))); 
				break;
			default:
				cerr << "gspplsq: we can't deal with this WeightObs method!";
				return false;
			}
		}

		// get the combine obs factor
		auto b1 = gobs1.band();
		if (type == TYPE_C || type == TYPE_P)
		{
			sigRange = sigRange /** 2*/;
		}
		else if (type == TYPE_L)
		{
			sigPhase = sigPhase * _crt_obs.wavelength(b1) /** 2*/;
		}
		else
		{
			sigRange = 0.0;
			sigPhase = 0.0;
		}

		if (sigPhase == 0.0 && sigRange == 0.0) return false;

		switch (type)
		{
		case TYPE_L:
			if (sigPhase == 0.0) return false;
			wgt = 1.0 / pow(factor * sigPhase, 2);
			break;
		case TYPE_C:
		case TYPE_P:
			if (sigRange == 0.0) return false;
			wgt = 1.0 / pow(factor * sigRange, 2);
			if (_crt_obs.getoutliers(gobs1.gobs()) >= 1) wgt *= 0.2;
			break;
		default:
			return false;
		}
		return true;
	}



	bool t_gprecisemodel::_prt_obs(const t_gtime & epoch, t_gallpar& pars, t_gobs & gobs, vector<pair<int, double>>& coeff)
	{
		t_gtriple groundEll;
		xyz2ell(_trs_rec_crd, groundEll, false);
		int parNum = pars.parNumber();

		auto par_list = pars.getPartialIndex(_crt_rec, _crt_sat);
		for (int ipar: par_list)
		{
			double coeff_value = _Partial(epoch, _crt_obs, gobs, pars.getPar(ipar));
			if (coeff_value != 0.0)
			{
				coeff.push_back(make_pair(ipar+1, coeff_value));
			}
		}
		return true;
	}




	int t_gprecisemodel::outlierDetect(vector<t_gsatdata>& data)
	{
		vector<t_gsatdata>::iterator itMaxVcodeNORM = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseNORM = data.end();

		vector<t_gsatdata>::iterator itDataErase = data.end();

		double maxVcodeNORM = 0.0;
		double maxVphaseNORM = 0.0;

		double maxVcodeORIG = 0.0;
		double maxVphaseORIG = 0.0;

		// find maximal code/phase residuals
		maxVcodeNORM = _maxres(false, data, itMaxVcodeNORM, RESIDTYPE::RES_NORM);
		maxVphaseNORM = _maxres(true, data, itMaxVphaseNORM, RESIDTYPE::RES_NORM);

		int maxVphaseNORM_loc = itMaxVphaseNORM - data.begin() + 1;
		int maxVcodeNORM_loc = itMaxVcodeNORM - data.begin() + 1;


		if (_check_outl(true, maxVphaseNORM, itMaxVphaseNORM, maxVphaseORIG, itMaxVphaseNORM, itDataErase, data)) { data.erase(itDataErase); return maxVphaseNORM_loc; }
		if (_check_outl(false, maxVcodeNORM, itMaxVcodeNORM, maxVcodeORIG, itMaxVcodeNORM, itDataErase, data)) { data.erase(itDataErase);  return maxVcodeNORM_loc; }

		return 0;

	}

	double t_gprecisemodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs)
	{
		gtrace("t_gpppmodel::cmpObs");
		if (_gallobj == nullptr)
		{
			throw logic_error("gallobj is nullptr in precisemodel cmpObs");
		}

		double spp_model = t_gsppmodel::cmpObs(epoch, param, obsdata, gobs);
		if (spp_model < 0) return -1;

		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple ell;
		xyz2ell(_trs_rec_crd, ell, false);

		// band used
		GOBSBAND band = gobs.band();

		// Wind up correction 
		double wind = 0.0;
		if (gobs.is_phase())
		{
			double wavelength = obsdata.wavelength(gobs.band());
			if (fabs(obsdata.wind()) > 0)
			{
				wind = obsdata.wind() * wavelength;
			}
			else
			{
				wind = windUp(obsdata, _trs_rec_crd.crd_cvect()) * wavelength;
			}
		}

		// ion correction
		double ion = 0.0;
		auto band_1 = _band_index[obsdata.gsys()][FREQ_1];
		ion = ionDelay(epoch, param, obsdata, _ion_model, band_1, gobs);

		// Phase center variation correction
		double pcv_R = 0.0;
		double pcv_S = 0.0;
		double pco_R = 0.0;
		double pco_S = 0.0;

		if (_isCalPCO) {

			shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
			shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);

			shared_ptr<t_gpcv>  sat_pcv = (sat_obj != nullptr) ? sat_obj->pcv(epoch) : nullptr;
			shared_ptr<t_gpcv>  rec_pcv = (rec_obj != nullptr) ? rec_obj->pcv(epoch) : nullptr;

			if (sat_pcv != nullptr)
			{
				// Satellite phase center variation
				// -- Satellite phase center offset
				t_gtriple pco(0, 0, 0);
				if (sat_pcv->pcoS_raw(obsdata, pco, band) > 0)
				{
					string antenna = sat_pcv->anten();
					t_gtriple dx(0, 0, 0);
					ColumnVector i(3);
					ColumnVector j(3);
					ColumnVector k(3);
					if (_attitudes == ATTITUDES::YAW_NOMI)	   attitude(obsdata, "", i, j, k);
					else if (_attitudes == ATTITUDES::YAW_RTCM) attitude(obsdata, obsdata.yaw(), i, j, k);
					else							   attitude(obsdata, antenna, i, j, k);
					dx[0] = pco[0] * i(1) + pco[1] * j(1) + pco[2] * k(1);
					dx[1] = pco[0] * i(2) + pco[1] * j(2) + pco[2] * k(2);
					dx[2] = pco[0] * i(3) + pco[1] * j(3) + pco[2] * k(3);
					sat_pcv->pco_proj(pco_S, obsdata, _trs_rec_crd, dx);
					obsdata.addpco(dx);
				}
			}

			if (rec_pcv != nullptr)
			{
				// Receiver phase center offset
				t_gtriple pco(0.0, 0.0, 0.0);
				if (rec_pcv->pcoR_raw(obsdata, pco, band) > 0)
				{
					Matrix _rot_matrix = _RotMatrix_Ant(obsdata, _rec_epo, rec_obj, false);
					t_gtriple dx(_rot_matrix * (pco.crd_cvect()));
					rec_pcv->pco_proj(pco_R, obsdata, _trs_rec_crd, dx);
				}
				pco_R *= -1;
			}

			if (gobs.is_phase())
			{
				if (sat_pcv != 0)
				{
					// Satellite phase center variation
					sat_pcv->pcvS_raw(pcv_S, obsdata, band, _trs_rec_crd);
				}
				if (rec_pcv != 0)
				{
					// Receiver phase center variation
					rec_pcv->pcvR_raw(pcv_R, obsdata, band);
				}
			}
		}

		return spp_model +
			wind +
			pco_R +
			pco_S +
			pcv_R +
			pcv_S + 
			ion;
	}

	double t_gprecisemodel::getZHD(const string& site, const t_gtime& epo)
	{
		t_gtriple xyz, ell;
		_grec = _gallobj->obj(site);
		xyz = _grec->crd(epo);
		xyz2ell(xyz, ell, false);
		double ZHD = this->tropoModel()->getZHD(ell, epo);
		return ZHD;

	}
	double t_gprecisemodel::getZWD(const string& site, const t_gtime& epo)
	{
		t_gtriple xyz, ell;
		_grec = _gallobj->obj(site);
		xyz = _grec->crd(epo);
		xyz2ell(xyz, ell, false);
		double ZWD = this->tropoModel()->getZWD(ell, epo);
		return ZWD;
	}

	double t_gprecisemodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata)
	{
		gtrace("t_gprecisemodel::tropoDelay");

		if (site_ell[2] > 1E4)
		{
			return 0.0;
		}
		else
		{
			if (_tropoModel == 0)
			{
				if (_log)
				{
					_log->comment(0, "gppp", "Tropo Model setting is not correct. Default used! Check config.");
				}
				else
				{
					cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
				}
				_tropoModel = make_shared<t_saast>();
			}

			double ele = satdata.ele();
			double azi = satdata.azi();

			double delay = 0.0;
			double zwd = 0.0;
			double zhd = 0.0;
			int i, j, k;
			t_gtriple ell;

			xyz2ell(_trs_rec_crd, ell, false);

			if (ell[2] > 1E4) {
				return 100.0;
			}

			i = param.getParam(_site, par_type::TRP, "");
			j = param.getParam(_site, par_type::GRD_N, "");
			k = param.getParam(_site, par_type::GRD_E, "");

			if (i >= 0)
			{
				zwd = param[i].value();
				
				if (param[i].apriori()>1E-4 && (zwd == 0.0 || epoch == param[i].beg))
				{
					zwd = _tropoModel->getZWD(ell, epoch);
					param[i].value(zwd);
					
				}
				zhd = _tropoModel->getZHD(ell, epoch);
					param[i].zhd = zhd;
			}
			else
			{
				if (_tropoModel != 0)
				{
					zwd = _tropoModel->getZWD(ell, epoch);
					zhd = _tropoModel->getZHD(ell, epoch);
				}
			}

			double mfh, mfw, dmfh, dmfw;
			mfh = mfw = dmfh = dmfw = 0.0;
			if (_tropo_mf == ZTDMPFUNC::GMF)
			{
				t_gmf mf;
				mf.gmf(epoch.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
					mfh, mfw, dmfh, dmfw);
			}
			else if (_tropo_mf == ZTDMPFUNC::COSZ)
			{
				mfh = mfw = 1 / sin(ele);
				dmfh = dmfw = -(cos(ele)) / (sin(ele) * sin(ele));
			}
			else return 0.0;

			satdata.addmfH(mfh);
			satdata.addmfW(mfw);

			delay = mfh * zhd + mfw * zwd;

			double grdN, grdE;
			grdN = grdE = 0.0;

			if (j >= 0 && k >= 0) {
				if (_grad_mf == GRDMPFUNC::TILTING) {
					grdN = param[j].value() * dmfw * cos(azi);
					grdE = param[k].value() * dmfw * sin(azi);
					satdata.addmfG(dmfw);
				}
				else if (_grad_mf == GRDMPFUNC::CHEN_HERRING) {
					double mfg = 1.0 / (sin(ele)*tan(ele) + 0.0032);
					grdN = param[j].value() * 1000.0 * mfg * cos(azi);  grdN /= 1000.0;
					grdE = param[k].value() * 1000.0 * mfg * sin(azi);  grdE /= 1000.0;
					satdata.addmfG(mfg);
				}
				else if (_grad_mf == GRDMPFUNC::BAR_SEVER) {
					double mfg = mfw * (1 / tan(ele));
					grdN = param[j].value() * mfg * cos(azi);
					grdE = param[k].value() * mfg * sin(azi);
					satdata.addmfG(mfg);
				}

				delay += grdN + grdE;
			}
			return delay;
		}
	}

	double t_gprecisemodel::ionDelay(t_gtime & epoch, t_gallpar & param, t_gsatdata& satdata, IONMODEL& ion_model, GOBSBAND& band_1, t_gobs& gobs)
	{
		gtrace("t_gsppmodel::ionoDelay");
		if (band_1 == BAND || band_1 == BAND_A || band_1 == BAND_B || band_1 == BAND_C || band_1 == BAND_D)
		{
			return 0.0;
		}

		if (gobs.band() == BAND || gobs.band() == BAND_A || gobs.band() == BAND_B || gobs.band() == BAND_C || gobs.band() == BAND_D)
		{
			return 0.0;
		}

		if (ion_model == IONMODEL::VION)
		{
			throw logic_error(" not support VION model, use SION in XML");
		}

		double iono_delay = 0.0;

		double f1 = satdata.frequency(band_1);;
		double fk = satdata.frequency(gobs.band());
		double alfa = 0.0;

		if (gobs.is_phase())
		{
			alfa = -(f1 * f1) / (fk * fk);
		}
		if (gobs.is_code())
		{
			alfa = (f1 * f1) / (fk * fk);
		}

		// ionosphere slant delay parameter
		int i = param.getParam(_site, par_type::SION, satdata.sat());
		if (i >= 0)
		{
			iono_delay = alfa * param[i].value();
		}

		return iono_delay;
	}

	void t_gprecisemodel::_update_rot_matrix(const t_gtime& epoch)
	{
		if (_gdata_erp->isEmpty()) return;
		if (!_gdata_erp)
		{
			if (_log)
				_log->comment(1, "t_gorbmodel::_update_rotmatrix ", "have no poleut1 data,cant' compute trs2crs!!");
			throw exception();
		}

		t_gtime tdt = epoch;
		tdt.tsys(t_gtime::TT);

		auto find_iter = _trs2crs_list.find(tdt);
		if (find_iter == _trs2crs_list.end())
		{
			_trs2crs_2000 = make_shared<t_gtrs2crs>(false, _gdata_erp);
			_trs2crs_2000->calcRotMat(tdt, true, true, true);
			_trs2crs_list.insert(make_pair(tdt, _trs2crs_2000));

			auto before_iter = _trs2crs_list.lower_bound(tdt - 300.0);
			if (before_iter != _trs2crs_list.begin())
			{
				_trs2crs_list.erase(_trs2crs_list.begin(), --before_iter);
			}
		}
		else 
		{
			_trs2crs_2000 = find_iter->second;
		}

		return;

	}

	bool t_gprecisemodel::_apply_rec(const t_gtime& crt_epo, const t_gtime& rec_epo, t_gallpar& pars)
	{
		if (_rec_obj_flag == make_pair(_crt_rec, rec_epo))
		{
			return true;
		}
		else
		{
			_rec_obj_flag = make_pair(_crt_rec, rec_epo);
		}

		int ix = pars.getParam(_crt_rec, par_type::CRD_X, "");
		int iy = pars.getParam(_crt_rec, par_type::CRD_Y, "");
		int iz = pars.getParam(_crt_rec, par_type::CRD_Z, "");

		t_gtriple trs_rec_xyz(0.0, 0.0, 0.0);
		if (ix >= 0 && iy >= 0 && iz >= 0)
		{
			trs_rec_xyz = t_gtriple(pars.getParValue(ix), pars.getParValue(iy), pars.getParValue(iz));
			if (trs_rec_xyz.zero())
			{
				trs_rec_xyz = _crt_obj->crd(rec_epo);
				pars.setParValue(ix, trs_rec_xyz[0]);
				pars.setParValue(iy, trs_rec_xyz[1]);
				pars.setParValue(iz, trs_rec_xyz[2]);
			}
		}
		else
		{
			if (_crd_est != CONSTRPAR::FIX)
			{
				write_log_info(_log, 1, "NOTE", "can't get the crd of" + _site + " in time" + rec_epo.str_ymdhms());
				return false;
			}
			trs_rec_xyz = _crt_obj->crd(rec_epo);
		}

		if (trs_rec_xyz.zero())
		{
			write_log_info(_log, 1, "NOTE", "can't get the crd of" + _site + " in time" + rec_epo.str_ymdhms());
			return false;
		}

		bool tide_valid = _apply_rec_tides(rec_epo, trs_rec_xyz);
		if (!tide_valid)
		{
			write_log_info(_log, 1, "NOTE", "apply tide failed for " + _site + " in time" + rec_epo.str_ymdhms());
			return false;
		}

		string strEst = dynamic_cast<t_gsetgen*>(_settings)->estimator();
		bool isFLT = (strEst == "FLT");
		// ARP Correction
		if(!isFLT || _crd_est == CONSTRPAR::FIX) trs_rec_xyz = trs_rec_xyz + _crt_obj->eccxyz(rec_epo);
		_trs_rec_crd = trs_rec_xyz;

		// TRS2CRS 
		_update_rot_matrix(rec_epo);
		Matrix trs2crs = _trs2crs_2000->getRotMat();
		ColumnVector crs_rec_xyz = trs2crs * trs_rec_xyz.crd_cvect();
		_crs_rec_crd = t_gtriple(crs_rec_xyz);

		Matrix dtrs2crs = _trs2crs_2000->getMatDu() * OMGE_DOT;
		ColumnVector trs_rec_vel = dtrs2crs * trs_rec_xyz.crd_cvect();
		_crs_rec_vel = t_gtriple(trs_rec_vel);

		return true;
	}

	bool t_gprecisemodel::_apply_sat(const t_gtime& rec_epo, t_gtime& sat_epo)
	{
		// ITERATION
		// compute sat coord(CRS)  clk(estimated) (rewrite in orb model)
		double delay = 0.0;
		while (true)
		{
			sat_epo = rec_epo - delay;

			// Get CRS
			bool sat_pos_valid = _get_crs_sat_crd(sat_epo, _crt_sat, _crs_sat_crd);
			if (!sat_pos_valid)
			{
				write_log_info(_log, 1, "NOTE", "can not get sat pos for " + _crt_sat);
				return false;
			}
			if (double_eq(_crs_sat_crd[0] * _crs_sat_crd[1] * _crs_sat_crd[2], 0.0) || abs(_crs_sat_crd[0]) >= 1E18)
			{
				write_log_info(_log, 1, "NOTE", "can not get sat pos for " + _crt_sat);
				return false;
			}
			// SET TRS in epoch TR [include earth rotation]
			_update_rot_matrix(rec_epo);
			_trs_sat_crd = t_gtriple(_trs2crs_2000->getRotMat().t() * _crs_sat_crd.crd_cvect());

			bool sat_vel_valid = _get_crs_sat_vel(sat_epo, _crt_sat, _crs_sat_vel);
			if (!sat_vel_valid)
			{
				write_log_info(_log, 1, "NOTE", "can not get sat vel for " + _crt_sat);
				return false;
			}

			// PCO corr sat
			double pco_R = 0.0, pco_S = 0.0;
			if (_gallobj != 0)
			{
				_update_rot_matrix(rec_epo);
				_crt_obs.addcrd(_trs_sat_crd);
				ColumnVector x_earth = (_crs_sat_crd.crd_cvect().t()*_trs2crs_2000->getMatDu() / RAD2TSEC).t();
				_crt_obs.addvel(t_gtriple(_trs2crs_2000->getRotMat().t()*_crs_sat_vel.crd_cvect() - x_earth));
				shared_ptr<t_gobj>  sat_obj = _gallobj->obj(_crt_sat);
				shared_ptr<t_gobj>  rec_obj = _gallobj->obj(_crt_rec);

				shared_ptr<t_gpcv>  sat_pcv;
				shared_ptr<t_gpcv>  rec_pcv;

				if (sat_obj != 0)  sat_pcv = sat_obj->pcv(_crt_epo);
				if (rec_obj != 0)  rec_pcv = rec_obj->pcv(_crt_epo);

				GOBS_LC lc = LC_L1;
				if (sat_pcv)
				{
					// Satellite phase center offset
					t_gtriple pco(0, 0, 0);

					if (sat_pcv->pcoS(_crt_obs, pco, lc, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
					{
						Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, sat_epo, sat_obj, false);
						t_gtriple dx(_rot_matrix*(pco.crd_cvect()));
						sat_pcv->pco_proj(pco_S, _crt_obs, _trs_rec_crd, dx);
					}

				}

				if (rec_pcv)
				{
					// Receiver phase center offset
					t_gtriple pco(0.0, 0.0, 0.0);
					if (rec_pcv->pcoR(_crt_obs, pco, lc, _band_index[_crt_sys][FREQ_1], _band_index[_crt_sys][FREQ_2]) > 0)
					{
						Matrix _rot_matrix = _RotMatrix_Ant(_crt_obs, rec_epo, rec_obj, false);
						t_gtriple dx(_rot_matrix*(pco.crd_cvect()));
						rec_pcv->pco_proj(pco_R, _crt_obs, _trs_rec_crd, dx);
					}
					pco_R *= -1;
				}
			}

			double delay_temp;
			delay_temp = (_crs_sat_crd - _crs_rec_crd).norm();
			delay_temp += pco_R + pco_S;

			if (abs(delay_temp / CLIGHT - delay) < 1E-9) break;
			delay = delay_temp / CLIGHT;
		}

		_crt_obs.addcrd(_trs_sat_crd);
		_crt_obs.addcrdcrs(_crs_sat_crd);
		_crt_obs.addvel_crs(_crs_sat_vel);

		ColumnVector x_earth = (_crs_sat_crd.crd_cvect().t()*_trs2crs_2000->getMatDu() / RAD2TSEC).t();
		_crt_obs.addvel(t_gtriple(_trs2crs_2000->getRotMat().t()*_crs_sat_vel.crd_cvect() - x_earth));
		return true;
	}


	bool t_gprecisemodel::_get_crs_sat_crd(const t_gtime & sat_epoch, const string & sat, t_gtriple & crs_sat_crd)
	{
		if (_gall_nav)	return _get_crs_sat_crd(sat_epoch, sat, _gall_nav, crs_sat_crd);
		write_log_info(_log, 2, "NOTE", "there are no orb and nav");
		return false;
	}

	bool t_gprecisemodel::_get_crs_sat_crd(const t_gtime & sat_epoch, const string & sat, t_gallnav * nav, t_gtriple & crs_sat_crd)
	{
		bool pos_valid = false;
		if (nav)
		{
			double xyz_sat[3];
			pos_valid = (nav->pos(sat, sat_epoch, xyz_sat) >= 0) ? true : false;

			// TRS2CRS
			_update_rot_matrix(sat_epoch);
			crs_sat_crd = t_gtriple(xyz_sat);
			crs_sat_crd = t_gtriple(_trs2crs_2000->getRotMat()*crs_sat_crd.crd_cvect());
		}

		if (pos_valid) return pos_valid;
		else write_log_info(_log, 2, "NOTE", sat_epoch.str_ymdhms("can not get sat crd in epoch",false));
		return false;
	}

	bool t_gprecisemodel::_get_crs_sat_vel(const t_gtime& sat_epoch, string sat, t_gtriple& sat_vel)
	{
		if (!_gall_nav)
		{
			write_log_info(_log, 1, "NOTE", "cannot get vel " + sat);
			return false;
		}

		// get vel in CRS
		double xyz0[3], xyz1[3];
		if (!_gall_nav->pos(sat, sat_epoch, xyz0)) return false;
		if (!_gall_nav->pos(sat, sat_epoch + 1, xyz1)) return false;
		t_gtriple vel_trs(xyz1[0]- xyz0[0], xyz1[1] - xyz0[1], xyz1[2] - xyz0[2]);
		ColumnVector x_earth = (_crs_sat_crd.crd_cvect().t() * _trs2crs_2000->getMatDu() / RAD2TSEC).t();
		sat_vel = t_gtriple(_trs2crs_2000->getRotMat() * (vel_trs.crd_cvect() + x_earth));

		return true;

		return true;
	}

	double t_gprecisemodel::_Partial(const t_gtime& epoch, t_gsatdata& obsdata, const t_gobs & gobs, t_gpar& par)
	{
		double mfw, dmfw, mfh, dmfh;
		mfw = dmfw = mfh = dmfh = 0.0;
		t_gtriple ell(0.0, 0.0, 0.0);

		switch (par.parType)
		{
		case par_type::CRD_X: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(0)) / obsdata.rho();
		}
							else return 0.0;
		case par_type::CRD_Y: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(1)) / obsdata.rho();
		}
							else return 0.0;
		case par_type::CRD_Z: if (obsdata.site() == par.site)
		{
			return (par.value() - obsdata.satcrd().crd(2)) / obsdata.rho();
		}
							else return 0.0;
		case par_type::CLK:
		{
			if (obsdata.site() == par.site) return 1.0 - obsdata.drate();
			else return 0.0;
		}

		case par_type::CLK_SAT:    if (obsdata.sat() == par.prn)              return -1.0; else return 0.0;
		case par_type::TRP:
			if (obsdata.site() == par.site) {
				xyz2ell(_trs_rec_crd, ell, false);
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				return mfw;
			}
			else return 0.0;
		case par_type::SION:
		{
			if (obsdata.site() == par.site)
			{
				auto gsys = obsdata.gsys();
				double f1 = obsdata.frequency(_band_index[gsys][FREQ_1]);
				double fk = obsdata.frequency(gobs.band());
				double alfa = 0.0;
				if (gobs.is_phase() && par.prn == obsdata.sat())
				{
					alfa = -(f1 * f1) / (fk * fk);
				}
				if (gobs.is_code() && par.prn == obsdata.sat())
				{
					alfa = (f1 * f1) / (fk * fk);
				}
				return alfa;
			}
			else return 0.0;
		}
		case par_type::VION:
		{
			if (obsdata.site() == par.site)
			{
				auto gsys = obsdata.gsys();
				double f1 = obsdata.frequency(_band_index[gsys][FREQ_1]);
				double fk = obsdata.frequency(gobs.band());
				double mf = 1.0 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - obsdata.ele()), 2));
				double alfa = 0.0;
				if (gobs.is_phase() && par.prn == obsdata.sat()) { alfa = -(f1 * f1) / (fk * fk); }
				if (gobs.is_code() && par.prn == obsdata.sat()) { alfa = (f1 * f1) / (fk * fk); }
				return alfa * mf;
			}
			else return 0.0;
		}
		case par_type::P1P2G_REC:
		{
			if (obsdata.site() == par.site)
			{
				double f1 = G01_F;
				double fk = obsdata.frequency(gobs.band());
				double alfa = (f1 * f1) / (fk * fk);
				double beta = (G02_F * G02_F) / (G01_F * G01_F - G02_F * G02_F);
				FREQ_SEQ freq = t_gsys::band2freq(obsdata.gsys(), gobs.band());
				if (obsdata.gsys() == GPS && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
				else return 0.0;
			}
			else return 0.0;
		}
		case par_type::P1P2E_REC:
		{
			if (obsdata.site() == par.site) {
				double f1 = E01_F;
				double fk = obsdata.frequency(gobs.band());
				double alfa = (f1 * f1) / (fk * fk);
				double beta = (E05_F * E05_F) / (E01_F * E01_F - E05_F * E05_F);
				FREQ_SEQ freq = t_gsys::band2freq(obsdata.gsys(), gobs.band());
				if (obsdata.gsys() == GAL && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2)) { return -alfa * beta; }
				else return 0.0;
			}
			else return 0.0;
		}
		case par_type::GRD_N:
			if (obsdata.site() == par.site) {
				xyz2ell(_trs_rec_crd, ell, false);
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				if (_mf_grd == GRDMPFUNC::CHEN_HERRING) {
					double sinel = sin(obsdata.ele());
					double tanel = tan(obsdata.ele());
					double cosaz = cos(obsdata.azi());
					return (1.0 / (sinel * tanel + 0.0032)) * cosaz;
				}
				else if (_mf_grd == GRDMPFUNC::TILTING) {
					double cosaz = cos(obsdata.azi());
					return dmfw * cosaz;
				}
				else if (_mf_grd == GRDMPFUNC::BAR_SEVER) {
					double tanel = tan(obsdata.ele());
					double cosaz = cos(obsdata.azi());
					return mfw * (1.0 / tanel) * cosaz;
				}
				else cerr << "Grad N mapping function is not set up correctly!!!" << endl;
			}
			else return 0.0;
		case par_type::GRD_E:
			if (obsdata.site() == par.site) {
				_getmf(par, obsdata, ell, epoch, mfw, mfh, dmfw, dmfh);
				if (_mf_grd == GRDMPFUNC::CHEN_HERRING) {
					double sinel = sin(obsdata.ele());
					double tanel = tan(obsdata.ele());
					double sinaz = sin(obsdata.azi());
					return (1.0 / (sinel * tanel + 0.0032)) * sinaz;
				}
				else if (_mf_grd == GRDMPFUNC::TILTING) {
					double sinaz = sin(obsdata.azi());
					return dmfw * sinaz;
				}
				else if (_mf_grd == GRDMPFUNC::BAR_SEVER) {
					double tanel = tan(obsdata.ele());
					double sinaz = sin(obsdata.azi());
					return mfw * (1 / tanel) * sinaz;
				}
				else cerr << "Grad E mapping function is not set up correctly!!!" << endl;
			}
			else return 0.0;

		case par_type::GLO_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == GLO/*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::GLO_IFCB:
			if (obsdata.site() == par.site && !gobs.is_phase() && obsdata.gsys() == GLO && par.prn == obsdata.sat()) return 1.0;
			else return 0.0;
		case par_type::GLO_IFPB:
			if (obsdata.site() == par.site && gobs.is_phase() && obsdata.gsys() == GLO && par.prn == obsdata.sat()) return 1.0;
			else return 0.0;
		case par_type::GLO_IFB:
			if (obsdata.site() == par.site && obsdata.gsys() == GLO && par.channel == obsdata.channel() /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::GAL_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::BDS_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && (!_bds2_isb || obsdata.sat() > "C16") /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::BD2_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && (!_bds2_isb || obsdata.sat() < "C17") /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::QZS_ISB:
			if (obsdata.site() == par.site && obsdata.gsys() == QZS /*&& gobs.is_code()*/) return 1.0;
			else return 0.0;
		case par_type::IFB_QZS:
			if (obsdata.site() == par.site && obsdata.gsys() == QZS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case par_type::IFB_GPS:
			if (obsdata.site() == par.site && obsdata.gsys() == GPS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case par_type::IFB_GAL:
			if (obsdata.site() == par.site && obsdata.gsys() == GAL && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		case par_type::IFB_BDS:
			if (obsdata.site() == par.site && obsdata.gsys() == BDS && _freq_index[obsdata.gsys()][gobs.band()] == FREQ_3 && gobs.is_code()) return 1.0;
			else return 0.0;
		}
		
		return 0.0;
	}

	void t_gprecisemodel::_getmf(t_gpar& par, t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh)
	{
		if (par.parType != par_type::TRP && par.parType != par_type::GRD_N && par.parType != par_type::GRD_E) return;

		double ele = satData.ele();

		if (_mf_ztd == ZTDMPFUNC::COSZ)
		{
			mfw = mfh = 1.0 / sin(ele);
		}
		else if (_mf_ztd == ZTDMPFUNC::GMF)
		{
			t_gmf mf;
			mf.gmf(epoch.mjd(), crd[0], crd[1], crd[2], G_PI / 2.0 - ele, mfh, mfw, dmfh, dmfw);
		}
		else cerr << "ZTD mapping function is not set up correctly!!!" << endl;
	}

	Matrix t_gprecisemodel::_RotMatrix_Ant(t_gsatdata& obsdata, const t_gtime& epoch, shared_ptr<t_gobj> obj, bool isCRS)
	{
		t_gdata::ID_TYPE type = obj->id_type();
		Matrix rotmatrix(3, 3);
		t_gtriple ell(0.0, 0.0, 0.0);
		double sinPhi, cosPhi, sinLam, cosLam;
		if (type == t_gdata::TRN) {
			string antenna = (obj->pcv(_crt_epo))->anten();
			t_gtriple dx(0, 0, 0);
			ColumnVector i(3);
			ColumnVector j(3);
			ColumnVector k(3);
			if (_attitudes == ATTITUDES::YAW_NOMI)	   attitude(obsdata, "", i, j, k);
			else if (_attitudes == ATTITUDES::YAW_RTCM) attitude(obsdata, obsdata.yaw(), i, j, k);
			else							   attitude(obsdata, antenna, i, j, k);
			rotmatrix.Column(1) = i;
			rotmatrix.Column(2) = j;
			rotmatrix.Column(3) = k;
		}
		else if (type == t_gdata::REC) {
			xyz2ell(_trs_rec_crd, ell, false);
			sinPhi = sin(ell[0]);
			cosPhi = cos(ell[0]);
			sinLam = sin(ell[1]);
			cosLam = cos(ell[1]);
			rotmatrix << -sinPhi * cosLam << -sinLam << +cosPhi * cosLam
				<< -sinPhi * sinLam << +cosLam << +cosPhi * sinLam
				<< +cosPhi << 0.0 << +sinPhi;
		}
		else {
			throw exception();
		}

		if (isCRS) {
			_update_rot_matrix(epoch);
			rotmatrix = _trs2crs_2000->getRotMat()*rotmatrix;
		}

		return rotmatrix;
	}


	bool t_gprecisemodel::_apply_rec_tides(const t_gtime& epoch, t_gtriple& rec)
	{
		if (!_tide)
		{
			_log->comment(0, "ERROR", "The tide ptr is nullptr.");
			return false;
		}

		_update_rot_matrix(epoch);
		double xpole = _trs2crs_2000->getXpole();
		double ypole = _trs2crs_2000->getYpole();
		double gast = _trs2crs_2000->getGmst();
		Matrix rot_trs2crs = _trs2crs_2000->getRotMat();

		t_gtriple tide(0.0, 0.0, 0.0);
		try
		{
			// solid tide
			t_gtideIERS* solid_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!solid_ptr)
			{
				_log->comment(0, "ERROR", "can not get the solid tide ptr.");
				return false;
			}
			t_gtriple solid_earth = solid_ptr->tide_solid(epoch, rec, rot_trs2crs, _gdata_navde);

			// ocean load
			t_gtriple load_ocean = _tide->load_ocean(epoch, _site, rec);

			//pole tide
			t_gtideIERS* pole_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!pole_ptr)
			{
				_log->comment(0, "ERROR", "can not get the pole tide ptr.");
				return false;
			}
			t_gtriple tide_pole = pole_ptr->tide_pole_pod(epoch, xpole, ypole, rec);

			//atmosph
			t_gtriple load_atmosph = _tide->load_atmosph();

			// freq tide
			t_gtideIERS* freq_ptr = dynamic_cast<t_gtideIERS*>(_tide.get());
			if (!freq_ptr)
			{
				_log->comment(0, "ERROR", "can not get the freq tide ptr.");
				return false;
			}
			t_gtriple tide_freq = freq_ptr->tide_freq(_site, rec, gast);

			tide = solid_earth +
				load_ocean +
				tide_pole +
				load_atmosph +
				tide_freq;

		}
		catch (...)
		{
			_log->comment(0, "ERROR", "can not get the tide ptr.");
			return false;
		}

		//unit to m
		rec = rec + tide * 1.e3;
		return true;

	}

	bool t_gprecisemodel::_check_obs_valid(t_gobsgnss& obsdata, const OBSCOMBIN& obscom, const GOBSBAND& b1, const GOBSBAND& b2)
	{
		if (!_gall_nav || !_gallobj || !_gallbias)
		{
			write_log_info(_log, 1, "NOTE", "no nav or obj or bias");
			return false;
		}

		GSYS gs = obsdata.gsys();
		t_gobs gobs1;
		t_gobs gobs2;

		double ObsP = 0.0, ObsL = 0.0;
		switch (obscom)
		{
		case OBSCOMBIN::IONO_FREE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			gobs2 = t_gobs(obsdata.select_range(b2));
			ObsP = obsdata.P3(gobs1, gobs2);
			gobs1 = t_gobs(obsdata.select_phase(b1));
			gobs2 = t_gobs(obsdata.select_phase(b2));
			ObsL = obsdata.L3(gobs1, gobs2);
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case OBSCOMBIN::RAW_SINGLE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			ObsP = obsdata.getobs(gobs1.gobs());
			gobs1 = t_gobs(obsdata.select_phase(b1));
			ObsL = obsdata.getobs(gobs1.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case OBSCOMBIN::RAW_DOUBLE:
		{
			gobs1 = t_gobs(obsdata.select_range(b1));
			gobs2 = t_gobs(obsdata.select_phase(b1));
			ObsP = obsdata.getobs(gobs1.gobs());
			ObsL = obsdata.getobs(gobs2.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			gobs1 = t_gobs(obsdata.select_range(b2));
			gobs2 = t_gobs(obsdata.select_phase(b2));
			ObsP = obsdata.getobs(gobs1.gobs());
			ObsL = obsdata.getobs(gobs2.gobs());
			if (double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0))
			{
				write_log_info(_log, 1, "NOTE", "no range or phase data");
				return false;
			}
			break;
		}
		case OBSCOMBIN::RAW_ALL:
		{
			std::map<FREQ_SEQ, GOBSBAND>  sys_bands = _band_index[gs];
			for (const auto& band : sys_bands)
			{
				FREQ_SEQ freqX = band.first;
				GOBSBAND bandX = band.second;
				gobs1 = t_gobs(obsdata.select_range(bandX, true));
				gobs2 = t_gobs(obsdata.select_phase(bandX, true));
				ObsP = 0.0;
				ObsL = 0.0;
				ObsP = obsdata.getobs(gobs1.gobs());
				ObsL = obsdata.getobs(gobs2.gobs());
				if ((double_eq(ObsP, 0.0) || double_eq(ObsL, 0.0)))
				{
					write_log_info(_log, 1, "NOTE", "Missing observations");
					return false;
				}
			}
			break;
		}
		default:
		{
			write_log_info(_log, 1, "NOTE", "no obs mode valid");
			return false;
		}
		}

		return true;
	}


	bool apply_reldelay(t_gtriple crd_site, t_gtriple vel_site, t_gtriple crd_sat, t_gtriple vel_sat,
		double& reldelay)
	{
		reldelay = 2.0*(DotProduct(crd_sat.crd_cvect(), vel_sat.crd_cvect()))/ CLIGHT;
		ColumnVector xsat = crd_sat.crd_cvect();
		ColumnVector xsite = crd_site.crd_cvect();

		double r = xsite.norm_Frobenius() + xsat.norm_Frobenius();
		ColumnVector xsat2site = (xsite - xsat);
		double r_site2sat = xsat2site.norm_Frobenius();

		reldelay += 2.0*GM_CGCS / CLIGHT / CLIGHT * log((r + r_site2sat) / (r - r_site2sat));

		return true;

	}

	bool correct_DCB(string ac, t_gallbias* allbias, t_gsatdata& satdata, double& P, OBSCOMBIN obscombin, t_gobs* gobs1, t_gobs* gobs2, t_gobs* gobs3)
	{
		if (!allbias || !gobs1 || obscombin != OBSCOMBIN::RAW_ALL && !gobs2)
		{
			return false;
		}

		GOBSBAND band1, band2, band3;
		t_gobs fobs1, fobs2, fobs3;            // local gobs
		double wavelen1, wavelen2, wavelen3;
		double coef1, coef2;                   // coefficient

		// constant
		t_gtime epoch = satdata.epoch();
		string sat = satdata.sat();
		GSYS gsys = satdata.gsys();

		band1 = t_gsys::band_priority(gsys, FREQ_1);
		band2 = t_gsys::band_priority(gsys, FREQ_2);
		band3 = t_gsys::band_priority(gsys, FREQ_3);

		// Judge gobs
		if (gobs1->attr() == ATTR) fobs1 = satdata.id_range(gobs1->band());      // automatic GOBS selection
		else fobs1 = gobs1->gobs();                                              // specific GOBS
		fobs1.gobs2to3(gsys);                                                    // 2.xx to 3.xx

		if (gobs2 != 0) {
			if (gobs2->attr() == ATTR) fobs2 = satdata.id_range(gobs2->band());  // automatic GOBS selection
			else fobs2 = gobs2->gobs();                                          // specific GOBS
			fobs2.gobs2to3(gsys);                                                // 2.xx to 3.xx
		}
		if (gobs3 != 0) {
			if (gobs3->attr() == ATTR) fobs3 = satdata.id_range(gobs3->band());  // automatic GOBS selection
			else fobs3 = gobs3->gobs();                                          // specific GOBS
			fobs3.gobs2to3(gsys);                                                // 2.xx to 3.xx
		}

		// coefficients
		wavelen1 = satdata.wavelength(band1);
		wavelen2 = satdata.wavelength(band2);
		coef2 = wavelen1 / wavelen2;
		coef1 = 1.0 / (1.0 - coef2 * coef2);
		coef2 = coef2 * coef2 * coef1;

		if (obscombin == OBSCOMBIN::RAW_ALL)
		{
			FREQ_SEQ freq = t_gsys::band2freq(gsys, fobs1.band());
			if (freq == FREQ_1) {
				fobs2 = t_gobs(satdata.select_range(band2));
				fobs2.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				P -= -coef2 * corr;
			}
			else if (freq == FREQ_2) {
				fobs2 = fobs1;
				fobs1 = t_gobs(satdata.select_range(band1));
				fobs1.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				P -= -coef1 * corr;
			}
			else if (freq == FREQ_3) {
				fobs3 = fobs1;
				fobs1 = t_gobs(satdata.select_range(band1));
				fobs2 = t_gobs(satdata.select_range(band2));
				fobs1.gobs2to3(gsys);   // 2.xx to 3.xx
				fobs2.gobs2to3(gsys);   // 2.xx to 3.xx
				double corr1 = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				double corr2 = allbias->get(epoch, sat, fobs1.gobs(), fobs3.gobs(), ac);
				P -= -coef1 * corr2 + coef2 * (corr2 - corr1);
			}
		}
		else
		{
			// Triple IF DCB correction
			if (gobs3 != 0) {
				if (obscombin == OBSCOMBIN::IONO_FREE && gsys != GLO) {
					double corr1 = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
					double corr2 = allbias->get(epoch, sat, fobs1.gobs(), fobs3.gobs(), ac);
					wavelen3 = satdata.wavelength(band3);
					double coef3, coef4;
					coef4 = wavelen1 / wavelen3;
					coef3 = 1.0 / (1.0 - coef4 * coef4);
					coef4 = coef4 * coef4 * coef3;
					P -= -coef2 * corr1 + coef4 * corr2;
				}
			}
			else
			{
				double corr = allbias->get(epoch, sat, fobs1.gobs(), fobs2.gobs(), ac);
				double corr1 = coef2 * corr;
				double corr2 = coef1 * corr;

				if (obscombin == OBSCOMBIN::MW_COMBIN)
				{
					double fact = wavelen1 / wavelen2;
					P += -(corr1 / wavelen1 + corr2 / wavelen2) * (1.0 - fact) / (1.0 + fact);
				}
				else if (obscombin == OBSCOMBIN::EWL_COMBIN)
				{
					wavelen3 = satdata.wavelength(band3);
					double alpha = (wavelen2 * wavelen3 - wavelen1 * wavelen1) / (wavelen2 * wavelen2 - wavelen1 * wavelen1);
					P += -(corr1 * (1 - alpha) + corr2 * alpha) * (wavelen3 - wavelen2) / (wavelen2 * wavelen3);
				}
				else if (obscombin == OBSCOMBIN::IONO_FREE)
				{
					double alfa, beta;
					satdata.coef_ionofree(band1, alfa, band2, beta);
					P += alfa * corr1 + beta * corr2;
				}

			}
		}
		return true;
	}


	int Check_Range(vector<double>& omc, double& recclk, double& sigma)
	{
		double mean = 0.0;
		double sig = 0.0;
		vector<double> res;

		for (auto it = omc.begin(); it != omc.end();) {
			if (fabs(*it) < 10e-8) {
				it = omc.erase(it);
				continue;
			}
			mean += *it;
			it++;
		}

		if (omc.size() == 1 || (omc.size() == 2 && fabs(omc[0] - omc[1]) > 1000)) {
			recclk = 0.0;
			return 1;
		}

		mean /= omc.size();

		for (auto it = omc.begin(); it != omc.end(); it++) {
			sig += (*it - mean) * (*it - mean);
			res.push_back(fabs(*it - mean));
		}
		sig = sqrt(sig / (omc.size() - 1));

		double sig1 = 2000;
		double mean1 = 0.0;
		bool isFound = true;

		while (isFound && sig1 > 1000 && omc.size() >= 3)
		{
			auto maxRes = max_element(res.begin(), res.end());
			auto maxPos = omc.begin() + distance(res.begin(), maxRes);
			double MaxValue = *maxPos;

			mean1 = 0.0;
			sig1 = 0.0;
			res.clear();
			for (auto it = omc.begin(); it != omc.end(); it++) {
				if (it == maxPos) continue;
				mean1 += *it / (omc.size() - 1);
			}

			for (auto it = omc.begin(); it != omc.end(); it++) {
				if (it == maxPos) continue;
				sig1 += (*it - mean1) * (*it - mean1);
				res.push_back(fabs(*it - mean1));
			}
			sig1 = sqrt(sig1 / (omc.size() - 2));

			isFound = fabs(MaxValue - mean1) > 3 * (sig1 > 1000 / 3.0 ? sig1 : 1000 / 3.0);
			if (!isFound) continue;
			mean = mean1;
			sig = sig1;
			omc.erase(maxPos);
		}

		recclk = mean;
		sigma = sig;
		return omc.size();
	}
}
