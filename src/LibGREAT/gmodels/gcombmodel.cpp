/**
 * @file         gcombmodel.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        base combine biase model
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/gcombmodel.h"
#include "gutils/gstring.h"
#include <cmath>
#include <utility>


namespace great
{
	t_gcombmodel::t_gcombmodel(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		_bias_model(std::move(bias_model)),
		_glog(log)
	{
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

		// =======================================================================================
		// sigma
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

		log->logInfo("t_gcombmodel", "sigCodeGPS ", format("%16.8f", _sigCodeGPS));
		log->logInfo("t_gcombmodel", "sigCodeGLO ", format("%16.8f", _sigCodeGLO));
		log->logInfo("t_gcombmodel", "sigCodeGAL ", format("%16.8f", _sigCodeGAL));
		log->logInfo("t_gcombmodel", "sigCodeBDS ", format("%16.8f", _sigCodeBDS));
		log->logInfo("t_gcombmodel", "sigCodeQZS ", format("%16.8f", _sigCodeQZS));
		log->logInfo("t_gcombmodel", "sigPhaseGPS", format("%16.8f", _sigPhaseGPS));
		log->logInfo("t_gcombmodel", "sigPhaseGLO", format("%16.8f", _sigPhaseGLO));
		log->logInfo("t_gcombmodel", "sigPhaseGAL", format("%16.8f", _sigPhaseGAL));
		log->logInfo("t_gcombmodel", "sigPhaseBDS", format("%16.8f", _sigPhaseBDS));
		log->logInfo("t_gcombmodel", "sigPhaseQZS", format("%16.8f", _sigPhaseQZS));

		_frequency = dynamic_cast<t_gsetproc*>(setting)->frequency();

		_gallbias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gprecisemodel = dynamic_cast<t_gprecisemodel*>(_bias_model->precisemodel());

		_update_amb_lite = dynamic_cast<t_gsetturboedit*>(setting)->liteMode();

	}

	t_gcombmodel::~t_gcombmodel() = default;


	map<GSYS, map<FREQ_SEQ, GOBSBAND>> t_gcombmodel::get_band_index()
	{
		return _band_index;
	}

	map<GSYS, map<GOBSBAND, FREQ_SEQ>> t_gcombmodel::get_freq_index()
	{
		return _freq_index;
	}

	bool t_gcombmodel::_wgt_raw_obs(const t_gobs& gobs, const t_gsatdata& satdata, const double& factorP, const OBSWEIGHT& wgt_type, double& wgt)
	{
		auto obs_type = gobs.type();
		if (obs_type != TYPE_C && obs_type != TYPE_P && obs_type != TYPE_L)
		{
			_glog->logError("t_gcombmodel", "_wgt_raw_obs", " type should be TYPE_C/TYPE_P/TYPE_L");
			return false;
		}

		double factor = factorP;

		auto gsat = satdata.sat();
		auto gsys = satdata.gsys();
		if (gsys == GSYS::BDS) factor = 2.0;
		if (gsys == GSYS::BDS && t_gsys::bds_geo(gsat)) factor = 5.0;

		auto obj_type = satdata.id_type();
		double sigRange = 0.0;
		if (obs_type == TYPE_C || obs_type == TYPE_P)
		{
			if (obj_type == t_gdata::REC)
			{
				switch (gsys)
				{

				case GPS: sigRange = _sigCodeGPS; break;
				case GLO: sigRange = _sigCodeGLO; break;
				case GAL: sigRange = _sigCodeGAL; break;
				case BDS: sigRange = _sigCodeBDS; break;
				case QZS: sigRange = _sigCodeQZS; break;
				default:
					_glog->logError("t_gcombmodel", "_wgt_raw_obs", " sys should be GPS/GAL/GLO/BDS/QZS");
					return false;
				}
			}
			else
			{
				_glog->logError("t_gcombmodel", "_wgt_raw_obs", " obj_type should be REC/REC_LEO");
				return false;
			}
		}

		// case of Phase
		double sigPhase = 0.0;
		if (obs_type == TYPE_L)
		{
			if (obj_type == t_gdata::REC)
			{
				switch (gsys)
				{
				case GPS: sigPhase = _sigPhaseGPS; break;
				case GLO: sigPhase = _sigPhaseGLO; break;
				case GAL: sigPhase = _sigPhaseGAL; break;
				case BDS: sigPhase = _sigPhaseBDS; break;
				case QZS: sigPhase = _sigPhaseQZS; break;
				default:
					_glog->logError("t_gcombmodel", "_wgt_raw_obs", " sys should be GPS/GAL/GLO/BDS/QZS");
					return false;
				}
			}
			else
			{
				_glog->logError("t_gcombmodel", "_wgt_raw_obs", " obj_type should be REC/REC_LEO");
				return false;
			}
		}
		else
		{
			return	false;
		}

		if (obj_type == t_gdata::REC)
		{
			double leo_ele = satdata.ele_leo();
			double leo_ele_deg = satdata.ele_leo_deg();
			double sin_leo_ele = sin(leo_ele);
			// get the _weight factor
			switch (wgt_type)
			{
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
			case OBSWEIGHT::DEF_OBSWEIGHT:
			default:
				_glog->logError("t_gcombmodel", "_wgt_raw_obs", " WeightObs (default) should not happened!\n");
				return false;
			}
		}

		// get the combine obs factor
		switch (obs_type)
		{
		case TYPE_L:
			sigPhase = sigPhase * satdata.wavelength(gobs.band());
			if (sigPhase == 0.0) return false;
			wgt = 1.0 / pow(factor * sigPhase, 2);
			break;
		case TYPE_C:
		case TYPE_P:
			sigRange = sigRange;
			if (sigRange == 0.0) return false;
			wgt = 1.0 / pow(factor * sigRange, 2);
			break;
		default:
			return false;
		}
		return true;
	}

	t_gcombIF::t_gcombIF(t_gsetbase* setting, t_glog* log, shared_ptr<t_gbiasmodel> bias_model, t_gallproc* data) :
		t_gcombmodel(setting, log, std::move(bias_model), data)
	{

	}

	t_gcombIF::~t_gcombIF() = default;

	bool t_gcombIF::cmb_equ(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, t_gbaseEquation& result)
	{
		// ========================================================================================================================================
		// check Obs type
		GOBSBAND b1 = _band_index[obsdata.gsys()][FREQ_1];
		GOBSBAND b2 = _band_index[obsdata.gsys()][FREQ_2];

		obsdata.apply_bias(*_gallbias);

		return this->cmb_equ_IF(epoch, params, obsdata, b1, b2, result);
	}

	bool t_gcombIF::cmb_equ_IF(t_gtime& epoch, t_gallpar& params, t_gsatdata& obsdata, GOBSBAND b1, GOBSBAND b2, t_gbaseEquation& result)
	{
		// ========================================================================================================
		if (b1 == BAND || b2 == BAND) return false;

		// ========================================================================================================
		t_gobs gobsP1 = obsdata.select_range(b1);
		t_gobs gobsP2 = obsdata.select_range(b2);
		t_gobs gobsL1 = obsdata.select_phase(b1);
		t_gobs gobsL2 = obsdata.select_phase(b2);

		vector<pair<t_gobs, t_gobs> > type_list;
		type_list.emplace_back(gobsP1, gobsP2);
		type_list.emplace_back(gobsL1, gobsL2);

		auto gsys = obsdata.gsys();
		// ========================================================================================================
		// IF coef
		double coef1, coef2;
		obsdata.coef_ionofree(b1, coef1, b2, coef2);

		// ========================================================================================================
		t_glsqEquationMatrix equ_IF;
		for (const auto &item : type_list)
		{
			t_gobs gobs1 = item.first;
			t_gobs gobs2 = item.second;

			//combine f1 and f2
			t_gbaseEquation temp_equ;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs1, temp_equ)) return false;
			if (!_bias_model->cmb_equ(epoch, params, obsdata, gobs2, temp_equ)) return false;

			// ========================================================================================================
			if (temp_equ.B[0].size() != temp_equ.B[1].size()) throw logic_error("coeff size is not equal in f1 and f2");

			vector<pair<int, double> > coef_IF;
			double P_IF = 0.0, l_IF = 0.0;
			for (int i = 0; i < temp_equ.B[0].size(); i++)
			{
				if (temp_equ.B[0][i].first != temp_equ.B[1][i].first) 
				{
					throw logic_error("coeff par is not the same in f1 and f2");
				}
				// combine coeff
				coef_IF.emplace_back(temp_equ.B[0][i].first, coef1 * temp_equ.B[0][i].second + coef2 * temp_equ.B[1][i].second);
			}
			// combine P and l
			P_IF = 1.0 / (pow(coef1, 2) / temp_equ.P[0] + pow(coef2, 2) / temp_equ.P[1]);
			l_IF = coef1 * temp_equ.l[0] + coef2 * temp_equ.l[1];

			// ========================================================================================================
			// CORRECT AMB
			if (gobs1.is_phase())
			{
				par_type amb_type = par_type::NO_DEF;
				if (_freq_index[gsys][b2] == FREQ_2)  amb_type = par_type::AMB_IF;

				int idx = params.getParam(obsdata.site(), amb_type, obsdata.sat());
				if (idx < 0) return false;

				// if amb is new  then init value with Obs - Modelobs				
				if (double_eq(params[idx].value(), 0.0) || params[idx].beg == epoch)
				{
					double obs_P3 = obsdata.P3(gobsP1, gobsP2);
					double obs_L3 = obsdata.L3(gobsL1, gobsL2);

					params[idx].value(obs_L3 - obs_P3);
				}
				// update B l
				coef_IF.emplace_back(idx + 1, 1.0);
				l_IF -= params[idx].value();
			}

			t_gobscombtype type(gobs1, b1, b2, _freq_index[gsys][b1], _freq_index[gsys][b2], OBSCOMBIN::IONO_FREE);
            equ_IF.add_equ(coef_IF, P_IF,l_IF,obsdata.site(),obsdata.sat(), type, false);
        }

        if (dynamic_cast<t_glsqEquationMatrix*>(&result))
        {
            auto *lsq_result = dynamic_cast<t_glsqEquationMatrix *>(&result);
            lsq_result->add_equ(equ_IF);
        }
        else
        {
            result = result + equ_IF;
        }
	
        return true;
    }
}