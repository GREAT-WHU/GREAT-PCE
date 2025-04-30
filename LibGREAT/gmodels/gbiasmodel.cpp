/**
 * @file         gbiasmodel.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        mainly about how to cacultae B P l single
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/gbiasmodel.h"
#include "gutils/ginfolog.h"


namespace great
{
    t_gbiasmodel::t_gbiasmodel()
    {

    }

    t_gbiasmodel::~t_gbiasmodel()
    {

    }

    t_gprecisebias::t_gprecisebias(t_gallproc *data, t_glog *log, t_gsetbase *setting):
		_multi_thread_flag(false),
		gmodel(data,log,setting),
        _log(log)
    {

    }

    t_gprecisebias::~t_gprecisebias()
    {

    }

	void t_gprecisebias::set_multi_thread(const set<string>& sites)
	{
		_multi_thread_flag = true;
		for (auto rec_item : sites)
		{
			_map_site_model[rec_item] = shared_ptr<t_gprecisemodel>(new t_gprecisemodel(gmodel));
			_map_flag[rec_item] = make_tuple("", "", t_gtime());
		}
	}

    bool t_gprecisebias::cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result)
    {
        // check obs_type valid
        double Obs_value = obsdata.getobs(gobs.gobs());
        double snr_value = obsdata.getobs(pl2snr( gobs.gobs())  );
        if (double_eq(Obs_value,0.0))
        {
            write_log_info(_log, 1, "NOTE", "Obs_value is 0.0");
            return false;
        }

		t_gprecisemodel* gmodel_ptr = &gmodel;
		tuple<string, string, t_gtime>* flag = &_rec_sat_before;
		if (_multi_thread_flag) {
			gmodel_ptr = _map_site_model.at(obsdata.site()).get();
			flag = &_map_flag.at(obsdata.site());
		}

        // skip the same obsdata caculate common bias
        if (make_tuple(obsdata.site(), obsdata.sat(),epoch) != *flag)
        {
            bool update_valid = gmodel_ptr->_update_obs_info(epoch, obsdata, params);
            if (!update_valid)
            {
                write_log_info(_log, 1, "NOTE", "update obs information failed" + epoch.str_ymdhms("", false));
                return false;
            }

            // prepare Caculate some common bias[sat_pos,rec_pos,relative,rho]
            bool pre_valid = gmodel_ptr->_prepare_obs(epoch, params);
            if (!pre_valid)
            {
                write_log_info(_log, 1, "NOTE", "prepare obs information failed" + epoch.str_ymdhms("", false));
                return false;
            }

            *flag = make_tuple(obsdata.site(), obsdata.sat(),epoch);
        }

        // combine equ
        t_glsqEquationMatrix equ;
        double omc=0.0, wgt=0.0;
        vector<pair<int, double> > coef;

        if (!gmodel_ptr->_omc_obs_ALL(epoch,params,gobs , omc))
        { 
			write_log_info(_log, 1, "NOTE", "omc obs failed"); 
			return false; 
		};

		if (!gmodel_ptr->_wgt_obs_ALL(gobs,gmodel_ptr->_factorP, wgt))
		{ 
			write_log_info(_log, 1, "NOTE", "weight obs failed"); 
			return false; 
		};

        if (!gmodel_ptr->_prt_obs(epoch, params, gobs, coef))
		{ 
			write_log_info(_log, 1, "NOTE", "partialrange obs failed"); 
			return false; 
	    }

        result.B.push_back(coef);
        result.P.push_back(wgt);
        result.l.push_back(omc);

        gmodel_ptr->_update_obs_info(obsdata);
        return true;
    }
    
    void t_gprecisebias::update_obj_clk(const string& obj, const t_gtime& epo, double clk)
    {
		t_gprecisemodel* gmodel_ptr = &gmodel;

		if (_multi_thread_flag) {
			gmodel_ptr = _map_site_model.at(obj).get();
		}

        gmodel_ptr->_obj_clk[obj].first = epo;
        gmodel_ptr->_obj_clk[obj].second = clk;
        gmodel_ptr->_rec_clk[obj] = clk;
    }
}