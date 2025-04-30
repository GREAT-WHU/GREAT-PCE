/**
 * @file         gupdateparIF.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gproc/gupdateparIF.h"
#include "gutils/gturboedit.h"

namespace great
{
	void t_gupdateparIF::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info)
	{
		if (_lite_update_amb)
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, par_type::AMB_IF, update_info);
		}
		else
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, par_type::AMB_IF, _cycleslip.get(), update_info);
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,par_type ambtype, t_gcycleslip* _slip, t_gupdateparinfo& update_info)
	{

		map<string, set<string> > mapPrn;
		int idx_max = allpars.maxIndex();
		// int num = 0;
		for (auto& obs_data : obsdata)
		{
			if (ambtype == par_type::AMB_IF   && !obs_data.tb12()) continue;

			string rec = obs_data.site();
			string sat = obs_data.sat();
			mapPrn[rec].insert(sat);

			int idx_amb = allpars.getParam(rec, ambtype, sat);

			t_gtime crt_epoch = epoch;
			// first add new amb
			// second update old amb if cycleslip
			if (idx_amb < 0 && _slip->num_of_amb_arc(rec, sat, crt_epoch) ||
				idx_amb >= 0 && _slip->cycle_slip(obs_data, crt_epoch))
			{
				t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1 + _new_par, sat);

				int amb_flag = _slip->num_of_amb_arc(rec, sat, epoch);
				_slip->set_amb_flag(rec, sat, amb_flag);
				auto tb_slip = dynamic_cast<t_gturboedit*>(_slip);

				t_gtime end_time = LAST_TIME;
				if (tb_slip) {
					end_time = tb_slip->get_crt_amb_end(rec, sat);
				}
				_new_par++;
				par_newamb.value(0.0);
				par_newamb.setTime(epoch, end_time);
				par_newamb.apriori(10000.0);

				// add new amb
				if (idx_amb < 0)
				{
					update_info.add(par_newamb);
				
				}
				// reset the new amb
				else
				{
					if (allpars[idx_amb].end > epoch - _intv) allpars[idx_amb].end = epoch - _intv;
					update_info.add(idx_amb + 1, par_newamb);
				}
			}
		}

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
			{
				if (allpars[i].end < epoch)
				{
					if (!update_info.exist(i + 1))
					{
						if (allpars[i].end > epoch - _intv) allpars[i].end = epoch - _intv;
						update_info.add(i + 1);
					}
				}
				else
				{
					continue;
				}
			}
		}

	}

	void t_gupdateparIF::_update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, par_type ambtype, t_gupdateparinfo& update_info)
	{		
		int idx_max = allpars.maxIndex();
		for (auto obs_data : obsdata)
		{
			string rec = obs_data.site();
			string sat = obs_data.sat();
			GSYS   sys = obs_data.gsys();

			int idx_amb = allpars.getParam(rec, ambtype, sat);
			GOBS obsL1, obsL2;
			if (ambtype == par_type::AMB_IF)
			{
				obsL1 = obs_data.select_phase(_band_index[sys][FREQ_1]);
				obsL2 = obs_data.select_phase(_band_index[sys][FREQ_2]);
			}
			else
			{
				return;
			}

			if (idx_amb < 0)
			{
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1 + _new_par, sat);

					par_newamb.value(0.0);
					par_newamb.setTime(epoch, epoch);
					par_newamb.apriori(9000.0);
					update_info.add(par_newamb);
					_new_par++;
				}
			}
			else
			{
				if (obsL1 != GOBS::X && obsL2 != GOBS::X)
				{
					if (obs_data.getlli(obsL1) < 1 && obs_data.getlli(obsL2) < 1)
					{
						allpars[idx_amb].end = epoch;
					}
					else
					{
						allpars[idx_amb].end = epoch - _intv;

						int idx_max = allpars[allpars.parNumber() - 1].index;
						t_gpar par_newamb = t_gpar(rec, ambtype, idx_max + 1, sat);

						par_newamb.value(0.0);
						par_newamb.setTime(epoch, epoch);
						par_newamb.apriori(9000.0);

						update_info.add(idx_amb + 1, par_newamb);
					}
				}
				else  // obsL1==X || obsL2==X
				{
					update_info.add(idx_amb + 1);
				}
			}
		}

		// delete old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == ambtype)
			{
				if (allpars[i].end < epoch - 600)
				{
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1);
					}
				}
				else if (allpars[i].amb_ini == true && allpars[i].beg + -1 * _intv < epoch)
				{
					t_gpar newamb = allpars[i];
					newamb.beg = epoch;
					if (!update_info.exist(i + 1))
					{
						update_info.add(i + 1, newamb);
					}
				}
			}
		}

	}

	t_gupdateparIF_1X::t_gupdateparIF_1X(shared_ptr<t_gcycleslip> cycleslip12, shared_ptr<t_gcycleslip> cycleslip13, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index):
		t_gupdateparIF(cycleslip12, band_index)
	{
		_cycleslip13 = cycleslip13;
	}

	t_gupdateparIF_1X::~t_gupdateparIF_1X()
	{
	}

	void t_gupdateparIF_1X::set_cycleslip13(shared_ptr<t_gcycleslip> cycleslip13)
	{
		_cycleslip13 = cycleslip13;
	}

	void t_gupdateparIF_1X::set_cycleslip14(shared_ptr<t_gcycleslip> cycleslip14)
	{
		_cycleslip14 = cycleslip14;
	}

	void t_gupdateparIF_1X::set_cycleslip15(shared_ptr<t_gcycleslip> cycleslip15)
	{
		_cycleslip15 = cycleslip15;
	}

	t_gupdateparinfo t_gupdateparIF_1X::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_info;
		// update AMB
		t_gupdateparIF_1X::_update_amb_pars(epoch, allpars, obsdata, update_info);

		// update gps rec ifb pars
		if (allpars.orbParNumber() > 0)
		{
			_udpate_gps_rec_ifb_pars(epoch, allpars, obsdata, update_info);
		}

		// update Other pars
		_update_process_pars(epoch, allpars, update_info);

		return update_info;
	}

	void t_gupdateparIF_1X::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info)
	{
		t_gupdateparIF::_update_amb_pars(epoch, allpars, obsdata, update_info);
		if (_lite_update_amb)
		{
			_update_one_type_amb_pars(epoch, allpars, obsdata, par_type::AMB_IF, update_info);
		}
	}

	t_updateparWL::t_updateparWL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index):
		t_gupdateparIF(cycleslip, band_index)
	{
	}
	t_updateparWL::~t_updateparWL()
	{
	}
	void t_updateparWL::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		_update_one_type_amb_pars(epoch, allpars, obsdata, par_type::AMB_WL, _cycleslip.get(), update_info);
	}
}