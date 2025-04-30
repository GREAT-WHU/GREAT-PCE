/**
 * @file         gupdateparALL.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gproc/gupdateparALL.h"
#include "gutils/gturboedit.h"

namespace great
{
	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> > & band_index, const IONO_ORDER& iono_order)
	{
		_cycleslip = cycleslip;
		_band_index = band_index;
		_iono_order = iono_order;
		_freq = 2;
	}

	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index, const IONO_ORDER& iono_order)
	{
		_cycleslip = cycleslip_12;
		_cycleslip_13 = cycleslip_13;
		_band_index = band_index;
		_iono_order = iono_order;
		_freq = 3;
	}
	
	t_gupdateparALL::t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, shared_ptr<t_gcycleslip> cycleslip_14, shared_ptr<t_gcycleslip> cycleslip_15, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index, const IONO_ORDER& iono_order)
	{
		_cycleslip = cycleslip_12;
		_cycleslip_13 = cycleslip_13;
		_cycleslip_14 = cycleslip_14;
		_cycleslip_15 = cycleslip_15;
		_band_index = band_index;
		_iono_order = iono_order;

		if (_cycleslip) _freq = 2;
		if (_cycleslip && _cycleslip_13) _freq = 3;
		if (_cycleslip && _cycleslip_13 && _cycleslip_14 && _cycleslip_15) _freq = 5;

	}

	t_gupdateparALL::~t_gupdateparALL()
	{
	}

	
	t_gupdateparinfo t_gupdateparALL::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_par_info;

		_new_par = 0;
		// update ISB
		update_isb_pars(epoch, allpars, obsdata, update_par_info);

		// update_ion_info
		_udpate_ion_pars(epoch, allpars, obsdata, update_par_info);

		// update gps rec ifb pars
		if(_cycleslip_13 != nullptr && allpars.orbParNumber() > 0) _udpate_gps_rec_ifb_pars(epoch, allpars, obsdata, update_par_info);

		// update_amb_info
		if (_lite_update_amb)
		{
			_lite_update_amb_pars(epoch, allpars, obsdata, update_par_info);
		}
		else
		{
			_update_amb_pars(epoch, allpars, obsdata, update_par_info);
		}
		

		// update_process_info
		_update_process_pars(epoch, allpars, update_par_info);

		return update_par_info;
	}
	
    void t_gupdateparALL::_update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		// add new amb
		auto tb_slip_12 = dynamic_cast<t_gturboedit*>(_cycleslip.get());
		auto tb_slip_13 = dynamic_cast<t_gturboedit*>(_cycleslip_13.get());
		auto tb_slip_14 = dynamic_cast<t_gturboedit*>(_cycleslip_14.get());
		auto tb_slip_15 = dynamic_cast<t_gturboedit*>(_cycleslip_15.get());

		map<string, set<string> > mapPrn;

		if (!tb_slip_12) throw runtime_error("tb is empty!");
		for (auto& data : obsdata)
		{
			string rec = data.site();
			string sat = data.sat();
			GSYS   sys = data.gsys();

			mapPrn[rec].insert(sat);
			bool isUpL1 = false; //update L1
			bool isUpL2 = false; //update L2
			bool isUpL3 = false; //update L3
			bool isUpL4 = false;
			bool isUpL5 = false;

			int idx_amb_L2 = allpars.getParam(rec, par_type::AMB_L2, sat);
			int idx_amb_L3 = allpars.getParam(rec, par_type::AMB_L3, sat);
			int idx_amb_L4 = allpars.getParam(rec, par_type::AMB_L4, sat);
			int idx_amb_L5 = allpars.getParam(rec, par_type::AMB_L5, sat);


			int Flag12 = 0, Flag13 = 0, Flag14 = 0, Flag15 = 0;
			Flag12 = tb_slip_12->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_13) Flag13 = tb_slip_13->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_14) Flag14 = tb_slip_14->num_of_amb_arc(rec, sat, epoch);
			if (tb_slip_15) Flag15 = tb_slip_15->num_of_amb_arc(rec, sat, epoch);

			t_gtime end_time    = LAST_TIME;
			t_gtime end_time_12 = LAST_TIME;
			t_gtime end_time_13 = LAST_TIME;
			t_gtime end_time_14 = LAST_TIME;
			t_gtime end_time_15 = LAST_TIME;

			if (idx_amb_L2 < 0 || allpars[idx_amb_L2].end < epoch)
			{
				if (Flag12 > 0)
				{
					tb_slip_12->set_amb_flag(rec, sat, Flag12);
					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
					isUpL1 = true;
					isUpL2 = true;
				}
			}
			else
			{
				bool Slip_12 = tb_slip_12->cycle_slip(data, epoch);
				if (Slip_12)
				{
					tb_slip_12->set_amb_flag(rec, sat, Flag12);
					end_time_12 = tb_slip_12 ? tb_slip_12->get_crt_amb_end(rec, sat) : LAST_TIME;
					isUpL1 = true;
					isUpL2 = true;
					if (Flag12 == 0) throw runtime_error("Flag12 is zero while idx_L2 >= 0");
				}
			}

			if (tb_slip_13 || tb_slip_14 || tb_slip_15)
			{
				if (idx_amb_L3 < 0)
				{
					if (Flag13 > 0)
					{
						tb_slip_13->set_amb_flag(rec, sat, Flag13);
						end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL3 = true;
					}
				}
				else
				{
					bool Slip_13 = tb_slip_13->cycle_slip(data, epoch);
					if (Slip_13)
					{
						tb_slip_13->set_amb_flag(rec, sat, Flag13);
						end_time_13 = tb_slip_13 ? tb_slip_13->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL3 = true;
						if (Flag13 == 0) throw runtime_error("Flag13 is zero while idx_L3 >= 0");
					}
				}

				if (idx_amb_L4 < 0)
				{
					if (Flag14 > 0)
					{
						tb_slip_14->set_amb_flag(rec, sat, Flag14);
						end_time_14 = tb_slip_14 ? tb_slip_14->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL4 = true;
					}
				}
				else
				{
					bool Slip_14 = tb_slip_14->cycle_slip(data, epoch);
					if (Slip_14)
					{

						tb_slip_14->set_amb_flag(rec, sat, Flag14);
						end_time_14 = tb_slip_14 ? tb_slip_14->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL4 = true;
						if (Flag14 == 0) throw runtime_error("Flag14 is zero while idx_L3 >= 0");
					}

				}				
				
				if (idx_amb_L5 < 0)
				{
					if (Flag15 > 0)
					{
						tb_slip_15->set_amb_flag(rec, sat, Flag15);
						end_time_15 = tb_slip_15 ? tb_slip_15->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL5 = true;
					}
				}
				else
				{
					bool Slip_15 = tb_slip_15->cycle_slip(data, epoch);
					if (Slip_15)
					{
						tb_slip_15->set_amb_flag(rec, sat, Flag15);
						end_time_15 = tb_slip_15 ? tb_slip_15->get_crt_amb_end(rec, sat) : LAST_TIME;
						isUpL1 = true;
						isUpL5 = true;
						if (Flag15 == 0) throw runtime_error("Flag12 is zero while idx_L5 >= 0");
					}

				}
				if (!isUpL2 && !isUpL3 && !isUpL4 && !isUpL5) continue;
				end_time = end_time.diff(end_time_12) < 0 ? end_time : end_time_12;
				end_time = end_time.diff(end_time_13) < 0 ? end_time : end_time_13;
				end_time = end_time.diff(end_time_14) < 0 ? end_time : end_time_14;
				end_time = end_time.diff(end_time_15) < 0 ? end_time : end_time_15;
			}
			else
			{
				if (!isUpL2) continue;
				end_time = end_time_12;
			}

			int num = 0;
			for (const auto& it : _band_index[sys])
			{
				par_type amb_type = par_type::AMB_IF;

				if (it.first == FREQ_1) amb_type = par_type::AMB_L1;
				else if (it.first == FREQ_2) amb_type = par_type::AMB_L2;
				else if (it.first == FREQ_3) amb_type = par_type::AMB_L3;
				else if (it.first == FREQ_4) amb_type = par_type::AMB_L4;
				else if (it.first == FREQ_5) amb_type = par_type::AMB_L5;
				else  throw std::logic_error("Unknown band!!");

				int idx_amb = allpars.getParam(rec, amb_type, sat);


				if (amb_type == par_type::AMB_L2 && Flag12 == 0) continue;
				if (amb_type == par_type::AMB_L3 && Flag13 == 0) continue;
				if (amb_type == par_type::AMB_L4 && Flag14 == 0) continue;
				if (amb_type == par_type::AMB_L5 && Flag15 == 0) continue;


				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1 + num, sat);


				par_newamb.value(0.0);
				par_newamb.setTime(epoch, end_time);
				par_newamb.apriori(9000.0);
				num++;

				if (idx_amb < 0)
				{
					update_info.add(par_newamb);
				}
				// reset the new amb
				else
				{
					allpars[idx_amb].end = epoch - _intv;
					update_info.add(idx_amb + 1, par_newamb);
				}
			}
		}

		// remove old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == par_type::AMB_L1 ||
				allpars[i].parType == par_type::AMB_L2 ||
				allpars[i].parType == par_type::AMB_L3 ||
				allpars[i].parType == par_type::AMB_L4 ||
				allpars[i].parType == par_type::AMB_L5)
			{
				if (allpars[i].end < epoch)
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

	void t_gupdateparALL::_lite_update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		for (auto data : obsdata)
		{
			int num_band = 0;
			string rec = data.site();
			string sat = data.sat();
			GSYS   sys = data.gsys();

			for (const auto& it : _band_index[sys])
			{
				par_type amb_type = par_type::AMB_IF;
				if (it.first == FREQ_1)      amb_type = par_type::AMB_L1;
				else if (it.first == FREQ_2) amb_type = par_type::AMB_L2;
				else if (it.first == FREQ_3) amb_type = par_type::AMB_L3;
				else if (it.first == FREQ_4) amb_type = par_type::AMB_L4;
				else if (it.first == FREQ_5) amb_type = par_type::AMB_L5;
				else {
					throw std::logic_error("Unknown band!!");
				};
				num_band++;
				if (num_band > _freq) continue;

				GOBS obsL = data.select_phase(it.second, true);

				int idx_amb = allpars.getParam(rec, amb_type, sat);
				if (idx_amb < 0)
				{
					if (obsL != GOBS::X)
					{
						int idx_max = allpars[allpars.parNumber() - 1].index;
						t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

						par_newamb.value(0.0);
						par_newamb.setTime(epoch, epoch);
						par_newamb.apriori(9000.0);
						update_info.add(par_newamb);
					}
				}
				else
				{
					if (obsL != GOBS::X)
					{
						if (data.getlli(obsL) < 1)
						{
							allpars[idx_amb].end = epoch;
						}
						else
						{
							allpars[idx_amb].end = epoch - _intv;

							int idx_max = allpars[allpars.parNumber() - 1].index;
							t_gpar par_newamb = t_gpar(rec, amb_type, idx_max + 1, sat);

							par_newamb.value(0.0);
							par_newamb.setTime(epoch, epoch);
							par_newamb.apriori(9000.0);
							update_info.add(idx_amb + 1, par_newamb);
						}
					}
					else  // obsL==X
					{
						update_info.add(idx_amb + 1);
					}
				}
			}
		}

		// remove old amb
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == par_type::AMB_L1 || allpars[i].parType == par_type::AMB_L2 ||
				allpars[i].parType == par_type::AMB_L3 || allpars[i].parType == par_type::AMB_L4 || allpars[i].parType == par_type::AMB_L5)
			{
				if (allpars[i].end < epoch)
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

	void t_gupdateparALL::_udpate_ion_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		map<string, set<string> > mapPrn;
		for (auto& data : obsdata)
		{
			mapPrn[data.site()].insert(data.sat());

			// add new sion par
			int idx_ion = allpars.getParam(data.site(), par_type::SION, data.sat());
			if (idx_ion < 0)
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_ion = t_gpar(data.site(), par_type::SION, idx_max + 1, data.sat());
				par_ion.value(0.0);
				par_ion.setTime(epoch, epoch);
				par_ion.apriori(_sig_ion);
				update_info.add(par_ion);
			}

			if (_iono_order == IONO_ORDER::THIRD_ORDER)
			{
				int idx_ion_3 = allpars.getParam(data.site(), par_type::SION_3, data.sat());
				if (idx_ion_3 < 0)
				{
					int idx_max = allpars[allpars.parNumber() - 1].index;
					t_gpar par_ion = t_gpar(data.site(), par_type::SION_3, idx_max + 1, data.sat());
					par_ion.value(0.0);
					par_ion.setTime(epoch, epoch);
					par_ion.apriori(9000);
					update_info.add(par_ion);
				}
			}
		}

		// remove old sion par
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == par_type::SION && mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
			{
				update_info.add(i + 1);
			}

			if (_iono_order == IONO_ORDER::THIRD_ORDER)
			{
				if (allpars[i].parType == par_type::SION_3 && mapPrn[allpars[i].site].count(allpars[i].prn) == 0)
				{
					update_info.add(i + 1);
				}
			}
		}

	}

	void t_gupdateparALL::_check_log(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13)
	{
	}
}
