/**
 * @file         grecoverdata.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Storage all recover observ equation/ parameter information
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#include  <gdata/grecoverdata.h>
#include <sstream>


namespace great
{

	t_grcover_head::t_grcover_head():
		interval(0.0)
	{
	}

	t_grcover_head::~t_grcover_head()
	{
	}

	t_gtime t_grcover_head::get_beg_time() const
	{
		return *time_list.begin();
	}

	t_gtime t_grcover_head::get_end_time() const
	{
		return *time_list.end();
	}

	t_grecover_data::t_grecover_data()
	{
	}

	t_grecover_data::~t_grecover_data()
	{
		
	}

	t_grecover_equation::t_grecover_equation(const t_gtime & time, const string site, string sat) :
		time(time),
		sat_name(sat),
		site_name(site),
		weight(0.0),
		resuidal(0.0),
		is_newamb(0)
	{
		_type = t_gdata::RESOBS;
	}

	t_grecover_equation::t_grecover_equation(const t_grecover_equation & other) :
		site_name(other.site_name),
		sat_name(other.sat_name),
		time(other.time),
		obstype(other.obstype),
		weight(other.weight),
		resuidal(other.resuidal),
		is_newamb(other.is_newamb)
	{
		time = other.time;
	}

	t_grecover_equation::~t_grecover_equation()
	{
	}

	void t_grecover_equation::operator=(const t_grecover_equation & other)
	{
		site_name = other.site_name;
		sat_name = other.sat_name;
		time = other.time;
		obstype = other.obstype;
		weight = other.weight;
		resuidal = other.resuidal;
		is_newamb = other.is_newamb;
	}

	t_gtime t_grecover_equation::get_recover_time() const
	{
		return time;
	}

	string t_grecover_equation::convert2strline() const
	{
		stringstream strline("");

		// format control
		strline << setiosflags(ios::right) << setprecision(4) << setiosflags(ios::fixed)
			<< setw(5) << "RES:="
			<< setw(25) << time.str()
			<< setw(5) << is_newamb
			<< setw(8) << site_name
			<< setw(8) << sat_name
			<< setw(8) << obstype.convert2str()
			<< setw(15) << weight
			<< setw(15) << resuidal << endl;

		return strline.str();
	}

	void t_grecover_equation::set_recover_equation(const t_gobscombtype & obstype, const pair<double, double>& resinfo,int is_newamb) 
	{
		this->obstype = obstype;
		this->weight = resinfo.first;
		this->resuidal = resinfo.second;
		this->is_newamb = is_newamb;
	}

	t_grecover_par::t_grecover_par(const t_grecover_par & other):
		par(other.par),
		correct_value(other.correct_value)
	{
	}

	t_grecover_par::t_grecover_par(const t_gpar & par, double correct_value):
		par(par),
		correct_value(correct_value)
	{
		_type = t_gdata::RESPAR;
	}

	t_grecover_par::~t_grecover_par()
	{
	}

	t_gtime t_grecover_par::get_recover_time() const
	{
		return par.end;
	}

	string t_grecover_par::convert2strline() const
	{
		stringstream strline("");

		if (double_eq(correct_value, 0.0)) return strline.str();

		strline << setiosflags(ios::right) << setiosflags(ios::fixed) 
			<< setw(5) << "PAR:="
			<< setw(25) << gpar2str(par)
			<< setw(25) << par.beg.str()
			<< setw(25) << par.end.str()
			<< setw(25) << right << setprecision(7) << fixed << par.value()
			<< setw(25) << right << setprecision(7) << scientific << uppercase << correct_value
			<< setw(25) << right << setprecision(7) << fixed << par.value() + correct_value
			<< endl;
		return strline.str();
	}

	bool t_grecover_par::operator<(const t_grecover_par& data) const
	{
		if (this->par > data.par)
			return true;
		return false;
	}

}
