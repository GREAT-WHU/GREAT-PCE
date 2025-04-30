/**
 * @file         gsetturboedit.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetturboedit.h"

using namespace std;
using namespace pugi;


namespace great
{

	t_gsetturboedit::t_gsetturboedit()
	{
		_set.insert(XMLKEY_TURBOEDIT);
	    _defaultPCLimit = 250.0;       // unit: meter
	    _defaultMWLimit = 4.0;         // unit: cycle
	    _defaultGFLimit = 1.0;         // unit: cycle
	    _defaultGFRmsLimit = 2.0;      // unit: cycle
	    _defaultSingleFreqLimit = 1.0; // unit: meter
		_defaultGapArcLimit = 20;      // unit: epoch
		_defaultShortArcLimit = 10;    // unit: epoch
	    _defaultMinPercent = 60.0;     // unit: %
	    _defaultMinMeanNprn = 4;       // unit: 
	    _defaultMaxMeanNamb = 3;       // unit: 
	}

	bool t_gsetturboedit::liteMode()
	{
		bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).attribute("lite_mode").as_bool();
		return tmp;
	}

	bool t_gsetturboedit::isAmbOutput()
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("amb_output");
		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::isEphemeris()
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("ephemeris");
		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkPC(double& pc_limit)
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_pc");
		
		string tmp = tmp_set.attribute("pc_limit").value();
		if (tmp.empty())  pc_limit = _defaultPCLimit;
		else              pc_limit=str2dbl(tmp);
		
		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkMW(double& mw_limit)
	{
		bool is_litemodel = this->liteMode();
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_mw");
		string tmp = tmp_set.attribute("mw_limit").value();

		if (is_litemodel)
		{
			if (tmp.empty())  mw_limit = 0.0;
			else              mw_limit = str2dbl(tmp);
		}
		else
		{
			if (tmp.empty())  mw_limit = _defaultMWLimit;
			else              mw_limit = str2dbl(tmp);
		}

		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkGF(double& gf_limit, double& gf_rms_limit)
	{
		bool is_litemodel = this->liteMode();
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_gf");
		string tmp_gf = tmp_set.attribute("gf_limit").value();
		string tmp_rms = tmp_set.attribute("gf_rms_limit").value();
		if (is_litemodel)
		{
			if (tmp_gf.empty()) gf_limit = 0.0;
			else             gf_limit = str2dbl(tmp_gf);

			if (tmp_rms.empty()) gf_rms_limit = 0.0;
			else             gf_rms_limit = str2dbl(tmp_rms);
		}
		else
		{
			if (tmp_gf.empty()) gf_limit = _defaultGFLimit;
			else             gf_limit = str2dbl(tmp_gf);

			if (tmp_rms.empty()) gf_rms_limit = _defaultGFRmsLimit;
			else             gf_rms_limit = str2dbl(tmp_rms);
		}

		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkSingleFreq(double& sf_limit)
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_sf");

		string tmp = tmp_set.attribute("sf_limit").value();
		if (tmp.empty())  sf_limit = _defaultSingleFreqLimit;
		else              sf_limit = str2dbl(tmp);

		return  tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkGap(int& gap_limit)
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_gap");

		string tmp = tmp_set.attribute("gap_limit").value(); 

		if (tmp.empty())  gap_limit = _defaultGapArcLimit;
		else              gap_limit = str2int(tmp);

		return tmp_set.attribute("valid").as_bool();
	}

	int t_gsetturboedit::smoothWindows()
	{
		double value;
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("smooth_win");

		string tmp = tmp_set.attribute("value").value();

		if (tmp.empty())  value = 25;
		else              value = str2int(tmp);

		return value;
	}

	bool t_gsetturboedit::checkShort(int& short_limit)
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_short");

		string tmp = tmp_set.attribute("short_limit").value();

		if (tmp.empty())  short_limit = _defaultShortArcLimit;
		else              short_limit = str2int(tmp);

		return tmp_set.attribute("valid").as_bool();
	}

	bool t_gsetturboedit::checkStatistics(double& minPercent, int& minMeanNprn, int& maxMeanNamb)
	{
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_TURBOEDIT).child("check_statistics");

		string tmp = tmp_set.attribute("min_percent").value();
		if (tmp.empty())  minPercent = _defaultMinPercent;
		else              minPercent = str2dbl(tmp);

		tmp = tmp_set.attribute("min_mean_nprn").value();
		if (tmp.empty())  minMeanNprn = _defaultMinMeanNprn;
		else              minMeanNprn = str2int(tmp);

		tmp = tmp_set.attribute("max_mean_namb").value();
		if (tmp.empty())  maxMeanNamb = _defaultMaxMeanNamb;
		else              maxMeanNamb = str2int(tmp);

		return  tmp_set.attribute("valid").as_bool();
	}

	void t_gsetturboedit::help()
	{	
		cerr << "<turboedit  lite_mode=\"false\" >" << endl
			<< "<amb_output  valid=\"true\"  />  " << endl
			<< "<simulation  valid=\"false\" />  " << endl
			<< "<ephemeris   valid=\"true\"  />  " << endl
			<< "<check_pc   pc_limit=\"250\" valid=\"true\"  />" << endl
			<< "<check_mw   mw_limit=\"4\"   valid=\"true\"  />" << endl
			<< "<check_gf   gf_limit=\"1\"   gf_rms_limit=\"2\"   valid=\"true\" />" << endl
			<< "<check_sf   sf_limit=\"1\"   valid=\"false\" />         " << endl
			<< "<check_gap    gap_limit=\"20\"    valid=\"true\" />     " << endl
			<< "<check_short  short_limit=\"10\"  valid=\"true\" />     " << endl
			<< "<check_statistics  min_percent=\"60\"  min_mean_nprn=\"4\"  max_mean_namb=\"3\"  valid=\"true\" /> " << endl
			<< "</turboedit>" << endl
			<< endl;
		return;
	}

}