/**
*
* @verbatim
History
-1.0 BoWong  2019-02-20  creat the file.
-1.1 BoWong  2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallambflag.cpp
* @brief	Storage the XXXXddd0.yyo.log files'(ambiguity flag) data(more than one site)
*				XXXX  ---- SITE name
*				ddd  ---- Doy of the file
*				yy  ---- year
*
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#include "gall/gallambflag.h"
#include "gutils/gtime.h"
#include <algorithm>

using namespace std;

namespace great {

	t_gallambflag::t_gallambflag(ID_TYPE t)
	{
		if (t != t_gdata::AMBFLAG && t != t_gdata::AMBFLAG13)
		{
			id_type(t_gdata::AMBFLAG);
		}
		else
		{		
			id_type(t);
		}

	}

	t_gallambflag::~t_gallambflag()
	{
		_map_allambflag.clear();
	}

	void t_gallambflag::addAmbFlagHead(string sta, const ambflag_hd& hd)
	{
		_map_allambflag[sta].setAmbFlagHead(hd);
		_map_allambflag[sta].setSite(sta);
	}

	void t_gallambflag::addAmbFlag(string sta, const t_gambflag& ambflag)
	{
		_map_allambflag[sta] = ambflag;
	}

	void t_gallambflag::addAmbFlagData(string sta, string sat, const ambflag_data& data)
	{
		_map_allambflag[sta].addAmbFlagData(sat, data);
	}

	t_gambflag& t_gallambflag::getOneAmbFlag(string sta)
	{
		transform(sta.begin(), sta.end(), sta.begin(), ::tolower);
		if (_map_allambflag.find(sta) != _map_allambflag.end())
		{
			return _map_allambflag[sta];
		}
		transform(sta.begin(), sta.end(), sta.begin(), ::toupper);
		if (_map_allambflag.find(sta) != _map_allambflag.end())
		{
			return _map_allambflag[sta];
		}
		return _map_allambflag[sta];
	}

	t_map_allambflag& t_gallambflag::getAllAmbFlag()
	{
		return _map_allambflag;
	}

	set<string> t_gallambflag::getSatList()
	{
		set<string> satlist;
		set<string> tmp;
		for (auto it_site = _map_allambflag.begin(); it_site != _map_allambflag.end(); it_site++)
		{
			tmp = it_site->second.getAllSatSet();
			satlist.insert(tmp.begin(),tmp.end());
		}
		return satlist;
	}

	set<string> t_gallambflag::getSiteList()
	{
		set<string> sitelist;
		string tmp;
		for (auto it_site = _map_allambflag.begin(); it_site != _map_allambflag.end(); it_site++)
		{
			tmp = it_site->first;
			transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
			sitelist.insert(tmp);
		}
		return sitelist;
	}

} // namespace

