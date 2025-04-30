/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_base.cpp
* @brief	main decode function for all files in base part
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#include "gcoders/gallcoder_base.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gcoders/poleut1.h"
#include "gcoders/dvpteph405.h"
#include "gio/gio.h"
#include "gio/glog.h"
#include "gio/gfile.h"

#include "gall/gallprec.h"
#include "gdata/gleapsecond.h"


using namespace gnut;
namespace great
{
	t_gallcoder_base::t_gallcoder_base()
	{
		_base_support.insert(SP3_INP);
		_base_support.insert(RINEXO_INP);
		_base_support.insert(RINEXC_INP);
		_base_support.insert(RINEXN_INP);
		_base_support.insert(ATX_INP);
		_base_support.insert(BLQ_INP);
		_base_support.insert(BIASINEX_INP);
		_base_support.insert(BIABERN_INP);
		_base_support.insert(SINEX_INP);
	}

	t_gallcoder_base::~t_gallcoder_base()
	{
	}

	bool t_gallcoder_base::support(const IFMT& ifmt)
	{
		if (find(_base_support.begin(), _base_support.end(), ifmt) == std::end(_base_support)) return false;
		return true;
	}

	void t_gallcoder_base::coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_base::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case SP3_INP: { gcoder = new t_sp3(this->_gset, "", 8172);         break; }
		case RINEXO_INP: 
		{ 
			gcoder = new t_rinexo(this->_gset, "", 4096);      
			break; 
		}
		case RINEXC_INP: 
		{ 
			gcoder = new t_rinexc(this->_gset, "", 4096);     
			break; 
		}
		case RINEXN_INP: { gcoder = new t_rinexn(this->_gset, "", 4096);      break; }
		case ATX_INP: { gcoder = new t_atx(this->_gset, "", 4096);         break; }
		case BLQ_INP: { gcoder = new t_blq(this->_gset, "", 4096);         break; }
		case BIASINEX_INP: { gcoder = new t_biasinex(this->_gset, "", 20480);   break; }
		case BIABERN_INP: { gcoder = new t_biabernese(this->_gset, "", 20480); break; }
		case SINEX_INP: { gcoder = new t_sinex(this->_gset, "", 20480);      break; }
		default:
		{
			if (_glog)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, _class_id, "coder", "not support for this format " + t_gsetinp::ifmt2str(ifmt));
			}
			return;
		}
		}

		if (gcoder)
		{
			// prepare gio for the file
			t_gio* gio = new t_gfile;
			gio->glog(this->_glog);
			gio->path(path);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path(path);
			gcoder->glog(this->_glog);

			// Put the data container into gcoder
			_data_id++;
			gcoder->add_data("ID"+int2str(_data_id), data);
			gcoder->add_data("OBJ", obj);

			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			gio->coder(gcoder);

			// Read the data from file here
			t_gtime runepoch = t_gtime::current_time(t_gtime::GPS);
			gio->run_read();
			t_gtime lstepoch = t_gtime::current_time(t_gtime::GPS);
			this->_glog->comment(0, "main", "READ: " + path + " time: "
				+ dbl2str(lstepoch.diff(runepoch)) + " sec");

			if (gio)    { delete gio;    gio = nullptr; }
			if (gcoder) { delete gcoder; gcoder = nullptr; }
		}
	}

	void t_gallcoder_base::coder(t_gtime beg, t_gtime end, const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj)
	{
		if (!t_gallcoder_base::support(ifmt)) return;

		t_gcoder* gcoder = nullptr;
		switch (ifmt)
		{
		case SP3_INP: { gcoder = new t_sp3(this->_gset, "", 8172);         break; }
		case RINEXO_INP:
		{
			gcoder = new t_rinexo(beg, end, this->_gset, "", 4096);
			break;
		}
		case RINEXC_INP:
		{
			gcoder = new t_rinexc(this->_gset, "", 4096);
			break;
		}
		case RINEXN_INP: { gcoder = new t_rinexn(this->_gset, "", 4096);      break; }
		case ATX_INP: { gcoder = new t_atx(this->_gset, "", 4096);         break; }
		case BLQ_INP: { gcoder = new t_blq(this->_gset, "", 4096);         break; }
		case BIASINEX_INP: { gcoder = new t_biasinex(this->_gset, "", 20480);   break; }
		case BIABERN_INP: { gcoder = new t_biabernese(this->_gset, "", 20480); break; }
		case SINEX_INP: { gcoder = new t_sinex(this->_gset, "", 20480);      break; }
		default:
		{
			if (_glog)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, _class_id, "coder", "not support for this format " + t_gsetinp::ifmt2str(ifmt));
			}
			return;
		}
		}

		if (gcoder)
		{
			// prepare gio for the file
			t_gio* gio = new t_gfile;
			gio->glog(this->_glog);
			gio->path(path);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path(path);
			gcoder->glog(this->_glog);

			// Put the data container into gcoder
			_data_id++;
			gcoder->add_data("ID" + int2str(_data_id), data);
			gcoder->add_data("OBJ", obj);

			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			gio->coder(gcoder);

			// Read the data from file here
			t_gtime runepoch = t_gtime::current_time(t_gtime::GPS);
			gio->run_read();
			t_gtime lstepoch = t_gtime::current_time(t_gtime::GPS);
			this->_glog->comment(t_glog::LOG_LV::LOG_INFO, "coder", "main", "READ: " + path + " time: "
				+ dbl2str(lstepoch.diff(runepoch)) + " sec");

			if (gio) { delete gio;       gio = nullptr; }
			if (gcoder) { delete gcoder; gcoder = nullptr; }
		}
	}

	void t_gallcoder_base::creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		if (!t_gallcoder_base::support(ifmt)) return;

		auto type = ifmt2data(ifmt);
		if (map_gdata.count(type) > 0) return;

		switch (ifmt)
		{
		case RINEXC_INP:
		case RINEXN_INP:
		case SP3_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallprec();
			dynamic_cast<t_gallnav*>(map_gdata[type])->glog(_glog);
			dynamic_cast<t_gallprec*>(map_gdata[type])->use_clknav(true);
			dynamic_cast<t_gallprec*>(map_gdata[type])->use_clkrnx(true);
			dynamic_cast<t_gallprec*>(map_gdata[type])->use_clksp3(true);
			dynamic_cast<t_gallprec*>(map_gdata[type])->use_posnav(true);
			break;
		}
		case RINEXO_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallobs();
			dynamic_cast<t_gallobs*>(map_gdata[type])->glog(_glog);
			dynamic_cast<t_gallobs*>(map_gdata[type])->gset(_gset);
			break;
		}
		case ATX_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallpcv();
			dynamic_cast<t_gallpcv*>(map_gdata[type])->glog(_glog);
			break;
		}
		case BLQ_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallotl();
			dynamic_cast<t_gallotl*>(map_gdata[type])->glog(_glog);
			break;
		}
		case BIASINEX_INP:
		case BIABERN_INP:
		{
			map_gdata[type] = nullptr;
			map_gdata[type] = new t_gallbias();
			dynamic_cast<t_gallbias*>(map_gdata[type])->glog(_glog);
			break;
		}
		case SINEX_INP:
		{
			// no need save in t_gdata, save in t_gallobj
			break;
		}
		default:
		{
			map_gdata[type] = nullptr;
			if (_glog)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, _class_id, "creat", "not support for this format : " + t_gsetinp::ifmt2str(ifmt));
			}
			return;
		}
		}
	}
}