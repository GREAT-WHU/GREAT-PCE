/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder.cpp
* @brief	main decode function for all files
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#include "gcoders/gallcoder.h"
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

namespace great
{
	great::t_gallcoder::~t_gallcoder()
	{

	}

	void t_gallcoder::log(t_glog* log)
	{
		this->_glog = log;
	}

	void t_gallcoder::set(t_gsetbase* set)
	{
		this->_gset = set;
	}

	IFMT             t_gallcoder::data2ifmt(t_gdata::ID_TYPE type)
	{
		switch (type)
		{
		case UNDEF:
		default:
			break;
		}
		return IFMT::UNDEF;
	};

	t_gdata::ID_TYPE t_gallcoder::ifmt2data(IFMT             ifmt)
	{
		switch (ifmt)
		{
		case SP3_INP:
		case RINEXC_INP:
		case RINEXN_INP:
			return t_gdata::ID_TYPE::ALLNAV;
		case ATX_INP:
			return t_gdata::ID_TYPE::ALLPCV;
		case BLQ_INP:
			return t_gdata::ID_TYPE::ALLOTL;
		case OCEANTIDE_INP:
			return t_gdata::ID_TYPE::OCEANTIDE;
		case RINEXO_INP:
			return t_gdata::ID_TYPE::ALLOBS;
		case BIASINEX_INP:
		case BIABERN_INP:
			return t_gdata::ID_TYPE::ALLBIAS;
		case POLEUT1_INP:
			return t_gdata::ID_TYPE::ALLPOLEUT1;
		case LEAPSECOND_INP:
			return t_gdata::ID_TYPE::LEAPSECOND;
		case SATPARS_INP:
			return t_gdata::ID_TYPE::SATPARS;
		case DE_INP:
			return t_gdata::ID_TYPE::ALLDE;
		case SINEX_INP:
		default:
			return t_gdata::ID_TYPE::LAST;
		}
	}

	void t_gallcoder::destroy(map<t_gdata::ID_TYPE, t_gdata*>& map_gdata)
	{
		for (auto& iter : map_gdata)
		{
			if (iter.second) {
				delete iter.second; iter.second = nullptr;
			}
		}
	}
}
