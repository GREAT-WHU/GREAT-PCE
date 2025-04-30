/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder.h
* @brief	main decode function for all files
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/
#ifndef GALLCODER_H
#define GALLCODER_H

#include "gexport/ExportLibGREAT.h"
#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"

using namespace gnut;
namespace great
{
	/** @brief class for t_gallcoder. */
	class LibGREAT_LIBRARY_EXPORT t_gallcoder
	{
	public:
		/** @brief default constructor. */
		t_gallcoder() {};
		virtual ~t_gallcoder();

		virtual void log(t_glog* log);
		virtual void set(t_gsetbase* set);

		IFMT                     data2ifmt(t_gdata::ID_TYPE type);
		t_gdata::ID_TYPE ifmt2data(IFMT             ifmt);

		virtual bool support(const IFMT& ifmt) = 0;
		virtual void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) = 0;
		virtual void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) = 0;
		void destroy(map<t_gdata::ID_TYPE, t_gdata*>& map_gdata);
	protected:
		t_gsetbase*              _gset = nullptr;
		t_glog*                  _glog = nullptr;
	};
}

#endif