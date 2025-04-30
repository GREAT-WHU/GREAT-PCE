/**
*
* @verbatim
History
	-1.0 Jiande Huang  2019-10-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gallcoder_base.h
* @brief	main decode function for all files in base part
*
*
* @author   Jiande Huang, Wuhan University
* @version	1.0.0
* @date		2019-10-25
*
*/

#ifndef GALLCODER_BASE_H
#define GALLCODER_BASE_H

#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"
#include "gall/gallobj.h"
#include "gcoders/gallcoder.h"

using namespace gnut;
namespace great
{
	/** @brief class for t_gallcoder_base. */
	class LibGREAT_LIBRARY_EXPORT t_gallcoder_base : public virtual t_gallcoder
	{
	public:
		/** @brief default constructor. */
		t_gallcoder_base();
		virtual ~t_gallcoder_base();

		bool support(const IFMT& ifmt) override;
		void coder(const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj) override;
		void coder(t_gtime beg, t_gtime end, const IFMT& ifmt, string path, t_gdata* data, t_gallobj* obj);
		void creat(const IFMT& ifmt, map<t_gdata::ID_TYPE, t_gdata*>& map_gdata) override;


	protected:
		std::set<IFMT>   _base_support;
		string  _class_id = "t_gallcoder_base";
		int  _data_id = 0;
	};
}

#endif