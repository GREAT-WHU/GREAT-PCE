/**
 * @file         gsetlsq.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        The setting of lsq estimator
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GSETFLT_H
#define GSETFLT_H

#define XMLKEY_FLT "filter"

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"
#include "gexport/ExportLibGnut.h"
using namespace std;

namespace gnut
{
	/// The class of settings for lsq estimator
	class LibGnut_LIBRARY_EXPORT t_gsetlsq : public virtual t_gsetbase
	{
	public:
		t_gsetlsq();
		virtual ~t_gsetlsq();

		void check();                                  // settings check
		void help();                                   // settings help

	};

} // namespace

#endif
