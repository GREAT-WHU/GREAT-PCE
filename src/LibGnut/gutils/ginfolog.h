
/**
*
* @verbatim
	History
	2016-06-26  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        ginfolog.h
* @brief       Purpose: file conversion utilities
* @author      JD
* @version     1.0.0
* @date        2016-06-26
*
*/

#ifndef GINFOLOG_H
#define GINFOLOG_H

#include "gio/glog.h"
#include "gexport/ExportLibGnut.h"
/*
	ERROR 0 : error information
	NOTE  0 : important note information
	WARN  0 : important warnning information
	WARN  1 : warnning information that not important
	NOTE  2 :
*/

#define LV1 "ERROR  "
#define LV2 "WANNING"
#define LV3 "NOTE   "
#define LV4 "OUTPUT "

namespace gnut
{
	enum class LOG_LV { LOG_ERROR, LOG_WARN, LOG_NOTE, LOG_OUTPUT, LOG_DEF};

	void   LibGnut_LIBRARY_EXPORT write_log(t_glog* log, const LOG_LV& LV, const string& class_id, const string& function_id, const string& description);

	void   LibGnut_LIBRARY_EXPORT write_log(t_glog* log, int l, string description);
	void   LibGnut_LIBRARY_EXPORT write_log(t_glog* log, int l, string ID, string description);

	void   LibGnut_LIBRARY_EXPORT writeLogInfo(t_glog* log, int l, string description);
	void   LibGnut_LIBRARY_EXPORT writeLogInfo(t_glog* log, int l, string ID, string description);

	void   LibGnut_LIBRARY_EXPORT write_log_info(t_glog* log, int l, string description);
	void   LibGnut_LIBRARY_EXPORT write_log_info(t_glog* log, int l, string ID, string description);



}
#endif // !GINFOLOG_H
