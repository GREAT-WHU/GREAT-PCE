
/**
*
* @verbatim
	History
	2016-06-26  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        ginfolog.cpp
* @brief       Purpose: file conversion utilities
* @author      JD
* @version     1.0.0
* @date        2016-06-26
*
*/

#include "ginfolog.h"
#include <iostream>
#include <iomanip>

void LibGnut_LIBRARY_EXPORT gnut::write_log(t_glog* log, const LOG_LV& LV, const string& class_id, const string& function_id, const string& description)
{

}

void gnut::write_log(t_glog * log, int l, string description)
{
	write_log_info(log, l, description);
}

void gnut::write_log(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}

void gnut::writeLogInfo(t_glog * log, int l, string description)
{
	if (log)
	{
		log->comment(l, description);
	}
	else
	{
		std::cout << setw(88) << description << endl;
	}
}

void gnut::writeLogInfo(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}

void gnut::write_log_info(t_glog * log, int l, string description)
{
	if (log)
	{
		log->comment(l, description);
	}
	else
	{
		std::cout << setw(88) << description << endl;
	}
}

void gnut::write_log_info(t_glog * log, int l, string ID, string description)
{
	if (log)
	{
		log->comment(l, ID, description);
	}
	else
	{
		std::cout << std::setw(8) << ID + " : " << std::setw(88) << description << endl;
	}
}
