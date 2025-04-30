/**
 * @file         gcycleslip.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        cycleslip
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gcycleslip.h"
using namespace gnut;
great::t_gcycleslip::t_gcycleslip()
{
}

great::t_gcycleslip::t_gcycleslip(t_gsetbase * set, t_glog * log)
{
	_slip_model = dynamic_cast<t_gsetproc*>(set)->slip_model();
	_set = set;
	_log = log;
}

great::t_gcycleslip::~t_gcycleslip()
{
}
