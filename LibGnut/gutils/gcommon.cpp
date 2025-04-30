
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>

#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

	size_t t_gpair_string_hash::operator()(const pair<string, string>& a) const
	{
		return std::hash<string>()(a.first + a.second);
	}
} // namespace
