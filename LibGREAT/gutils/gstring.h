/**
 * @file         gstring.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        for string
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include <string>
#include <vector>
#include <set>
#include "gexport/ExportLibGREAT.h"
using namespace std;

namespace great
{
	LibGREAT_LIBRARY_EXPORT std::string trim(std::string s);
	LibGREAT_LIBRARY_EXPORT std::string clean(std::string s);
	LibGREAT_LIBRARY_EXPORT void split(const std::string& s, std::string delim, std::vector< std::string >& ret);
	
	LibGREAT_LIBRARY_EXPORT std::string set2str(std::set<std::string> &s);
	LibGREAT_LIBRARY_EXPORT std::string str2upper(const std::string& str);
	LibGREAT_LIBRARY_EXPORT std::string str2lower(const std::string& str);
	LibGREAT_LIBRARY_EXPORT std::string do_fraction(long double value, int digits);
	LibGREAT_LIBRARY_EXPORT std::string format(const char* fmt, ...);

	LibGREAT_LIBRARY_EXPORT std::string remove_surplus_spaces(const std::string& src);

}
