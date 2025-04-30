/**
 * @file         gstring.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        for string
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gutils/gstring.h"
#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <math.h>
namespace great
{
	std::string trim(std::string s)
	{
		if (s.empty())
		{
			return s;
		}

		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);
		return s;
	}

	std::string clean(std::string s)
	{
		if (s.empty())
		{
			return s;
		}

		s.erase(0, s.find_first_not_of(" \r\t\n"));
		s.erase(s.find_last_not_of(" \r\t\n") + 1);

		s.erase(0, s.find_first_not_of(" \r\t\n"));
		s.erase(s.find_last_not_of(" \r\t\n") + 1);
		return s;
	}

	void split(const std::string& s, std::string delim, std::vector< std::string >& ret)
	{
		size_t last = 0;
		size_t index = s.find_first_of(delim, last);

		while (index != std::string::npos)
		{
			ret.push_back(s.substr(last, index - last));
			last = index + 1;
			index = s.find_first_of(delim, last);
		}
		if (index - last > 0)
		{
			ret.push_back(s.substr(last, index - last));
		}
	}

	std::string set2str(std::set<std::string> &s)
	{
		bool isChangeLine = false;
		std::string tmp;
		std::set<std::string>::const_iterator beg = s.begin();
		std::set<std::string>::const_iterator end = s.end();
		std::set<std::string>::const_iterator crt = beg;

		if (s.size() > 10)  isChangeLine = true;
		int num = 0;
		for (crt; crt != end; crt++)
		{
			if (num % 10 == 0)
			{
				if (isChangeLine) tmp = tmp + "\n\t\t\t";
			}

			num++;
			tmp = tmp + *crt + " ";
		}
		if (isChangeLine) tmp = tmp + "\n\t\t";
		return tmp;
	}

	string str2upper(const string & str)
	{
		string tmp = str;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		return tmp;
	}

	string str2lower(const string & str)
	{
		string tmp = str;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
		return tmp;
	}

	string do_fraction(long double value, int digits)
	{
		long double val1 = value;
		long double  resVal;

		val1 = value * pow(10, digits);

		resVal = int(val1);
		resVal = resVal / (pow(10, digits));

		ostringstream os;
		os << resVal;
		string result;
		istringstream is(os.str());
		is >> result;

		return result;
	}

	std::string format(const char* fmt, ...)
	{
		va_list ap;
		va_start(ap, fmt);
		int len = vsnprintf(nullptr, 0, fmt, ap);
		va_end(ap);
		std::string buf(len + 1, '\0');
		va_start(ap, fmt);
		vsnprintf(&buf[0], buf.size(), fmt, ap);
		va_end(ap);
		buf.pop_back();
		return buf;
	}

	std::string remove_surplus_spaces(const std::string& src)
	{
		std::string result = "";
		for (int i = 0; src[i] != '\0'; i++)
		{
			if (src[i] != ' ')
				result.append(1, src[i]);
			else
				if (src[i + 1] != ' ')
					result.append(1, src[i]);
		}
		return result;
	}
}