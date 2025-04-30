/**
 * @file         gsetamb.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        control set from XML
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gset/gsetamb.h"
#include <sstream>

using namespace std;
using namespace pugi;

namespace great
{
	/** @brief default constructor. */
	t_gsetamb::t_gsetamb()
	 : t_gsetbase()
	{
		_gmutex.lock();

		_set.insert(XMLKEY_AMBIGUITY);

		_gmutex.unlock();
	}

	/** @brief default destructor. */
	 t_gsetamb::~t_gsetamb()
	{

	}

	 /**
	 * @brief settings check.
	 * @return	  void
	 */
	 void t_gsetamb::check()
	 {
		 _gmutex.lock();
		 _gmutex.unlock();
		 return;
	 }

	 /**
	 * @brief settings help.
	 * @return	  void
	 */
	 void t_gsetamb::help()
	 {
		 _gmutex.lock();
		 cerr << "<ambiguity>\n"
			 << "<upd> EWL/WL/NL/NONE </upd>\n"
			 << "<fix_mode> ROUND/SEARCH/NO </fix_mode>\n"
			 << "<ratio> 3.0 </ratio>\n"
			 << "<all_baselines> NO </all_baselines>\n"
			 << "<min_common_time> 30 </min_common_time>\n"
			 << "<widelane_decision     maxdev = \"0.15\" maxsig = \"0.10\" alpha = \"1000\"/>\n"
			 << "<narrowlane_decision   maxdev = \"0.15\" maxsig = \"0.10\" alpha = \"1000\"/>\n";

		 _gmutex.unlock();
		 return;
	 }

	 FIX_MODE t_gsetamb::fix_mode()
	 {
		 if (!_doc) str2fixmode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("fix_mode");
		 return str2fixmode(trim(tmp));
	 }

	 UPD_MODE t_gsetamb::upd_mode() {
		 if (!_doc) str2upd_mode(string());
		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("upd_mode");
		 return str2upd_mode(trim(tmp));
	 }

	 double t_gsetamb::lambda_ratio()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("ratio");

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::bootstrapping()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("boot");

		 if (tmp.empty())
		 {
			 _gmutex.unlock();
			 return -0.001;
		 }

		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 double t_gsetamb::min_common_time()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("min_common_time");

		 _gmutex.unlock(); return tmp.empty() ? 0 : str2dbl(tmp);
	 }

	 map<string, double> t_gsetamb::get_amb_decision(string type)
	 {
		 _gmutex.lock();
		 map<string, string> type2child = { 
			 {"EWL", "extra_widelane_decision"}, {"WL", "widelane_decision"}, {"NL", "narrowlane_decision"}
		 };
		 map<string, double> amb_decision;
		 if (_default_decision.find(type) == _default_decision.end()) return amb_decision;
		 amb_decision = _default_decision[type];
		 for (auto iter = amb_decision.begin(); iter != amb_decision.end(); ++iter) {
			 double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child(type2child[type].c_str()).attribute(iter->first.c_str()).as_double();
			 if (!double_eq(tmp, 0)) iter->second = tmp;
		}
		 _gmutex.unlock(); return amb_decision;
	 }

	 bool t_gsetamb::part_ambfix()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("part_fix"));
		 string tmp;
		 is >> tmp;
		 bool is_part_fix = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return is_part_fix;
	 }

	 int t_gsetamb::part_ambfix_num()
	 {
		 _gmutex.lock();

		 string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("part_fix_num");

		 if (tmp.empty())
		 {
			 _gmutex.unlock();
			 return 2;
		 }
		 _gmutex.unlock(); return str2dbl(tmp);
	 }

	 /**
	  * @brief change string to FIXMODE.
	  * @param[in]  str   fix mode string
	  * @return	  FIXMODE
	  */
	 FIX_MODE t_gsetamb::str2fixmode(string str)
	 {
		 if (str == "NO") 
		 {
			 return FIX_MODE::NO;
		 }
		 else if (str == "ROUND") 
		 {
			 return FIX_MODE::ROUND;
		 }
		 else if (str == "SEARCH") 
		 {
			 return FIX_MODE::SEARCH;
		 }
		 else if (str == "HOLD")
		 {
			 return FIX_MODE::HOLD;
		 }
		 else
		 {
			 cout << "*** warning: not defined ambiguity fixing mode [" << str << "]\n"; cout.flush();
		 }
		 return FIX_MODE::NO;
	 }

	 /**
	 * @brief change string to UPD_MODE.
	 * @param[in]  str   upd mode string
	 * @return	  UPD_MODE
	 */
	 UPD_MODE t_gsetamb::str2upd_mode(string str) {
		 if (str == "UPD" || str == "upd") {
			 return UPD_MODE::UPD;
		 }
		 else {
			 cout << "*** warning: not defined upd mode [" << str << "]\n"; cout.flush();
			 return UPD_MODE::UPD;
		 }
	 }

	 bool t_gsetamb::isSetRefSat()
	 {
		 _gmutex.lock();

		 istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_AMBIGUITY).child_value("set_refsat"));
		 string tmp;
		 is >> tmp;
		 bool isSetRefsat = (tmp == "YES" || tmp == "yes");

		 _gmutex.unlock(); return isSetRefsat;
	 }


}
