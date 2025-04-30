/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file		gsetgen.cpp
* @brief	implements common general settings
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gutils/gnss.h"
#include "gutils/gsys.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetgen.h"

using namespace std;
using namespace pugi;

namespace gnut
{

	// Constructor
	// ----------
	t_gsetgen::t_gsetgen(bool gnss)
		: t_gsetbase(),
		_gnss(gnss),
		_dec(0)
	{
		_set.insert(XMLKEY_GEN);

		// initiate GNSS string
		if (_gnss) {
			t_map_sats gnss_sats = GNSS_SATS();
			for (auto itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS) {
				_sys += " " + t_gsys::gsys2str(itGNS->first);
			}
		}
	}

	// Destructor
	// ----------
	t_gsetgen::~t_gsetgen()
	{}


	// Return gtime
	// ----------
	t_gtime t_gsetgen::beg(bool conv)
	{
		_gmutex.lock();

		string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value("beg");
		substitute(str, "\n", "");
		substitute(str, "\"", "");

		t_gtime gt(t_gtime::GPS);

		if (str.empty()) {
			_gmutex.unlock();
			t_gtime tmp(FIRST_TIME);
			return tmp;
		}

		gt.from_str("%Y-%m-%d %H:%M:%S", trim(str), conv);
		_gmutex.unlock(); return gt;
	}


	// Return gtime
	// ----------
	t_gtime t_gsetgen::end(bool conv)
	{
		_gmutex.lock();

		string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value("end");
		substitute(str, "\n", "");
		substitute(str, "\"", "");

		t_gtime gt(t_gtime::GPS);

		if (str.empty()) 
		{
			_gmutex.unlock(); t_gtime tmp(LAST_TIME); return tmp;
		}

		gt.from_str("%Y-%m-%d %H:%M:%S", trim(str), conv);
		_gmutex.unlock(); return gt;
	}


	// Return sampling
	// ----------
	double t_gsetgen::sampling()
	{
		_gmutex.lock();

		string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value("int");

		// delete spaces
		str.erase(remove(str.begin(), str.end(), ' '), str.end());

		double tmp = str2dbl(str);

		if (str.find(".") != string::npos) {
			_dec = str.substr(str.find(".") + 1).length(); // decimal digits resolution
		}

		_gmutex.unlock(); return tmp;
	}


	// Return set
	// ----------
	set<string> t_gsetgen::sys()
	{
		_gmutex.lock();

		set<string> xcl, tmp = t_gsetbase::_setval(XMLKEY_GEN, "sys");

		// exclude starting with '-'
		set<string>::iterator itSYS, itTMP;
		for (itTMP = tmp.begin(); itTMP != tmp.end(); ) {
			if ((*itTMP)[0] == '-') { xcl.insert(*itTMP); itSYS = itTMP; ++itTMP; tmp.erase(itSYS); }//"-"mean without this system?
			else ++itTMP;
		}

		// if empty, complete, i.e. if only exclusions listed (and gnss requested!)
		if (tmp.size() == 0 && _gnss) {
			t_map_sats gnss_sats = GNSS_SATS();
			t_map_sats::const_iterator itGNS;
			// loop over all systems
			for (itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS) {
				string gs = t_gsys::gsys2str(itGNS->first);
				if (xcl.find("-" + gs) == xcl.end()) tmp.insert(gs);
			}
		}
		_gmutex.unlock(); return tmp;
	}


	// Return crd of site
	// ----------
	vector<double> t_gsetgen::crd(string site)
	{
		_gmutex.lock();
		vector<double> tmp;
		xml_node tmp_set = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child("crd");
		string tpm_str = tmp_set.attribute(site.c_str()).value();
		if (trim(tpm_str) != "")
		{
			stringstream temp(tpm_str);
			double x, y, z;
			temp >> x >> y >> z;
			tmp.push_back(x);
			tmp.push_back(y);
			tmp.push_back(z);
		}
		_gmutex.unlock(); return tmp;
	}

	// Return reference clk set
	// ------------------------
	string t_gsetgen::refsat()
	{
		_gmutex.lock();
	
		set<string> src = t_gsetbase::_setval(XMLKEY_GEN, "refsat");
		if (src.empty()) {
			string tmp("");
	        _gmutex.unlock();
			return tmp; 
		}
		else { _gmutex.unlock(); return *src.begin(); }
	}

	set<string> t_gsetgen::recs()
	{
		_gmutex.lock();

		set<string> tmp = t_gsetbase::_setvals(XMLKEY_GEN, "rec");


		_gmutex.unlock(); return tmp;
	}

    set<string> t_gsetgen::rec_all()
	{
		_gmutex.lock();

		set<string> tmp = t_gsetbase::_setvals(XMLKEY_GEN, "rec");

		_gmutex.unlock(); return tmp;
	}

	vector<string> t_gsetgen::list_base()
	{
		_gmutex.lock();
		vector<string> vals;
		string word;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value("base"));
		_gmutex.unlock();

		while (is >> word) {
			transform(word.begin(), word.end(), word.begin(), ::toupper);
			vals.push_back(word);
		}
		return vals;
	}

	vector<string> t_gsetgen::list_rover()
	{
		_gmutex.lock();
		vector<string> vals;
		string word;
		istringstream is(_doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value("rover"));

		while (is >> word) {
			transform(word.begin(), word.end(), word.begin(), ::toupper);
			vals.push_back(word);
		}
		_gmutex.unlock(); return vals;
	}

	// Return reference clk set
	// ------------------------
	string t_gsetgen::refsite()
	{
		_gmutex.lock();
	
		set<string> src = t_gsetbase::_setval(XMLKEY_GEN, "refsite");
		if (src.empty()) {
			string tmp("");
	        _gmutex.unlock();
			return tmp;
		}
		else { _gmutex.unlock(); return *src.begin(); }
	}

	// Return reference clk sigma
	double t_gsetgen::sig_refclk()
	{
	
		string child = "refsite";
		string refsite = this->refsite();
		if (refsite.empty()) {
			child = "refsat";
		}
	
		double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child(child.c_str()).attribute("sigclk").as_double(0.0);
		return tmp;
	}

	// Return estimator model
	// ------------------------
	string t_gsetgen::estimator()
	{
		_gmutex.lock();

		set<string> src = t_gsetbase::_setval(XMLKEY_GEN, "est");
		_gmutex.unlock();

		if (src.empty() || ((*src.begin() != "LSQ" && *src.begin() != "FLT")))//lvhb changed in 20200724
		{
			string tmp("LSQ");
			return tmp;
		}
		else
		{

			return *src.begin();
		}
	}

	// Return satellites which no used in calculation
    // ------------------------
	set<string> t_gsetgen::sat_rm()
	{
		_gmutex.lock();

		set<string> tmp = t_gsetbase::_setval(XMLKEY_GEN, "sat_rm");
		_gmutex.unlock(); return tmp;
	}

	// settings check
	// ----------
	void t_gsetgen::check()
	{
		_gmutex.lock();

		// check existence of nodes/attributes
		xml_node parent = _doc.child(XMLKEY_ROOT);
		xml_node node = _default_node(parent, XMLKEY_GEN);

		_default_node(node, "BEG", "");                         // all!
		_default_node(node, "END", "");	                      // all!
		_default_node(node, "SYS", "");	                      // all!
		_default_node(node, "REC", "");	                      // none 
		_default_node(node, "INT", int2str(DEF_SAMPLING).c_str());  // default

	  //  xml_node attr_int = _default_node(node, "INT", int2str(DEF_SAMPLING).c_str());  // default
	  //  _default_attr(attr_int, "unit", bool(DEF_SAMPUNIT));             // default

		// TO CHECK USER-SETUP CONSISTENCY
		_gmutex.unlock();
		if (floor(sampling()) < 1 || _dec > 0) {

			if (sampling() < 0.0) {
				_default_node(node, "INT", "0.0", true); // reset!
				if (_log) _log->comment(0, "gsetgen", "Warning: sampling rate settings negative: reset to 0");
				cerr << "gsetgen - Warning: sampling rate settings negative: reset to 0\n";
			}
			else {
				if (_log) _log->comment(1, "gsetgen: sampling rate settings above 1Hz recognized");
			}
		}

		return;
	}


	// help body
	// ----------
	void t_gsetgen::help()
	{
		_gmutex.lock();
		t_gtime beg(t_gtime::GPS); beg = beg - beg.sod();
		t_gtime end(t_gtime::GPS); end = beg + 86399;

		cerr << "\n <gen>\n"
			<< "   <beg> \"" << beg.str_ymdhms() << "\" </beg>\n"  // FIRST_TIME.str("\"%Y-%m-%d %H:%M:%S\"")
			<< "   <end> \"" << end.str_ymdhms() << "\" </end>\n"; // LAST_TIME.str("\"%Y-%m-%d %H:%M:%S\"")

		if (_gnss)
			cerr << "   <sys> " << _sys << " </sys>\n"; // GNSS systems

		cerr << "   <rec> GOPE WTZR POTS                </rec>\n"  // list of site identificators
			<< "   <int>" + int2str(DEF_SAMPLING) + "</int>\n"
			<< " </gen>\n";

		cerr << "\t<!-- general description:\n"
			<< "\t beg    .. beg time          (default: all)\n"
			<< "\t end    .. end time          (default: all)\n"
			<< "\t int    .. data sampling     (default: 30s)\n";

		if (_gnss)
			cerr << "\t sys    .. GNSS system(s)    (default: all)\n";

		cerr << "\t rec    .. GNSS receiver(s)  (rec active list, e.g.: GOPE ONSA WTZR ... )\n"
			<< "\t -->\n\n";

		_gmutex.unlock(); return;
	}

} // namespace
