/**
*
* @verbatim
(c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
Ondrejov 244, 251 65, Czech Republic
@endverbatim
*
* @file		gsetinp.cpp
* @brief	implements input setting class
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut 
{

	IFMT t_gsetinp::str2ifmt(const string& s)
	{
		string tmp = s;
		transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
		if (tmp == "RINEXC") return RINEXC_INP;
		
		if (tmp == "RINEXO") return RINEXO_INP;
		if (tmp == "RNO") return RINEXO_INP;

		if (tmp == "RINEXN") return RINEXN_INP;
		if (tmp == "RNN") return RINEXN_INP;

		if (tmp == "SINEX")     return SINEX_INP;
		if (tmp == "SNX")       return SINEX_INP;
		if (tmp == "SP3") return SP3_INP;
		if (tmp == "ATX") return ATX_INP;

		if (tmp == "OTL") return BLQ_INP;
		if (tmp == "BLQ") return BLQ_INP;
		if (tmp == "STA") return STA_INP;

		if (tmp == "BIASINEX") return BIASINEX_INP;
		if (tmp == "BIABERN") return BIABERN_INP;
		if (tmp == "DCB")     return BIABERN_INP;

		if (tmp == "DE")  return DE_INP;

		if (tmp == "OCEANTIDE") return OCEANTIDE_INP;
		if (tmp == "OTD")       return OCEANTIDE_INP;


		if (tmp == "POLEUT1")return POLEUT1_INP;
		if (tmp == "ERP")    return POLEUT1_INP;

		if (tmp == "LEP")    return LEAPSECOND_INP;
 		if (tmp == "LEAPSECOND")return LEAPSECOND_INP;

		if (tmp == "SATPARS")return SATPARS_INP;
		if (tmp == "SAT")    return SATPARS_INP;
		
		if (tmp == "EPODIR") return EPODIR_INP;
		if (tmp == "UPD") return UPD_INP;

		return IFMT(-1);
	}


	// Convertor for INP formats
	// ----------
	string t_gsetinp::ifmt2str(const IFMT& f)
	{
		switch (f) {
		case RINEXO_INP:  return "RINEXO";
		case RINEXC_INP:  return "RINEXC";

		case RINEXN_INP:  return "RINEXN";
		case SINEX_INP:   return "SINEX";
		case SP3_INP:     return "SP3";

		case STA_INP:      return "STA";
		case ATX_INP:     return "ATX";
		case BLQ_INP:     return "BLQ";
		case POLEUT1_INP: return "ERP";
		case BIASINEX_INP:   return "BIASINEX";
		case BIABERN_INP:    return "BIABERN";
		case UPD_INP:        return "UPD";
		case LEAPSECOND_INP:  return "LEAPSECOND";
		case EPODIR_INP: return "EPODIR";
		case DE_INP:     return "DE";
		default:             return "UNDEF";
		}
		return "UNDEF";
	}

	// Constructor
	// ----------
	t_gsetinp::t_gsetinp()
		: t_gsetbase()
	{
		_set.insert(XMLKEY_INP);
		_chkNavig = true;
		_chkHealth = true;
		_corrStream = "";
	}


	// Destructor
	// ----------
	t_gsetinp::~t_gsetinp()
	{}


	// Get formats input size
	// ----------
	int t_gsetinp::input_size(const string& fmt)
	{
		_gmutex.lock();

		int tmp = _inputs(fmt).size();

		_gmutex.unlock(); return tmp;
	}

	bool t_gsetinp::check_input(const string& fmt)
	{
		_gmutex.lock();
		int tmp = _inputs(fmt).size();
		_gmutex.unlock();
		return (tmp > 0);
	}


	// Get formats inputs (all in multimap)
	// ----------
	multimap<IFMT, string> t_gsetinp::inputs_all()
	{
		_gmutex.lock();

		multimap<IFMT, string> map;

		set<string> ifmt = _iformats();
		set<string>::const_iterator itFMT = ifmt.begin();

		while (itFMT != ifmt.end())
		{
			string fmt = *itFMT;
			IFMT  ifmt = str2ifmt(fmt);
			vector<string> inputs = _inputs(fmt);//get file name in input node
			vector<string>::const_iterator itINP = inputs.begin();
			while (itINP != inputs.end()) 
			{
				map.insert(map.end(), pair<IFMT, string>(ifmt, *itINP));
				itINP++;
			}
			itFMT++;
		}
		_gmutex.unlock(); return map;
	}


	// Get formats inputs
	// ----------
	vector<string> t_gsetinp::inputs(const string& fmt)
	{
		IFMT ifmt = str2ifmt(fmt);
		return _inputs(ifmt);
	}

	vector<string> t_gsetinp::inputs(const IFMT & ifmt)
	{
		return _inputs(ifmt);
	}

	// Get formats inputs
	// ----------
	vector<string> t_gsetinp::_inputs(const string& fmt)
	{
		vector<string> tmp;
		set<string> list;
		string str;

		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			if (node.name() == fmt) 
			{
				istringstream is(node.child_value());
				while (is >> str && !is.fail())
				{
					if (str.find("://") == string::npos) str = GFILE_PREFIX + str;
					if (list.find(str) == list.end()) 
					{
						tmp.push_back(str);
						list.insert(str);
					}
					else 
					{
						if (_log) _log->comment(1, "gsetinp", "READ: " + str + " multiple request ignored");
					}
				}
			}
		}
		return tmp;
	}

	vector<string> t_gsetinp::_inputs(const IFMT &fmt)
	{
		vector<string> tmp;
		set<string> list;
		string str;

		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			IFMT ifmt = str2ifmt(node.name());
			if (ifmt == IFMT::UNDEF) continue;
			if (ifmt == fmt)
			{
				istringstream is(node.child_value());
				while (is >> str && !is.fail())
				{
					if (str.find("://") == string::npos)
						str = GFILE_PREFIX + str;
					if (list.find(str) == list.end())
					{
						tmp.push_back(str);
						list.insert(str);
					}
					else
					{
						if (_log)
							_log->comment(1, "gsetinp", "READ: " + str + " multiple request ignored");
					}
				}
			}
		}
		return tmp;
	}

	// Get input formats
	// ----------
	set<string> t_gsetinp::_iformats()
	{
		set<string> tmp;
		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling())
		{
			tmp.insert(node.name());
		}
		return tmp;
	}

	// settings check
	// ----------
	void t_gsetinp::check()
	{
		_gmutex.lock();

		// check existence of nodes/attributes
		xml_node parent = _doc.child(XMLKEY_ROOT);
		xml_node node = _default_node(parent, XMLKEY_INP);

		// check supported input formats (see IFMT enum !)
		set<string> ifmt = _iformats();
		set<string>::const_iterator itFMT = ifmt.begin();
		while (itFMT != ifmt.end()) {
			string fmt = *itFMT;
			IFMT  ifmt = str2ifmt(fmt);
			if (ifmt < 0) {
				_doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
				if (_log) _log->comment(0, "Warning: " + fmt + " inp format not implemented [gsetinp::check()]!");
				itFMT++;
				continue;
			}
			// check application-specific output format
			if (_IFMT_supported.find(ifmt) == _IFMT_supported.end()) {
				_doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child(node.child(fmt.c_str()));
				if (_log) _log->comment(0, "Warning: " + fmt + " inp format not supported by this application!");
				else                cerr << "Warning: " + fmt + " inp format not supported by this application!\n";
			}
			itFMT++;
		}

		_default_attr(node, "chk_nav", _chkNavig);
		_default_attr(node, "chk_health", _chkHealth);

		xml_node nodeBNCRTCM = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child("bncrtcm");
		_default_attr(nodeBNCRTCM, "_corrStream", _corrStream);

		_gmutex.unlock(); return;
	}


	// settings help
	// ----------
	void t_gsetinp::help()
	{
		_gmutex.lock();


		cerr << " <inputs>\n"
			<< "   <rinexo> file://dir/name </rinexo> \t\t <!-- obs RINEX decoder -->\n"
			<< "   <rinexn> file://dir/name </rinexn> \t\t <!-- nav RINEX decoder -->\n"
			<< " </inputs>\n";

		cerr << "\t<!-- inputs description:\n"
			<< "\t <decoder> path1 path2 path3  </decoder>\n"
			<< "\t ... \n"
			<< "\t where path(i) contains [file,tcp,ntrip]:// depending on the application\n"
			<< "\t -->\n\n";

		_gmutex.unlock(); return;
	}

} // namespace
