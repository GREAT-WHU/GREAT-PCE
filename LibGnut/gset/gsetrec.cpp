/*
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file        gsetrec.cpp
* @brief       implements receiver object setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetrec.h"
#include "gutils/gsysconv.h"

using namespace std;
using namespace pugi;

#define XMLKEY_REC_REC "rec"
namespace gnut {

	// Constructor
	// ----------
	t_gsetrec::t_gsetrec()
		: t_gsetbase()
	{
		_set.insert(XMLKEY_REC);
		_id = "";
		_name = "";
		_desc = "";
		_domes = "";
		_X = 0.0;
		_Y = 0.0;
		_Z = 0.0;
		_dX = 0.0;
		_dY = 0.0;
		_dZ = 0.0;
		_dE = 0.0;
		_dN = 0.0;
		_dU = 0.0;
		_ZTD = 0.0;
		_ZTD_sig = 0.0;
		_CRD_sig = "";
		_rec = "";
		_ant = "";
		_beg = FIRST_TIME;
		_end = LAST_TIME;
		_overwrite = false;
	}


	// Destructor
	// ----------
	t_gsetrec::~t_gsetrec()
	{}

	// Return value
	// ----------
	string t_gsetrec::rec(string s)
	{
		_gmutex.lock();

		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());
		string tmp = site.attribute("rec").value();

		_gmutex.unlock();
		return tmp;
	}


	// Return value
	// ----------
	string t_gsetrec::ant(string s)
	{
		_gmutex.lock();

		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());
		string tmp = site.attribute("ant").value();

		_gmutex.unlock(); return tmp;
	}

	// Return value
	// ----------
	shared_ptr<t_grec> t_gsetrec::grec(string s, t_glog* glog)
	{
		_gmutex.lock();

		string  str;
		shared_ptr<t_grec>  tmp = make_shared<t_grec>();

		// return if no rec id found
		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());
		if (!site)
		{
			_gmutex.unlock();
			return 0;
		}

		t_gtime begdef(FIRST_TIME);
		t_gtime enddef(LAST_TIME);
		t_gtime beg(begdef), end(enddef);
		t_gtriple xyz = _get_crd_xyz(s);
		t_gtriple blh = _get_crd_blh(s);
		t_gtriple eccneu = _get_ecc_neu(s);

		if (double_eq(xyz.crd(0), 0.0) &&
			double_eq(xyz.crd(1), 0.0) &&
			double_eq(xyz.crd(2), 0.0))
		{
			if (double_eq(blh.crd(2), HSL_UNKNOWN))
			{
				if (_log) _log->comment(0, "gsetrec", "warning - invalid coordinates for rec: " + s);
				else               cerr << "gsetrec:  warning - invalid coordinates for rec: " + s << endl;
			}
			else
			{
				ell2xyz(blh.crd_array(), xyz.crd_array(), true);
			}
		}

		string str_beg_gen = _doc.child(XMLKEY_ROOT).child("gen").child_value("beg");
		substitute(str_beg_gen, "\"", "");
		t_gtime gen_beg(t_gtime::GPS);
		if (!str_beg_gen.empty()) gen_beg.from_str("%Y-%m-%d %H:%M:%S", trim(str_beg_gen));
		else gen_beg = FIRST_TIME;

		string str_end_gen = _doc.child(XMLKEY_ROOT).child("gen").child_value("end");
		substitute(str_end_gen, "\"", "");
		t_gtime gen_end(t_gtime::GPS);
		if (!str_end_gen.empty()) gen_end.from_str("%Y-%m-%d %H:%M:%S", trim(str_end_gen));
		else gen_end = FIRST_TIME;

		tmp->glog(glog);
		tmp->overwrite(site.attribute("overwrite").as_bool());
		tmp->id(site.attribute("id").value());

		t_gtriple std(10, 10, 10);
		tmp->crd(xyz, std, beg, end);
		tmp->eccneu(eccneu, beg, end);
		tmp->rec(site.attribute("rec").value(), beg, end);
		tmp->ant(site.attribute("ant").value(), beg, end);

		string x = site.attribute("X").value();
		string y = site.attribute("Y").value();
		string z = site.attribute("Z").value();
		if (!x.empty() && !y.empty() && !z.empty()) {
			xyz[0] = str2dbl(x);
			xyz[1] = str2dbl(y);
			xyz[2] = str2dbl(z);
			t_gtriple std(10, 10, 10);
			tmp->crd(xyz, std, beg, end);
		}

		string dx = site.attribute("DX").value();
		string dy = site.attribute("DY").value();
		string dz = site.attribute("DZ").value();
		if (!dx.empty() && !dy.empty() && !dz.empty()) {
			xyz[0] = str2dbl(dx);
			xyz[1] = str2dbl(dy);
			xyz[2] = str2dbl(dz);
			tmp->eccxyz(xyz, beg, end);
		}

		string dn = site.attribute("DN").value();
		string de = site.attribute("DE").value();
		string du = site.attribute("DU").value();
		if (!dx.empty() && !dy.empty() && !dz.empty()) {
			xyz[0] = str2dbl(dn);
			xyz[1] = str2dbl(de);
			xyz[2] = str2dbl(du);
			tmp->eccneu(xyz, beg, end);
		}

		_gmutex.unlock(); return tmp;
	}

	t_gobj* t_gsetrec::obj(string s)
	{

		_gmutex.lock();
		t_gobj* obj = nullptr;
		_gmutex.unlock();
		return obj;
	}

	// Get formats inputs
	// ----------
	set<string> t_gsetrec::objects()
	{
		return _objects();
	}

	// Get formats inputs
	// ----------
	set<string> t_gsetrec::_objects()
	{
		set<string> tmp;
		string str;

		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).first_child(); node; node = node.next_sibling()) {
			string name = node.name();
			if (name.compare(XMLKEY_REC_REC) == 0)
			{
				string str = node.attribute("id").as_string();
				if (!str.empty())
				{
					tmp.insert(str);
				}
			}
		}
		return tmp;
	}



	// Return value --> IN FUTURE OBSOLETE, BUT STILL USED IN GPPPFILTER (CHANGE TO _get_crd() below )
	// ----------
	int t_gsetrec::get_crd_xyz(t_gtriple& xyz, string s)
	{
		_gmutex.lock();

		xyz = _get_crd_xyz(s);

		if (double_eq(xyz.crd(0), 0.0) &&
			double_eq(xyz.crd(1), 0.0) &&
			double_eq(xyz.crd(2), 0.0)) {

			_gmutex.unlock(); return 0;
		}

		_gmutex.unlock(); return 1;
	}

	t_gtriple t_gsetrec::get_crd_xyz(string s)
	{
		return _get_crd_xyz(s);
	}
	t_gtriple t_gsetrec::get_std_xyz(string s)
	{
		return _get_std_xyz(s);
	}


	// Return set
	// ----------
	set<string> t_gsetrec::recs()
	{
		return _objects();
	}

	set<string> t_gsetrec::all_rec()
	{
		set<string> tmp;
		string temp;
		auto xml_node = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).children("rec");
		for (auto rec_node = xml_node.begin(); rec_node != xml_node.end(); rec_node++)
		{
			temp += rec_node->child_value();
		}
		istringstream is(temp);
		string word;
		while (is >> word) {
			transform(word.begin(), word.end(), word.begin(), ::toupper);
			tmp.insert(word);
		}

		return tmp;
	}

	// Return value
	// ----------
	t_gtriple t_gsetrec::_get_crd_xyz(string s)
	{
		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());

		t_gtriple xyz(site.attribute("X").as_double(),
			site.attribute("Y").as_double(),
			site.attribute("Z").as_double());

		return xyz;
	}

	t_gtriple t_gsetrec::_get_std_xyz(string s)
	{
		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());

		t_gtriple xyz(site.attribute("dX").as_double(),
			site.attribute("dY").as_double(),
			site.attribute("dZ").as_double());

		return xyz;
	}


	// Return value
	// ----------
	t_gtriple t_gsetrec::_get_ecc_neu(string s)
	{
		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());

		t_gtriple neu(site.attribute("DN").as_double(),
			site.attribute("DE").as_double(),
			site.attribute("DU").as_double());
		return neu;
	}

	// Return value
	// ----------
	t_gtriple t_gsetrec::_get_ecc_xyz(string s)
	{
		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());

		t_gtriple xyz(site.attribute("DX").as_double(),
			site.attribute("DY").as_double(),
			site.attribute("DZ").as_double());
		return xyz;
	}

	// Return value
	// ----------
	t_gtriple t_gsetrec::_get_crd_blh(string s)
	{
		t_gtriple zero(0.0, 0.0, HSL_UNKNOWN);

		xml_node site = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).find_child_by_attribute(XMLKEY_REC_REC, "id", s.c_str());
		xml_attribute attr;

		if ((attr = site.attribute("LAT")) &&
			(attr = site.attribute("LON")) &&
			(attr = site.attribute("HSL"))
			) {
			t_gtriple blh(site.attribute("LAT").as_double(),
				site.attribute("LON").as_double(),
				site.attribute("HSL").as_double());
			double pres, temp, undu = 0.0;
			_ggpt.gpt_v1(51544.0, blh[0] * D2R, blh[1] * D2R, blh[2], pres, temp, undu); // RAD
			blh[2] += undu; // ADD UNDULATION BECAUSE ONLY HSL (ABOVE MEAN SEE LEVEL) SUPPORTED
			return blh;
		}

		t_gtriple xyz = _get_crd_xyz(s);

		if (!double_eq(xyz.crd(0), 0.0) &&
			!double_eq(xyz.crd(1), 0.0) &&
			!double_eq(xyz.crd(2), 0.0)
			) {
			t_gtriple blh;
			if (xyz2ell(xyz, blh, true) == 1) return blh;
			else                             return zero;
		}

		return zero;
	}

	// settings check
	// ----------
	void t_gsetrec::check()
	{
		_gmutex.lock();

		// check existence of nodes/attributes
		xml_node parent = _doc.child(XMLKEY_ROOT);
		xml_node node = _default_node(parent, XMLKEY_REC);

		// check existence of attributes
		_default_attr(node, "X", _X);
		_default_attr(node, "Y", _Y);
		_default_attr(node, "Z", _Z);
		_default_attr(node, "DX", _dX);
		_default_attr(node, "DY", _dY);
		_default_attr(node, "DZ", _dZ);
		_default_attr(node, "DN", _dN);
		_default_attr(node, "DE", _dE);
		_default_attr(node, "DU", _dU);
		_default_attr(node, "overwrite", _overwrite);

		_default_attr(node, "ZTD", _ZTD);
		_default_attr(node, "ZTD_sig", _ZTD_sig);
		_default_attr(node, "CRD_sig", _CRD_sig);

		// check duplicities
		map<string, xml_node> mchild;
		for (xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC); node; )
		{
			xml_node rec = node;
			node = node.next_sibling(XMLKEY_REC);

			// check ID & NAME & DOMES
			string id = rec.attribute("id").value();
			string name = rec.attribute("name").value();
			string domes = rec.attribute("domes").value();

			if (name.empty()) {  // complete NAME
				rec.insert_attribute_after("name", rec.attribute("id")) = id.c_str();
			}
			if (domes.empty()) {  // complete NAME
				if (name.length() < 14) { domes = "XXXXXXXXX"; }
				if (name.length() >= 14) { domes = tail(name, 9); name = name.substr(0, name.size() - 10); }
				rec.insert_attribute_after("domes", rec.attribute("name")) = domes.c_str();
				rec.insert_attribute_after("name", rec.attribute("id")) = name.c_str();
			}

			if (mchild.find(id) == mchild.end()) {
				mchild[id] = rec;
			}
			else {
				_doc.child(XMLKEY_ROOT).remove_child(rec);
				if (_log) _log->comment(0, "gsetrec", "warning - removed duplicated record [id] for rec:" + id);
				else               cerr << "gsetrec:  warning - removed duplicated record [id] for rec:" + id << endl;
			}

			map<string, xml_node> mattr;
			for (xml_node set = rec.child("set"); set; ) {
				xml_node tmp = set;
				set = set.next_sibling("set");  

				string beg = tmp.attribute("beg").value();
				if (mattr.find(beg) == mattr.end()) {
					mattr[beg] = set;
				}
				else {
					rec.remove_child(tmp);
					if (_log) _log->comment(0, "gsetrec", "warning - removed duplicated record [set] for rec:" + id + ", beg:" + beg);
					else               cerr << "gsetrec:  warning - removed duplicated record [set] for rec:" + id + ", beg:" + beg << endl;
				}
			}
		}

		_gmutex.unlock(); return;
	}


	// help body
	// ----------
	void t_gsetrec::help()
	{
		_gmutex.lock();

		cerr << " <rec"
			<< " id=\"GOPE\""
			<< " name=\"GOPE\""
			<< " domes=\"11502M002\""
			<< " desc=\"Geodetic Observatory Pecny, Czech Republic\" >\n"
			<< "  <set"
			<< " beg=\"1995 05 13 00 00 00\""
			<< " end=\"1997 06 11 00 00 00\""
			<< " rec=\"TRIMBLE 4000SSE\""
			<< " ant=\"TRM14532.00     NONE\""
			<< "\n\t"
			<< " X=\"3979316.0\""
			<< " Y=\"1050312.0\""
			<< " Z=\"4857067.0\""
			<< " dX=\"0.0\""
			<< " dY=\"0.0\""
			<< " dZ=\"0.0\""
			<< " dN=\"0.0\""
			<< " dE=\"0.0\""
			<< " dU=\"0.0\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"1997 06 11 00 00 00\""
			<< " end=\"1997 06 20 00 00 00\""
			<< " rec=\"SPP GEOTRACER100\""
			<< " ant=\"TRM14532.00     NONE\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"1997 06 20 00 00 00\""
			<< " end=\"1999 11 04 00 00 00\""
			<< " rec=\"TRIMBLE 4000SSE\""
			<< " ant=\"TRM14532.00     NONE\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"1999 11 04 00 00 00\""
			<< " end=\"2000 07 24 00 00 00\""
			<< " rec=\"ASHTECH Z18\""
			<< " ant=\"ASH701073.3     SNOW\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"2000 07 24 00 00 00\""
			<< " end=\"2000 10 04 00 00 00\""
			<< " rec=\"TRIMBLE 4000SSE\""
			<< " ant=\"TRM14532.00     NONE\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"2000 10 04 00 00 00\""
			<< " end=\"2001 07 18 00 00 00\""
			<< " rec=\"ASHTECH Z18\""
			<< " ant=\"ASH701073.3     SNOW\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"2006 07 14 00 00 00\""
			<< " end=\"2009 12 14 00 00 00\""
			<< " rec=\"ASHTECH Z18\""
			<< " ant=\"TPSCR3_GGD      CONE\""
			<< "  />\n"
			<< "  <set"
			<< " beg=\"2009 12 14 00 00 00\""
			<< " end=\"2013 03 19 00 00 00\""
			<< " rec=\"TPS NETG3\""
			<< " ant=\"TPSCR.G3        TPSH\""
			<< "  />\n"
			<< " </rec>\n";

		cerr << "\t<!-- receiver description:\n"
			<< "\t id     .. ID name\n"
			<< "\t name   .. site name\n"
			<< "\t domes  .. site domes\n"
			<< "\t desc   .. description\n"
			<< "\t X[YZ]  .. X,Y,Z-coordinate [m]\n"
			<< "\t dX[YZ] .. X,Y,Z-eccentricity [m]\n"
			<< "\t -->\n\n";

		_gmutex.unlock(); return;
	}

} // namespace
