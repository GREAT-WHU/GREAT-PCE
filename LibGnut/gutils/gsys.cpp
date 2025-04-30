
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gutils/gsys.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut
{

	// ---------------------------------------------------------------------------------
	// class GSYS
	// ---------------------------------------------------------------------------------

	// constructor
	// ---------
	t_gsys::t_gsys(GSYS sys) : _gsys(sys)
	{}

	// constructor
	// ---------
	t_gsys::t_gsys(string sys)
		: _gsys(GNS)
	{
		if (sys.size() > 0) _gsys = str2gsys(sys);
	}

	// constructor
	// ---------
	t_gsys::t_gsys(char c)
	{
		_gsys = char2gsys(c);
	}


	// --------------------------------------------------------
	// STATIC FUNCTIONS
	// --------------------------------------------------------
	// get system band
	// ----------
	GOBSBAND t_gsys::gfrq2band(GSYS gs, GFRQ gfreq)
	{
		//t_map_freq m = GNSS_FREQ_PRIORITY;
		if (gs == GNS || gfreq == LAST_GFRQ) return BAND;
		try
		{
			t_vec_freq sys_freq = GNSS_FREQ_PRIORITY.at(gs);
			int size = sys_freq.size();
			for (size_t i = 0; i < size; ++i)
			{
				if (sys_freq.at(i) == gfreq) return band_priority(gs, FREQ_SEQ(i));
			}
			return BAND;
		}
		catch (...)
		{
			return BAND;
		}
	}

	// get system freq
	// ----------
	GFRQ t_gsys::band2gfrq(GSYS gs, GOBSBAND gband)
	{

		if (gs == GNS || gband == BAND) return LAST_GFRQ;
		try
		{
			t_vec_band sys_band = GNSS_BAND_PRIORITY.at(gs);
			int size = sys_band.size();
			for (size_t i = 0; i < size; ++i)
			{
				if (sys_band.at(i) == gband) return freq_priority(gs, FREQ_SEQ(i));
			}
			return LAST_GFRQ;
		}
		catch (...)
		{
			return LAST_GFRQ;
		}	
	}

	// get system sequence
	// ----------
	FREQ_SEQ t_gsys::gfrq2freq(GSYS gs, GFRQ gfreq)
	{

		if (gs == GNS || gfreq == LAST_GFRQ) return FREQ_X;
		try
		{
			t_vec_freq sys_freq = GNSS_FREQ_PRIORITY.at(gs);
			int size = sys_freq.size();
			for (size_t i = 0; i < size; ++i) 
			{
				if (sys_freq.at(i) == gfreq) return FREQ_SEQ(i);
			}
			return FREQ_X;
		}
		catch (...)
		{
			return FREQ_X;
		}
	}

	// get system sequence
	// ----------
	GSYS t_gsys::gfrq2gsys(GFRQ freq)
	{
		string s = t_gfreq::gfreq2str(freq);
		return t_gsys::char2gsys(s[0]);
	}

	FREQ_SEQ t_gsys::str2freq(string freq)
	{
		if (freq.empty()) return FREQ_SEQ::FREQ_X;

		char cfreq = freq[0];
		switch (cfreq)
		{
		case '0': return  FREQ_SEQ::FREQ_X;
		case '1': return FREQ_SEQ::FREQ_1;
		case '2':return FREQ_SEQ::FREQ_2;
		case '3':return FREQ_SEQ::FREQ_3;
		case '4':return FREQ_SEQ::FREQ_4;
		case '5':return FREQ_SEQ::FREQ_5;
		case '6':return FREQ_SEQ::FREQ_6;
		case '7':return FREQ_SEQ::FREQ_7;
		default:return FREQ_SEQ::FREQ_X;
		}
	}

	string t_gsys::freq2str(FREQ_SEQ freq)
	{
		switch (freq)
		{
		case FREQ_SEQ::FREQ_X: return string("0");
		case FREQ_SEQ::FREQ_1: return string("1");
		case FREQ_SEQ::FREQ_2: return string("2");
		case FREQ_SEQ::FREQ_3: return string("3");
		case FREQ_SEQ::FREQ_4: return string("4");
		case FREQ_SEQ::FREQ_5: return string("5");
		case FREQ_SEQ::FREQ_6: return string("6");
		case FREQ_SEQ::FREQ_7: return string("7");

		default:return string("0");
		}
	}


	// get system freq sequence
	// ----------
	FREQ_SEQ t_gsys::band2freq(GSYS gs, GOBSBAND gband)
	{

		if (gs == GNS || gband == BAND) return FREQ_X;
		try
		{
			t_vec_band sys_band = GNSS_BAND_PRIORITY.at(gs);
			int size = sys_band.size();
			for (size_t i = 0; i < size; ++i)
			{
				if (sys_band.at(i) == gband) return FREQ_SEQ(i);
			}
			return FREQ_X;
		}
		catch (...)
		{
			return FREQ_X;
		}
	}


	// get band selection
	// ----------
	GOBSBAND t_gsys::band_priority(GSYS gs, FREQ_SEQ iseq)  // FREQ_SEQ priority band
	{

		try
		{
			if (gs == GNS || iseq >= GNSS_BAND_PRIORITY.at(gs).size()) return BAND;
			return GNSS_BAND_PRIORITY.at(gs).at(iseq);
		}
		catch (...)
		{
			return BAND;
		}
	}

	GOBS t_gsys::gobs_priority(GSYS gs, GOBS gobs)
	{
		if (gs == GPS)
		{
			switch (gobs)
			{
			case C1C:
			case C1P:
			case C1W:
			{
				return C1W;
			}
			case C2C:
			case C2L:
			case C2X:
			case C2P:
			case C2W:
			{
				return C2W;
			}
			case C5Q:
			case C5X:
			{
				return C5X;
			}
			default:  return GOBS::X;
			}
		}
		else if (gs == GLO)
		{
			switch (gobs)
			{
			case C1C:
			case C1P: return C1P;
			case C2C:
			case C2P: return C2P;
			}
		}
		else
		{
			return gobs;
		}
	}

	GOBS t_gsys::gobs_priority(GSYS gs, GOBSBAND band)
	{
		if (gs == GPS)
		{
			switch (band)
			{
			case BAND_1:
			{
				return C1W;
			}
			case BAND_2:
			{
				return C2W;
			}
			case BAND_5:
			{
				return C5X;
			}
			default:  return GOBS::X;
			}
		}
		else if (gs == GLO)
		{
			switch (band)
			{
			case BAND_1: return C1P;
			case BAND_2: return C2P;
			default:  return GOBS::X;
			}
		}
		else if (gs = GAL)
		{
			switch (band)
			{
			case BAND_1: return C1X;
			case BAND_5: return C5X;
			case BAND_7: return C7X;
			default:  return GOBS::X;
			}
		}
		else
		{
			return GOBS::X;
		}
	}

	GOBS t_gsys::gobs_defaults(GSYS gs, GOBS  gobs1, GOBSBAND band)
	{
		if (gs == GPS)
		{
			switch (band)
			{
			case BAND_1:
			{
				return C1W;
			}
			case BAND_2:
			{
				return C2W;
			}
			case BAND_5:
			{
				return C5X;
			}
			default:  return GOBS::X;
			}
		}
		else if (gs == GLO)
		{
			switch (band)
			{
			case BAND_1:  return C1P;
			case BAND_2:  return C2P;
			default:  return GOBS::X;
			}
		}
		else if (gs == GAL)
		{
			string strGOBS = gobs2str(gobs1);
			switch (band)
			{
			case BAND_1:
			{
				if (strGOBS[strGOBS.size() - 1] == 'X')
				{
					return C1X;
				}
				else
				{
					return C1C;
				}
			}
			case BAND_5:
			{
				if (strGOBS[strGOBS.size() - 1] == 'X')
				{
					return C5X;
				}
				else
				{
					return C5Q;
				}
			}
			case BAND_7:
			{
				if (strGOBS[strGOBS.size() - 1] == 'X')
				{
					return C7X;
				}
				else
				{
					return C7Q;
				}
			}
			default:  return GOBS::X;
			}
		}
		else
		{
			return GOBS::X;
		}
	}

	// get freq selection
	// ----------
	GFRQ t_gsys::freq_priority(GSYS gs, FREQ_SEQ iseq)  // FREQ_SEQ priority frequency
	{

		try
		{
			if (gs == GNS || iseq >= GNSS_FREQ_PRIORITY.at(gs).size()) return LAST_GFRQ;

			return GNSS_FREQ_PRIORITY.at(gs).at(iseq);
		}
		catch (...)
		{
			return LAST_GFRQ;
		}
	}

	// get attr selection
	// ----------
	GOBSATTR t_gsys::attr_priority(GSYS gs, GOBSBAND gb, GOBSTYPE gt, unsigned int iseq)  // iseq priority sequence
	{
		t_map_gnss m = GNSS_DATA_PRIORITY();

		if (gs == GNS ||
			gt == TYPE ||
			gb == BAND || iseq > m[gs][gb][gt].size()) return ATTR;

		return m[gs][gb][gt][iseq];
	}


	// convert GSYS enum to GSYS string
	// ---------- 
	string t_gsys::gsys2str(GSYS sys)
	{
		gtrace("t_gsys::gsys2str");

		// PROBLEM WITH STATIC FUNCTION
		//  boost::mutex::scoped_lock lock(_mutex);

		switch (sys)
		{
		case GPS: return "GPS";
		case GLO: return "GLO";
		case GAL: return "GAL";
		case BDS: return "BDS";
		case SBS: return "SBS";
		case QZS: return "QZS";
		case IRN: return "IRN";
		case LEO: return "LEO";
		case GNS: return "GNS";

		default: { cout << "*** warning: unknown GNSS system!\n"; cout.flush(); }
		}

		return "GNS";
	}


	// convert GSYS enum to GSYS char
	// ---------- 
	char t_gsys::gsys2char(GSYS sys)
	{
		gtrace("t_gsys::gsys2char");

		// PROBLEM WITH STATIC FUNCTION
		//  boost::mutex::scoped_lock lock(_mutex);

		switch (sys) {
		case GPS: return 'G';
		case GLO: return 'R';
		case GAL: return 'E';
		case BDS: return 'C';
		case SBS: return 'S';
		case QZS: return 'J';
		case IRN: return 'I';
		case GNS: return 'X';

		default: { cout << "*** warning: unknown GNSS system \n"; cout.flush(); }
		}

		return 'X';
	}
	// convert GSYS string to GSYS enum
	// ---------- 
	string t_gsys::str2strfirst(string s)
	{
		gtrace("t_gsys::str2gsys [" + s + "]");

		// PROBLEM WITH STATIC FUNCTION
		//  boost::mutex::scoped_lock lock(_mutex);

		if (s.size() == 0) {
			cout << "*** warning: not defined GNSS system code [NULL]\n"; cout.flush();
			return "X";
		}

		transform(s.begin(), s.end(), s.begin(), ::toupper);

		if (s == "G" || s == "GPS" || s == "NAVSTAR") return "G";
		else if (s == "R" || s == "GLO" || s == "GLONASS") return "R";
		else if (s == "E" || s == "GAL" || s == "GALILEO") return "E";
		else if (s == "C" || s == "COMP" || s == "COMPASS") return "C";
		else if (s == "C" || s == "BDS" || s == "BEIDOU") return "C";
		else if (s == "S" || s == "SBS" || s == "EGNOS") return "S";
		else if (s == "S" || s == "SBAS") return "S";
		else if (s == "J" || s == "QZS" || s == "JAXA") return "J";
		else if (s == "J" || s == "QZSS") return "J";
		else if (s == "I" || s == "IRN" || s == "IRNSS") return "I";
		else { cout << "*** warning: not defined GNSS system code [" << s[0] << "]\n"; cout.flush(); }

		return "X";
	}

	// convert GSYS string to GSYS enum
	// ---------- 
	GSYS t_gsys::str2gsys(string s)
	{
		gtrace("t_gsys::str2gsys [" + s + "]");

		// PROBLEM WITH STATIC FUNCTION
		//  boost::mutex::scoped_lock lock(_mutex);

		if (s.size() == 0) {
			cout << "*** warning: not defined GNSS system code [NULL]\n"; cout.flush();
			return GNS;
		}

		transform(s.begin(), s.end(), s.begin(), ::toupper);

		if (s == "G" || s == "GPS" || s == "NAVSTAR") return GPS;
		else if (s.substr(0, 1) == "2" || s.substr(0, 1) == "3" || s.substr(0, 1) == "4") return GPS;
		else if (s == "R" || s == "GLO" || s == "GLONASS") return GLO;
		else if (s == "E" || s == "GAL" || s == "GALILEO") return GAL;
		else if (s == "C" || s == "COMP" || s == "COMPASS") return BDS;
		else if (s == "C" || s == "BDS" || s == "BEIDOU") return BDS;
		else if (s == "S" || s == "SBS" || s == "EGNOS") return SBS;
		else if (s == "S" || s == "SBAS") return SBS;
		else if (s == "J" || s == "QZS" || s == "JAXA") return QZS;
		else if (s == "J" || s == "QZSS") return QZS;
		else if (s == "I" || s == "IRN" || s == "IRNSS") return IRN;
		else if (s == "L" || s == "LEO") return LEO;
		else if (s == "X" || s == "GNS" || s == "GNSS") return GNS;
		else if (s == "M") return GNS;
		else { cout << "*** warning: not defined GNSS system code [" << s[0] << "]\n"; cout.flush(); }

		return GNS;
	}

	GSYS t_gsys::sat2gsys(const string& s)
	{
		string tmp = s.substr(0, 1);
		if (s.empty()) return GNS;
		else return str2gsys(tmp);
	}

	string t_gsys::sat2cgsys(const string& s)
	{
		char sys = s[0];
		switch (sys) {
		case 'G': return "GPS";
		case 'E':return "GAL";
		case 'C':return "BDS";
		case 'R':return "GLO";
		case 'J':return "QZS";
		case 'S':return "SBS";
		case 'I':return "IRN";
		case 'L':return "LEO";
		case 'X':
		case 'M': return "GNS";
		default: return "";
		}
		return "";
	}


	// convert GSYS string to GSYS enum
	// ---------- 
	char t_gsys::str2char(string s)
	{
		if (s.size() > 0) return gsys2char(str2gsys(s));
		return 'X';
	}


	// convert GSYS char to GSYS enum
	// ---------- 
	GSYS t_gsys::char2gsys(char c)
	{
		gtrace("t_gsys::char2gsys");

		if (c == 'G') return GPS;
		else if (c == 'R') return GLO;
		else if (c == 'E') return GAL;
		else if (c == 'C') return BDS;
		else if (c == 'S') return SBS;
		else if (c == 'J') return QZS;
		else if (c == 'I') return IRN;
		else if (c == 'M') return GNS;
		else if (c == 'X') return GNS;
		else
		{
			cout << "*** warning: not defined GNSS system char [" << c << "]\n"; cout.flush();
		}

		return GNS;
	}


	// convert GSYS char to GSYS enum
	// ---------- 
	string t_gsys::char2str(char c)
	{
		return gsys2str(char2gsys(c));
	}


	// convert satellite name
	// ---------- 
	string t_gsys::eval_sat(string sat, GSYS sys)
	{
		istringstream is(sat); // .substr(0,3) NE!
		GSYS gnss = sys;
		int  svn = 0;
		char chr = 'G';        // if empty, use GPS
		size_t l = is.str().length();

		if (l == 0) return "X00";

		// Slow!!
		//if( l < 3 || sat[0] == ' ' )  is        >> svn;
		//else                          is >> chr >> svn;

		// Change by ZhengHJ
		if (l < 3 || sat[0] == ' ')  sscanf(sat.c_str(), "%d", &svn);
		else                         sscanf(sat.c_str(), "%c%d", &chr, &svn);


		if (is.fail()) { return "X00"; }
		if (chr != 'G') gnss = char2gsys(chr);

		switch (gnss) {
		case GPS: chr = 'G'; break;
		case GLO: chr = 'R'; break;
		case GAL: chr = 'E'; break;
		case BDS: chr = 'C'; break;
		case SBS: chr = 'S'; break;
		case QZS: chr = 'J'; break;
		case IRN: chr = 'I'; break;
		case GNS: chr = 'X'; break;
		default: chr = 'X'; break;
		}

		if (svn > QZS_OFFSET) { svn -= QZS_OFFSET; }
		if (svn > SBS_OFFSET) { svn -= SBS_OFFSET; }

		char tmp[4];
		sprintf(tmp, "%c%02d", chr, svn);

		//return os.str();
		string tmp_svn(tmp);
		return tmp_svn;
	}


	// convert satellite name
	// ---------- 
	string t_gsys::eval_sat(int svn, GSYS sys)
	{
		char chr;
		switch (sys) {
		case GPS: chr = 'G'; break;
		case GLO: chr = 'R'; break;
		case GAL: chr = 'E'; break;
		case BDS: chr = 'C'; break;
		case SBS: chr = 'S'; break;
		case QZS: chr = 'J'; break;
		case IRN: chr = 'I'; break;
		case GNS: chr = 'X'; break;
		default: chr = 'X'; break;
		}
		if (svn > QZS_OFFSET) { svn -= QZS_OFFSET; }
		if (svn > SBS_OFFSET) { svn -= SBS_OFFSET; }

		ostringstream os;
		os << setw(1) << chr << setfill('0') << setw(2) << svn;

		return os.str();
	}

	// set GSYS from string
	// ---------- 
	void t_gsys::from_string(string sys)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gsys = str2gsys(sys);
	}


	// set GSYS from enum
	// ---------- 
	void t_gsys::from_gsys(GSYS sys)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gsys = sys;
	}


	// overloaded equivalence operator
	// ----------
	bool t_gsys::operator==(const string& sys) const
	{
		if (_gsys == str2gsys(sys)) return true;
		return false;
	}


	// overloaded equivalence operator
	// ----------
	bool t_gsys::operator==(const GSYS& sys) const
	{
		//  if( _gsys == sys.gsys() ) return true;
		return false;
	}

	// ---------------------------------------------------------------------------------
	// class FREQ
	// ---------------------------------------------------------------------------------

	// get GFRQ enum from string
	// ----------
	GFRQ t_gfreq::str2gfreq(string freq)
	{
		gtrace("t_gsys::str2gfreq");

		transform(freq.begin(), freq.end(), freq.begin(), ::toupper);

		if (freq == "G01") return G01;
		else if (freq == "G02") return G02;
		else if (freq == "G05") return G05;
		else if (freq == "R01") return R01;
		else if (freq == "R02") return R02;
		else if (freq == "R03") return R03_CDMA;
		else if (freq == "R05") return R05_CDMA;
		else if (freq == "E01") return E01;
		else if (freq == "E05") return E05;
		else if (freq == "E07") return E07;
		else if (freq == "E08") return E08;
		else if (freq == "E06") return E06;
		else if (freq == "C02") return C02;
		else if (freq == "C07") return C07;
		else if (freq == "C06") return C06;
		else if (freq == "C01") return C01;
		else if (freq == "C08") return C08;
		else if (freq == "C05") return C05;
		else if (freq == "C09") return C09;
		else if (freq == "J01") return J01;
		else if (freq == "J02") return J02;
		else if (freq == "J05") return J05;
		else if (freq == "J06") return J06;
		else if (freq == "S01") return S01;
		else if (freq == "S05") return S05;
		else if (freq == "I05") return I05;
		else if (freq == "LAST_GFRQ") return LAST_GFRQ;

		return LAST_GFRQ;
	}


	// get string from GFRQ enum
	// ---------- 
	string t_gfreq::gfreq2str(GFRQ freq)
	{
		gtrace("t_gsys::gfreq2str");

		switch (freq) {
		case G01: return "G01";
		case G02: return "G02";
		case G05: return "G05";
		case R01: return "R01";
		case R02: return "R02";
		case R01_CDMA: return "R01";
		case R02_CDMA: return "R02";
		case R03_CDMA: return "R03";
		case R05_CDMA: return "R05";
		case E01: return "E01";
		case E05: return "E05";
		case E07: return "E07";
		case E08: return "E08";
		case E06: return "E06";
		case C02: return "C02";
		case C06: return "C06";
		case C07: return "C07";
		case C01: return "C01";
		case C08: return "C08";
		case C05: return "C05";
		case C09: return "C09";
		case J01: return "J01";
		case J02: return "J02";
		case J05: return "J05";
		case J06: return "J06";
		case S01: return "S01";
		case S05: return "S05";
		case I05: return "I05";      
		case LAST_GFRQ: return "LAST_GFRQ";

		default: {  cout << "*** warning: unknown frequency code \n"; cout.flush(); }
		}

		return "LAST_GFRQ";
	}

	// get true if BDS GEO satellite
	// ------------------------------
	bool t_gsys::bds_geo(const string& sat)
	{
		set<string> geo;                                  // prns of geost. satellites
		geo.insert("C01");
		geo.insert("C02");
		geo.insert("C03");
		geo.insert("C04");
		geo.insert("C05");
		geo.insert("C59");
		geo.insert("C60");
		if (geo.find(sat) != geo.end()) return true;

		return false;
	}

	// get BDS satellite following CAST CYS attitude law
	// ------------------------------
	bool t_gsys::bds_cast(const string& sat)
	{
		set<string> cys;
		cys.insert("C13");
		cys.insert("C14");
		cys.insert("C06");
		cys.insert("C16");
		cys.insert("C19");
		cys.insert("C20");
		cys.insert("C21");
		cys.insert("C22");
		cys.insert("C23");
		cys.insert("C24");
		cys.insert("C32");
		cys.insert("C33");
		cys.insert("C36");
		cys.insert("C37");
		cys.insert("C41");
		cys.insert("C42");
		cys.insert("C45");
		cys.insert("C46");
		cys.insert("C38");
		cys.insert("C39");
		cys.insert("C40");
		if (cys.find(sat) != cys.end()) return true;
		return false;
	}

	// get BDS satellite following SECM CYS attitude law
	// ------------------------------
	bool t_gsys::bds_secm(const string& sat)
	{
		set<string> cys;
		cys.insert("C25");
		cys.insert("C26");
		cys.insert("C27");
		cys.insert("C28");
		cys.insert("C29");
		cys.insert("C30");
		cys.insert("C34");
		cys.insert("C35");
		cys.insert("C43");
		cys.insert("C44");
		if (cys.find(sat) != cys.end()) return true;

		return false;
	}
} // namespace
