
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "gall/gallnav.h" 
#include "gutils/gsys.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"
#include "gutils/gstat.h"



using namespace std;

namespace gnut
{

	// constructor
	// ----------
	t_gallnav::t_gallnav()
		: t_gdata(),
		_com(false),
		_offset(0),
		_nepoch(t_gtime::GPS),
		_multimap(false),
		_overwrite(false),
		_chk_health(false),
		_chk_navig(false),
		_chk_tot(false)
	{
		gtrace("t_gallnav::constructor");
		id_type(t_gdata::ALLNAV);
		id_group(t_gdata::GRP_EPHEM);
	}


	// destructor
	// ----------
	t_gallnav::~t_gallnav()
	{
		gtrace("t_gallnav::destructor");

		_mapsat.clear();
	}


	// return gnav element
	// ----------
	shared_ptr<t_geph> t_gallnav::find(string sat, const t_gtime& t, bool chk_mask)
	{
		gtrace("t_gallnav::find");

		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		_gmutex.unlock(); return tmp;
	};

	// return gnav elements
	// ----------
	vector<shared_ptr<t_geph> > t_gallnav::find_mult(string sat, const t_gtime& t)
	{
		gtrace("t_gallnav::find_mult");

		_gmutex.lock();

		vector<shared_ptr<t_geph> > vec = this->_find_mult(sat, t);

		_gmutex.unlock();
		return vec;
	}


	// return position
	// ----------
	int t_gallnav::pos(string sat, const t_gtime& t, double  xyz[],
		double  var[], double  vel[], bool chk_mask) // [m]
	{
		gtrace("t_gallnav::pos");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);
		if (tmp == _null) 
		{
			tmp = t_gallnav::_find(sat, t, false);
		}
		// end adding
		if (tmp == _null)
		{
			for (int i = 0; i < 3; i++)
			{
				xyz[i] = 0.0;
				if (var) var[i] = 0.0;
				if (vel) vel[i] = 0.0;
			}
			_gmutex.unlock(); return -1;
		}

		int irc = tmp->pos(t, xyz, var, vel, _chk_health && chk_mask);
		if (irc == -1) 
		{
			irc = tmp->pos(t, xyz, var, vel, false);
		}
		// end adding		

		_gmutex.unlock(); return irc;

		//  return find( sat, t )->pos( t, xyz, var, vel, _chk_health && chk_mask );
	}


	// return satellite health
	// ----------
	bool t_gallnav::health(string sat, const t_gtime& t)
	{
		gtrace("t_gallnav::health");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, false);

		if (tmp == _null) { _gmutex.unlock(); return false; }

		bool status = tmp->healthy();

		_gmutex.unlock();
		return status;
	}


	// return position aka navigation
	// ----------
	int t_gallnav::nav(string sat, const t_gtime& t, double  xyz[],
		double  var[], double  vel[], bool chk_mask) // [m]
	{
		gtrace("t_gallnav::nav");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		if (tmp == _null)
		{
			for (int i = 0; i < 3; i++)
			{
				xyz[i] = 0.0;
				if (var) var[i] = 0.0;
				if (vel) vel[i] = 0.0;
			}
			_gmutex.unlock(); return -1;
		}
		int irc = tmp->nav(t, xyz, var, vel, _chk_health && chk_mask);

		_gmutex.unlock();
		return irc;

		//  return find( sat, t )->pos( t, xyz, var, vel, _chk_health && chk_mask );
	}


	// return clock corrections
	// ----------
	int t_gallnav::clk(string sat, const t_gtime& t, double*  clk,
		double*  var, double* dclk, bool chk_mask) // [s]
	{
		gtrace("t_gallnav::clk");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t, _chk_health && chk_mask);

		if (tmp == _null) {
			*clk = 0.0;
			if (var) *var = 0.0;
			if (dclk) *dclk = 0.0;
			_gmutex.unlock(); return -1;
		}

		int irc = tmp->clk(t, clk, var, dclk, _chk_health && chk_mask);
		_gmutex.unlock(); return irc;

		//  return this->find( sat, t )->clk( t, clk, var, dclk, _chk_health && chk_mask );
	}

	// print function
	// -------------------
	void t_gallnav::print(string sat, const t_gtime& t)
	{
		gtrace("t_gallnav::print");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		shared_ptr<t_geph> tmp = t_gallnav::_find(sat, t);
		tmp->print();

		_gmutex.unlock(); return;
	}


	// return list of available constellation
	// ----------
	set<GSYS> t_gallnav::systems()
	{
		gtrace("t_gallnav::systems");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		set<GSYS> all_sys;
		auto itPRN = _mapsat.begin();

		while (itPRN != _mapsat.end()) {
			GSYS gsys = t_gsys::str2gsys(itPRN->first.substr(0, 1));
			if (all_sys.find(gsys) == all_sys.end()) all_sys.insert(gsys);
			itPRN++;
		}
		_gmutex.unlock(); return all_sys;
	}

	// return list of available satellites
	// ----------
	set<string> t_gallnav::satellites()
	{
		gtrace("t_gallnav::satellites");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		set<string> all_sat;
		auto itPRN = _mapsat.begin();

		while (itPRN != _mapsat.end()) {
			if (all_sat.find(itPRN->first) == all_sat.end()) all_sat.insert(itPRN->first);
			itPRN++;
		}

		_gmutex.unlock(); return all_sat;
	}


	// return first position for satellite
	// ----------
	t_gtime t_gallnav::beg_gnav(string prn)
	{
		gtrace("t_gallnav::beg_gnav");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		t_gtime tmp = LAST_TIME;
		if (!prn.empty()) {
			if (_mapsat.find(prn) != _mapsat.end()) tmp = _mapsat[prn].begin()->first;
		}
		else {
			for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT) {
				for (auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it) {
					if (_mapsat[itSAT->first].begin()->first < tmp) tmp = _mapsat[itSAT->first].begin()->first;
				}
			}
		}

		_gmutex.unlock(); return tmp;
	}


	// return last position for satellite
	// ----------
	t_gtime t_gallnav::end_gnav(string prn)
	{
		gtrace("t_gallnav::end_gnav");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		t_gtime tmp = FIRST_TIME;

		if (!prn.empty()) {
			if (_mapsat[prn].size() > 0 && _mapsat.find(prn) != _mapsat.end()) {
				tmp = _mapsat[prn].rbegin()->first;
			}

		}
		else {
			for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT) {
				for (auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it) {
					if (_mapsat[itSAT->first].rbegin()->first > tmp) {
						tmp = _mapsat[itSAT->first].rbegin()->first;
					}
				}
			}
		}

		_gmutex.unlock(); return tmp;
	}


	// add navigation message
	// ----------
	int t_gallnav::add(shared_ptr<t_gnav> nav)
	{
		gtrace("t_gallnav::add");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		t_gtime ep(nav->epoch());
		string sat = nav->sat();

		if (_chk_navig) {
			set<string> msg; nav->chk(msg);
			if (_log) {
				for (auto it = msg.begin(); it != msg.end(); it++)
					_log->comment(2, "gallnav", *it);
			}
		}

		// test navigation type
		bool add = false;
		if (_multimap) { add = true; }    // multimap
		else if (_mapsat[sat].find(ep) == _mapsat[sat].end()) { add = true; }    // non-existent
		else if (nav->id_type() == t_gdata::EPHGAL) {                  // check nav-type
			auto itB = _mapsat[sat].lower_bound(ep);
			auto itE = _mapsat[sat].upper_bound(ep);
			while (itB != itE) {
				if (dynamic_pointer_cast<t_gnav>(itB->second)->gnavtype() == nav->gnavtype()) {
					add = false; break;
				} // exclude the message and skip!
				else { add = true; } // ok
				++itB;
			}
		}
		if (!nav->valid()) add = false;    // validity test  

		if (add) {
			if (_log && _log->verb() >= 3) {
				ostringstream lg;
				lg << "add sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
					<< " iod: " << fixed << setw(3) << nav->iod()
					<< " flg: " << fixed << setw(3) << nav->healthy();
				if (nav->gio()) lg << " " << base_name(nav->gio()->path());
				_log->comment(3, "gallnav", lg.str());
			}
			_mapsat[sat].insert(make_pair(ep, nav));
			if (sat.substr(0, 1) == "R") {
				_glo_freq_num[sat] = nav->freq_num();
			}
			

		}
		else if (_log && _log->verb() >= 3) {
			ostringstream lg;
			lg << "skip sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
				<< " iod: " << fixed << setw(3) << nav->iod()
				<< " flg: " << fixed << setw(3) << nav->healthy();
			if (nav->gio()) lg << " " << base_name(nav->gio()->path());
			_log->comment(3, "gallnav", lg.str());
		}

		_gmutex.unlock(); return 0;
	}


	// return number of epochs
	// ----------
	unsigned int t_gallnav::nepochs(const string& prn)
	{
		gtrace("t_gallnav::nepochs");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		unsigned int tmp = 0;
		if (_mapsat.find(prn) != _mapsat.end()) tmp = _mapsat[prn].size();

		_gmutex.unlock(); return tmp;
	}


	// list of epochs  ( "" = all, "G" = system, "G01" = satellite )
	// ----------
	set<t_gtime> t_gallnav::epochs(string prn)
	{
		gtrace("t_gallnav::epochs");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		set<t_gtime> epochs;

		auto itFIRST = _mapsat.begin();
		auto itLAST = _mapsat.end();

		if (!prn.empty() && prn.size() == 3) itLAST++ = itFIRST = _mapsat.find(prn);

#ifdef DEBUG
		if (itFIRST == _mapsat.begin()) cout << prn << " FIRST = begin\n";
		if (itLAST == _mapsat.end()) cout << prn << "  LAST = end  \n";
		if (itFIRST == itLAST) cout << prn << "  LAST = FIRST\n";
#endif

		for (auto itSAT = itFIRST; itSAT != itLAST; ++itSAT)
		{
			if (itSAT->first == prn ||                    // individual satellite
				(itSAT->first[0] == prn[0] && prn.size() == 1)) // all satellites of system
			{
				for (auto itEPO = itSAT->second.begin(); itEPO != itSAT->second.end(); ++itEPO) {
					epochs.insert(itEPO->first);
				}
			}
		}

		_gmutex.unlock(); return epochs;
	}


	// list of nav messages
	// ----------
	vector<shared_ptr<t_geph> > t_gallnav::vec_nav(string prn)
	{
		gtrace("t_gallnav::vec_nav");

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

		vector<shared_ptr<t_geph> > v_nav;

		auto itFIRST = _mapsat.begin();
		auto itLAST = _mapsat.end();

		if (!prn.empty()) itLAST++ = itFIRST = _mapsat.find(prn);

		for (auto itSAT = itFIRST; itSAT != itLAST; ++itSAT) {
			for (auto itEPO = itSAT->second.begin(); itEPO != itSAT->second.end(); ++itEPO) {
				v_nav.push_back(itEPO->second);
			}
		}

		_gmutex.unlock(); return v_nav;
	}


	// list of calculated crd
	// ----------
	map<string, t_gtriple> t_gallnav::map_xyz(set<string> prns, const t_gtime& epo)
	{
		map<string, t_gtriple> m_xyz;

		t_gtriple xyz;

		for (auto itPRN = prns.begin(); itPRN != prns.end(); itPRN++) {
			string prn = *itPRN;

			double sat_xyz[3] = { 0.0, 0.0, 0.0 };
			double sat_var[3] = { 0.0, 0.0, 0.0 };
			double sat_vel[3] = { 0.0, 0.0, 0.0 };

			if (pos(prn, epo, sat_xyz, sat_var, sat_vel) >= 0) {
				xyz.set(0, sat_xyz[0]);
				xyz.set(1, sat_xyz[1]);
				xyz.set(2, sat_xyz[2]);

				m_xyz[prn] = xyz;

			}
			else {
				continue;
			}
		}

		return m_xyz;
	}

	// list of calculated crd and clk   
	map<string, double> t_gallnav::map_clk(set<string> prns, const t_gtime& epo)
	{
		map<string, double> m_clk;

		double sat_clk = 0;
		double var = 0;
		double dclk = 0;

		for (auto itPRN = prns.begin(); itPRN != prns.end(); itPRN++) {
			if (clk(*itPRN, epo, &sat_clk, &var, &dclk) >= 0) {
				m_clk[*itPRN] = sat_clk;
			}
		}

		return m_clk;
	}

	// clean invalid messages
	// ----------
	void t_gallnav::clean_invalid()
	{
		gtrace("t_gallnav::clean_invalid");

		_gmutex.lock();

		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT) {
			string prn = itSAT->first;

			auto itEPO = itSAT->second.begin();
			while (itEPO != itSAT->second.end()) {
				if (!itEPO->second->valid()) {
					if (_log) _log->comment(2, "gallnav", "Del NAV invalid: " + prn
						+ " " + itEPO->second->epoch().str_ymdhms()
						+ " " + itEPO->second->gio()->path());
					itEPO = _mapsat[prn].erase(itEPO);
				}
				else { itEPO++; }
			}
		}

		_gmutex.unlock();
	}

	// clean redundant messages
	// ----------
	void t_gallnav::clean_duplicit()
	{
		gtrace("t_gallnav::clean_duplicit");

		_gmutex.lock();

		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string prn = itSAT->first;
			map<t_gtime, set<int> > list_epo;

			for (auto itEPO = itSAT->second.begin();
				itEPO != itSAT->second.end(); )
			{
				t_gtime epo = itEPO->first;
				int src = itEPO->second->src(false); // clean all duplicites (false: distinguish INAV/FNAV type only)

				if (list_epo.find(epo) == list_epo.end() ||
					list_epo[epo].find(src) == list_epo[epo].end()
					) {
					list_epo[epo].insert(src);
					if (_log) _log->comment(2, "gallnav", "OK  NAV unique  : " + prn
						+ " " + itEPO->second->epoch().str_ymdhms()
						+ " " + gnavtype2str(itEPO->second->gnavtype()));
					++itEPO;
				}
				else { // remove redundant
					if (_log) _log->comment(2, "gallnav", "Del NAV multiple: " + prn
						+ " " + itEPO->second->epoch().str_ymdhms()
						+ " " + gnavtype2str(itEPO->second->gnavtype()));
					itEPO = _mapsat[prn].erase(itEPO);
				}
			}
		}

		_gmutex.unlock();
	}

	void t_gallnav::erase_time(const t_gtime& t) {
		gtrace("t_gallnav::erase_time");
		_gmutex.lock();

		for (auto& obj : _mapsat) {
			auto& obj_sat = obj.second;
			auto iter = obj_sat.lower_bound(t);
			obj_sat.erase(obj_sat.begin(), iter);
		}
		_gmutex.unlock();
	}

	// clean function
	// ----------
	void t_gallnav::clean_outer(const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::clean_outer");

		if (end < beg) return;
		if (beg == FIRST_TIME) return;
		if (end == LAST_TIME) return;

#ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();

#ifdef DEBUG   
		if (_log) _log->comment(1, "gallnav", "NAV cleaned interval: " + beg.str_ymdhms() + " - " + end.str_ymdhms());
#endif

		// loop over all satellites
		auto itPRN = _mapsat.begin();
		while (itPRN != _mapsat.end()) {
			string prn = itPRN->first;

			// find and CLEAN all data (epochs) out of the specified period !
			auto itFirst = _mapsat[prn].begin();
			auto itLast = _mapsat[prn].end();
			auto itBeg = _mapsat[prn].lower_bound(beg);  // greater only   old: // greater|equal
			auto itEnd = _mapsat[prn].upper_bound(end);  // greater only!

		//  cout << "distance = " << distance(itBeg,itEnd) << " prn " << prn << endl;

			// remove before BEGIN request
			if (itBeg != itFirst) {
				if (_log && _log->verb() >= 2) {
					_log->comment(2, "gallnav", "Del NAV before: " + itFirst->first.str_ymdhms() + " " + prn);
				}
				auto it = itFirst;

				while (it != itBeg && it != itLast) {
					if ((it->second).use_count() == 1) { it = _mapsat[prn].erase(it); }
					else { it++; }
				}
			}

			// remove after END request
			if (itEnd != itLast) {
				if (_log && _log->verb() >= 2) {
					_log->comment(2, "gallnav", "Del NAV after : " + itEnd->first.str_ymdhms() + " " + prn);
				}
				auto it = itEnd;

				while (it != itLast) {
					if ((it->second).use_count() == 1) { it = _mapsat[prn].erase(it); }
					else { it++; }
				}
			}
			itPRN++;
		}
#ifdef BMUTEX   
		lock.unlock();
#endif
		_gmutex.unlock();
	}


	// get number of satellites
	// ----------
	int t_gallnav::nsat(GSYS gs)
	{
		gtrace("t_gallnav::nsat(GSYS)");

		int nsatell = 0;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first;
			if (t_gsys::char2gsys(sat[0]) == gs || gs == GNS) nsatell++;
		}
		return nsatell;
	}


	// get interval between messages   
	// ----------
	int t_gallnav::intv(GSYS gs)
	{
		gtrace("t_gallnav::intv(GSYS)");

		int interval = 0;  if (gs == GNS) return interval;

		// collect sampling intervals
		map<int, int> m_intv;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first; if (t_gsys::char2gsys(sat[0]) != gs) continue;

			t_gtime lstEPO;
			for (auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO)
			{
				if (itEPO != _mapsat[sat].begin()) {
					int diff = int(floor(itEPO->first - lstEPO));
					m_intv[diff]++;
				}
				lstEPO = itEPO->first;
			}
		}

		// auto-detect sampling interval
		int count = 0;
		for (auto it = m_intv.begin(); it != m_intv.end(); ++it) {

			//    cout << " " << it->first << ":" << it->second;
			if (it->second > count) {
				interval = it->first;
				count = it->second;
			}
			//    cout << endl;
		}
		//  cout << t_gsys::gsys2str(gs) << " : " << interval << endl;

		return interval;
	}


	// get existing number of messages
	// ----------
	int t_gallnav::have(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::have(GSYS,beg,end)");

		int existing = 0;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string sat = itSAT->first; if (t_gsys::char2gsys(sat[0]) != gs) continue;

			for (auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO) {
				if (itEPO->first < beg - 900 || itEPO->first > end + 900) continue;
				else existing++;
			}
		}
		return existing;
	}


	// get expected number of messages
	// ----------
	int t_gallnav::expt(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::expt(GSYS,beg,end)");

		if (end < beg) return 0;
		int xint = intv(gs);
		int xsat = nsat(gs);
		int diff = int(floor(end - beg));
		if (xint == 0 || xsat == 0 || diff == 0.0) return 0;
		return int(floor(diff / xint * xsat));
	}


	// get excluded number of messages
	// ----------
	int t_gallnav::excl(GSYS gs, const t_gtime& beg, const t_gtime& end)
	{
		gtrace("t_gallnav::excl(GSYS,beg,end)");
		int exclude = 0;
		return exclude;
	}


	// consolidate NAV messsage
	// ----------
	int t_gallnav::consolidate(double cfdi)
	{
		gtrace("t_gallnav::consolidate");

		if (cfdi <= 0.0) cfdi = 10.0; // DEFAULT

		map<shared_ptr<t_geph>, double> m_penalty;

		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string prn = itSAT->first;

			NAVDATA str_KPL[] = { NAV_A, NAV_E, NAV_M, NAV_I, NAV_IDOT, NAV_OMEGA, NAV_OMG, NAV_OMGDOT, NAV_DN, NAV_F0 };
			NAVDATA str_XYZ[] = { NAV_X, NAV_XD, NAV_XDD, NAV_Y, NAV_YD, NAV_YDD, NAV_Z, NAV_ZD, NAV_ZDD };

			vector<NAVDATA> param;
			switch (prn[0]) {
			case 'R': case 'S': param.assign(str_XYZ, str_XYZ + 9); break;
			default: param.assign(str_KPL, str_KPL + 10); break;
			}

			for (auto itPAR = param.begin(); itPAR != param.end(); ++itPAR)
			{
				vector<shared_ptr<t_geph> > v_ptr;
				vector<shared_ptr<t_geph> >::iterator itPTR;

				vector<t_gtime> v_tim;
				vector<t_gtime>::const_iterator itTIM;

				vector<double> v_val, v_timdif, v_valdif;
				vector<double>::const_iterator itVAL, itDIF;

				// prepare values & differences vectors for statistics
				for (auto itEPH = _mapsat[prn].begin(); itEPH != _mapsat[prn].end(); ++itEPH)
				{
					t_timdbl tmp = itEPH->second->param(*itPAR);
					double mindif = 60.0; // minimum timedif applied for differences

					if (!v_tim.empty()) { // itEPH != _mapsat[prn].begin() ){
					  // check selected parameters (differences) for periodic switch
						if (itEPH->second->param_cyclic(*itPAR)) {
							if (tmp.second - v_val.back() > G_PI) tmp.second -= 2 * G_PI;
							if (tmp.second - v_val.back() < -G_PI) tmp.second += 2 * G_PI;
							mindif = 1.0; // periodic
						}
						v_timdif.push_back(tmp.first - v_tim.back());
						v_valdif.push_back(fabs(v_timdif.back()) > mindif // relax a short time difference to a minute
							? (tmp.second - v_val.back()) / (tmp.first - v_tim.back())
							: (tmp.second - v_val.back()) / mindif
						);
					}

					v_ptr.push_back(itEPH->second); // save pointer
					v_tim.push_back(tmp.first);
					v_val.push_back(tmp.second);

				}

				if (v_val.size() <= 1) { continue; } // no values for v_val or v_valdif

				// calculate values & differences statistics
				t_gstat stt_val(v_val);
				t_gstat stt_dif(v_valdif);

				if (!stt_val.valid() && _log) { _log->comment(1, "gallnav", "Invalid statistics for prn[" + prn + "] values"); }
				if (!stt_dif.valid() && _log) { _log->comment(1, "gallnav", "Invalid statistics for prn[" + prn + "] differences"); }

				ostringstream oSTT;
				oSTT << fixed << setprecision(5)
					<< prn << " stat: " << setw(12) << "parameter" << setw(12) << "difference" << endl
					<< prn << " medi: " << setw(12) << stt_val.get_median() << setw(12) << stt_dif.get_median() << endl
					<< prn << " mean: " << setw(12) << stt_val.get_mean() << setw(12) << stt_dif.get_mean() << endl
					<< prn << " sdev: " << setw(12) << stt_val.get_sdev() << setw(12) << stt_dif.get_sdev() << endl
					<< prn << " rms : " << setw(12) << stt_val.get_rms() << setw(12) << stt_dif.get_rms() << endl
					<< prn << " min : " << setw(12) << stt_val.get_min() << setw(12) << stt_dif.get_min() << endl
					<< prn << " max : " << setw(12) << stt_val.get_max() << setw(12) << stt_dif.get_max() << endl;

				if (_log) _log->comment(-3, oSTT.str());

				for (itPTR = v_ptr.begin(), itVAL = v_val.begin(), itTIM = v_tim.begin(), itDIF = v_valdif.begin();
					itPTR != v_ptr.end(), itVAL != v_val.end();
					++itPTR, ++itVAL, ++itTIM, ++itDIF
					) {
					bool badV = (fabs(*itVAL - stt_val.get_median()) > fabs(cfdi * stt_val.get_sdev()) && stt_val.get_sdev() > 0);
					bool badD = false;
					if (itDIF != v_valdif.end()) {
						badD = (fabs(*itDIF - stt_dif.get_median()) > fabs(cfdi * stt_dif.get_sdev()) && stt_dif.get_sdev() > 0);
					}
					else {
						itDIF--;
					}
					if (itPTR == v_ptr.begin()) {

						ostringstream os;
						os << fixed << setprecision(9) << prn << itTIM->str_ymdhms(" ")
							<< "   rmsV: " << setw(20) << cfdi * stt_val.get_sdev()
							<< "   rmsD: " << setw(20) << cfdi * stt_dif.get_sdev();
						if (_log) _log->comment(-3, os.str());
					}

					ostringstream os;
					os << fixed << setprecision(6)
						<< prn << itTIM->str_ymdhms(" ")
						<< " " << setw(2) << *itPAR
						<< " val: " << setw(18) << *itVAL
						<< " " << setw(18) << *itVAL - stt_val.get_median()
						<< " badO: " << setw(1) << badV;

					if (itDIF != v_valdif.end()) {
						os << setprecision(12)
							<< " dif: " << setw(20) << *itDIF
							<< " " << setw(20) << *itDIF - stt_dif.get_median()
							<< " badD: " << setw(1) << badD;
					}
					else { os << endl; }
					if (_log) _log->comment(-3, os.str());

					if (badV) { m_penalty[*itPTR] += 1.00; }
					if (badD) {
						m_penalty[*itPTR] += 0.25;
						if ((itPTR + 1) != v_ptr.end())  m_penalty[*(itPTR + 1)] += 0.25;
					}
				}
			}
		}

		for (auto it = m_penalty.begin(); it != m_penalty.end(); ++it) {
			if (it->second > 3) { // criteria
				shared_ptr<t_geph> itP = it->first; itP->valid(false); // SET INVALID

				if (_log) _log->comment(1, "gallnav", "Set NAV invalid: " + itP->sat() + " " + itP->epoch().str_ymdhms()
					+ " [penalty" + dbl2str(it->second, 1) + "] " + itP->gio()->path());
			}
		}

		return 0;
	}


	// consolidate healthy status & biases
	// ----------
	int t_gallnav::consolidate_others()
	{
		gtrace("t_gallnav::consolidate_others");

		// identify special issue
		map<string, set<NAVDATA> > excl;
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string prn = itSAT->first;
			GSYS gs = t_gsys::char2gsys(itSAT->first[0]);

			for (auto itEPH = _mapsat[prn].begin(); itEPH != _mapsat[prn].end(); ++itEPH)
			{
				shared_ptr<t_geph> gnav = itEPH->second;

				GNAVTYPE nav = gnav->gnavtype(false); // distinguish INAV/FNAV only (no source)
				GNAVTYPE ful = gnav->gnavtype(true);  // distinguish all
				string   hlt = gnav->health_str();
				string   pth = gnav->gio()->path();

				if (gs == GAL) {  // identify TRIMBLE NETR9 issue (GALILEO)
					if ((nav == INAV && hlt.compare("110000000") == 0) ||
						(ful == INAV_E01 && hlt.compare("000000111") == 0) ||
						(ful == INAV_E01 && hlt.compare("000111111") == 0))
					{
						if (_log) { _log->comment(1, "gallnav", "TRIMBLE NETR9 identified [" + prn + "] " + gnavtype2str(ful) + " " + hlt + " " + pth); }
						excl[pth].insert(NAV_HEALTH); // mark for exclusion
						excl[pth].insert(NAV_TGD0);   // mark for exclusion
						excl[pth].insert(NAV_TGD1);   // mark for exclusion
					}
				}
			}
		}

		// consolidate
		for (auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT)
		{
			string prn = itSAT->first;
			GSYS gs = t_gsys::char2gsys(itSAT->first[0]);

			vector<NAVDATA> param;
			//                  param.push_back(NAV_IOD);    // don't do this
			param.push_back(NAV_HEALTH); // should be before TGD to identify exclusions

			switch (prn[0]) {
			case 'R': break;
			case 'S': break;
			case 'E': param.push_back(NAV_TGD0); param.push_back(NAV_TGD1); break;
			case 'C': param.push_back(NAV_TGD0); param.push_back(NAV_TGD1); break;
			case 'G': param.push_back(NAV_TGD0); param.push_back(NAV_TGD1); param.push_back(NAV_TGD2); param.push_back(NAV_TGD3); break;
			case 'Q': param.push_back(NAV_TGD0); param.push_back(NAV_TGD1); param.push_back(NAV_TGD2); param.push_back(NAV_TGD3); break;
			case 'I': param.push_back(NAV_TGD0); param.push_back(NAV_TGD1); param.push_back(NAV_TGD2); param.push_back(NAV_TGD3); break;
			}

			for (auto itPAR = param.begin(); itPAR != param.end(); ++itPAR)
			{
				NAVDATA navd = *itPAR;
				bool tgd_param = (navd == NAV_TGD0 || navd == NAV_TGD1 || navd == NAV_TGD2 || navd == NAV_TGD3);
				bool tgd_combine = true; // use all NAV types for TGDs!!! (and see below)

				map<GNAVTYPE, map<double, int> >    day_sort;
				map<t_gtime, map<GNAVTYPE, map<double, int> > >  data_sort;
				map<GNAVTYPE, map<double, int>::const_iterator> itDAY;
				map<GNAVTYPE, map<double, int>::const_iterator> itMAX;

				for (auto itEPH = _mapsat[prn].begin(); itEPH != _mapsat[prn].end(); ++itEPH)
				{
					shared_ptr<t_geph> gnav = itEPH->second;

					t_gtime  epo = itEPH->first;
					t_timdbl tmp = gnav->param(*itPAR);
					GNAVTYPE nav = gnav->gnavtype(false); // distinguish INAV/FNAV only (no source)
					GNAVTYPE ful = gnav->gnavtype(true);  // full
					string   hlt = gnav->health_str();
					string   pth = gnav->gio()->path();
					double   val = tmp.second;

					if (gs == GAL) { // skip issue

						bool xxx = (excl.find(pth) != excl.end() && excl[pth].find(navd) != excl[pth].end());

						if (ful == INAV && (hlt.compare("111000000") == 0 || hlt.compare("000000111") == 0)) { // INAV (E1B+E5b with different HEALTH)
							if (_log) { _log->comment(1, "gallnav", "INAV - inconsistent E1B/E5b HEALTHY flags [" + prn + "] " + gnavtype2str(nav) + " " + hlt + " " + pth); }
							continue;
						} // use excl for HEALTHY INAV
						if (xxx && nav == FNAV && hlt.compare("000000000") == 0) { continue; } // use excl for HEALTHY FNAV
						if (xxx && nav == FNAV && tgd_param) { continue; } // use excl for GALILEO TGD
			  //        if(        val == 0.0  && tgd_param                     ){ continue; } // use excl zero values
						if (xxx &&navd == NAV_HEALTH) { continue; }
					}

					if (tgd_param && tgd_combine) { nav = NAV; } // use all NAV types for TGDs!!! (and see below)

			//        double xval = val; if(tgd_param) xval *= 1e9;
			//        if( _log ){ _log->comment(3,"gallnav",prn+"["+int2str(navd)+"] add to statistics "+epo.str_ymdhms()+" "+dbl2str(xval)
			//                                                                   +" "+gnavtype2str(nav)+" "+hlt+" "+pth); }
					data_sort[epo][nav][val]++;
					day_sort[nav][val]++;
				}

				// find DAY maximum
				for (auto itNAV = day_sort.begin(); itNAV != day_sort.end(); ++itNAV) {
					GNAVTYPE xnav = itNAV->first; if (tgd_param && tgd_combine) { xnav = NAV; }
					for (auto it = itNAV->second.begin(); it != itNAV->second.end(); ++it) {
						if (itDAY.find(xnav) == itDAY.end() || double_eq(itDAY[xnav]->first, 0.0) ||
							(it->second > itDAY[xnav]->second && !double_eq(it->first, 0.0))) {
							itDAY[xnav] = it;
						}
					}
				}

				// loop over nav records
				for (auto itEPH = _mapsat[prn].begin(); itEPH != _mapsat[prn].end(); ++itEPH)
				{
					shared_ptr<t_geph> gnav = itEPH->second;

					t_gtime  epo = itEPH->first;
					t_timdbl tmp = gnav->param(*itPAR);
					GNAVTYPE nav = gnav->gnavtype(false); // distinguish INAV/FNAV only (no source)
					string   hlt = gnav->health_str();
					string   pth = gnav->gio()->path();
					string  name = base_name(pth);

					auto itLOW = _mapsat[prn].lower_bound(epo);
					auto itUPP = _mapsat[prn].upper_bound(epo);

					if (itEPH == itLOW) // once find epoch-wise frequency maximum
					{
						if (data_sort.find(epo) != data_sort.end()) {
							for (auto itNAV = data_sort[epo].begin(); itNAV != data_sort[epo].end(); ++itNAV) {
								GNAVTYPE xnav = itNAV->first; if (tgd_param && tgd_combine) { xnav = NAV; }
								for (auto it = itNAV->second.begin(); it != itNAV->second.end(); ++it) {
									if (itMAX.find(xnav) == itMAX.end() || double_eq(itMAX[xnav]->first, 0.0) ||
										(it->second > itMAX[xnav]->second && !double_eq(it->first, 0.0))) {
										itMAX[xnav] = it;
									}
								}
								if (_log && _log->verb() >= 3) {
									ostringstream os;
									double xval = itMAX[xnav]->first; if (tgd_param) xval *= 1e9;
									os << setw(8) << gnavtype2str(xnav) << " Epoch maximum [" << navd << "]: " << prn << epo.str_ymdhms(" ")
										<< fixed << setprecision(3)
										<< setw(9) << xval
										<< setw(9) << itMAX[xnav]->second
										<< setw(9) << distance(itLOW, itUPP);
									if (_log) _log->comment(3, "gallnav", os.str());
								}
							}
						}
						else {
							for (auto itNAV = day_sort.begin(); itNAV != day_sort.end(); ++itNAV) {
								GNAVTYPE xnav = itNAV->first; if (tgd_param && tgd_combine) { xnav = NAV; }
								for (auto it = itNAV->second.begin(); it != itNAV->second.end(); ++it) {
									itMAX[xnav] = it;
									if (_log && _log->verb() >= 3) {
										ostringstream os;
										double xval = itDAY[xnav]->first; if (tgd_param) xval *= 1e9;
										os << setw(8) << gnavtype2str(xnav) << " Daily maximum used [" << navd << "]: " << prn << epo.str_ymdhms(" ")
											<< fixed << setprecision(3)
											<< setw(9) << xval
											<< setw(9) << itDAY[xnav]->second;
										if (_log) _log->comment(3, "gallnav", os.str());
									}
								}
							}
						}
					}

					if (tgd_param)  // consolidate TGDs
					{
						GNAVTYPE navfix = nav; if (tgd_combine) { navfix = NAV; } // fix all nav to TGD

						if (_log && _log->verb() >= 3) {
							ostringstream os;
							string xxx = (excl.find(pth) != excl.end() && excl[pth].find(navd) != excl[pth].end()) ? "TRIMBLE" : "ok";
							os << setw(8) << gnavtype2str(nav)
								<< setw(5) << gnavtype2str(navfix) << " testing epoch TGD*1e9 [" << navd << "]: " << prn << epo.str_ymdhms(" ")
								<< fixed << setprecision(3) << setw(9) << tmp.second * 1e9
								<< " " << itMAX[navfix]->first * 1e9
								<< ":" << itMAX[navfix]->second
								<< " " << itDAY[navfix]->first * 1e9
								<< ":" << itDAY[navfix]->second
								<< " " << xxx
								<< " " << name;
							if (_log) _log->comment(3, "gallnav", os.str());
						}

						if ((int((tmp.second - itMAX[navfix]->first)*1e12) != 0) &&
							(int((itMAX[navfix]->first)*1e12) != 0) &&
							((itMAX[navfix]->second)) > 0)
						{
							gnav->param(navd, itMAX[navfix]->first); // change data to EPOCH MAX
							if (_log) {
								ostringstream os;
								os << setw(8) << gnavtype2str(nav) << " epoch TGD*1e9 changed [" << navd << "]: " << prn << epo.str_ymdhms(" ")
									<< fixed << setprecision(3)
									<< setw(8) << tmp.second*1e9 << " --> "
									<< setw(8) << itMAX[navfix]->first*1e9 << setprecision(0)
									<< " [" << itMAX[navfix]->second << "] " << name;
								_log->comment(2, "gallnav", os.str());
							}

						}
						else if ((int(tmp.second*1e12) == 0) &&
							(int(itDAY[navfix]->first*1e12) != 0))
						{
							gnav->param(navd, itDAY[navfix]->first); // change data to DAILY MAX
							if (_log) {
								ostringstream os;
								os << setw(8) << gnavtype2str(nav) << " daily TGD*1e9 changed [" << navd << "]: " << prn << epo.str_ymdhms(" ")
									<< fixed << setprecision(3)
									<< setw(8) << tmp.second * 1e9 << " --> "
									<< setw(8) << itDAY[navfix]->first * 1e9 << setprecision(0)
									<< " [" << itDAY[navfix]->second << "] " << name;
								_log->comment(2, "gallnav", os.str());
							}
						}

					}
					else {   // consolidate HEALTHY+IOD

						if (_log && _log->verb() >= 3) {
							if (navd == NAV_HEALTH) {
								ostringstream os;
								os << gnavtype2str(nav) << " testing epoch HEALTHY [" << navd << "]: " << prn << epo.str_ymdhms(" ")
									<< fixed << setprecision(3) << setw(9) << tmp.second << " " << name;
								if (_log) _log->comment(3, "gallnav", os.str());
							}
						}

						if ((int(tmp.second - itMAX[nav]->first) != 0) &&
							(int(itMAX[nav]->first) != 0)) {
							gnav->param(navd, itMAX[nav]->first); // change data

							if (_log) {
								ostringstream os;
								os << setw(8) << gnavtype2str(nav) << " epoch HEALTHY changed [" << navd << "]: " << prn << epo.str_ymdhms(" ")
									<< fixed << setprecision(0)
									<< setw(8) << tmp.second << " --> "
									<< setw(8) << itMAX[nav]->first << " " << "[" << itMAX[nav]->second << "] " << name;
								_log->comment(2, "gallnav", os.str());
							}
						}
					}
				}
			}
		}

		return 0;
	}


	// get IONO CORR
    // -----------------
	t_iono_corr t_gallnav::get_iono_corr(const IONO_CORR c) const
	{
		if (_brdc_iono_cor.find(c) != _brdc_iono_cor.end()) {
			return _brdc_iono_cor.at(c);
		}
		t_iono_corr io;
		// default settings , according rtklib->rtkcmn.c->ionmodel
		if (c == IO_GPSA)
		{
			io.x0 = 0.1118e-07;
			io.x1 =-0.7451e-08;
			io.x2 =-0.5961e-07;
			io.x3 = 0.1192e-06;
		}
		else if (c == IO_GPSB)
		{
			io.x0 = 0.1167e+06;
			io.x1 =-0.2294e+06;
			io.x2 =-0.1311e+06;
			io.x3 = 0.1049e+07;
		}
		
		return io;
	}

	// set IONO CORR
	// -----------------
	void t_gallnav::add_iono_corr(const IONO_CORR c, t_iono_corr io)
	{
		if (_brdc_iono_cor.find(c) == _brdc_iono_cor.end()) {
			_brdc_iono_cor[c] = io;
		}
	}

	// return list of available satellites
	// ----------
	shared_ptr<t_geph> t_gallnav::_find(string sat, const t_gtime& t, bool chk_mask)
	{
		gtrace("t_gallnav::_find sat/time");

		if (_mapsat.find(sat) == _mapsat.end()) return _null; // make_shared<t_geph>();

		GSYS gnss = t_gsys::str2gsys(sat.substr(0, 1));

		if (_mapsat[sat].size() == 0) {
			if (_log) _log->comment(2, "gallnav", sat + " no gnav elements found");
			return _null;
		}

		auto it = _mapsat[sat].lower_bound(t);  // greater|equal  (can be still end())
		if (it == _mapsat[sat].end()) it--;                   // size() > 0 already checked above
		else if (it != _mapsat[sat].begin()) {             // use the nearest record
			auto it1 = it; --it1;
			if (abs(it1->first.diff(t)) < abs(it->first.diff(t))) it = it1;
		}

		double maxdiff = t_gnav::nav_validity(gnss)*1.1;

		if (gnss == GAL) {
			// Galileo ephemerides are valid for toc+10min -> toc+20...180min
			for (int bck = 1; bck <= 5; bck++) {     // max eph for going back is 5
				t_gtime toc = it->second->epoch();
				if ((t < toc + 10 * 60 || t > toc + maxdiff) ||
					(_chk_tot && !it->second->chktot(t))) {        // check ToT for using past messages only
					if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
					else break;
				}
				if ((t >= toc + 10 * 60 && t <= toc + maxdiff) &&
					(_chk_tot && it->second->chktot(t))) break;
			}
		}
		else if (gnss == BDS) {
			// BeiDou ephemerides are valid for toc -> toc+60min
			for (int bck = 1; bck <= 5; bck++) {     // max eph for going back is 5
				t_gtime toc = it->second->epoch();
				if ((t < toc || t > toc + maxdiff) ||
					(_chk_tot && !it->second->chktot(t))) {        // check ToT for using past messages only
					if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
					else break;
				}
				if ((t >= toc && t <= toc + maxdiff) &&
					(_chk_tot && it->second->chktot(t))) break;
			}
		}
		else if (gnss == GLO) {
			t_gtime tt = t - maxdiff;  // time span of Glonass ephemerides is 15 min

			auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal
			if (itt == _mapsat[sat].end()) {                         // size() > 0 already checked above
				if (_log) _log->comment(2, "gallnav", sat + " gnav element not found: " + t.str_ymdhms());
				return _null;
			}

			double dt1 = itt->first - t;
			double dt2 = it->first - t;
			if (fabs(dt1) < fabs(dt2)) it = itt;
		}
		else {
			auto itLAST = _mapsat[sat].end(); itLAST--;
			auto itSAVE = it;
			bool found = false;
			while (_chk_tot && !itSAVE->second->chktot(t) && itSAVE->second->epoch().sod() % 3600 != 0 && itSAVE != itLAST) {
				itSAVE++; // one more step forward (special case for non-regular ephemerides (e.g. XX:59:44))
			}
			if (_chk_tot && itSAVE->second->chktot(t)) {   // check found non-regular ephemeris
				found = true;
				it = itSAVE;
			}

			// 2 conditions for all cases when NAV is in fugure and the iterator should be moved back (it--)
			if (fabs(t - it->second->epoch()) > maxdiff ||  // too far navigation message in future!
				(_chk_tot && !it->second->chktot(t) && !found)        // check ToT for using past messages only
				) {
				if (_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
			}

		}

		// tested found ephemeris
		if (fabs(t - it->second->epoch()) > maxdiff ||
			(_chk_tot && !it->second->chktot(t))) {   // simulation real-time epoch search 
			if (_log) {
				string lg(sat + " gnav element not found: " + t.str_ymdhms());
				_log->comment(2, "gallnav", lg);
			}
			return _null;
		}

		if (_chk_health && chk_mask) {
			if (it->second->healthy()) return it->second;
			else                      return _null;
		}
		return it->second;

	}

	// return list of available satellites
	// ----------
	vector<shared_ptr<t_geph> > t_gallnav::_find_mult(string sat, const t_gtime& t)
	{
		gtrace("t_gallnav::_find vec sat/time");

		vector<shared_ptr<t_geph> > vec_geph;

		if (_mapsat.find(sat) == _mapsat.end()) {
			vec_geph.push_back(_null);
			return vec_geph;
		}

		GSYS gnss = t_gsys::str2gsys(sat.substr(0, 1));

		auto it = _mapsat[sat].lower_bound(t);  // greater|equal

		double maxdiff = t_gnav::nav_validity(gnss)*1.1;

		if (gnss != GLO) {
			if (it == _mapsat[sat].end()) {
				t_gtime tt = t - maxdiff;
				auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal
				if (itt == _mapsat[sat].end()) {
					if (_log) {
						string lg(sat + " gnav element not found: " + t.str_ymdhms());
						_log->comment(3, "gallnav", lg);
					}
					vec_geph.push_back(_null); return vec_geph;
				}
				else it = itt;
			}
		}
		else {
			t_gtime tt = t - maxdiff;  // time span of Glonass ephemerides is 15 min
			auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal     
			if (it == _mapsat[sat].end()) {
				if (itt == _mapsat[sat].end()) {
					vec_geph.push_back(_null); return vec_geph;
				}
				else it = itt;
			}
			else {
				if (itt == _mapsat[sat].end()) {
					vec_geph.push_back(_null); return vec_geph;
				}
				double dt1 = itt->first - t;
				double dt2 = it->first - t;
				if (fabs(dt1) < fabs(dt2)) it = itt;
			}
		}

		t_gtime epo = it->first;
		unsigned int cnt = _mapsat[sat].count(epo);

		if (cnt >= 1) {
			for (auto itMULT = _mapsat[sat].equal_range(epo).first; itMULT != _mapsat[sat].equal_range(epo).second; ++itMULT) {
				vec_geph.push_back(itMULT->second);
			}
		}

		return vec_geph;
	}


	// return list of available satellites
	// ----------
	shared_ptr<t_geph> t_gallnav::_find(string sat, int iod, const t_gtime& t)
	{
		gtrace("t_gallnav::_find sat/iod");

		if (_mapsat.find(sat) == _mapsat.end()) return _null; // make_shared<t_geph>();

		auto it = _mapsat[sat].begin();
		while (it != _mapsat[sat].end()) {
			shared_ptr<t_gnav> pt_nav;
			if ((pt_nav = dynamic_pointer_cast<t_gnav>(it->second)) != nullptr) {
				//       cout << "nav->epoch() " << pt_nav->epoch().str_ymdhms() << " nav->iod() = " << pt_nav->iod()  << " rtcm iod: " << iod << endl;
				if (pt_nav->iod() == iod && abs(pt_nav->epoch().diff(t)) <= MAX_GPS_TIMEDIFF) {
					break; // found
				}
			}
			it++;
		}

		// not found !
		if (it == _mapsat[sat].end()) {
			if (_log && _log->verb() >= 3) {
				ostringstream lg;
				lg << sat + " gnav element not found "
					<< " [iod: " << iod << "]";
				_log->comment(3, "gallnav", lg.str());
			}
			return _null; // make_shared<t_geph>();
		}

		return it->second;
	}

} // namespace
