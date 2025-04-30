
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <assert.h>

#include "gall/gallpar.h"

using namespace std;

namespace gnut {

	// Add t_gpar to t_gallpar
	// Parameter is stored at the end of the vector
	// -----------------------------------------------
	void t_gallpar::addParam(t_gpar newPar)
	{
		gtrace("t_gallpar::addParam");
		this->_vParam.push_back(newPar);
		this->_point_par.push_back(_max_point++);
		this->_index_par[newPar.get_head()][newPar.get_timearc()] = this->_point_par[this->_point_par.size() - 1];
	}

	// Delete paremeter
	// Parameter is deleted according to index value
	// ----------------------------------------------------
	void t_gallpar::delParam(int i)
	{
		gtrace("t_gallpar::delParam");

		auto& all = this->_index_par.at(_vParam[i].get_head());
		for (auto iter = all.begin(); iter != all.end(); iter++) 
		{
			if (iter->second == _point_par[i]) 
			{
				all.erase(iter);
				break;
			}
		}
		if (this->_index_par[_vParam[i].get_head()].size() == 0) 
		{
			this->_index_par.erase(_vParam[i].get_head());
		}
		_point_par.erase(_point_par.begin() + i);

		_vParam.erase(_vParam.begin() + i);
	}

	// get position of item according: station name, par type, PRN, begin time, end time
	// -----------------------------------------------------
	int t_gallpar::getParam(const string& site, const par_type& par, const string& prn,
		const t_gtime& beg, const t_gtime& end, const  int& channel) const
	{
		t_gparhead parhead(par, site, prn, channel);
		if (_index_par.find(parhead) == _index_par.end()) 
		{
			return -1;
		}
		else 
		{
			const auto& all = _index_par.at(parhead);
			t_gtimearc dst_timearc(beg, end);

			auto all_end = all.end();
			auto par_beg = _point_par.begin();
			auto par_end = _point_par.end();
			for (auto iter = all.begin(); iter != all_end; iter++)
			{
				if (iter->first.inside(dst_timearc)) 
				{
					auto ans = lower_bound(par_beg, par_end, iter->second);
					assert(ans != par_end && *ans == iter->second);
					return ans - par_beg;
				}
			}
			return -1;
		}
		return -1;
	}

	// get position of item according: index (position in covariance matrix)
	// -----------------------------------------------------
	int t_gallpar::getParam(int index)
	{
		gtrace("t_gallpar::getParam(int)");

		for (unsigned int i = 0; i <= _vParam.size() - 1; i++)
		{
			if (_vParam[i].index == index)
			{
				return i;
				break;
			}
		}
		return -1;
	}

	int t_gallpar::getParIndex(int idx)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			return _vParam[idx].index;
		}
		else
		{
			return -1;
		}
	}

	int t_gallpar::getAmbParam(string site, string prn, par_type  type, const t_gtime& beg, const t_gtime& end) const
	{
		gtrace("t_gallpar::getAmbParam");

		if (this->_index_par.count(t_gparhead(type, site, prn)) == 0) {
			return -1;
		}
		else {
			auto all = this->_index_par.find(t_gparhead(type, site, prn))->second;
			t_gtimearc dst_timearc(beg, end);
			for (auto iter = all.begin(); iter != all.end(); iter++) {
				if (dst_timearc.inside(iter->first)) {
					auto ans = lower_bound(_point_par.begin(), _point_par.end(), iter->second);
					assert(ans != _point_par.end() && *ans == iter->second);
					return ans - _point_par.begin();
				}
			}
			return -1;
		}
		return -1;

	}

	t_gpar& t_gallpar::getAmbParam(int idx)
	{
		return _vParam[idx];
	}

	double t_gallpar::getParValue(int idx)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			return _vParam[idx].value();
		}
		else
		{
			return 0.0;
		}
	}

	double t_gallpar::getParSig(int idx)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			return _vParam[idx].apriori();
		}
		else
		{
			return 0.0;
		}
	}

	void t_gallpar::setParValue(int idx, double value)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			_vParam[idx].value(value);
		}
	}

	void t_gallpar::setParSig(int idx, double value)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			_vParam[idx].apriori(value);
		}
	}

	void t_gallpar::setParRemove(int idx, bool judge)
	{
		if (idx >= 0 && idx < _vParam.size())
		{
			_vParam[idx].lremove = judge;
		}
	}

	// Reindexing parametres.
	// New indexes are reordered form 1 to n
	// -----------------------------------------------
	void t_gallpar::reIndex()
	{
		gtrace("t_gallpar::reIndex");

		int index_new = 1;
		for (unsigned int iPar = 0; iPar <= _vParam.size() - 1; iPar++)
		{
			_vParam[iPar].index = index_new;
			index_new++;
		}
	}

	// Reindexing parametres.
	// All indexes larger than "i" is decresed by 1
	// -----------------------------------------------
	void t_gallpar::decIndex(int i)
	{
		gtrace("t_gallpar::decIndex(int)");

		for (unsigned int iPar = 0; iPar <= _vParam.size() - 1; iPar++)
		{
			if (_vParam[iPar].index > i) _vParam[iPar].index -= 1;
		}
	}

	// Reindexing parametres.
	// All indexes larger than "i" is incresed by 1
	// -----------------------------------------------
	void t_gallpar::incIndex(int i)
	{
		gtrace("t_gallpar::incIndex(int)");

		for (unsigned int iPar = 0; iPar <= _vParam.size() - 1; iPar++)
		{
			if (_vParam[iPar].index >= i) _vParam[iPar].index += 1;
		}
	}

	// Get number of parametres
	// ---------------------------------------------
	unsigned int t_gallpar::parNumber() const
	{
		gtrace("t_gallpar::parNumber");

		return _vParam.size();
	}

	// Get number of ambiguities
	// ---------------------------------------------
	unsigned int t_gallpar::ambNumber() const
	{
		gtrace("t_gallpar::ambNumber");

		int i = 0;
		for (const auto& par : _vParam)
			if (t_gpar::is_amb(par.parType)) ++i;

		return i;
	}

	// Get number of orbit parametres
   // ---------------------------------------------
	unsigned int t_gallpar::orbParNumber() const
	{
		gtrace("t_gallpar::orbParNumber");

		unsigned int orbparnum = 0;
		for (auto sat : _vOrbParam)
		{
			orbparnum += sat.second.size();
		}

		return orbparnum;
	}

	int t_gallpar::maxIndex() const
	{
		int val = 0;
		for (const auto& par : _vParam) {
			if (par.index > val) val = par.index;
		}
		return val;
	}


	// Strore Coordinates parametres to t_triple crd
	// ------------------------------------------------

	int t_gallpar::getCrdParam(string station, t_gtriple& crd, t_gtime Tbeg, t_gtime Tend) const
	{
		gtrace("t_gallpar::getCrdParam");

		int found = 0;

		int idx = this->getParam(station, par_type::CRD_X, "", Tbeg, Tend);
		int idy = this->getParam(station, par_type::CRD_Y, "", Tbeg, Tend);
		int idz = this->getParam(station, par_type::CRD_Z, "", Tbeg, Tend);

		if (idx != -1) {
			crd.set(0, _vParam[idx].value());
			found++;
		}
		if (idy != -1) {
			crd.set(1, _vParam[idy].value());
			found++;
		}
		if (idz != -1) {
			crd.set(2, _vParam[idz].value());
			found++;
		}

		if (found == 3)  return  1;    // all three crd were found
		if (found == 1)  return -1;    // just one crd was found
		if (found == 2)  return -2;    // just two crd was found

		if (found == 0) return  -3;    // crd not found

		return -1;
	}

	int t_gallpar::getVelParam(string station, t_gtriple & vel, t_gtime Tbeg, t_gtime Tend) const
	{
		gtrace("t_gallpar::getCrdParam");

		int found = 0;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {
			if (iter->parType == par_type::VEL_X && iter->site.compare(station) == 0 &&
				iter->beg == Tbeg && iter->end == Tend) {
				vel.set(0, iter->value());
				found++;
			}
			else if (iter->parType == par_type::VEL_Y && iter->site.compare(station) == 0 &&
				iter->beg == Tbeg && iter->end == Tend) {
				vel.set(1, iter->value());
				found++;
			}
			else if (iter->parType == par_type::VEL_Z && iter->site.compare(station) == 0 &&
				iter->beg == Tbeg && iter->end == Tend) {
				vel.set(2, iter->value());
				found++;
			}
		}
		if (found == 3)  return  1;    // all three crd were found
		if (found == 1)  return -1;    // just one crd was found
		if (found == 2)  return -2;    // just two crd was found

		if (found == 0) return  -3;    // crd not found

		return -1;
	}

	vector<int> t_gallpar::getPartialIndex(string site, string sat) 
	{
		_update_partial_index();

		vector<int> ans;
		pair<string, string> type_list[4] =
		{
			make_pair(site,sat),
			make_pair(site,""),
			make_pair("",sat),
			make_pair("","")
		};

		_allpar_mtx.lock();
		for (int i = 0; i < 4; i++)
		{
			ans.insert(ans.end(), _index_for_parital[type_list[i]].begin(), _index_for_parital[type_list[i]].end());
		}
		_allpar_mtx.unlock();

		return ans;
	}

	t_gpar& t_gallpar::getPar(int idx)
	{
		return _vParam[idx];
	};

	// get single t_gpar element from container
	// ----------------------------------------------
	t_gpar& t_gallpar::operator[](const size_t idx)
	{
		return _vParam[idx];
	}

	// Operator -
	// ----------------------------------------
	t_gallpar t_gallpar::operator-(t_gallpar& gallpar)
	{
		t_gallpar diff;
		if (this->parNumber() != gallpar.parNumber()) {
			cerr << "t_gallpar::operator-: Incompatible dimension ("
				<< this->parNumber() << ", " << gallpar.parNumber() << ")"
				<< endl;
			return diff;
		}

		diff = (*this);
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {
			int i = gallpar.getParam(iter->site, iter->parType, iter->prn, iter->beg, iter->end);
			if (i >= 0) {
				diff[i] = (*iter) - gallpar[i];
			}
			//      else cerr << "NENASEL JSEM " << iter->str_type() << endl;
		}
		return diff;
	}

	// Operator +
	// ----------------------------------------
	t_gallpar t_gallpar::operator+(t_gallpar& gallpar)
	{
		t_gallpar diff;

		if (this->parNumber() != gallpar.parNumber()) {
			cerr << "t_gallpar::operator+: Incopatible dimension ("
				<< this->parNumber() << ", " << gallpar.parNumber() << ")"
				<< endl;
			return diff;
		}

		diff = (*this);
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {
			int i = gallpar.getParam(iter->site, iter->parType, iter->prn, iter->beg, iter->end);
			if (i >= 0) {
				diff[i] = (*iter) + gallpar[i];
			}

		}
		return diff;
	}

	//Doplnil Gabo
	void t_gallpar::delAllParam()
	{
		gtrace("t_gallpar::delAllParam");

		_vParam.clear();
		this->_index_par.clear();
		this->_point_par.clear();
		this->_max_point = 0;
		this->_last_point = make_pair(0, 0);
	}

	// Delete all ambiguity params
	// ------------------------------------
	vector<int> t_gallpar::delAmb()
	{
		gtrace("t_gallpar::addAmb");

		vector<int> ind;
		vector<t_gpar>::iterator iter;
		iter = _vParam.begin();
		while (iter != _vParam.end())
		{
			if (iter->parType == par_type::AMB_IF
				|| iter->parType == par_type::AMB_L1
				|| iter->parType == par_type::AMB_L2
				|| iter->parType == par_type::AMB_L3
				|| iter->parType == par_type::AMB_L4
				|| iter->parType == par_type::AMB_L5
				|| iter->parType == par_type::AMB_WL
				) {
				ind.push_back(iter->index);
				int i = iter - _vParam.begin();

				auto& all = this->_index_par[_vParam[i].get_head()];
				for (auto iter = all.begin(); iter != all.end(); iter++) {
					if (iter->second == _point_par[i]) {
						all.erase(iter);
						break;
					}
				}
				if (this->_index_par[_vParam[i].get_head()].size() == 0) {
					this->_index_par.erase(_vParam[i].get_head());
				}
				_point_par.erase(_point_par.begin() + i);

				iter = _vParam.erase(iter);
			}
			else iter++;
		}
		return ind;
	}

	// Reset all parems value
	// -------------------------------
	void t_gallpar::resetAllParam()
	{
		gtrace("t_gallpar::resetAllParam");

		vector<t_gpar>::iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++) iter->value(0);
	}


	// set site name for all pars in gallpar
	// -----------------------------------------
	void t_gallpar::setSite(string site)
	{
		gtrace("t_gallpar::setSite");

		vector<t_gpar>::iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)  iter->site = site;

	}

	set<string> t_gallpar::amb_prns(par_type type)
	{
		gtrace("t_gallpar::amb_prns");

		set<string> prns;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->parType == type)
			{
				prns.insert(iter->prn);
			}
		}
		return prns;

	}

	// Multiple matrix and gallpar
	// -----------------
	int t_gallpar::mult(const Matrix& K, t_gallpar& mult, t_gallpar& res)
	{
		gtrace("t_gallpar::mult");

		size_t parN = mult.parNumber();
		if (parN != K.Ncols()) {
			cerr << "t_gallpar::mult - incorrect dimension " << parN << " " << K.Ncols() << endl;
			return -1;
		}

		res = mult;

		double c = 0;
		int pos = -1;
		for (size_t i = 1; i <= K.Nrows(); i++) {
			for (size_t j = 1; j <= K.Ncols(); j++) {
				pos = mult.getParam(j);
				if (pos < 0) { return -1; }
				c += K(i, j) * mult[j - 1].value();
			}
			pos = mult.getParam(i);
			res[i - 1].value(c);
			c = 0;
		}

		return 1;
	}

	set<string> t_gallpar::amb_prns()
	{
		gtrace("t_gallpar::amb_prns");

		set<string> prns;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if(iter->str_type().find("AMB") != string::npos)
			{
				prns.insert(iter->prn);
			}
		}
		return prns;
	}

	set<string> t_gallpar::amb_recs(par_type type)
	{
		gtrace("t_gallpar::amb_recs");

		set<string> recs;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->parType == type)
			{
				recs.insert(iter->site);
			}
		}
		return recs;
	}

	set<string> t_gallpar::amb_recs()
	{
		gtrace("t_gallpar::amb_recs");

		set<string> recs;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->str_type().find("AMB") != string::npos)
			{
				recs.insert(iter->site);
			}
		}
		return recs;
	}

	map<string, string> t_gallpar::rec_sys_map()
	{		
		gtrace("t_gallpar::amb_sys");
		set<string> rec_list = amb_recs();
		map<string, string> sys_map;

		string rec;
		string sat;
		string sys;

		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->parType == par_type::AMB_IF ||
				iter->parType == par_type::AMB_L1 ||
				iter->parType == par_type::AMB_L2 ||
				iter->parType == par_type::AMB_L3 ||
				iter->parType == par_type::AMB_L4 ||
				iter->parType == par_type::AMB_L5)
			{
				sat = iter->prn;
				sys = iter->prn.substr(0, 1);
				rec = iter->site;
				if (sys_map.count(rec) <= 0)
				{
					sys_map[rec] = sys;
				}
				else
				{
					if (sys_map[rec].find(sys) == string::npos) sys_map[rec] += sys;
					sort(sys_map[rec].begin(), sys_map[rec].end());
				}

			}
		}
		return sys_map;
	}

	// -----------------------------
	ostream& operator<<(ostream& os, t_gallpar& x)
	{
		for (unsigned int i = 0; i < x.parNumber(); i++)
		{
			if (x[i].parType == par_type::AMB_IF ||
				x[i].parType == par_type::AMB_L1 ||
				x[i].parType == par_type::AMB_L2 ||
				x[i].parType == par_type::AMB_L3 ||
				x[i].parType == par_type::AMB_L4 ||
				x[i].parType == par_type::AMB_L5 ||
				x[i].parType == par_type::SION ||
				x[i].parType == par_type::VION ||
				x[i].parType == par_type::GPS_REC_IFB_C3||
				x[i].parType == par_type::IFCB_F3 ||
				x[i].parType == par_type::IFCB_F4 ||
				x[i].parType == par_type::IFCB_F5)
			{
				os << x[i].str_type() << "_" << x[i].prn << " ";
			}
			else
			{
				os << x[i].str_type() << " ";
			}

			if (x[i].parType == par_type::GRD_N || x[i].parType == par_type::GRD_E)
			{
				os << "value: " << x[i].value() * 1000 << " " << "index:" << x[i].index;
			}
			else
			{
				os << "value: " << x[i].value() << " " << "index:" << x[i].index;
			}
			os << endl;
		}
		return os;
	}

	// X1 + X2
	//---------------------------------
	int t_gallpar::sum(t_gallpar& X1, t_gallpar& X2)
	{
		if (X1.parNumber() != X2.parNumber()) return -1;

		for (unsigned int i = 0; i < _vParam.size(); i++) {
			int id1 = X1.getParam(_vParam[i].site, _vParam[i].parType, _vParam[i].prn);
			int id2 = X2.getParam(_vParam[i].site, _vParam[i].parType, _vParam[i].prn);
			if (id1 >= 0 && id2 >= 0) _vParam[i].value(X1[id1].value() + X2[id2].value());
			else return -1;
		}
		return 1;
	}

	// get ColumnVector w.r.t par indexes
	// ---------------------------
	ColumnVector t_gallpar::get_cvect(t_gallpar& par)
	{
		int n = par.parNumber();
		ColumnVector vec(n);

		for (int i = 1; i <= n; i++) {
			int idx = this->getParam(i);
			vec(i) = _vParam[idx].value();
		}
		return vec;
	}

	void t_gallpar::setOrbPar(string sat, vector<t_gpar> par)
	{
		this->_vOrbParam[sat] = par;
	}

	vector<t_gpar> t_gallpar::getOrbPar(string sat)
	{
		return _vOrbParam[sat];
	}

	vector<t_gpar> t_gallpar::getAllPar()
	{
		return _vParam;
	}

	vector<t_gpar> t_gallpar::getAmbPar(string sat, string rec, par_type type)
	{
		vector<t_gpar> tmp;

		for (const auto& par : _vParam) {
			if (par.prn == sat && par.site == rec && par.parType == type) tmp.push_back(par);
		}

		return tmp;
	}
	map<string, vector<t_gpar>> t_gallpar::getAllOrbPar()
	{
		return _vOrbParam;
	}
	void t_gallpar::_update_partial_index()
	{
		_allpar_mtx.lock();
		auto point_now = make_pair(_point_par[_point_par.size() - 1], int(_point_par.size() - 1));

		if (point_now!= _last_point)
		{
			_last_point = point_now;
			_index_for_parital.clear();
			set<string> site_list, sat_list;

			for (int i = 0; i < _vParam.size(); i++)
			{
				_index_for_parital[make_pair(_vParam[i].site, _vParam[i].prn)].push_back(i);
			}
		}
		_allpar_mtx.unlock();
		return;
	}


	// Get Single/double/triple frequency satellites, separately
	// --------------
	map<string, int> t_gallpar::freq_sats_num(int freq)
	{
		gtrace("t_gallpar::freq_sats_num");

		set<string> prns1, prns2, prnsif, prns13if, prns14if, prns15if;
		map<string, int>  Nsats;
		vector<t_gpar>::const_iterator iter;
		for (iter = _vParam.begin(); iter != _vParam.end(); iter++)
		{
			if (iter->parType == par_type::AMB_L1)
			{
				prns1.insert(iter->prn);
			}
			else if (iter->parType == par_type::AMB_L2)
			{
				prns2.insert(iter->prn);
			}
			else if (iter->parType == par_type::AMB_IF)
			{
				prnsif.insert(iter->prn);
			}
		}

		if (prnsif.size() != 0 || prns13if.size() != 0)
		{
			set<string> tmp;
			set_intersection(prnsif.begin(), prnsif.end(), prns13if.begin(), prns13if.end(), inserter(tmp, tmp.begin()));
			Nsats["Triple"] = tmp.size();
			Nsats["Double"] = prnsif.size() + prns13if.size() - 2 * tmp.size();
		}
		else if (prns1.size() != 0 || prns2.size() != 0)
		{
			set<string> tmp;
			set_intersection(prns1.begin(), prns1.end(), prns2.begin(), prns2.end(), inserter(tmp, tmp.begin()));
			Nsats["Double"] = tmp.size();
			Nsats["Single"] = prns1.size() + prns2.size() - 2 * tmp.size();
		}

		return Nsats;
	}


} // namespace
