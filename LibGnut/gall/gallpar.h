/**
*
* @verbatim
    History
*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     gallpar.h
* @brief    Purpose: parametres container
*
* @author   PV
* @version  1.0.0
* @date     2011-04-18
*
*/

#ifndef GALLPAR_H
#define GALLPAR_H 


#include <string>
#include <set>

#include "gexport/ExportLibGnut.h"
#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gmodels/gpar.h"
#include <unordered_map>
#include <algorithm>

using namespace std;

namespace gnut
{
	/**
	 *@brief Class for t_gallpar
	 */
	class LibGnut_LIBRARY_EXPORT  t_gallpar
	{
	public:
		/**
		*@brief add parameter
		*/
		void addParam(t_gpar);
		/**
		*@brief delete parameter
		*/
		void delParam(int i);
		/**
		*@brief delete all parameter
		*/
		void delAllParam();
		/**
		*@brief reset all parameter
		*/
		void resetAllParam();

		set<int> delEmpty();
		vector<int> delAmb();
		/**
		*@brief get parameter
		*/
		int getParam(const string& site, const par_type& par, const string& prn,
			const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME, const int& channal = DEF_CHANNEL) const;
		int getParam(int index);
		/**
		*@brief get par index
		*/
		int getParIndex(int idx);
		/**
		*@brief get par value
		*/
		double getParValue(int idx);
		/**
		*@brief get par sigma
		*/
		double getParSig(int idx);
		/**
		*@brief set par value
		*/
		void setParValue(int idx, double value);
		/**
		*@brief set par sigma
		*/
		void setParSig(int idx, double value);
		void setParRemove(int idx, bool judge);
		/**
		*@brief get par/amb/orbPar Number
		*/
		unsigned int parNumber() const;
		unsigned int ambNumber() const;
		/**
		 * @brief
		 *
		 * @return unsigned int
		 */
		unsigned int orbParNumber() const;
		int maxIndex() const;
		/**
		*@brief get Amb Param
		*/
		int getAmbParam(string site, string prn, par_type  type, 
			const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME) const;
		t_gpar& getAmbParam(int idx);
		/**
		*@brief get Crd Param
		*/
		int getCrdParam(string site, t_gtriple& crd,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
		/**
		*@brief get Vel Param
		*/
		int getVelParam(string site, t_gtriple& crd,
			t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
		/**
		*@brief get partial Index
		*/
		vector<int> getPartialIndex(string site, string sat);

		ColumnVector get_cvect(t_gallpar& par);
		/**
		 * @brief Get the Par object
		 *
		 * @param idx
		 * @return const t_gpar&
		 */
		t_gpar&   getPar(int idx);
		t_gpar&   operator[](const size_t idx);
		t_gallpar operator-(t_gallpar& par);
		t_gallpar operator+(t_gallpar& par);
		static int mult(const Matrix& K, t_gallpar& par, t_gallpar& res);
		void   reIndex();
		void   decIndex(int i);
		void   incIndex(int i);
		/**
		 * @brief Set the Site object
		 *
		 * @param site
		 */
		void   setSite(string site);
		set<string> amb_prns(par_type type);
		set<string> amb_prns();
		set<string> amb_recs(par_type type);
		set<string> amb_recs();
		map<string, string> rec_sys_map();
		/**
		 * @brief
		 *
		 * @param X1
		 * @param X2
		 * @return int
		 */
		int sum(t_gallpar& X1, t_gallpar& X2);

		void setOrbPar(string sat, vector<t_gpar> pars);
		/**
		 * @brief Get the All Orb object
		 *
		 * @return vector<t_gpar>
		 */
		vector<t_gpar>  getOrbPar(string sat);
		/**
		 * @brief Get the All Par object
		 *
		 * @return vector<t_gpar>
		 */
		vector<t_gpar>  getAllPar();
		vector<t_gpar>  getAmbPar(string sat, string rec, par_type type);
		map<string, vector<t_gpar>> getAllOrbPar();
		/**
		 * @brief
		 *
		 * @param os
		 * @param x
		 * @return ostream&
		 */
		friend ostream& operator<<(ostream& os, t_gallpar& x);
		/**
		 * @brief
		 *
		 * @param freq
		 * @return map<string, int>
		 */
		map<string, int> freq_sats_num(int freq);

	private:
		vector<t_gpar> _vParam;
		map<string, vector<t_gpar>> _vOrbParam;

		// Fast get param index
		unordered_map<t_gparhead, map<t_gtimearc, long >, t_gparhead_hash > _index_par;
		long _max_point=0;
		vector<long> _point_par;

		// Fast get parital  
		unordered_map<pair<string, string>, vector<int>, t_gpair_string_hash > _index_for_parital;
		pair<long, int> _last_point;
		/** @brief update partial index. */
		void _update_partial_index();

		t_gmutex _allpar_mtx;
	};

} // namespace

#endif
