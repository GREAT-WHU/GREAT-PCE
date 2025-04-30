/**
  *
  * @verbatim
    History
    -1.0    Jiande Huang : Optimizing the efficiency of "get" function
  *
  @endverbatim
  * Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  *
  * @file     gallbias.h
  * @brief    container for all biases
  *
  * @author   JD
  * @version  1.0.0
  * @date     2012-11-05
  *
  */

#ifndef GALLBIAS_H
#define GALLBIAS_H


#include "gdata/gdata.h"
#include "gmodels/gbias.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut
{

	typedef shared_ptr<t_gbias> t_spt_bias;

	/**
	 *@brief Class for bias setting derive from t_gdata
	 */
	class LibGnut_LIBRARY_EXPORT t_gallbias : public t_gdata
	{

		typedef map<GOBS, t_spt_bias>   t_map_gobs;     ///< first : GNSS Observations, second : bias
		typedef map<string, t_map_gobs> t_map_sat;    ///< first : sat name, second : observations
		typedef map<t_gtime, t_map_sat> t_map_epo;    ///< first : time, second : sat data
		typedef map<string, t_map_epo>  t_map_ac;      ///< first : ac name, second : epoch

	public:
		/** @brief default constructor. */
		t_gallbias();
		/** @brief default destructor. */
		virtual ~t_gallbias();

		/**
		* @brief set single bias element value.
		*
		* @param[in]  ac        ac name of the data
		* @param[in]  epo       epoch of the data
		* @param[in]  obj        object of the data
		* @param[in]  pt_bias    pt bias of the data
		* @return void
		*/
		void   add(const string& ac, const t_gtime& epo, const string& obj, t_spt_bias pt_bias);
		/**
		* @brief get DCB.
		*
		* @param[in]  epo       epoch of the data
		* @param[in]  obj        object of the data
		* @param[in]  gobs1        first observation of the data
		* @param[in]  gobs2        second observation bias of the data
		* @param[in]  ac            ac of the data
		* @return DCB
		*/
		double get(const t_gtime& epo, const string& obj, const GOBS& gobs1, const GOBS& gobs2, string ac = "");
		/**
		* @brief get single bias element.
		*
		* @param[in]  prd
		* @param[in]  epo       epoch of the data
		* @param[in]  obj        object of the data
		* @param[in]  gobs1        observation of the data
		* @param[in]  meter        unit
		* @return single bias
		*/
		double get(const string prd, const t_gtime& epo, const string& obj, const GOBS& gobs1, const bool meter = true);
		/**
		* @brief get ac list.
		* @return ac list
		*/
		vector<string> get_ac();
		/**
		* @brief set ac priority.
		* @return priority of ac
		*/
		string ac_priority();
		/**
		* @brief set used av.
		* @return
		*/
		void   set_used_ac(string ac);
		/**
		* @brief get used av.
		* @return used ac
		*/
		string get_used_ac();
		/**
		* @brief judge is osb.
		* @return
		*/
		bool is_osb();
		/**
		* @brief set overwrite mode
		* @return
		*/
		void overwrite(bool b) { _overwrite = b; }
		/**
		* @brief get overwrite mode
		* @return
		*/
		bool overwrite() { return _overwrite; }
		/**
		* @brief clean outtime
		* @return
		*/
		void clean_outer(const t_gtime& beg, const t_gtime& end);
		/**
		* @brief clean all data
		* @return
		*/
		void clean_all() { _mapbias.clear(); }
	protected:
		/**
		* @brief get single bias element pointer.
		*
		* @param[in]  ac            ac of the data
		* @param[in]  epo       epoch of the data
		* @param[in]  obj        object of the data
		* @param[in]  gobs        observation of the data
		* @return    pt bias
		*/
		t_spt_bias _find(const string& ac, const t_gtime& epo, const string& obj, const GOBS& gobs);   // get single bias element pointer
		/**
		* @brief get single bias element pointer.
		*
		* @param[in]  ac            ac of the data
		* @param[in]  epo       epoch of the data
		* @param[in]  obj        object of the data
		* @param[in]  ref
		* @return    vec bias
		*/
		vector<t_spt_bias> _find_ref(const string& ac, const t_gtime& epo, const string& obj, const GOBS& ref);
		/**
		* @brief convert type of observations.
		*
		* @param[in]  ac            ac of the data
		* @param[in]  obj        object of the data
		* @param[in]  obstype    observation type
		* @return    void
		*/
		void	_convert_obstype(const string& ac, const string& obj, GOBS& obstype);
		/**
		* @brief connect DCB pt_cb2 with first GOBS.
		*
		* @param[in]  pt_cb1
		* @param[in]  pt_cb2
		* @return    void
		*/
		void       _connect_first(const t_spt_bias pt_cb1, t_spt_bias pt_cb2);
		/**
		* @brief connect DCB pt_cb2 with second GOBS.
		*
		* @param[in]  pt_cb1
		* @param[in]  pt_cb2
		* @return    void
		*/
		void       _connect_second(const t_spt_bias pt_cb1, t_spt_bias pt_cb2);
		/**
		* @brief consolidate all biases with reference signal of pt_cb2.
		*
		* @param[in]  ac            ac of the data
		* @param[in]  obj        object of the data
		* @param[in]  pt_cb1
		* @param[in]  pt_cb2
		* @return    void
		*/
		void       _consolidate(const string& ac, const string& obj, const t_spt_bias pt_cb1, t_spt_bias pt_cb2);

		t_map_ac   _mapbias;			// map of all satellite biases (all ACs & all period & all objects)   
		string     _used_ac;			//
		bool       _overwrite = false;	// flag of overwrite

		map<string, int> _ac_order;		// map of all ACs
		bool             _isOrdered = false;	// if AC is ordered
		string           _pri_ac = "DLR_R";		// primary AC
	};

} // namespace

#endif
