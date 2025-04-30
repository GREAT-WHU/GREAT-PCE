/**
 * @file         gupdateparIF.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GUPDATEPARIF_H
#define GUPDATEPARIF_H

#include <map>
#include "gproc/gupdatepar.h"

namespace great {
	/**
	* @brief  class for update IF parameters information
	*/
	class LibGREAT_LIBRARY_EXPORT t_gupdateparIF :public t_gupdatepar
	{
	public:
		t_gupdateparIF(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index) {
			_cycleslip = cycleslip;
			_band_index = band_index;
		}
		virtual ~t_gupdateparIF() = default;

		void set_break_amb_by_obs(bool isSet) { isBreakObs = isSet; };
	protected:

		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) override;

		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata,par_type amb_type,t_gcycleslip* _cycleslip,t_gupdateparinfo& update_info);
		void _update_one_type_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, par_type amb_type, t_gupdateparinfo& update_info);
	
	protected:
		map<pair<string, string>, int> _amb_id;
		bool isBreakObs                = true;

	};

	class LibGREAT_LIBRARY_EXPORT t_gupdateparIF_1X :public t_gupdateparIF
	{
	public:
		t_gupdateparIF_1X(shared_ptr<t_gcycleslip> cycleslip12,shared_ptr<t_gcycleslip> cycleslip13, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index);
		virtual ~t_gupdateparIF_1X();

		void set_cycleslip13(shared_ptr<t_gcycleslip> cycleslip13);
		void set_cycleslip14(shared_ptr<t_gcycleslip> cycleslip14);
		void set_cycleslip15(shared_ptr<t_gcycleslip> cycleslip15);
		t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata) override;

	protected:

		void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info) override;


		shared_ptr<t_gcycleslip> _cycleslip13 = nullptr;
		shared_ptr<t_gcycleslip> _cycleslip14 = nullptr;
		shared_ptr<t_gcycleslip> _cycleslip15 = nullptr;
	};

	class LibGREAT_LIBRARY_EXPORT t_updateparWL : public t_gupdateparIF
	{
	public:
		t_updateparWL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index);
		virtual ~t_updateparWL();

	protected:
		virtual void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info) override;

	};

}
#endif // !GUPDATEPARIF_H
