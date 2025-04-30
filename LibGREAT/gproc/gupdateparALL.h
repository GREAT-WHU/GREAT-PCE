/**
 * @file         gupdateparALL.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GUPDATEPARALL_H
#define GUPDATEPARALL_H

#include <map>
#include "gproc/gupdatepar.h"

using namespace gnut;

namespace great
{
	/**
	* @brief  class for update all parameters information
	*/
	class LibGREAT_LIBRARY_EXPORT t_gupdateparALL :public t_gupdatepar
	{
	public:
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip, const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index, const IONO_ORDER& iono_order);
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index, const IONO_ORDER& iono_order);
		t_gupdateparALL(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13, shared_ptr<t_gcycleslip> cycleslip_14, shared_ptr<t_gcycleslip> cycleslip_15, map< GSYS, map<FREQ_SEQ, GOBSBAND> > band_index, const IONO_ORDER& iono_order);
		virtual ~t_gupdateparALL();

		t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata) override;

	protected:

		void _update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)override;
		void _lite_update_amb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		void _udpate_ion_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);

		void _check_log(shared_ptr<t_gcycleslip> cycleslip_12, shared_ptr<t_gcycleslip> cycleslip_13);

		shared_ptr<t_gcycleslip>              _cycleslip_13 = nullptr;
		shared_ptr<t_gcycleslip>              _cycleslip_14 = nullptr;
		shared_ptr<t_gcycleslip>              _cycleslip_15 = nullptr;

		IONO_ORDER                            _iono_order = IONO_ORDER::NONE;
		int                                   _freq;
	};


}
#endif // !GUPDATEPARALL_H
