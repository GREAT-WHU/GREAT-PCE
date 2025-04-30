/**
 * @file         gcombmodel.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        base combine biase model
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GCOMBMODEL_H
#define GCOMBMODEL_H

#include "gexport/ExportLibGREAT.h"
#include "gmodels/gbasemodel.h"
#include "gmodels/gbiasmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gall/gallbias.h"
#include "gall/gallproc.h"
#include "gutils/gcommon.h"
#include <map>
#include <memory>

using namespace gnut;
using namespace std;

namespace great
{
    /** @brief Class for combine biase model */
    class LibGREAT_LIBRARY_EXPORT t_gcombmodel : virtual public t_gbasemodel
    {
    public:
        /** @brief default constructor.
        *
        *param[in]    log                log control
        *param[in]    settings        setbase control
        *param[in]    bias_model        model of bias
        *param[in]    data            all data
        */
        t_gcombmodel(t_gsetbase* setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        virtual ~t_gcombmodel();
        /** @brief get index of band */
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > get_band_index();
        /** @brief get index of frequency */
        map< GSYS, map<GOBSBAND, FREQ_SEQ> > get_freq_index();
    protected:
		bool _wgt_raw_obs(const t_gobs& gobs, const t_gsatdata& satdata, const double& factorP, const OBSWEIGHT& wgt_type, double& wgt);

        map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;

        shared_ptr<t_gbiasmodel> _bias_model; ///< baise model

        t_gprecisemodel* _gprecisemodel = nullptr;
        t_gallbias*      _gallbias      = nullptr;
        t_glog*          _glog           = nullptr;

        int             _frequency;         ///< frequency
		bool            _update_amb_lite;
		

	protected:
		double _sigCodeGPS;		///< code bias of GPS
		double _sigCodeGLO;		///< code bias of GLO
		double _sigCodeGAL;		///< code bias of GAL
		double _sigCodeBDS;		///< code bias of BDS
		double _sigCodeQZS;		///< code bias of QZS

		double _sigPhaseGPS;	///< phase bias of GPS
		double _sigPhaseGLO;	///< phase bias of GLO
		double _sigPhaseGAL;	///< phase bias of GAL
		double _sigPhaseBDS;	///< phase bias of BDS
		double _sigPhaseQZS;	///< phase bias of QZS

    };

    class LibGREAT_LIBRARY_EXPORT t_gcombIF:virtual public t_gcombmodel
    {
    public:
        /** @brief constructor. */
        t_gcombIF(t_gsetbase *setting,t_glog* log,shared_ptr<t_gbiasmodel> bias_model,t_gallproc* data);
        ~t_gcombIF();
        /** @brief Combine equation */
        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gbaseEquation& result) override;
        bool cmb_equ_IF(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,GOBSBAND b1,GOBSBAND b2, t_gbaseEquation& result);

    };
	
}
#endif