/**
 * @file         gbiasmodel.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        mainly about how to cacultae B P l single
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GBIASEMODEL_H
#define GBIASEMODEL_H

#include "gexport/ExportLibGREAT.h"
#include "gutils/gtime.h"
#include "gall/gallpar.h"
#include "gdata/gsatdata.h"
#include "gdata/gobsgnss.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/gbasemodel.h"

using namespace std;
using namespace gnut;

namespace great
{
    /**
    *@brief t_gbiasmodel Class for bias model
    */
    class LibGREAT_LIBRARY_EXPORT t_gbiasmodel
    {
    public:
        /** @brief default constructor.*/
        t_gbiasmodel();
        virtual ~t_gbiasmodel();
        /** @brief Combined equation */
        virtual bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result)=0;
        /** @brief empty function */
        virtual void update_obj_clk(const string& obj, const t_gtime& epo, double clk) = 0;

        virtual t_gmodel* precisemodel() { return nullptr; }
    };
    /**
    *@brief t_gbiasmodel Class for precise bias model
    */
    class LibGREAT_LIBRARY_EXPORT t_gprecisebias:public t_gbiasmodel
    {
    public:
        /** @brief default constructor.*/
        t_gprecisebias(t_gallproc *data, t_glog *log, t_gsetbase *setting);
        ~t_gprecisebias();

		void set_multi_thread(const set<string>& sites);
        /** @brief Combined equation */
        bool cmb_equ(t_gtime& epoch,t_gallpar& params,t_gsatdata& obsdata,t_gobs& gobs,t_gbaseEquation& result) override;
        void update_obj_clk(const string& obj, const t_gtime& epo, double clk) override;

        t_gmodel* precisemodel() { return &gmodel; };
    protected:
        t_gprecisemodel gmodel;
		tuple<string, string, t_gtime> _rec_sat_before;

		bool _multi_thread_flag;
		map<string, shared_ptr<t_gprecisemodel> > _map_site_model;
		map < string, tuple<string, string, t_gtime> > _map_flag;
		t_gmutex _map_mtx;

    private:
        t_glog *_log;
    };

    class LibGREAT_LIBRARY_EXPORT t_gsppbias:public t_gbiasmodel
    {

    };

}

#endif