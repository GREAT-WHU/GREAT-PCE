/**
 * @file         gprecisemodel.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        mainly about how to cacultae B P l single
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GPRECISEMODEL_H
#define	GPRECISEMODEL_H

#include "gexport/ExportLibGREAT.h"
#include "gmodels/gpppmodel.h"
#include "gmodels/gtideIERS.h"
#include "gproc/glsqmatrix.h"
#include "gproc/glsq.h"
#include "gall/gallproc.h"
#include "gdata/gpoleut1.h"
#include "gutils/gtrs2crs.h"
#include "gdata/gnavde.h"
#include "gset/gsetproc.h"

namespace great
{
	/**
	* @brief precise correction model
	*/
	class LibGREAT_LIBRARY_EXPORT t_gprecisemodel : public t_gpppmodel
	{
	public:

		/**
		* @brief constructors
		* @param[in] data all proc data
		* @param[in] log log file for precise model
		* @param[in] set xml setting for precise model
		*/
		t_gprecisemodel(t_gallproc* data, t_glog* log, t_gsetbase* setting);

		/**
		* @brief default destructor
		*/
		virtual ~t_gprecisemodel();

		/**
		* @brief set otl data for tide correction
		* @param[in] all ocean tide data
		*/
		void setOTL(t_gallotl* otl);

		void set_multi_debug_output(string filename);

		bool _prepare_obs(const t_gtime& epoch, t_gallpar& pars);
		bool _omc_obs_ALL(const t_gtime& crt_epo, t_gallpar& pars, t_gobs& gobs,  double& omc);
		bool _wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt);
		bool _wgt_obs_ALL(t_gobs& gobs1, double factorP, double& wgt, double snr_value); // add wh
		bool _prt_obs(const t_gtime& epoch, t_gallpar& pars, t_gobs & gobs, vector<pair<int, double> >&coeff);

		bool _update_obs_info(const t_gtime& epoch, t_gsatdata& obsdata, t_gallpar& pars);
		bool _update_obs_info(t_gsatdata& obsdata);

		void  update_obj_clk(const string& obj, const t_gtime& epo, double clk);
		bool _update_obj_clk(const string& obj, const t_gtime& epo, t_gallpar& par, double& clk, map<string, pair<t_gtime, double> >& obj_clk);

		int outlierDetect(vector<t_gsatdata>& obsdata);


		/**
		* @brief compute windup correciton
		* @note override for fixed bugs in computing satellite attitude
		* @param[in] satdata observ data info
		* @param[in] rRec coord of receiver in TRS
		* @return correction of windup
		*/
		double windUp(t_gsatdata& satdata, const ColumnVector& rRec) override;


		/**
		* @brief combine obs equations
		* @note  override for setting site
		* @param[in] epoch proc epoch
		* @param[in] param all estimate param
		* @param[in] param observ data info
		* @param[in] gobs observ type
		* @param[in] com flag for computing pcv pco
		* @return correction of model obs
		*/
		//double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs, bool com) override;

		double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& obsdata, t_gobs& gobs);
		/**
		* @brief get dry tropo delay
		* @note overide for setting site
		* @param[in] site specified site
		* @param[in] epo  specified epoch
		* @return correction of dry delay
		*/
		double getZHD(const string& site, const t_gtime& epo) override;


		/**
		* @brief get wet torpo delay
		* @note overide for setting site
		* @param[in] site specified site
		* @param[in] epo  specified epoch
		* @return correciton of wet delay
		*/
		double getZWD(const string& site, const t_gtime& epo) override;


		/**
		* @brief get torpo delay
		* @note override ofr leo sattlite
		* @param[in] epoch specified epoch
		* @param[in] param all estimate parameters
		* @param[in] site_ell ell coord of reciver
		* @param[in] satdata observ data info
		* @return delay of tropo
		*/
		double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata)override;
		double ionDelay(t_gtime& epoch, t_gallpar& param, t_gsatdata& satdata, IONMODEL& ion_model, GOBSBAND& band_1st, t_gobs& gobs);
		
	public:
		/**
		 * @brief judge the validity of observation
		 * @param[in] obsdata observ data info
		 * @param[in] obscom combine type of observ
		 * @return true for sucess false for fail
		 */
		bool _check_obs_valid(t_gobsgnss& obsdata,const OBSCOMBIN& obscom, const GOBSBAND& band1, const GOBSBAND& band2);

		/**
		* @brief get rotmatrix from ant to TRS or CRS
		* @param[in] obsdata observ data info
		* @param[in] epoch specified epoch
		* @param[in] obj pcv info for site
		* @param[in] isCRS ture to CRS false to TRS
		* @return Rotmatrix
		*/
		Matrix _RotMatrix_Ant(t_gsatdata& obsdata, const t_gtime& epoch, shared_ptr<t_gobj> obj, bool isCRS);

		/**
		* @brief update TRS2CRS rotmatrix
		* @param[in] epoch specified epoch
		*/
		void _update_rot_matrix(const t_gtime& epoch);

		/**
		* @brief get the site crd after arp,tide correction in TRS and CRS
		* @param[in] epoch_Rec epoch of site receive signal
		* @param[in] param all estimate parameters
		* @return true for sucess false for fail
		*/
		bool _apply_rec(const t_gtime& crt_epo, const t_gtime& rec_epo, t_gallpar& pars);

		/**
		* @brief get the sat crd from rinexn,sp3,orb file
		* @param[in] satdata observ data info
		* @param[in] epoch_Rec give the epoch of site receive signal
		* @param[out] epoch_Sat get the epoch of sat send signal
		* return true for sucess  false for fail
		*/
		bool _apply_sat(const t_gtime& rec_epo, t_gtime& sat_epo);


		/**
		* @brief compute tides correction for site
		* @param[in] epoch specifie epoch
		* @param[out] rec coord of site after adding correction
		* @return ture for correciont sucess false for fail
		*/
		bool _apply_rec_tides(const t_gtime& epoch, t_gtriple& rec);

		/**
		* @brief iterator for computing sat pos
		* @param[in] epoch_Sat epoch of sat send signal
		* @param[in] sat_name  name of sat
		* @param[out] sat_crd get coord of sat in CRS
		* @return true for sucess false for fail
		*/
		bool _get_crs_sat_crd(const t_gtime& sat_epoch, const string& sat, t_gtriple& crs_sat_crd);
		bool _get_crs_sat_crd(const t_gtime& sat_epoch, const string& sat, t_gallnav*  nav, t_gtriple& crs_sat_crd);

		/**
		* @brief compute sat vel
		* @param[in] epoch_Sat epoch of sat send signal
		* @param[in] sat_name  name of sat
		* @param[out] sat_crd get vel of sat in CRS
		* @return true for sucess false for fail
		*/
		bool _get_crs_sat_vel(const t_gtime& epoch_Sat, string sat_name, t_gtriple& sat_vel);

		/**
		* @brief compute coeff of different par
		* @param[in] par get the specified type of par
		* @param[in] epoch the pecified epoch
		* @param[in] obsdata observ data info
		* @param[in] gobs observ type
		* return get the partial value of par
		*/	
		//double _Partial(const t_gtime& epoch, t_gsatdata& obsdata, const t_gobs & gobs1, const t_gobs& gobs2, t_gpar& par);
		double _Partial(const t_gtime& epoch, t_gsatdata& obsdata, const t_gobs & gobs, t_gpar& par);

		void   _getmf(t_gpar & par, t_gsatdata & satData, const t_gtriple & crd, const t_gtime & epoch, double & mfw, double & mfh, double & dmfw, double & dmfh);

		t_gallnav*      _gall_nav = nullptr; 		  ///< all nav data include rinexn,sp3,clk
		t_gpoleut1*     _gdata_erp    = nullptr;  	  ///< all poleut1 data
		t_gnavde *      _gdata_navde  = nullptr;	  ///< all panetnav info
		shared_ptr<t_gtide> _tide ;				      ///< tide correction model

		map<t_gtime, shared_ptr<t_gtrs2crs> > _trs2crs_list;
		shared_ptr<t_gtrs2crs>    _trs2crs_2000;	  ///< trs2crs matrix

		double _minElev;  		///< min ele for prepare 
		OBSWEIGHT _weight;		///< weight calculation method
		ZTDMPFUNC _mf_ztd;      ///< mapping function for ZTD
		GRDMPFUNC _mf_grd;      ///< mapping function for GRD 

		map<string, map<string, map<t_gtime, double> > > _phase_windup; ///< recording for calculating windup

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

		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;
		IONMODEL                             _ion_model;
		GSYS                _crt_sys;
		t_gsatdata			_crt_obs;
		string				_crt_rec;
		string				_crt_sat;


		CONSTRPAR           _crd_est = CONSTRPAR::FIX;

		int       base_size;
		bool      _bds2_isb = false;

	public:
		double              _crt_rec_clk;
		double              _crt_sat_clk;
		
		t_gtime				_crt_epo;
		
		shared_ptr<t_gobj>	_crt_obj;
		t_gdata::ID_TYPE    _crt_obj_type;
		map<string, pair<t_gtime, double> > _obj_clk;
		map<string, double> _rec_clk;

		t_gtime _rec_epo;
		t_gtime _sat_epo;
		t_gtriple _trs_rec_crd;
		t_gtriple _crs_rec_crd;
		t_gtriple _trs_sat_crd;
		t_gtriple _crs_sat_crd;
		t_gtriple _crs_rec_pco;
		t_gtriple _crs_sat_pco;
		t_gtriple _crs_rec_vel;
		t_gtriple _crs_sat_vel;

		Matrix _rec_rotmat;
		Matrix _dr_dxpole;
		Matrix _dr_dypole;
		Matrix _dr_dut1;
		Matrix _sat_partial;
		Matrix _rec_partial;

		int _sat_index;
		int _rec_index;

		pair<string, t_gtime> _rec_obj_flag;
		double _factorP;
		double _factorL;
		bool _isCalPCO = true;
	};

	/**
	* @brief compute DCB correction for specified obs
	* @param[in] allbias code bias info for caculation
	* @param[in] satdata observ data info
	* @param[in] P obs value
	* @param[in] obscombin combine type of observ
	* @param[in] gobs1 get the first observ type of obs
	* @param[in] gobs2 get the second observ type(only in IONO FREE model)
	* @return true for sucess false for fail
	*/
	bool LibGREAT_LIBRARY_EXPORT correct_DCB(string ac, t_gallbias* allbias, t_gsatdata& satdata, double& P, OBSCOMBIN obscombin, t_gobs* gobs1, t_gobs* gobs2 = 0, t_gobs* gobs3 = 0);
	
	/**
	* @brief compute relative delay
	* @note all coord and vel must be in the same coordinate system
	* @param[in] crd_site  coord of site
	* @param[in] crd_sat coord of sat
	* @param[in] vel_site  velocity of site
	* @param[in] vel_sat velocity of sat
	* @param[out] reldelay get the delay of relative
	* @return true for sucess false for fail
	*/
	bool apply_reldelay(t_gtriple crd_site, t_gtriple vel_site, t_gtriple crd_sat, t_gtriple vel_sat, double& reldelay);


	/**
	* @brief check code range res
	* @param[in] omc all compute res of equations
	* @param[out] recclk  add average of res to recclk
	* @param[out] sigma  sigma of res
	* @return number of use res
	*/
	int LibGREAT_LIBRARY_EXPORT Check_Range(vector<double>& omc, double& recclk, double& sigma);

};
#endif /* PGRECISEMODEL_H */