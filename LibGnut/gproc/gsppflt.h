/**
*
* @verbatim
	History
	2014-11-24  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*.
* @file       gsppflt.h
* @brief      Purpose: implements spp client
*.
* @author     PV
* @version    1.0.0
* @date       2014-11-24
*
*/

#ifndef GSPPFLT_H
#define GSPPFLT_H

#include "gproc/gspp.h"
#include "gmodels/gdop.h"
#include "gproc/gflt.h"
#include "gall/gallpar.h"
#include "gall/gallbias.h"
#include "gall/gallprod.h"
#include "gutils/gsysconv.h"
#include "gmodels/gstochasticmodel.h"
#include "gset/gsetflt.h"

namespace gnut
{

#define SPP_MINSAT  (6)    // minimum number of satellites
	/** @brief class for t_gsppflt derive from t_gspp. */
	class LibGnut_LIBRARY_EXPORT t_gsppflt : public virtual t_gspp
	{
	public:
		/** @brief constructor 1. */
		t_gsppflt(string mark, t_gsetbase* set);
		virtual ~t_gsppflt();
		/** @brief Process observation batch. */
		virtual int processBatch(const t_gtime& beg, const t_gtime& end);
		/** @brief min sat. */
		void minsat(size_t minsat) { _minsat = (minsat < 5) ? 5 : minsat; }
		/** @brief get crd. */
		t_gtriple getCrd(const t_gtime& time);
		/** @brief get outlier sat. */
		vector<string> get_outlier_sat() { return _outlier_sat; }

	protected:
		// Predict state vector and covariance matrix
		virtual void   _predict();
		// Restore state and covariance matrix
		virtual void   _restore(const SymmetricMatrix& Qsav, const t_gallpar& Xsav);
		// Satelite position
		virtual int    _satPos(t_gtime&, t_gsatdata&);
		// Prepare data: filter, bancroft, members in gsatdata
		virtual int    _prepareData();
		// Tides are not applied in SPP
		virtual int    _apply_tides(t_gtime& _epoch, t_gtriple& xRec);
		// Process one epoch
		virtual int    _processEpoch(const t_gtime& runEpoch);
		int    _addObsP(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
		// add data to A, l, P - carrier phase
		virtual int    _addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
		// add pseudo observations as a constrains
		virtual int    _addPseudoZTD(unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
		// weight coef for observations  
		double   _weightObs(t_gsatdata& satdata, t_gobs& go);
		// Update time for RDW processes
		virtual void   _timeUpdate(const t_gtime& epo);
		// Sync inter-GNSS systems bias
		void   _syncSys();
		// Add/Remove ionosphere delay
		void   _syncIono();
		// Add/Remove inter-freq. biases
		void   _syncIFB();
		// save observations residuals
		void   _save_residuals(ColumnVector& v, vector<t_gsatdata>& satdata, RESIDTYPE restype);
		// Apply SICB correction
		double _applySICB(string prn, double elev, GOBSBAND freq);
		// Apply DCB correction
		int    _applyDCB(t_gsatdata& satdata, double& P, t_gobs* gobs1, t_gobs* gobs2 = 0);

		vector<t_gsatdata> _data; ///< data
		map<string, int> _newAMB; ///< newAMB
		set<string> _slips;       ///< slips
		unsigned int _minsat;     ///< minsat
		double _sig_unit;         ///< sig unit
		int _frequency;           ///< frequency


		// Models   
		t_randomwalk*       _trpStoModel;
		t_randomwalk*       _ionStoModel;
		t_randomwalk*       _gpsStoModel;
		t_randomwalk*       _gloStoModel;
		t_randomwalk*       _galStoModel;
		t_randomwalk*       _bdsStoModel;
		t_randomwalk*       _qzsStoModel;
		t_whitenoise*       _clkStoModel;
		t_whitenoise*       _crdStoModel;


		// Parameters and covariance matrix   
		t_gallpar           _param;
		SymmetricMatrix     _Qx;

		// Noise matrix
		DiagonalMatrix      _Noise;

		// Epoch time
		t_gtime           _epoch;

		// Estimation objects   
		t_gflt*           _filter;

		int               _numSat(GSYS gsys);
		ColumnVector      _vBanc;
		t_gdop            _dop;
		int               _cntrep;
		int               _numcor;
		bool              _smooth;
		unsigned int      _n_NPD_flt; ///< NPD flt number
		unsigned int      _n_ALL_flt; ///< NPD smt number
		unsigned int      _n_NPD_smt; ///< NPD smt number
		unsigned int      _n_ALL_smt; ///< ALL smt number

		RESIDTYPE         _resid_type;

		int                  _getgobs(string prn, GOBSTYPE type, GOBSBAND band, t_gobs& gobs);
		map<string, int>     _frqNum;
		map<string, t_gtime> _lastEcl;

		bool                 _ifb3_init;
		bool                 _ifb4_init;
		bool                 _ifb5_init;

		bool                 _auto_band;
		CBIASCHAR            _cbiaschar;

		map<t_gtime, t_gtriple> _map_crd;

		vector<string> _outlier_sat;
	};

} // namespace

#endif // GSPPFLT_H
