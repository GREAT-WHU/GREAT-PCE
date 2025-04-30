/**
*
* @verbatim
    History
    2011-09-04  JD: created
    2012-09-24  JD: selections of code, phase etc.
    2020-10-11  Jiande Huang : Optimizing the efficiency of "select" function

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file      gobsgnss.h
* @brief     Purpose: implementation of GNSS observation element
*            Todo: multi-GNSS & various LCs !
*            gobj implementation
*.
* @author     JD
* @version    1.0.0
* @date       2011-09-04
*
*/

#ifndef GOBSGNSS_H
#define GOBSGNSS_H


#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>

#include "gdata/gdata.h"
#include "gutils/gnss.h"
#include "gutils/gobs.h"
#include "gutils/gtime.h"
#include "gutils/gsys.h"

#define DEF_CHANNEL 255

using namespace std;

namespace gnut
{

	const static double NULL_GOBS = 0.0;


	// priority tables for choice of available signals (code [m])
	const static GOBS code_choise[9][19] = {
	   {   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X }, //
	   { C1A, C1B, C1C,   X, C1I, C1L, C1M,   X, C1P, C1Q, C1S, C1W, C1X, C1Y, C1Z,  P1,  C1,  CA,  CB }, //  C1
	   {   X,   X, C2C, C2D, C2I, C2L, C2M,   X, C2P, C2Q, C2S, C2W, C2X, C2Y,   X,  P2,  C2,  CC,  CD }, //  C2
	   {   X,   X,   X,   X, C3I,   X,   X,   X,   X, C3Q,   X,   X, C3X,   X,   X,   X,   X,   X,   X }, //
	   {   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X }, //
	   {   X,   X,   X,   X, C5I,   X,   X,   X,   X, C5Q,   X,   X, C5X,   X,   X,  P5,  C5,   X,   X }, //  C5
	   { C6A, C6B, C6C,   X, C6I,   X,   X,   X,   X, C6Q,   X,   X, C6X,   X, C6Z,   X,  C6,   X,   X }, //  C6
	   {   X,   X,   X,   X, C7I,   X,   X,   X,   X, C7Q,   X,   X, C7X,   X,   X,   X,  C7,   X,   X }, //  C7
	   {   X,   X,   X,   X, C8I,   X,   X,   X,   X, C8Q,   X,   X, C8X,   X,   X,   X,  C8,   X,   X }  //  C8
	};

	// priority tables for choice of available signals (phase [full-cycles])
	const static GOBS phase_choise[9][19] = {
	   {   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X }, //
	   { L1A, L1B, L1C,   X, L1I, L1L, L1M, L1N, L1P, L1Q, L1S, L1W, L1X, L1Y, L1Z,   X,  L1,  LA,  LB }, //  L1
	   {   X,   X, L2C, L2D, L2I, L2L, L2M, L2N, L2P, L2Q, L2S, L2W, L2X, L2Y,   X,   X,  L2,  LC,  LD }, //  L2
	   {   X,   X,   X,   X, L3I,   X,   X,   X,   X, L3Q,   X,   X, L3X,   X,   X,   X,   X,   X,   X }, //
	   {   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X }, //
	   {   X,   X,   X,   X, L5I,   X,   X,   X,   X, L5Q,   X,   X, L5X,   X,   X,   X,  L5,   X,   X }, //  L5
	   { L6A, L6B, L6C,   X, L6I,   X,   X,   X,   X, L6Q,   X,   X, L6X,   X, L6Z,   X,  L6,   X,   X }, //  L6
	   {   X,   X,   X,   X, L7I,   X,   X,   X,   X, L7Q,   X,   X, L7X,   X,   X,   X,  L7,   X,   X }, //  L7
	   {   X,   X,   X,   X, L8I,   X,   X,   X,   X, L8Q,   X,   X, L8X,   X,   X,   X,  L8,   X,   X }  //  L8
	};

	const static map<GSYS, map<GOBSBAND, string> >	range_order_attr_raw = {
		{ GPS, { {BAND_1,"CPW" },{BAND_2,"CLXPW"}, {BAND_5,"QX"} } },
		{ GAL, { {BAND_1,"CX"  },{BAND_5,"IQX"  }, {BAND_7,"IQX"},{BAND_8,"IQX"},{BAND_6,"ABCXZ"}} },
		{ BDS, { {BAND_2,"IQX" },{BAND_7,"IQX"  }, {BAND_6,"IQX"},{BAND_5,"DPX"},{BAND_9,"DPZ"},{BAND_8,"DPX"},{BAND_1,"DPX"}} },
		{ GLO, { {BAND_1,"CP"  },{BAND_2,"CP"   }                       } },
		{ QZS, { {BAND_1,"CSLX"},{BAND_2,"LX"   },{BAND_5,"IQX"}} }
	};

	const static map<GSYS, map<GOBSBAND, string> >	range_order_attr_cmb = {
		{ GPS, { {BAND_1,"CPW"},{BAND_2,"CLXPW"},{BAND_5,"QX"} }},
		{ GAL,{{BAND_1,"CX"},{BAND_5,"IQX"}, {BAND_7,"IQX"},{BAND_8,"IQX"},{BAND_6,"ABCXZ"}} },
		{ BDS,{{BAND_2,"IQX"},{BAND_7,"IQX"},{BAND_6,"IQX"},{BAND_5,"DPX"},{BAND_9,"DPZ"},{BAND_8,"DPX"},{BAND_1,"DPX"}} },
		{ GLO,{{BAND_1,"CP"},{BAND_2,"CP"}}},
		{ QZS,{{BAND_1,"CSLX"},{BAND_2,"LX"},{BAND_5,"IQX"}} }
	};

	const static map<GSYS, map<GOBSBAND, string> >	phase_order_attr_raw = {
		{GPS,{{BAND_1,"CSLXPWYM"},{BAND_2,"CDLXPWYM"},{BAND_5,"IQX"}}},
		{GAL,{{BAND_1,"ABCXZ"},{BAND_5,"IQX"},{BAND_7,"IQX"},{BAND_8,"IQX"},{BAND_6,"ABCXZ"}}},
		{BDS,{{BAND_2,"XIQ"},{BAND_7,"IQX"},{BAND_6,"IQX"},{BAND_5,"DPX"},{BAND_9,"DPZ"},{BAND_8,"DPX"},{BAND_1,"DPX"}}},
		{GLO,{{BAND_1,"PC"},{BAND_2,"CP"}}},
		{QZS,{{BAND_1,"CSLX"},{BAND_2,"LX"},{BAND_5,"IQX"}}}
	};

	const static map<GSYS, map<GOBSBAND, string> >	phase_order_attr_cmb =
	{
		{GPS,{{BAND_1,"CSLXPWYM"},{BAND_2,"CDLXPWYM"},{BAND_5,"IQX"}}},
		{GAL,{{BAND_1,"ABCXZ"},{BAND_5,"IQX"},{BAND_7,"IQX"},{BAND_8,"IQX"},{BAND_6,"ABCXZ"}}},
		{BDS,{{BAND_2,"QXI"},{BAND_7,"IQX"},{BAND_6,"IQX"},{BAND_5,"DPX"},{BAND_9,"DPZ"},{BAND_8,"DPX"},{BAND_1,"DPX"}}},
		{GLO,{{BAND_1,"PC"},{BAND_2,"CP"}}},
		{QZS,{{BAND_1,"CSLX"},{BAND_2,"LX"},{BAND_5,"IQX"}}}

	};

	class t_obscmb;
	class t_gallbias;
	/** @brief class for t_gobsgnss derive from t_gdata. */
	class LibGnut_LIBRARY_EXPORT t_gobsgnss : public t_gdata {

	public:
		/** @brief default constructor. */
		t_gobsgnss();
		t_gobsgnss(const string& sat);
		t_gobsgnss(const string& site, const string& sat, const t_gtime& t);
		t_gobsgnss(const string& site, const string& sat, const t_gtime& t, const bool is_leo);
		virtual ~t_gobsgnss();
		void addobs(const GOBS& obs, const double& d);                                   // add a single observation
		void addlli(const GOBS& obs, const int& i);                                      // add a lost-of-lock indicator
		void addslip(const GOBS& obs, const int& i);                                     // add an estimated cycle slip   
		void addele(double d);                                                           // add approximate elevation
		void addoutliers(const GOBS& obs, const int& i);
		void setrangestate(string name, bool b);
		bool getrangestate(string name) const;

		size_t size() const;                                                             // get number of observations
		double frequency(GOBSBAND b) const;                                              // get system-specific frequency  for band i
		double frequency(int b) const;                                                   // get system-specific frequency  for band i
		double wavelength(GOBSBAND b) const;                                             // get system-specific wavelength for band i
		double wavelength_L3(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2);                // get wavelength for iono-free LC
		double wavelength_WL(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2);                // get wavelength for wild lane
		double wavelength_NL(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2);                // get wavelength for narrow lane
		GSYS gsys() const;                                                               // get GNSS system from satellite IDOC

		GOBS id_range(GOBSBAND b) { return _id_range(b); };                              //gppflt NEW!
		GOBS id_phase(GOBSBAND b) { return _id_phase(b); };                              //gppflt NEW!

		GOBS select_range(GOBSBAND b, const bool& isRawAll = true)const;				 // get code obs of selected band
		GOBS select_phase(GOBSBAND b, const bool& isRawAll = true)const;				 // get pahse obs of selected band

		vector<GOBS> obs();                                                              // get vector of available observations
		set<GOBS>    obs_phase(const int& band);                                         // get set of available phase observations for a band
														                                 
		double getobs(const GOBS& obs) const;                                                  // get a single observation (DIFFERENT UNITS!)
		double getobs(const string& obs);                                                // get a single observation (DIFFERENT UNITS!)
		double getele();                                                                 // get approximate elevation
		int    getlli(const GOBS& obs) const;                                                  // get a lost-of-lock indicator
		int    getslip(const GOBS& obs);                                                 // get an estimated cycle slip   
		int    getlli(const string& obs) const;                                                // get a lost-of-lock indicator
		void   nbands(pair<int, int>& nb);                                               // get number of code/phase available bands
		int    getoutliers(const GOBS& obs);
		void   resetobs(const GOBS& obs, double v);
		bool   apply_bias(t_gallbias& allbias);
		bool   apply_code_phase_bias(t_gallbias& allbias);
		bool   apply_bias_tmp(t_gallbias& allbias);
		bool   apply_bias_tmp2(t_gallbias& allbias);
		bool   apply_dcb(t_gallbias& allbias);
		bool   apply_osb(t_gallbias& allbias);
		double interFreqDcb(t_gallbias& allbias, GOBSBAND band1, GOBSBAND band2);


		void   channel(int canal);                                                       // set channel number for Glonass satellites
		int    channel() const;                                                          // get channel number for Glonass satellites   
															                             
		void rtcm_end(unsigned int end) { _rtcm_end = end; }                             // set RTCM Multiple Message bit
		unsigned int rtcm_end() { return _rtcm_end; }                                    // get RTCM Multiple Message bit
															                             
		double obs_range(const t_gband& gb);                                             // get  code observation [m] only requested type NEW + AUTO!
		double obs_phase(const t_gband& gb);                                             // get phase observation [m] only requested type NEW + AUTO !   

		double obs_C(const t_gobs& gobs);  // NEW                                        // get code observation [m] only requested type!
		double obs_L(const t_gobs& gobs);  // NEW                                        // get phase observation [m] only requested type
		double obs_D(const t_gobs& gobs);  // NEW                                        // get doppler observation [m/s] only requested type!
		double obs_S(const t_gobs& gobs);  // NEW                                        // get snr  observation [dbHz] only requested type! by wh
		double obs_C(const t_gband& gb);   // NEW                                        // get code observation [m] only requested type
		double obs_L(const t_gband& gb);   // NEW                                        // get phase observation [m] only requested type
		double obs_D(const t_gband& gb);   // NEW                                        // get doppler observation [m/s] only requested type
		double obs_S(const t_gband& gb);   // NEW                                        // get snr observation [dbHz] only requested type by wh

		int mod_L(const double& dL, const GOBS& gobs = X, const int i = 1);  // modify phase observation by adding dL: i = 0 [cycles], i = 1 [m]
		int mod_L(const double& dL, const GOBSBAND& band, const int i = 1);  // modify phase observation by adding dL: i = 0 [cycles], i = 1 [m], i = 2 [set val to 0]

		double frequency_lc(const    int& band1,  // return value [Hz] of phase linear combination frequency (c1*O1 + ... )
			const double& coef1,  // for 2 or 3 bands with given coefficients
			const    int& band2,
			const double& coef2,
			const    int& band3 = 0,
			const double& coef3 = 0); // -->> OLD INTERFACE !!

		double isf_lc(const int& band1,  // return value of phase linear combination ionospheric scale factor (c1*O1 + ... )
			const double& coef1,  // for 2 or 3 bands with given coefficients
			const    int& band2,
			const double& coef2);
		//     const    int& band3 = 0, 
		//     const double& coef3 = 0); // -->> OLD INTERFACE !!

		double pnf_lc(const int& band1,  // return value of linear combination phase noise factor
			const double& coef1,  // for 2 or 3 bands with given coefficients
			const    int& band2,
			const double& coef2,
			const    int& band3 = 0,
			const double& coef3 = 0); // -->> OLD INTERFACE !!
// POTENCIALNE NOVY (uvidime jak bude uzitecne, ale kdyz uz mame nastaveno t_gobs, melo by byt)
// --------------------------------------------------------------------------------------------
		double P3(const t_gobs& g1, const t_gobs& g2);  // get ionosphere-free combination for code  [m]
		double P4(const t_gobs& g1, const t_gobs& g2);  // get geometry-free combination for code    [m]

		// POTENCIALNE NOVY (uvidime jak bude uzitecne, ale kdyz uz mame nastaveno t_gobs, melo by byt)
		// --------------------------------------------------------------------------------------------
		double L3(const t_gobs& g1, const t_gobs& g2);  // get ionosphere-free combination for phase [m]
		double L3_cycle(const t_gobs& g1, const t_gobs& g2);  // get ionosphere-free combination for phase [CYCLE]
		double L4(const t_gobs& g1, const t_gobs& g2);  // get geometry-free combination for phase   [m]
		double L4_cycle(const t_gobs& g1, const t_gobs& g2);  // get geometry-free combination for phase   [cycle]
		double MP(const t_gobs& code,
			const t_gobs& L1, const t_gobs& L2);  // get multipath for C-obs + Li/Lj
		double GFIF_meter(const t_gobs& gL1, const t_gobs& gL2, const t_gobs& gL3);

		double PWL(const t_gobs& g1, const t_gobs& g2);
		double LWL(const t_gobs& g1, const t_gobs& g2); // get wide-lane combination for phase [m]!
		double LNL(const t_gobs& g1, const t_gobs& g2); // get narrow-lane combination for phase [m]!      
		double MW(const t_gobs& g1, const t_gobs& g2,  // get Melbourne-Wuebenna combination for phase & code [m]!
			bool phacode_consistent = false);    // --> select phase observation + (optionally) use for code as well
		double MW_cycle(const t_gobs& gL1, const t_gobs& gL2,
			const t_gobs& gC1, const t_gobs& gC2);  // get Melbourne-Wuebenna combination for phase & code [cycle]!);
		double EWL_cycle(const t_gobs& gL1, const t_gobs& gL2, const t_gobs& gL3,
			const t_gobs& gC1, const t_gobs& gC2);  // get extra-wide-lane Melbourne-Wuebenna combination for phase & code [cycle]!);
		double LW_meter(const t_gobs& gL1, const t_gobs& gL2);  // get LW combination for phase [meter]!);
		double LE_meter(const t_gobs& gL2, const t_gobs& gL3);  // get LE combination for phase [meter]!);
		double LWL_factor13(const t_gobs& gL1, const t_gobs& gL3);  // get LWL combination factor);
		// NOVE FUNKCE!
		// --------------------------------------------------------------------------------------------
		t_gtime epoch() const { return _epoch; }                      // get reference epoch
		string  site()  const { return _staid; }                      // get station id
		string  sat()   const { return _satid; }                      // get satellite id
		string  sys()   const { return _satid.substr(0, 1); }          // get satellite system id
		map<GOBS, int> lli() const { return _glli; }                 // get satellite lost-of-lock identifications
		map<GOBS, int> slip() const { return _gslip; }               // get maps of estimated cycle slips

		void   sat(string id) {
			_satid = id;
			_gsys = t_gsys::char2gsys(id[0]);
		} // set satellite id         
		void  site(string id) { _staid = id; }                       // set site id
		void   epo(t_gtime t) { _epoch = t; }                       // set epoch   

		t_gobsgnss operator-(t_gobsgnss& obs);
		void clear();
		bool valid();
		bool obs_empty();

		set<GFRQ>     freq_avail();                             // get available freq  ---> NOVA IMPLEMENTACE !!!!!
		set<GOBSBAND> band_avail(bool _phase = true);           // get available band  ---> NOVA IMPLEMENTACE !!!!!
		bool          contain_freq(FREQ_SEQ freq);

		double P3(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2); // get ionosphere-free combination for code  [m]
		double P4(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2); // get geometry-free combination for code    [m]

		double L3(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2); // get ionosphere-free combination for phase [m]
		double L4(GOBSBAND b1 = BAND_1, GOBSBAND b2 = BAND_2); // get geometry-free combination for phase   [m]   

		int coef_ionofree(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the ionosphere-free linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int coef_ionofree_phi(GOBSBAND b1, double& c1,
			GOBSBAND b2, double& c2);				// return coefficients for carrier phase observtion

		int coef_geomfree(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the geometry-free linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int coef_narrlane(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the narrow-lane linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int coef_widelane(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the wide-lane linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		bool health() const { return _health; }                // get sat health
		void health(bool health) { _health = health; } // set sat health
		void dcb_mark(bool h) { _dcb_correct_mark = h; }

	protected:
		double _obs_range(const t_gobs& go);        // get  code observation [m] only requested type NEW + AUTO!
		double _obs_range(const t_gband& gb);        // get  code observation [m] only requested type NEW + AUTO!
		double _obs_phase(const t_gband& gb);        // get phase observation [m] only requested type NEW + AUTO !
		double _obs_doppler(const t_gband& gb);      // get doppler observation [m] only requested type NEW + AUTO !
		double _obs_snr(const t_gband& gb);      // get snr observation [dbhz] only requested type NEW + AUTO !


		double _lcf_range(const t_gobs* g1, const int& cf1,        // return value [m] of general pseudo-range
			const t_gobs* g2, const int& cf2,        // LC = c1*f1*O1 + ... / c1*f1 +
			const t_gobs* g3 = 0, const int& cf3 = 0);     // for 2 or 3 bands with given coefficients  NEW !

		double _lcf_phase(const t_gobs* g1, const int& cf1,        // return value [m] of general carrier-phase
			const t_gobs* g2, const int& cf2,        // LC = c1*f1*O1 + ... / c1*f1 +           
			const t_gobs* g3 = 0, const int& cf3 = 0);     // for 2 or 3 bands with given coefficients   NEW !

		double _lc_range(const t_gobs* g1, const double& cf1,      // return value [m] of general pseudo-range
			const t_gobs* g2, const double& cf2,      // LC = c1*O1 + ... 
			const t_gobs* g3 = 0, const double& cf3 = 0);   // for 2 or 3 bands with given coefficients  NEW !   

		double _lc_phase(const t_gobs* g1, const double& cf1,      // return value [m] of general carrier-phase
			const t_gobs* g2, const double& cf2,   // LC = c1*O1 + ... 
			const t_gobs* g3 = 0, const double& cf3 = 0);   // for 2 or 3 bands with given coefficients   NEW !

		int _coef_multpath(GOBSBAND bC, double& cC,     // get c1,c2,c3 coefficients of code multipath LC
			GOBSBAND b1, double& c1,     // (2 bands, 1x code & 2x phase)
			GOBSBAND b2, double& c2);

		int _coef_ionofree(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the ionosphere-free linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int _coef_geomfree(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the geometry-free linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int _coef_narrlane(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the narrow-lane linear combination
			GOBSBAND b2, double& c2);  // (2 bands)

		int _coef_widelane(GOBSBAND b1, double& c1,   // return coefficients (c1,c2) of the wide-lane linear combination
			GOBSBAND b2, double& c2);  // (2 bands)


		GOBS _id_range(GOBSBAND b);                   // get  code ID of selected band (according to table)
		GOBS _id_phase(GOBSBAND b);                   // get phase ID of selected band (according to table)
		GOBS _id_doppler(GOBSBAND b);                 // get doppler ID of selected band (according to table)
		GOBS _id_snr(GOBSBAND b);                     // get snr ID of selected band (according to table)


		GOBS _cod_id(const int& band);                // get  code ID of selected band (according to table)
		GOBS _pha_id(const int& band);                // get phase ID of selected band (according to table)

		virtual void _clear();
		virtual bool _valid() const;
		virtual bool _valid_obs() const;

		set<GFRQ>     _freq_avail();                  // get available freq  for phase ---> NOVA IMPLEMENTACE !!!!!
		set<GOBSBAND> _band_avail();                  // get available bands for phase ---> NOVA IMPLEMENTACE !!!!!
		set<GOBSBAND> _band_avail_code();             // get available bands for code  ---> NOVA IMPLEMENTACE !!!!!   

		map<GOBS, double>   _gobs;   // maps of observations
		map<GOBS, int>      _glli;   // maps of lost-of-lock identifications
		map<GOBS, int>      _gslip;  // maps of estimated cycle slips
		map<GOBS, int>      _goutlier;  // maps of gross error of range/doppler observations

		string              _staid;  // station id
		string              _satid;  // satellite id ["G??", "R??", "E??" ...]
		GSYS                _gsys;   // system 
		t_gtime             _epoch;  // epoch of the observation
		double              _apr_ele;// approximate elevation
		int                 _channel;// satellite channel number
		unsigned int        _rtcm_end; // RTCM Multiple Message bit (0 = end, 1 = cont.)

		bool                _health;
	private:
		bool                _isDoubleIf = false;
		bool                _dcb_correct_mark;      // mark for correct dcb; true: dcb corrected; false: not corrected; 
		bool                _range_smooth_mark;     // true, smooth; false, not smooth
		bool                _bds_code_bias_mark;    // true, correct; false, not correct
	};

	class t_obscmb
	{
	public:
		t_obscmb() { num = 0.0; lam = 0.0; };
		double num;
		double lam;
		t_gobs first;
		t_gobs second;
		bool operator<(const t_obscmb& t) const;
	};

} // namespace

#endif
