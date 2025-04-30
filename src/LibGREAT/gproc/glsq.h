/**
 * @file         glsq.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GLSQ_H
#define GLSQ_H 


#include "gexport/ExportLibGREAT.h"

#include "newmat/newmat.h"
#include "newmat/newmatap.h"

#include "gall/gallpar.h"
#include "gmodels/gstochasticmodel.h"
#include "gmodels/gmodel.h"
#include "gutils/gsysconv.h"
#include "gutils/gcycleslip.h"
#include "gproc/glsqmatrix.h"
#include "gio/giobigf.h"
#include "gall/gallrecover.h"
#include "gproc/gupdatepar.h"
#include "gproc/ginverse.h"
#include <thread>

using namespace std;
namespace great
{

	/**
	*@brief	   Class for lsq estimator include upadte,slove,recover
	*/
	class LibGREAT_LIBRARY_EXPORT t_glsq
	{
	public:

	public:

		/** @brief	constructor */
		t_glsq(t_gsetbase* set = nullptr);

		/** @brief  copy constructor */
		t_glsq(const t_glsq& Other);

		/** @brief  default destructor */
		virtual ~t_glsq();

		/**
		* @brief  add now parameter
		* @param[in] par new par preapred to add
		* return -1 for fail 1 for success
		*/
		virtual int add_parameter(const t_gpar& par);

		/**
		* @brief  add the partype and the state mode
		* @param[in] par_type type of parameter
		* @param[in] order 0:white noise 1:random walk
		* @param[in] dt [unit:hour]
		* @param[in] noise process noise
		*/
		void add_par_state_equ(par_type par_type, int order, double dt, double noise);

		void set_update_par(shared_ptr<t_gupdatepar> updatepar);

		void set_solve_matrix(shared_ptr<t_ginverse<Matrix> > solve_matrix);

		/**
		* @brief update all lsq par with now obs data
		* @note according to obsdata update amb par , and time update outsate par
		* @		use matrix remove way may have little loss of accruacy
		* @param[in] epoch
		* @param[in] obsdata observ info data
		* @param[in] if use matrix remove way
		* @return -1 for fail 1 for success
		*/
		virtual bool update_parameter(const t_gtime& epoch, vector<t_gsatdata>& obsdata, bool matrix_remove, bool write_temp=true);

		bool remove_all_ambiguity(const t_gtime& epoch, par_type ambtype, bool matrix_remove, bool write_temp = true);

		/**
		* @brief set lsq par
		* @note this set will clean the par and other info before
		* @param[in] parameters all parmaters preapred to add
		*/
		virtual int set_parameter(const t_gallpar& parameters);


		/**
		* @brief add new observ equations[newmat matrix format]
		* @param[in] B coeff of equations
		* @param[in] P weight of equations
		* @param[in] l res of equations
		* @param[in] epoch time of equations
		*/
		int add_equation(const Matrix& B, const DiagonalMatrix& P, const ColumnVector& l, const t_gtime& epoch = t_gtime());

		/**
		* @brief add new observ equations[lsqmatrix format]
		* @param[in] equ equations info
		* @param[in] epoch time of equations
		*/
		int add_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch = t_gtime(),bool write_temp=true);

		int get_equ_obs_total_num();

		/** @brief del old observ equations[lsqmatrix format]
		* @param[in] equ equations info
		* @param[in] epoch time of equations
		*/
		
		int del_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch = t_gtime());
		

		/**
		* @brief write the equtions information to tempfile for recovering later[lsqmatrix format
		* @param[in] equ equations of lsqmatrix format
		* @param[in] epoch time of equations
		* @param[in] info information of equations
		* @return true for success false for fail
		*/
		bool  write_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch);

		/**
		* @brief remove the specified par in the lsq estimator
		* @note idx begin from 1
		*/
		virtual int remove_parameter(const int& idx, bool write_temp = true);

		virtual int remove_parameter(vector<int>& idx, bool write_temp = true);

		// sovle the least squares
		/**
		* @brief solve the equtions by least squares
		* @note result save into the dx
		*/
		virtual int solve_NEQ();
		
		virtual void print_Qx();
	
		/**
		* @brief solve the equtions by least squares without inverse the matrix,generally use for lsq_epo
		* @note result save into the dx
		*/
		virtual void solve_x();

		/**
		* @brief recover the removed par by information in tempfile 
		* @note recover should after solving equations and have temp file
		*		recover result will save into the t_gallrecover data struct
		*/
		virtual bool recover_parameters(t_gallrecover& allrecover);


		/** @brief get the final solution of parameters;save result to t_gallrecover */
		void get_result_parameter(t_gallrecover& allrecover);

		/** @brief add apriori for all parameter */
		void add_apriori_weight();

		/**
		* @brief add apriori for specified parameter
		* @note  idx begin from 1
		*/
		void add_apriori_weight(const int& idx);

		/**
		* @brief add specified apriori for specified parameter
		* @note  idx begin from 1
		* @param[in] idx specified location
		* @param[in] value given intial value
		* @param[in] weight given apriori weight
		*/
		void add_apriori_weight(const int& idx, const double& value, const double& weight);

		/** @brief reset filename of tempfile */
		void reset_tempfile(string filename);


		/** @brief add ISB/IFB constraint */
		int lsq_sysbias_constraint();

		/** @brief add CLK zero-mean constraint */
		int lsq_clk_constraint();

		/** @brief print matrix of NEQ(BTPB) and W(BTPL) */
		void print_matrx();

		/** @brief whether the parameter is alive in current epoch */
		bool par_alive(const int& idx);

		/** @brief Get NEQ Matrix*/
		SymmetricMatrix	NEQ() const;
		ColumnVector W();
		V_SymmetricMatrix v_NEQ()const;
		V_ColumnVector v_W()const;


		/** @brief get covariance matrix */
		const SymmetricMatrix& Qx() const;

		/** @brief get coerrection of parametes */
		ColumnVector 	dx() const;


		/** @brief get covariance of coerrection of parametes */
		ColumnVector	stdx() const;

		double Qx(const int& col, const int& row) const;

		/** @brief get coerrection of specified parametes */
		double dx(int idx) const;

		/** @brief get the covariance of specified parameter */
		double stdx(int idx) const;

		/** @brief get sgima0 */
		double sigma0() const;

		/** @brief get sum of omc */
		double res_obs() const { return _res_obs; }

		/** @brief get sum of res squares */
		double vtpv() const;

		/** @brief get total observation numbers*/
		int nobs_total()  const;

		/** @brief get nubmer of parameter*/
		int npar_number() const;


		/** @brief get current time  */
		t_gtime          epo() const;

		/** @brief set current time */
		void             epo(const t_gtime& t);

		/** @brief get lsq mode */
		LSQMODE mode();

		/** @brief set log file */
		void setlog(t_glog* log);

		/** @brief set NEQ Matrix specified location */
		void change_NEQ(int row, int col, double xx);

		/** @brief set Qx Matrix specified location */
		void change_Qx(int row, int col, double xx);
		void change_Qx(SymmetricMatrix Qx);

		/** @brief set dx  specified location */
		void change_dx(int n, double xx);
		void change_dx(ColumnVector dx);

		/** @brief set stdx specified location */
		void change_stdx(int n, double xx);
		void change_stdx(ColumnVector stdx);

		/** @brief set sigma0 by  specified sigma */
		void change_sigma0(double sigma);

		/** @brief set vtpv by  specified value */
		void change_vtpv(double xx);
		/** @brief change NEW and W */
		void set_new_NEQ(const V_ColumnVector& W, const V_SymmetricMatrix& NEQ);

		///** @brief get one site max amb in one epoch */
		int reset_npar_total(const t_gpar& par);

		t_gallpar  _x_solve;	///< all parameter in lsq estimator
		int _obs_total_num_epo = 0;       ///< totoal observation number of current epoch

	protected:
		/**
		* @brief  write the removed NEQ W to tempfile
		* @note idx from1
		*/
		int _write_coefficient(int idx);
		// remove_id from 1;
		int _write_parchage(const vector<int>& remove_id);


		/** @brief recover parameter */
		void _recover_par(t_greadtemp& tempfile_in, t_gallrecover& _allrecover);
		void _recover_par_part(t_greadtemp& tempfile_in, t_gallrecover& _allrecover);
		void _recover_par_swap(t_greadtemp& tempfile_in);

		/** @brief recover obs */
		void _recover_obs(t_greadtemp& tempfile_in, t_gallrecover& _allrecover);


		/** @brief remove zero element in NEQ matrix */
		void _remove_zero_element(SymmetricMatrix& B, ColumnVector& l, int idx);

		/** @brief solve equation*/
		void _solve_equation(const SymmetricMatrix& NEQ, const ColumnVector& W, ColumnVector& ans, ColumnVector& Q);
		void _solve_x(const SymmetricMatrix& NEQ, const ColumnVector& W, ColumnVector& ans);


		t_glog*	_log;						///< log file

		V_ColumnVector    _W;				///< BTPL Matrix
		V_SymmetricMatrix _NEQ;			    ///< BTPB Matrix

		SymmetricMatrix _Qx;					///< storage Qx after solve
		ColumnVector	_dx;					///< correction of all parameter
		ColumnVector    _dx_final;				///< last correction of all parameter
		ColumnVector	_stdx;					///< covarience of correction of all parameter

		double _res_obs;						///< sum of residuals of currecnt epoch
		double _vtpv;							///< sum of residuals of all epoch
		double _sigma0;							///< error in unit weight
		long long _obs_total_num;				///< totoal observation number
		int _npar_tot_num;						///< total par number include ll eliminated parameters 


		t_gtime _epo;						///< current time
		t_gtime	_beg;						///< begin time
		t_gtime	_end;						///< end time
		double  _interval;                     ///< intv
		t_giobigf*	_tempfile;					///< temp file
		LSQMODE  _mode;						///< lsq model
		int _buffer_size=1024*1000*10; 		///< tempfile buffer size

		
		shared_ptr<t_gupdatepar> _update_lsqpar;
		shared_ptr<t_ginverse<Matrix> > _solve_matrix;

		t_gmutex _lsq_mtx;

		int size_int = sizeof(int);
		int size_dbl = sizeof(double);

	};

	// for recover or exctract par
	LibGREAT_LIBRARY_EXPORT t_gtime str2gtime(string str_time);


} // namespace

#endif //GLSQ_H

