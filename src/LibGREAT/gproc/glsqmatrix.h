/**
 * @file         glsqmatrix.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef  GLSQMATRIX_H
#define  GLSQMATRIX_H

#include "gexport/ExportLibGREAT.h"
#include <gset/gsetproc.h>
#include "Eigen/Core"
#include "newmat/newmat.h"
#include "gutils/gtypeconv.h"
#include "gutils/gturboedit.h"
#include "gmodels/gbasemodel.h"

using namespace std;
using namespace gnut;
namespace great
{

	/**
	* @brief  Matrix for storage observe equations of B,P,L
	*/
	class LibGREAT_LIBRARY_EXPORT t_glsqEquationMatrix :public t_gbaseEquation
	{
	public:

		/** @brief default constructor */
		t_glsqEquationMatrix();

		/** @brief default destructor  */
		~t_glsqEquationMatrix();

		/**
		* @brief add equations
		* @param[in] B_value coeff of observ equation
		* @param[in] P_value weight of observ equation
		* @param[in] l_value res of observ equaion
		*/
		void add_equ(const vector<pair<int, double> >& B_value, const double& P_value, const double& l_value, const string& stie_name, const string& sat_name, const t_gobscombtype& obscombtype, const bool& is_newamb);
		
		//void add_equ(const vector<pair<int, double> >& B_value, const vector<pair<int, double> >& P_raw, const double& P_value, const double& l_value, const string& stie_name, const string& sat_name, const t_gobscombtype& obscombtype, const bool& is_newamb);

		void add_equ(const t_glsqEquationMatrix& Other);
		/**
		* @brief add equations by newmat format
		* @param[in] B_value coeff of observ equation newmat format
		* @param[in] P_value weight of observ equation newmat format
		* @param[in] l_value res of observ equaion newmat format
		*/
		void add_equ(const Matrix& B_value, const DiagonalMatrix& P_value, const ColumnVector& l_value, const vector<string>& site_name, const vector<string>&  sat_name, const vector<t_gobscombtype>& obstype);


		void set_newamb(const int& idx, const bool& is_newamb);
		void set_newamb(t_gturboedit* tb_slip);
		void set_newamb(const FREQ_SEQ& freq_1, const FREQ_SEQ& freq_2, t_gturboedit* tb_slip);
		/**
		* @brief remove specified equations
		* @note idx from 1
		* param[in] idx specified loc
		*/
		void remove(const int& idx);


		void remove_last_equation();

		/**
		* @brief change equations to newmat format
		* @param[out] B_value coeff of observ equation newmat format
		* @param[out] P_value weight of observ equation newmat format
		* @param[out] l_value res of observ equaion newmat format
		*/
		void chageNewMat(Matrix& B_value, DiagonalMatrix& P_value, ColumnVector& l_value, const int& par_num);
		void chageNewMat(Matrix& B_value, SymmetricMatrix& P_value, ColumnVector& l_value, const int& par_num);
	
		/** @brief print B P L Matrix*/
		void print();

		/**
		* @brief get numbers of equations
		* @return size of equations
		*/
		int num_equ()const;

		/**
		* @brief get resiuals of equations
		* @return LTPL of equations
		*/
		double res_equ() const;

		/**
		* @brief get resiuals of equations
		* @param[in] phase phase only or both
		* @return LTPL of equations
		*/
		double res_equ(bool phase) const;

		string get_sitename(int equ_idx) const;
		set<string> get_satlist(string rec) const;

		string get_satname(int equ_idx) const;
		vector<pair<string, string> > get_site_sat_pair() const;
		string get_obscombtype2str(int equ_idx) const;

		t_gobscombtype get_obscombtype(int equ_idx) const;

		vector<double> get_codeomc() const;
		vector<double> get_phaseomc() const;

		bool is_newamb(int equ_idx) const;

		friend LibGREAT_LIBRARY_EXPORT void  print_equ_debinfo(const t_glsqEquationMatrix& equ);


		vector<int> find_equ(const string& site);
		vector<int> find_equ(const string& site, const string & sat);
		int find_equ(const string& site, const string& sat, const t_gobscombtype& obscomtype);

		void clear_allequ();
		void  swap_allequ();

	protected:

		vector<pair<string, string> >  _site_sat_pairlist;
		vector<t_gobscombtype> _obstypelist;
		vector<bool> _newamb_list;
	};


	/**
	* @brief  Matrix for storage Matrix NEQ
	*/
	class LibGREAT_LIBRARY_EXPORT t_glsqSymmetricMatrix
	{
	public:


		/**
		* @brief default destructor
		* @return size of NEQ matrix
		*/
		virtual int num() const = 0;

		virtual double& num(int a, int b) = 0;
		//double  num(int a, int b) const ;

		/**
		* @brief resize the NEQ matrix according to specified size
		* @note clean the data before
		* @param[in] size specified size
		*/
		virtual void resize(int size) = 0;

		/**
		* @brief add observ equations
		* @parma[in] equ observ equations
		*/
		virtual void add(const t_glsqEquationMatrix& equ) = 0;

		/** @brief add one zero dimension at last*/
		virtual void addBackZero() = 0;

		/**
		* @brief removed specified dimension in the NEQ matrix
		* @note idx from 1
		* @param[in] idx specified dimension
		*/
		virtual void remove(int idx) = 0;

		/**
		* @brief print NEQ matrix in std out
		*/
		virtual void print() = 0;

		/**
		* @brief get diagonal value in NEQ matrix by specified loc
		* @note idx from 1
		* @param[in] idx specified loc
		*/
		virtual double center_value(int idx) const = 0;

		/**
		* @brief change NEQ matrix to NewMat format
		* @return NEQ matrix in newmat format
		*/
		virtual SymmetricMatrix changeNewMat() const = 0;

	};


	/**
	* @brief class for store BTPL(W) matrix
	*/
	class LibGREAT_LIBRARY_EXPORT t_glsqColumnVector
	{
	public:

		/**
		* @brief get size of column vector
		* @return size
		*/
		virtual int num() const = 0;

		virtual double& num(int a) = 0;
		//double num(int a) const;

		/**
		* @brief add observ equation
		* @param[in] equ observ equation
		*/
		virtual void add(const t_glsqEquationMatrix& equ) = 0;

		/**
		* @brief add zero dimension in the last loc
		*/
		virtual void addBackZero() = 0;

		/**
		* @brief resize the column vector
		* @note resize will clean data before
		* @param[in] size specified size of column vector
		*/
		virtual void resize(int size) = 0;

		/**
		* @brief remove specified dimension of W matrix
		* @param[in] idx specified loc
		*/
		virtual void remove(int idx) = 0;

		/** @brief print W matrix in stdout */
		virtual void print() = 0;

		virtual ColumnVector changeNewMat() = 0;

	};

	class LibGREAT_LIBRARY_EXPORT L_SymmetricMatrix :public t_glsqSymmetricMatrix
	{
	public:

		/** @brief default constructor */
		L_SymmetricMatrix();

		/** @brief default destructor */
		~L_SymmetricMatrix();

		/**
		* @brief default destructor
		* @return size of NEQ matrix
		*/
		int num() const;

		double& num(int a, int b);
		//double  num(int a, int b) const ;

		/**
		* @brief resize the NEQ matrix according to specified size
		* @note clean the data before
		* @param[in] size specified size
		*/
		void resize(int size);

		/**
		* @brief add observ equations
		* @parma[in] equ observ equations
		*/
		void add(const t_glsqEquationMatrix& equ);

		/** @brief add one zero dimension at last*/
		void addBackZero();

		/**
		* @brief removed specified dimension in the NEQ matrix
		* @note idx from 1
		* @param[in] idx specified dimension
		*/
		void remove(int idx);

		/**
		* @brief print NEQ matrix in std out
		*/
		void print();

		/**
		* @brief get diagonal value in NEQ matrix by specified loc
		* @note idx from 1
		* @param[in] idx specified loc
		*/
		double center_value(int idx) const;

		/**
		* @brief change NEQ matrix to NewMat format
		* @return NEQ matrix in newmat format
		*/
		SymmetricMatrix changeNewMat() const;

	private:

		list<list<double> > _element;		///< element in NEQ Matrix,store element in lower triangle
		list<list<double> >::iterator row_it;
		list<double>::iterator col_it;
		int row_record, col_record;
	};

	class LibGREAT_LIBRARY_EXPORT L_ColumnVector : public t_glsqColumnVector
	{
	public:
		/** @brief default constructor */
		L_ColumnVector();
		/** @brief default destructor */
		~L_ColumnVector();

		/**
		* @brief get size of column vector
		* @return size
		*/
		int num() const;

		double& num(int a);
		//double num(int a) const;

		/**
		* @brief add observ equation
		* @param[in] equ observ equation
		*/
		void add(const t_glsqEquationMatrix& equ);

		/**
		* @brief add zero dimension in the last loc
		*/
		void addBackZero();

		/**
		* @brief resize the column vector
		* @note resize will clean data before
		* @param[in] size specified size of column vector
		*/
		void resize(int size);

		/**
		* @brief remove specified dimension of W matrix
		* @param[in] idx specified loc
		*/
		void remove(int idx);

		/** @brief print W matrix in stdout */
		void print();

		ColumnVector changeNewMat();

	private:

		list<double> _element;		///< element in W matrix
		list<double>::iterator row_it;
		int row_record;
	};

	class LibGREAT_LIBRARY_EXPORT V_SymmetricMatrix;

	class LibGREAT_LIBRARY_EXPORT V_ColumnVector :public t_glsqColumnVector
	{
	public:

		int num() const;
		double& num(int a);
		void resize(int num);
		void add(const t_glsqEquationMatrix& equ);
		void del(const t_glsqEquationMatrix& equ);
		void add(const t_glsqEquationMatrix& equ, bool phase);
		void addBackZero();
		void remove(int idx);
		void print();
		ColumnVector changeNewMat();
		void swap(int a, int b);

		friend bool LibGREAT_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);
		friend bool LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero);
		friend int LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar);
		friend class t_glsq;
	private:
		vector<double> _element;
	};


	class LibGREAT_LIBRARY_EXPORT V_SymmetricMatrix :public t_glsqSymmetricMatrix
	{
	public:

		int num() const;
		inline double& num(int a, int b) { return (a < b) ? _element[b - 1][a - 1] : _element[a - 1][b - 1]; }
		void resize(int num);
		void add(const t_glsqEquationMatrix& equ);
		void add_related(const t_glsqEquationMatrix& equ);
		void del(const t_glsqEquationMatrix& equ);
		void add(const t_glsqEquationMatrix& equ, bool phase);
		void addBackZero();
		void remove(int idx);
		void print();
		double center_value(int idx) const;
		SymmetricMatrix changeNewMat() const;
		void swap(int a, int b);

		friend bool LibGREAT_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);
		friend bool LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero);
		friend int LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W, t_gallpar& allpar);
		friend class t_glsq;

	private:
		vector< vector<double> > _element;
	};




	bool LibGREAT_LIBRARY_EXPORT remove_lsqmatrix(int idx, V_SymmetricMatrix& NEQ, V_ColumnVector& W);


	// by matrix part
	bool LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(const vector<int>& remove_idx, t_gallpar& allpar, V_SymmetricMatrix& NEQ, V_ColumnVector& W, Eigen::MatrixXd& N11, Eigen::MatrixXd& N21, Eigen::MatrixXd& N22, Eigen::VectorXd& W1, Eigen::VectorXd& W2, bool idx_from_zero = true);
	int LibGREAT_LIBRARY_EXPORT rearrange_lsqmatrix(vector<int>& remove_idx, V_SymmetricMatrix & NEQ, V_ColumnVector & W,t_gallpar& allpar);
}


#endif /*  GLSQMAT_H  */