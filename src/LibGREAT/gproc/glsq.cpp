/**
 * @file         glsq.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gproc/glsq.h"
#include "gutils/ginfolog.h"
#include "gio/gfile.h"
#include "gcoders/recover.h"
#include "gutils/gturboedit.h"
#include <sstream>
#include <stdlib.h>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "gproc/ginverse_Eigen.h"
#include <stack>


using namespace std;

namespace great
{

	// Constructors

	t_glsq::t_glsq(t_gsetbase* set) :
		_mode(LSQMODE::EPO),
		_log(nullptr),
		_tempfile(0),
		_obs_total_num(0),
		_npar_tot_num(0),
		_sigma0(0.0), _res_obs(0.0)
	{


		_dx.ReSize(0);

		if (set != nullptr)
		{
			_beg = dynamic_cast<t_gsetgen*>(set)->beg();
			_end = dynamic_cast<t_gsetgen*>(set)->end();
			_interval = dynamic_cast<t_gsetgen*>(set)->sampling();
			_mode = dynamic_cast<t_gsetproc*>(set)->lsq_mode() ;
			_buffer_size = dynamic_cast<t_gsetproc*>(set)->lsq_buffer_size() * 1024 * 1000;
		}

		// get random_tempfile
		stringstream this_addr;
		this_addr << this;
		reset_tempfile("tempfile_" + this_addr.str());

		// default using Eigen Cholosky method
		_solve_matrix = shared_ptr<t_ginverse<Matrix> >(new t_ginverse_Eigen_Cholosky<Matrix>(false));
	}

	t_glsq::~t_glsq()
	{
		if (_tempfile)
		{
			if (_tempfile->is_open())
			{
				_tempfile->close();
			}
			remove(_tempfile->name().c_str());
			delete _tempfile;
			_tempfile = nullptr;
		}

	}

	t_glsq::t_glsq(const t_glsq& Other)
		:_epo(Other._epo),
		_beg(Other._beg),
		_end(Other._end),
		_tempfile(nullptr),
		_log(Other._log),
		_NEQ(Other._NEQ),
		_dx(Other.dx()),
		_dx_final(Other.dx()),
		_W(Other._W),
		_Qx(Other._Qx),
		_mode(Other._mode),
		_update_lsqpar(Other._update_lsqpar),
		_x_solve(Other._x_solve),
		_stdx(Other._stdx),
		_res_obs(Other._res_obs),
		_npar_tot_num(Other._npar_tot_num),
		_interval(Other._interval),
		_vtpv(Other._vtpv), _sigma0(Other._sigma0), _obs_total_num(Other._obs_total_num),
		_obs_total_num_epo(Other._obs_total_num_epo),
		_solve_matrix(Other._solve_matrix),
		_buffer_size(Other._buffer_size)
	{
		stringstream this_addr;
		this_addr << this;
		reset_tempfile("tempfile_" + this_addr.str());
	}

	int t_glsq::add_parameter(const t_gpar& par)
	{
		// add the par
		_x_solve.addParam(par);
		_NEQ.addBackZero();
		_W.addBackZero();
		_npar_tot_num++;
		return 1;
	}

	void t_glsq::add_par_state_equ(par_type par_type, int order, double dt, double noise)
	{
		_update_lsqpar->set_par_state_mode(par_type, order, dt, noise);
	}

	void t_glsq::set_update_par(shared_ptr<t_gupdatepar> updatepar)
	{
		_update_lsqpar = updatepar;
	}

	void t_glsq::set_solve_matrix(shared_ptr<t_ginverse<Matrix> > solve_matrix)
	{
		_solve_matrix = solve_matrix;
	}

	bool t_glsq::update_parameter(const t_gtime & epoch, vector<t_gsatdata>& obsdata, bool matrix_remove, bool write_temp)
	{
		if (!_update_lsqpar)
		{
			throw std::runtime_error("_update_lsqpar is empty!!!");
		}
		
		t_gupdateparinfo remove_info = _update_lsqpar->get_all_update_parameters(epoch, _x_solve, obsdata);

		vector<int> remove_id;
		vector<t_gpar> new_par_list;
		vector<t_gpar> equ_par_list;
		t_glsqEquationMatrix virtual_equ;
		remove_info.get(remove_id, new_par_list,equ_par_list, virtual_equ);

		// First. Add virtual equ
		for (const auto& par : equ_par_list) 
		{
			add_parameter(par);
		}
		add_equation(virtual_equ,epoch,false);

		// Second. Remove old pars
		remove_info.get(remove_id);
		if (matrix_remove) 
		{
			// use matrix remove way
			remove_parameter(remove_id, write_temp);
		}
		else 
		{
			// use one by one remove way
			//rearrange before remove
			sort(remove_id.begin(), remove_id.end());

			//write_change
			rearrange_lsqmatrix(remove_id, _NEQ, _W, _x_solve);
			_write_parchage(remove_id);

			int max_i = _NEQ.num();
			int beg_i = _NEQ.num()-remove_id.size();
			for (int i =max_i;i>beg_i; i--) 
			{
				remove_parameter(i, write_temp);
			}

		}
		
		// Third. Add new pars
		remove_info.get(new_par_list);
		for (const auto& par : new_par_list)
		{
			add_parameter(par);
		}
		_x_solve.reIndex();
		return true;
	}

	bool t_glsq::remove_all_ambiguity(const t_gtime& epoch, par_type ambtype, bool matrix_remove, bool write_temp)
	{
		t_gupdateparinfo remove_info;
		
		for (unsigned int i = 0; i < _x_solve.parNumber(); i++)
		{
			if (_x_solve[i].parType == ambtype)
			{
				if (_x_solve[i].end > epoch) _x_solve[i].end = epoch;
				remove_info.add(i + 1);
			}
		}
		vector<int> remove_id;
		remove_info.get(remove_id);
		if (matrix_remove)
		{
			// use matrix remove way
			remove_parameter(remove_id, write_temp);
		}
		else
		{
			// use one by one remove way
			sort(remove_id.begin(), remove_id.end());
			for (int i = 0; i < remove_id.size(); i++)
			{
				remove_parameter(remove_id[i] - i, write_temp);
			}
		}
		_x_solve.reIndex();
		return true;
	}

	int t_glsq::set_parameter(const t_gallpar& parameters)
	{
		_x_solve = parameters;

		int numpar = _x_solve.parNumber();

		_W.resize(numpar);
		_NEQ.resize(numpar);

		return 1;
	}
	int t_glsq::add_equation(const Matrix& B, const DiagonalMatrix& P, const ColumnVector& l, const t_gtime& epoch)
	{
		int numpar = _x_solve.parNumber();

		if (B.Ncols() != numpar ||
			B.Nrows() != P.Nrows() ||
			B.Nrows() != l.Nrows() ||
			P.Nrows() != l.Nrows())
		{
			if (_log) _log->comment(1, "t_glsq::set_parameter", "input matrix is wrong row/col!");
			throw exception();
		}

		_epo = epoch;

		t_glsqEquationMatrix equ;
		vector<string> empty_site(B.Nrows(), "");
		vector<string> empty_sat(B.Nrows(), "");
		vector<t_gobscombtype> empty_obstype(B.Nrows(), t_gobscombtype());
		equ.add_equ(B, P, l,empty_site,empty_sat,empty_obstype);

		_NEQ.add(equ);
		_W.add(equ);

		_obs_total_num += equ.num_equ();
		_res_obs += equ.res_equ();

		return 1;

	}
	int t_glsq::add_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch,bool write_temp)
	{
		if (write_temp) 
		{
			this->write_equation(equ, epoch);
		}

		_epo = epoch;

		_NEQ.add(equ);
		_W.add(equ);

		_obs_total_num += equ.num_equ();
		_res_obs       += equ.res_equ();

		return 1;
	}

	int t_glsq::get_equ_obs_total_num()
	{
		if (_obs_total_num_epo == 0)
		{
			_obs_total_num_epo = _obs_total_num;
			return _obs_total_num_epo;
		}
		else
		{
			int result = _obs_total_num - _obs_total_num_epo;
			_obs_total_num_epo = _obs_total_num;
			return result;
		}
		return -1;
	}
	
	int t_glsq::del_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch)
	{
		_epo = epoch;

		_NEQ.del(equ);
		_W.del(equ);

		_obs_total_num -= equ.num_equ();
		_res_obs -= equ.res_equ();

		return 1;
	}
	
	
	bool  t_glsq::write_equation(const t_glsqEquationMatrix& equ, const t_gtime& epoch)
	{
		if (_tempfile == nullptr)
		{
			if (_log) _log->comment(1, "t_glsq::write_equation", "the tempfile don't exist!");
			throw exception();
		}

		int len_size = 0;

		// write coefficinet Matrix
		for (int Row = 1; Row <= equ.num_equ(); Row++)
		{

			// length of record
			len_size = 0;

			// write identifier
			_tempfile->write("obs", 3);
			len_size += 3;

			// write parameter number
			int par_num = equ.B[Row - 1].size();
			_tempfile->write((char*)&par_num, SIZE_INT);
			len_size += SIZE_INT;

			// write time
			_tempfile->write((char*)&epoch, sizeof(t_gtime));
			len_size += sizeof(t_gtime);

			// write ambflag
			int is_newamb = equ.is_newamb(Row - 1) ? 1 : 0;
			_tempfile->write((char*)&is_newamb, SIZE_INT);
			len_size += SIZE_INT;

			// write info
			string out_info;
			out_info = equ.get_sitename(Row - 1) + "    "
				+ equ.get_satname(Row - 1) + "     "
				+ equ.get_obscombtype2str(Row - 1) + "     ";
			int len_out_info = out_info.size();
			_tempfile->write((char*)&len_out_info, SIZE_INT);
			len_size += SIZE_INT;
			_tempfile->write(out_info.c_str(),len_out_info);
			len_size += len_out_info;

			// write coefficient
			for (int Col = 1; Col <= par_num; Col++) {
				double coeff = equ.B[Row - 1][Col - 1].second;
				int loc = equ.B[Row - 1][Col - 1].first;

				_tempfile->write((char*)&loc, SIZE_INT);
				len_size += SIZE_INT;
				_tempfile->write((char*)&coeff, SIZE_DBL);
				len_size += SIZE_DBL;

			}

			// Write Max loc 
			int Max_loc = _x_solve.parNumber();
			_tempfile->write((char*)&Max_loc, SIZE_INT);
			len_size += SIZE_INT;

			// write P
			double p = equ.P[Row - 1];
			_tempfile->write((char*)&p, SIZE_DBL);
			len_size += SIZE_DBL;

			// write res
			double res = equ.l[Row - 1];
			_tempfile->write((char*)&res, SIZE_DBL);
			len_size += SIZE_DBL;

			_tempfile->write((char*)&len_size, SIZE_INT);
		}

		return true;
	}

	int t_glsq::solve_NEQ()
	{
		if (_NEQ.num() != _W.num())
		{
			if (_log) _log->comment(1, "t_glsq::solve_NEQ", "input Matrix NEQ or W is Wrong!");
			return 0;
		}

		// add apriori
		add_apriori_weight();

		SymmetricMatrix NEQ_Matrix = _NEQ.changeNewMat();
		ColumnVector W_Matrix = _W.changeNewMat();

		if (double_eq(NEQ_Matrix.maximum_absolute_value(), 0.0))
		{
			_dx.ReSize(_x_solve.parNumber()); _dx = 0.0;
			throw NPDException(Matrix(0.0,0,0));
		}
		else
		{
			vector<int> zero_idx;

			for (int Row = 1; Row <= NEQ_Matrix.Nrows(); Row++)
			{
				if (double_eq(NEQ_Matrix(Row, Row), 0.0))
				{
					zero_idx.push_back(Row);
				}
			}

			if (zero_idx.size() == 0)
			{
				_solve_equation(NEQ_Matrix, W_Matrix, _dx, _stdx);
			}
			else
			{

				V_SymmetricMatrix remove_NEQ = _NEQ;
				V_ColumnVector remove_W = _W;
				ColumnVector temp_dx, temp_dq;
				for (auto iter = zero_idx.rbegin(); iter != zero_idx.rend(); iter++) {
					remove_NEQ.remove(*iter);
					remove_W.remove(*iter);
				}

				SymmetricMatrix temp_NEQ = remove_NEQ.changeNewMat();
				ColumnVector temp_W = remove_W.changeNewMat();
				_solve_equation(temp_NEQ, temp_W, temp_dx, temp_dq);

				_dx.ReSize(NEQ_Matrix.Nrows()); _dx = 0.0;
				_stdx.ReSize(NEQ_Matrix.Nrows()); _stdx = 0.0;
				SymmetricMatrix temp_Qx = _Qx;
				_Qx.resize(NEQ_Matrix.Nrows()); _Qx = 0.0;

				int idx = 1;
				int idy = 1;
				for (int Row = 1; Row <= NEQ_Matrix.Nrows(); Row++) {
					// not zero
					if (find(zero_idx.begin(), zero_idx.end(), Row) == zero_idx.end()) {
						idy = 1;
						for (int Col = 1; Col <= NEQ_Matrix.Ncols(); Col++)
						{
							// not zero
							if (find(zero_idx.begin(), zero_idx.end(), Col) == zero_idx.end()) {
								_Qx(Row, Col) = temp_Qx(idx, idy);
								//cout << setw(12) << setprecision(8) << temp_Qx(idx, idy);
								idy++;
							}
						}
						_dx(Row) = temp_dx(idx);
						_stdx(Row) = temp_dq(idx);
						idx++;
						//cout << endl;
					}
				}

				if (idx != temp_dx.Nrows() + 1)
					throw exception();

			}


			// slove sigama0
			_vtpv = _res_obs;
			for (int i = 1; i <= _dx.Nrows(); i++) {
				_vtpv -= W_Matrix(i) * (_dx(i));
			}
			_sigma0 = sqrt(abs(_vtpv) / (_obs_total_num - _npar_tot_num));
			cout << " sigma0 = " << abs(_sigma0) << " ntot = " << _obs_total_num << " npar = " << _npar_tot_num << endl;

			if (_obs_total_num - _npar_tot_num < 0) _sigma0 = -1.0;

			for (int i = 1; i <= _stdx.Nrows(); i++) {
				_stdx(i) = sqrt(_stdx(i)) * _sigma0;
			}
		}

		_dx_final << _dx;
		return 1;
	}

	void t_glsq::print_Qx()
	{
		cout << "#PAR: " << endl;
		for (unsigned int i = 0; i < _x_solve.parNumber(); i++) {
			cout << " " << setw(14) << gpar2str(_x_solve[i]);
		}
		cout << endl;
		cout << "#STD:" << endl;
		for (int i = 1; i <= _stdx.Nrows(); i++) {
			cout << setw(15) << setprecision(3) << fixed << scientific << _stdx(i) ;
		}
		cout << endl;
		cout << "#DXX: " << endl;
		for (int i = 1; i <= _Qx.Nrows(); ++i) {
			for (int j = 1; j <= _Qx.Ncols(); ++j) {
				cout << setw(15) << setprecision(3) << fixed << scientific << _Qx(i, j)*_sigma0*_sigma0;
			}
			cout << endl;
		}
		cout << endl;
	}

	void t_glsq::solve_x()
	{
		if (_NEQ.num() != _W.num())
		{
			if (_log) _log->comment(1, "t_glsq::solve_x", "input Matrix NEQ or W is Wrong!");
			throw runtime_error("t_glsq::solve_x:input Matrix NEQ or W is Wrong!");
		}


		// add apriori
		add_apriori_weight();

		SymmetricMatrix NEQ_Matrix = _NEQ.changeNewMat();
		ColumnVector W_Matrix = _W.changeNewMat();

		if (double_eq(NEQ_Matrix.maximum_absolute_value(), 0.0))
		{
			_dx.ReSize(_x_solve.parNumber()); _dx = 0.0;
			throw NPDException(Matrix(0.0, 0, 0));
		}
		else
		{
			vector<int> zero_idx;

			for (int Row = 1; Row <= NEQ_Matrix.Nrows(); Row++)
			{
				if (double_eq(NEQ_Matrix(Row, Row), 0.0))
				{
					zero_idx.push_back(Row);
				}
			}

			if (zero_idx.size() == 0)
			{
				_solve_x(NEQ_Matrix, W_Matrix, _dx);
			}
			else
			{

				V_SymmetricMatrix remove_NEQ = _NEQ;
				V_ColumnVector remove_W = _W;
				ColumnVector temp_dx, temp_dq;
				for (auto iter = zero_idx.rbegin(); iter != zero_idx.rend(); iter++) {
					remove_NEQ.remove(*iter);
					remove_W.remove(*iter);
				}

				SymmetricMatrix temp_NEQ = remove_NEQ.changeNewMat();
				ColumnVector temp_W = remove_W.changeNewMat();
				_solve_x(temp_NEQ, temp_W, temp_dx);

				_dx.ReSize(NEQ_Matrix.Nrows()); _dx = 0.0;

				int idx = 1;
				int idy = 1;
				for (int Row = 1; Row <= NEQ_Matrix.Nrows(); Row++) {
					// not zero
					if (find(zero_idx.begin(), zero_idx.end(), Row) == zero_idx.end()) {
						_dx(Row) = temp_dx(idx);
						idx++;
					}
				}

				if (idx != temp_dx.Nrows() + 1)
					throw exception();
			}

		}


		// slove sigama0
		_vtpv = _res_obs;
		for (int i = 1; i <= _dx.Nrows(); i++) {
			_vtpv -= W_Matrix(i) * (_dx(i));
		}
		_sigma0 = sqrt(abs(_vtpv) / (_obs_total_num - _npar_tot_num));
		cout << " sigma0 = " << abs(_sigma0) << " ntot = " << _obs_total_num << " npar = " << _npar_tot_num << endl;

		if (_obs_total_num - _npar_tot_num < 0) _sigma0 = -1.0;

		_dx_final << _dx;
		return;
	}

	int t_glsq::remove_parameter(const int& idx, bool write_temp)
	{
		add_apriori_weight(idx);
		double center_value = _NEQ.center_value(idx);
		double w_remove = _W.num(idx);
		if (center_value!=0.0){
			// compute the ltpl
			_res_obs -= w_remove * w_remove / center_value;

		}
		// write the coeff to tempfile 
		if (write_temp) _write_coefficient(idx);
		// remove NEQ W Matrix
		if (!remove_lsqmatrix(idx, _NEQ, _W))
		{
			throw "remove lsq matrix error!";
		}

		// remove the par
		_x_solve.delParam(idx - 1);
		_npar_tot_num--;

		return 1;
	}

	int t_glsq::remove_parameter(vector<int>& idx, bool write_temp)
	{
		if (idx.empty()) 
		{
			return 1;
		}

		chrono::high_resolution_clock::time_point beg_t;
		chrono::high_resolution_clock::time_point end_t;

		beg_t = chrono::high_resolution_clock::now();

		//rearrange before remove
		sort(idx.begin(), idx.end());
		int zero_size = rearrange_lsqmatrix(idx,_NEQ,_W,_x_solve);
		_write_parchage(idx);
		end_t = chrono::high_resolution_clock::now();
		cout << "Rearrange Spend time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

		beg_t = chrono::high_resolution_clock::now();
		//process zero
		int rows_new = _NEQ.num()-idx.size();
		int rows_remove =  idx.size()-zero_size;
		for (int zero_idx= _NEQ.num();zero_idx > rows_new+rows_remove; zero_idx--) 
		{
			remove_parameter(zero_idx, write_temp);
		}
		end_t = chrono::high_resolution_clock::now();
		cout << "Process zero Spend time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

		int rows = _NEQ.num();
		// add apriori 
		for (int i =rows;i>rows_new;i--)
		{
			add_apriori_weight(i);
		}

		beg_t = chrono::high_resolution_clock::now();
		vector<double> temp_NEQ(rows*rows,0.0);

		for (int i =0;i<rows;i++){
			memcpy(&temp_NEQ[i*rows],&_NEQ._element[i][0],sizeof(double)*(i+1));
		}

		vector<double> temp_Eigen_NEQ(rows*rows,0.0);
		for (int i =0;i<rows;i++){
			for (int j=0;j<=i;j++){
				temp_Eigen_NEQ[i*rows+j] = _NEQ._element[i][j];
				temp_Eigen_NEQ[j*rows+i] = _NEQ._element[i][j];
			}
		}

		Eigen::VectorXd W1(Eigen::Map<Eigen::VectorXd> (&_W._element[0],rows_new));
		Eigen::VectorXd W2(Eigen::Map<Eigen::VectorXd> (&_W._element[rows_new],rows_remove));

		Eigen:: MatrixXd A(Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<> >(&temp_Eigen_NEQ[0],rows_new,rows_new,Eigen::OuterStride<>(rows)));
		Eigen:: MatrixXd B(Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<> > (&temp_Eigen_NEQ[rows_new],rows_remove,rows_new,Eigen::OuterStride<>(rows)));
		Eigen:: MatrixXd C(Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<> > (&temp_Eigen_NEQ[rows_new*rows+rows_new],rows_remove,rows_remove,Eigen::OuterStride<>(rows)));
		
		end_t = chrono::high_resolution_clock::now();
		cout << "Init Matrix time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

		Eigen::MatrixXd CL,CLi,Ci,Bt,CiB,CiW2,NEQ_new,W_new;
		vector<double> CiB_vec(rows_remove*rows_new,0.0);
		vector<double> CiW2_vec(rows_remove,0.0);
		if (rows_remove!=0){
			beg_t = chrono::high_resolution_clock::now();
			//Colosky 
			CL = Eigen::LLT<Eigen::MatrixXd>(C).matrixL();
			CLi.noalias() = CL.inverse();
			Ci.noalias() = CLi.transpose() * CLi;
			
			end_t = chrono::high_resolution_clock::now();
			cout << "Eigen Cholosky Spend time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

			beg_t = chrono::high_resolution_clock::now();
			// transpose
			Bt.noalias() = B.transpose();
			
			// multiply
			CiB.noalias() = Ci * B;
			CiW2.noalias() = Ci * W2;

			// multiply decmial  
			NEQ_new.noalias() = A-Bt*CiB;
			W_new.noalias() = W1 - Bt*CiW2;


			//get new res_obs
			_res_obs = _res_obs - (W2.transpose()*CiW2)(0,0);

			// ouptut matrix
			for (int i = 0; i < rows_new; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					_NEQ._element[i][j] = NEQ_new(i, j);
				}
				_W._element[i] = W_new(i);
			}
			end_t = chrono::high_resolution_clock::now();
			cout << "Eigen Multiply and Resize Matrix Spend time is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;
		}

		// Resize
		for (int i = rows - 1; i >= rows_new; i--)
		{
			_NEQ.remove(i + 1);
			_W.remove(i + 1);
		}

		beg_t = chrono::high_resolution_clock::now();
		// _write_coefficient
		if (_tempfile) 
		{
			for (int i = rows_remove-1; i>=0; i--)
			{
				int len_size = 0;
				// write identifier
				_tempfile->write("prt", 3);
				len_size += 3;

				// write parameter number
				vector<pair<int, double> > par_record;
				double w_record;

				w_record = CiW2(i);
				for (int j = 0; j < rows_new; j++)
				{
					if (CiB(i, j)!=0.0)
					{
						par_record.push_back(make_pair(j+1, CiB(i, j)));
					}
				}

				// write par num
				int par_num = par_record.size();
				_tempfile->write((char*)&par_num, SIZE_INT);
				len_size += SIZE_INT;

				// write parameter index
				int par_idx = rows_new+i+1;
				_tempfile->write((char*)&(par_idx), SIZE_INT);
				len_size += SIZE_INT;
				par_idx--;

				// write par name station satlltie
				string name = _x_solve[par_idx].str_type();
				name = _x_solve[par_idx].site + "_" + name;
				if (_x_solve[par_idx].prn == "" && _x_solve[par_idx].parType != par_type::GLO_IFB) name = name + "_";
				_tempfile->write(name.c_str(), 20);
				len_size += 20;


				// write par time
				_tempfile->write((char*)&_x_solve[par_idx].beg, sizeof(t_gtime));
				len_size += sizeof(t_gtime);
				_tempfile->write((char*)&_x_solve[par_idx].end, sizeof(t_gtime));
				len_size += sizeof(t_gtime);


				// write parameter value
				double value = _x_solve[par_idx].value();
				_tempfile->write((char*)&value, SIZE_DBL);
				len_size += SIZE_DBL;

				// write zhd for ztd retrieval
				double zhd = _x_solve[par_idx].zhd;
				_tempfile->write((char*)&zhd, SIZE_DBL);
				len_size += SIZE_DBL;

				// write coefficient BTPB
				for (int ipar = 1; ipar <= par_num; ipar++)
				{
					_tempfile->write((char*)&(par_record[ipar - 1].first), SIZE_INT);
					len_size += SIZE_INT;
					_tempfile->write((char*)&(par_record[ipar - 1].second), SIZE_DBL);
					len_size += SIZE_DBL;
				}

				// write BTPL
				_tempfile->write((char*)&w_record, SIZE_DBL);
				len_size += SIZE_DBL;

				_tempfile->write((char*)&len_size, SIZE_INT);

				// remove par
				_x_solve.delParam(par_idx);
				_npar_tot_num--;
			}
		}
		else{
			for (int i =rows_remove-1;i>=0;i--){
				// remove par
				_x_solve.delParam(rows_new+i);
				_npar_tot_num--;
			}
		}

		end_t = chrono::high_resolution_clock::now();
		cout << "Write tempfile is " << chrono::duration_cast<chrono::milliseconds>(end_t-beg_t).count() / 1000.0 << " sec " << endl;

		return 1;
	}
	bool t_glsq::recover_parameters(t_gallrecover& allrecover)
	{
		gtrace("_glsq::recover_parameters");


		if (!_tempfile)
		{
			if (_log) 
			{
				_log->comment(0, " t_glsq::recover_parameters", "have no tempfile Can't recover!");
			}
			return false;
		}

		_tempfile->flush();
		_tempfile->close();

		t_greadtemp tempfile_in(_tempfile->mask());

		int len_size = size_int;
		tempfile_in.seekg_from_cur(-len_size);
		tempfile_in.read((char*)&len_size, size_int);

		// loop for tempfile
		while (len_size != -1) 
		{
			len_size += size_int;
			//_tempfile->seekg(-len_size, ios::cur);
			tempfile_in.seekg_from_cur(-len_size);

			char identifiy[3];
			//_tempfile->read(identifiy, 3);
			tempfile_in.read(identifiy, 3);

			string temp(identifiy, 3);

			if (temp == "obs") 
			{
				_recover_obs(tempfile_in,allrecover);
			}
			else if (temp == "par") 
			{
				_recover_par(tempfile_in,allrecover);
			}
			else if (temp == "prt") 
			{
				_recover_par_part(tempfile_in, allrecover);
			}
			else if (temp == "swp") {
				_recover_par_swap(tempfile_in);
			}
			else {
				break;
			}

			// back 
			len_size -= size_int;
			tempfile_in.seekg_from_cur(-len_size);

			// read length
			len_size = size_int;
			tempfile_in.seekg_from_cur(-len_size);
			tempfile_in.read((char*)&len_size, size_int);
		}

		allrecover.set_interval(_interval);
		allrecover.set_sigma0(_sigma0);
		return true;
	}
	void t_glsq::get_result_parameter(t_gallrecover& allrecover)
	{

		string par_name;
		for (unsigned int i = 0; i < _x_solve.parNumber(); i++)
		{
			t_gpar par = _x_solve[i];
			t_grecover_par recover_par(_x_solve[i], _dx(i+1));
			allrecover.add_recover_par(recover_par);
		}
		
	}
	void t_glsq::add_apriori_weight()
	{
		// add apriori weight
		if (_x_solve.parNumber() != _NEQ.num()) 
		{
			throw exception();
		}

		if (_x_solve.parNumber() == 0) 
		{
			return;
		}

		int neq_num = _NEQ.num();
		for (int Row = 1; Row <= neq_num; Row++)
		{
			double apriori = _x_solve[Row - 1].apriori();

			if (_NEQ.num(Row, Row) == 0.0) continue;
			if (apriori == 0.0)            continue;
			_NEQ.num(Row, Row) = _NEQ.num(Row, Row) + 1 / (apriori * apriori);
		}
	}
	void t_glsq::add_apriori_weight(const int& idx)
	{
		// add apriori weight
		double apriori = _x_solve[idx - 1].apriori();
		if (!double_eq(apriori, 0.0)) 
		{
			_NEQ.num(idx, idx) = _NEQ.num(idx, idx) + 1 / (apriori*apriori);
		}
		return;
	}
	void t_glsq::add_apriori_weight(const int& idx, const double& value, const double& weight)
	{
		// add apriori weight

		// add in NEQ
		double P = weight;
		_NEQ.num(idx, idx) = _NEQ.num(idx, idx) + P;

		// add in W
		double l = value - _x_solve[idx - 1].value();
		_W.num(idx) = _W.num(idx) + P * l;

		return;
	}
	int t_glsq::lsq_sysbias_constraint()
	{
		//impose constraint on IFB/ISB parameters, the sum of all IFB / ISB to a satellite is zero
		// fix the ISB constraint and add IFB constraint
		int idx0 = -1;
		int idx1 = -1;
		int nx = 0;
		map<int, bool> isused;
		map<int, int> iptx;
		for (unsigned int ipar = 0; ipar < _x_solve.parNumber(); ipar++) {
			isused[ipar] = false;
			// ignore the new ISB pars of current epoch
			if (t_gpar::is_sysbias(_x_solve[ipar].parType) && _x_solve[ipar].end < _epo) {
				if (idx0 == -1) idx0 = ipar;
				idx1 = ipar;
			}
		}
		if (idx0 == -1 && idx1 == -1) return 0;
		for (int ipar = idx0; ipar <= idx1 - 1; ipar++)
		{
			if (!t_gpar::is_sysbias(_x_solve[ipar].parType)) continue;
			if (isused[ipar]) continue;
			if (double_eq(_NEQ.num(ipar + 1, ipar + 1), 0)) continue;
			nx = 0;
			iptx[nx] = ipar;
			isused[ipar] = true;
			for (int jpar = ipar + 1; jpar <= idx1; jpar++)
			{
				if (isused[jpar]) continue;
				if (_x_solve[ipar].parType != _x_solve[jpar].parType) continue;
				if (_x_solve[ipar].prn != _x_solve[jpar].prn) continue;
				if (double_eq(_NEQ.num(jpar + 1, jpar + 1), 0)) continue;
				nx++;
				iptx[nx] = jpar;
				isused[jpar] = true;
			}
			if (nx == 0) continue;

			for (int i = 0; i <= nx; i++)
			{
				int ip1 = iptx[i];
				for (int j = 0; j <= i; j++)
				{
					int jp1 = iptx[j];
					_NEQ.num(ip1 + 1, jp1 + 1) = _NEQ.num(ip1 + 1, jp1 + 1) + 1e5;
				}
			}
		}
		return 1;
	}

	int t_glsq::lsq_clk_constraint()
	{
		vector< pair<int, double> > B;
		for (size_t ipar = 0; ipar < _x_solve.parNumber(); ++ipar) {
			if (_x_solve[ipar].str_type().substr(0, 7) != "CLK_SAT") continue;
			if (double_eq(_x_solve.getAmbParam(ipar).value(), 0.0)) continue;
			B.push_back(make_pair(ipar + 1, 1));
		}
		t_glsqEquationMatrix virtual_equ;
		t_gobscombtype type;
		virtual_equ.add_equ(B, 1e5, 0.0, "", "", type, false);
		this->add_equation(virtual_equ, _epo, false);
		return 1;
	}

	void t_glsq::print_matrx()
	{
		/*_NEQ.print();
		_W.print();*/
		cout << "neq center is" << endl;
		for (unsigned int i = 0; i < _x_solve.parNumber(); i++) 
		{
			if (_NEQ.center_value(i + 1) == 0.0)  continue;
			cout << _x_solve[i].site + "_" + _x_solve[i].str_type() 
				 << " " 
				 << setw(20) << setprecision(15) << scientific << _NEQ.center_value(i + 1) 
				 << endl;
		}
	}

	bool t_glsq::par_alive(const int& idx)
	{
		return !double_eq(_NEQ.center_value(idx+1), 0);
	}

	void t_glsq::setlog(t_glog* log)
	{
		_log = log;
	}
	const SymmetricMatrix& t_glsq::Qx() const
	{
		return _Qx;
	}
	SymmetricMatrix	t_glsq::NEQ() const
	{
		return _NEQ.changeNewMat();
	}
	ColumnVector t_glsq::W()
	{
		return _W.changeNewMat();
	}
	V_SymmetricMatrix t_glsq::v_NEQ() const
	{
		return _NEQ;
	}

	V_ColumnVector t_glsq::v_W() const
	{
		return _W;
	}
	ColumnVector t_glsq::dx() const
	{
		return _dx_final;
	}
	ColumnVector t_glsq::stdx() const
	{
		return _stdx;
	}
	double t_glsq::Qx(const int& col, const int& row) const
	{
		return _Qx(col, row);
	}
	double t_glsq::dx(int idx) const
	{
		return _dx_final(idx);
	}
	double t_glsq::sigma0() const
	{
		return _sigma0;
	}
	double t_glsq::stdx(int idx) const
	{
		return _stdx(idx);
	}
	double t_glsq::vtpv() const
	{
		return _vtpv;
	}
	int t_glsq::nobs_total() const
	{
		return _obs_total_num;
	}
	t_gtime t_glsq::epo() const {
		return _epo;
	}
	void t_glsq::epo(const t_gtime& t)
	{
		_epo = t;
	}
	gnut::LSQMODE t_glsq::mode()
	{
		return _mode;
	}
	void t_glsq::reset_tempfile(string filename)
	{
		if (_tempfile)
		{
			if (_tempfile->is_open())
			{
				_tempfile->close();
			}
			remove(_tempfile->name().c_str());
			delete _tempfile;
		}

		_tempfile = new t_giobigf(filename,_buffer_size);
		_tempfile->tsys(t_gtime::GPS);
		_tempfile->mask(filename);
		_tempfile->append(false);
		_tempfile->open(filename, ios::out | ios::trunc | ios::binary);

		// for begin
		int identify = -1;
		_tempfile->write((char*)&identify, SIZE_INT);
	}

	int t_glsq::_write_parchage(const vector<int>& remove_id)
	{
		if (!_tempfile)
		{
			if (_log)
			{
				_log->comment(1, "t_glsq::_write_coefficient", "don't init the tempfile");
			}
			return -1;

		}

		int len_size = 0;
		_tempfile->write("swp", 3);
		len_size += 3;

		int total_size = _x_solve.parNumber();
		_tempfile->write((char*)&total_size, sizeof(int));
		len_size += sizeof(int);
		
		int remove_size = remove_id.size();
		_tempfile->write((char*)&remove_size, sizeof(int));
		len_size += sizeof(int);

		for (int i = 0; i < remove_size; i++) {
			_tempfile->write((char*)&remove_id[i], sizeof(int));
			len_size += sizeof(int);
		}

		_tempfile->write((char*)&len_size, sizeof(int));
		return 0;
	}
	int t_glsq::_write_coefficient(int idx)
	{
		ostringstream os;
		if (idx < 1 || idx > _NEQ.num()) 
		{
			if (_log)
			{
				_log->comment(1, "t_glsq::_write_coefficient", "the input idx is Wrong!");
			}
			return -1;
		}
		if (!_tempfile)
		{
			if (_log)
			{
				_log->comment(1, "t_glsq::_write_coefficient", "don't init the tempfile");
			}
			return -1;

		}
		int len_size = 0;
		// write identifier
		_tempfile->write("par", 3);
		len_size += 3;

		// write parameter number
		vector<pair<int, double> > par_record;

		double w_record = _W.num(idx);
		int NEQ_num = _NEQ.num();
		for (int col = 1; col <= NEQ_num; col++)
		{
			double value = _NEQ.num(idx, col);
			if (!double_eq(value, 0.0))
			{
				par_record.push_back(make_pair(col, value));
			}
		}

		// write par num
		int par_num = par_record.size();
		_tempfile->write((char*)&par_num, SIZE_INT);
		len_size += SIZE_INT;

		// write parameter index
		_tempfile->write((char*)&idx, SIZE_INT);
		len_size += SIZE_INT;

		// write par name station satlltie
		const t_gpar& par_tmp = _x_solve[idx - 1];
		string name = par_tmp.site + "_" + par_tmp.str_type();
		if (par_tmp.prn == "" && par_tmp.parType != par_type::GLO_IFB) name = name + "_";
		_tempfile->write(name.c_str(), 20);
		len_size += 20;


		// write par time
		_tempfile->write((char*)&par_tmp.beg, sizeof(t_gtime));
		len_size += sizeof(t_gtime);
		_tempfile->write((char*)&par_tmp.end, sizeof(t_gtime));
		len_size += sizeof(t_gtime);


		// write parameter value
		double value = par_tmp.value();
		_tempfile->write((char*)&value, SIZE_DBL);
		len_size += SIZE_DBL;

		double zhd = par_tmp.zhd;
		_tempfile->write((char*)&zhd, SIZE_DBL);
		len_size += SIZE_DBL;

		// write coefficient BTPB
		for (int ipar = 1; ipar <= par_num; ipar++)
		{
			_tempfile->write((char*)&(par_record[ipar - 1].first), SIZE_INT);
			len_size += SIZE_INT;
			_tempfile->write((char*)&(par_record[ipar - 1].second), SIZE_DBL);
			len_size += SIZE_DBL;
		}

		// write BTPL
		_tempfile->write((char*)&w_record, SIZE_DBL);
		len_size += SIZE_DBL;

		_tempfile->write((char*)&len_size, SIZE_INT);

		return 1;
	}

	void t_glsq::_recover_par(t_greadtemp& tempfile_in,t_gallrecover& _allrecover)
	{

		// read par number
		int par_num;
		tempfile_in.read((char*)&par_num, SIZE_INT);

		if (par_num > _dx.Nrows() + 1) 
		{
			if (_log) 
			{
				_log->comment(1, "t_glsq::_recover_obs", "par recover more than par now Recover Fail!");
			}
			throw exception();
		}


		// read idx
		int idx;
		tempfile_in.read((char*)&idx, SIZE_INT);

		// read par name sat
		char par_info[20];
		tempfile_in.read(par_info, 20);

		t_gtime par_begtime;
		tempfile_in.read((char*)&par_begtime, sizeof(t_gtime));
		t_gtime par_endtime;
		tempfile_in.read((char*)&par_endtime, sizeof(t_gtime));


		// read par value
		double par_value;
		tempfile_in.read((char*)&par_value, SIZE_DBL);

		double par_zhd;
		tempfile_in.read((char*)&par_zhd, SIZE_DBL);

		// N for recording coefficient data
		vector<int> loc;
		vector<double> N;
		for (int ipar = 0; ipar < par_num; ipar++) 
		{
			int i; 
			double j;
			tempfile_in.read((char*)&i, SIZE_INT);
			loc.push_back(i);
			tempfile_in.read((char*)&j, SIZE_DBL);
			N.push_back(j);
		}

		// read omc
		double omc = 0.0;
		tempfile_in.read((char*)&omc, SIZE_DBL);


		// Resize for dx !! Size is par_num+1
		ColumnVector dx_temp = _dx;
		_dx.ReSize(_dx.Nrows() + 1); 
		_dx = 0.0;

		if (idx != 1) 
		{
			_dx.Rows(1, idx - 1) << dx_temp.Rows(1, idx - 1);
		}
		if (idx != _dx.Nrows()) 
		{
			_dx.Rows(idx + 1, _dx.Nrows()) << dx_temp.Rows(idx, _dx.Nrows() - 1);
		}

		dx_temp << _dx;



		// compute the dx
		if (find(loc.begin(), loc.end(), idx) == loc.end()) 
		{
			_dx(idx) = 0;
		}
		else
		{
			_dx(idx) = omc;
			int flag_idx = 0;
			for (int ipar = 0; ipar < par_num; ipar++) 
			{
				if (loc[ipar] != idx) 
				{
					_dx(idx) = _dx(idx) - (N[ipar] * dx_temp(loc[ipar]));
				}
				else 
				{
					flag_idx = ipar;
				}
			}
			_dx(idx) = _dx(idx) / N[flag_idx];
		}


		t_gpar par = str2gpar(par_info);
		par.beg = par_begtime;
		par.end = par_endtime;
		par.value(par_value);
		par.zhd = par_zhd;
		t_grecover_par recover_par(par, _dx(idx));
		_allrecover.add_recover_par(recover_par);

	}
	void t_glsq::_recover_par_part(t_greadtemp & tempfile_in, t_gallrecover & _allrecover)
	{
		// read par number
		int par_num;
		tempfile_in.read((char*)&par_num, SIZE_INT);

		// read idx
		int idx;
		tempfile_in.read((char*)&idx, SIZE_INT);

		// read par name sat
		char par_info[20];
		tempfile_in.read(par_info, 20);

		// read par time
		t_gtime par_begtime;
		tempfile_in.read((char*)&par_begtime, sizeof(t_gtime));
		t_gtime par_endtime;
		tempfile_in.read((char*)&par_endtime, sizeof(t_gtime));


		// read par value
		double par_value;
		tempfile_in.read((char*)&par_value, SIZE_DBL);

		double par_zhd;
		tempfile_in.read((char*)&par_zhd, SIZE_DBL);

		// N for recording coefficient data
		vector<int> loc;
		vector<double> N;
		for (int ipar = 0; ipar < par_num; ipar++) {
			int i; double j;
			tempfile_in.read((char*)&i, SIZE_INT);
			loc.push_back(i);
			tempfile_in.read((char*)&j, SIZE_DBL);
			N.push_back(j);
		}

		// read omc
		double omc = 0.0;
		tempfile_in.read((char*)&omc, SIZE_DBL);


		// Resize for dx !! Size is par_num+1
		ColumnVector dx_temp = _dx;
		_dx.ReSize(_dx.Nrows() + 1); _dx = 0.0;

		if (idx != 1) {
			_dx.Rows(1, idx - 1) << dx_temp.Rows(1, idx - 1);
		}
		if (idx != _dx.Nrows()) {
			_dx.Rows(idx + 1, _dx.Nrows()) << dx_temp.Rows(idx, _dx.Nrows() - 1);
		}

		dx_temp << _dx;



		// compute the dx
		_dx(idx) = omc;
		for (int ipar = 0; ipar < par_num; ipar++) {
			_dx(idx) -= (N[ipar] * dx_temp(loc[ipar]));
		}

		for (int i = 0; i < loc.size(); i++)
		{
			assert(loc[i] != idx);
		}


		t_gpar par = str2gpar(par_info);
		par.beg = par_begtime;
		par.end = par_endtime;
		par.value(par_value);
		par.zhd = par_zhd;
		t_grecover_par recover_par(par, _dx(idx));
		_allrecover.add_recover_par(recover_par);
	}
	void t_glsq::_recover_par_swap(t_greadtemp & tempfile_in)
	{
		// read total size
		int total_size = 0;
		tempfile_in.read((char*)&total_size, sizeof(int));

		// read remove size
		int remove_size = 0;
		tempfile_in.read((char*)&remove_size, sizeof(int));


		// read remove_id;
		vector<int> remove_id(remove_size, 0);
		for (int i = 0; i < remove_size; i++) {
			tempfile_in.read((char*)&remove_id[i], sizeof(int));
		}

		vector<int> remove_id_sort(remove_id.begin(),remove_id.end());
		sort(remove_id_sort.begin(),remove_id_sort.end());

		ColumnVector dx_temp = _dx;
		int idx_old = 0, idx_new = 0;
		for (int i = 1; i <= total_size; i++) {
			if (idx_new < remove_size && i == remove_id_sort[idx_new]) {
				idx_new++;
			}
			else {
				_dx(i) = dx_temp(++idx_old);
			}
		}

		for (int i =0;i<remove_size;i++){
			_dx(remove_id[i]) = dx_temp(total_size-remove_size+1+i);
		}

	}
	void t_glsq::_recover_obs(t_greadtemp& tempfile_in,t_gallrecover& _allrecover)
	{
		// read parnumber
		int par_num;
		tempfile_in.read((char*)&par_num, SIZE_INT);

		// read time
		t_gtime time;
		tempfile_in.read((char*)&time, sizeof(t_gtime));

		// read ambflag;
		int is_newamb;
		tempfile_in.read((char*)&is_newamb, SIZE_INT);

		// read info
		int len_info;
		tempfile_in.read((char*)&len_info, SIZE_INT);
		char* info = new char[len_info];
		tempfile_in.read(info, len_info);


		// B for recording the coefficient of obs equation
		if (par_num > _dx.Nrows()) {
			if (_log) {
				_log->comment(1, "t_glsq::_recover_obs", "par recover more than par now Recover Fail!");
			}
			throw exception();
		}

		//shared_ptr<double> B(new double[par_num]);
		vector<int> loc;
		vector<double> B;

		for (int ipar = 1; ipar <= par_num; ipar++) {
			int i; double j;
			tempfile_in.read((char*)&i, SIZE_INT);
			loc.push_back(i);
			tempfile_in.read((char*)&j, SIZE_DBL);
			B.push_back(j);
		}

		int Max_loc;
		tempfile_in.read((char*)&Max_loc, SIZE_INT);
		if (_dx.Nrows() > Max_loc) {
			ColumnVector _dx_temp = _dx;
			_dx.ReSize(Max_loc);
			_dx << _dx_temp.Rows(1, Max_loc);
		}

		// read p
		double p;
		tempfile_in.read((char*)&p, SIZE_DBL);

		// read old res
		double res;
		tempfile_in.read((char*)&res, SIZE_DBL);

		// compute the new res
		for (int ipar = 0; ipar < par_num; ipar++) {
			res -= B[ipar] * _dx(loc[ipar]);
		}

		stringstream obsinfo(string(info,len_info));
		delete[] info; info = nullptr;
		string site, sat, str_obstype ;
		obsinfo >> site >> sat >> str_obstype;
		t_grecover_equation recover_equ(time, site, sat);
		recover_equ.set_recover_equation(t_gobscombtype(str_obstype), make_pair(p, res),is_newamb);
		_allrecover.add_recover_equation(recover_equ);

	}
	void t_glsq::_remove_zero_element(SymmetricMatrix& B, ColumnVector& l, int idx)
	{
		if (idx < 1 || idx > B.Nrows() || B.Nrows() != l.Nrows()) {
			if (_log) {
				_log->comment(1, "t_glsq::_remove_zero_element", "input is Wrong!");
			}
			throw exception();
		}

		SymmetricMatrix temp_B = B;
		ColumnVector temp_l = l;
		B.ReSize(temp_B.Nrows() - 1); B = 0.0;
		l.ReSize(temp_l.Nrows() - 1); l = 0.0;


		if (idx == 1) {
			B << temp_B.SymSubMatrix(2, temp_B.Nrows());
			l << temp_l.Rows(2, temp_l.Nrows());
		}
		else if (idx == temp_B.Nrows()) {
			B << temp_B.SymSubMatrix(1, temp_B.Nrows() - 1);
			l << temp_l.Rows(1, temp_l.Nrows() - 1);
		}
		else {
			B.SymSubMatrix(1, idx - 1) << temp_B.SymSubMatrix(1, idx - 1);
			B.SubMatrix(idx, B.Nrows(), 1, idx - 1) << temp_B.SubMatrix(idx + 1, temp_B.Nrows(), 1, idx - 1);
			B.SubMatrix(idx, B.Nrows(), idx, B.Ncols()) << temp_B.SubMatrix(idx + 1, temp_B.Nrows(), idx + 1, temp_B.Ncols());

			l.Rows(1, idx - 1) << temp_l.Rows(1, idx - 1);
			l.Rows(idx, l.Nrows()) << temp_l.Rows(idx + 1, temp_l.Nrows());
		}

	}
	void t_glsq::_solve_equation(const SymmetricMatrix& NEQ, const ColumnVector& W, ColumnVector& ans, ColumnVector& Q)
	{

		Matrix Inverse_N = NEQ;
		Matrix x(NEQ.nrows(), 1);

		_solve_matrix->sovle_NEQ(Inverse_N.nrows(), Inverse_N, W, x);

		ans << x;

		_Qx << Inverse_N;
		Q.ReSize(Inverse_N.Nrows());
		for (int i = 1; i <= Inverse_N.Nrows(); i++)
		{
			Q(i) = Inverse_N(i, i);
		}

	}
	void t_glsq::_solve_x(const SymmetricMatrix & NEQ, const ColumnVector & W, ColumnVector & ans)
	{
		Matrix x(NEQ.nrows(), 1);

		_solve_matrix->solve_x(NEQ.nrows(), NEQ, W, x);

		ans << x;
	}
	void t_glsq::change_NEQ(int row, int col, double xx)
	{
		_NEQ.num(row, col) = xx;
	}

	void t_glsq::change_Qx(int row, int col, double xx)
	{
		_Qx(row, col) = xx;
	}

	void t_glsq::change_Qx(SymmetricMatrix Qx)
	{
		_Qx << Qx;
	}

	void t_glsq::change_dx(int n, double xx)
	{
		if (n == 0)
		{
			_dx_final = xx;
		}
		else
		{
			_dx_final(n) = xx;
		}
	}

	void t_glsq::change_dx(ColumnVector dx)
	{
		_dx_final << dx;
	}

	void t_glsq::change_stdx(int n, double xx)
	{

		_stdx(n) = xx;
	}

	void t_glsq::change_stdx(ColumnVector stdx)
	{

		_stdx << stdx;
	}

	void t_glsq::change_sigma0(double sigma)
	{
		_sigma0 = sigma;
	}

	void t_glsq::change_vtpv(double xx)
	{
		_vtpv = xx;
	}

	void t_glsq::set_new_NEQ(const V_ColumnVector& W, const V_SymmetricMatrix& NEQ)
	{
		_W = W;
		_NEQ = NEQ;
	}


	int t_glsq::reset_npar_total(const t_gpar& par)
	{
		// Ionosphere/Clk_Sat/IFCB  estimate per epoch
		if (par.parType == par_type::SION || par.parType == par_type::VION
			|| par.parType == par_type::CLK_SAT
			|| par.str_type().find("IFCB") != string::npos)
		{
			return 1;
		}
		// Ambiguity per Arc
		if (par.str_type().find("AMB") != string::npos && _epo == par.beg)
		{
			return 1;
		}

		// IFB per Arc
		if (par.str_type().find("IFB") != string::npos && _epo == par.beg)
		{
			return 1;
		}

		return 0;
	}



	int t_glsq::npar_number()const
	{
		return _npar_tot_num;
	}
	t_gtime str2gtime(string str_time)
	{

		replace(str_time.begin(), str_time.end(), '-', ' ');
		replace(str_time.begin(), str_time.end(), ':', ' ');
		stringstream time(str_time);
		int year, month, day, hour, minute;
		double sec;
		time >> year >> month >> day >> hour >> minute >> sec;
		t_gtime time_t; time_t.from_ymdhms(year, month, day, hour, minute, sec);
		return time_t;
	}
}
