/**
 * @file         ginverse.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief	     header files of matrix inverse
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef GINVERSE_H
#define GINVERSE_H

namespace great
{
	template<class MyMatrix>
	class t_ginverse
	{
	public:
		t_ginverse(bool is_from_zero) {
			_beg = is_from_zero ? 0 : 1;
		}

		virtual void solve_x(int num, const MyMatrix& NEQ, const MyMatrix& W, MyMatrix& X) = 0;
		virtual void sovle_NEQ(int num, MyMatrix& NEQ, const MyMatrix& W, MyMatrix& X) = 0;
		virtual void inverse_NEQ(int num, MyMatrix& NEQ) = 0;

	protected:

		int _beg;
	};

}
#endif // !GINVERSE_H
