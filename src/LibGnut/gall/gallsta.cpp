/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2020-03-06	creat the file.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gallsta.h
 * @brief		The base class used to store range biase information.
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2020-03-06
 */

#include "gall/gallsta.h"

using namespace std;

namespace great
{
	t_gallsta::t_gallsta()
	{
		id_type(t_gdata::ALLSTA);
	}

	t_gallsta::~t_gallsta()
	{

	}

	void t_gallsta::addsta(string name, double dx, double dy, double dz)
	{
		sta_cor cor;
		cor.dx = dx;
		cor.dy = dy;
		cor.dz = dz;
		_Sta_Cor[name] = cor;
	}
	
	map<string, sta_cor>& t_gallsta::getsta()
	{
		return _Sta_Cor;
	}



}