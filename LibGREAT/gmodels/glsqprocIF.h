/**
 * @file         glsqprocIF.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GLSQPROC_H
#define GLSQPROC_H

#include "gdata/gdata.h"
#include "gall/gallproc.h"
#include "gall/gallprod.h"
#include "gproc/glsq.h"
#include "gmodels/gmodel.h"
#include "gmodels/gprecisemodel.h"
#include "gmodels/glsqproc.h"
#include "gproc/glsqmatrix.h"
#include "gdata/gsatdata.h"
#include "gall/gallobs.h"

#include "gexport/ExportLibGREAT.h"

namespace great
{

	/** @brief class for t_glsqprocIF based on t_glsqproc. */
	class LibGREAT_LIBRARY_EXPORT t_glsqprocIF : public t_glsqproc
	{
	public:
		/** @brief default constructor. */
		t_glsqprocIF();
		t_glsqprocIF(t_gsetbase* set,t_gallproc* data, t_glog* log);
		virtual ~t_glsqprocIF();
		
	protected:

	};
}


#endif