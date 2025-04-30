/**
 * @file         glsqprocIF.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/glsqprocIF.h"
#include "gproc/gupdateparIF.h"
#include "gutils/ginfolog.h"

using namespace std;

namespace great
{
	t_glsqprocIF::t_glsqprocIF()
	{

	}

	t_glsqprocIF::t_glsqprocIF(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqproc(set, data, log)
	{
		if (_frequency == 2)
		{
			shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12, _band_index));
			updateIF->set_amb_update_way(_lite_turboedit);
			updateIF->set_interval(_obs_intv);
			_lsq->set_update_par(updateIF);

			shared_ptr<t_gbiasmodel> precise_bias(new t_gprecisebias(data, log, set));
			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, precise_bias,data));
		}
		else
		{
			throw logic_error("you should use DOUBLE in xml");
		}


		gtrace("t_glsqprocIF::t_glsqprocIF");
	}

	t_glsqprocIF::~t_glsqprocIF()
	{
		
	}


}
