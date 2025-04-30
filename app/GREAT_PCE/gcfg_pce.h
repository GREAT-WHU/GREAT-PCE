
#ifndef GCFG_PCE_H
#define GCFG_PCE_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetpar.h"
#include "gset/gsetrec.h"
#include "gset/gsetamb.h"
#include "gset/gsetturboedit.h"
#include "gset/gsetflt.h"

#include "gall/gallproc.h"
#include "gall/gallprec.h"
#include "gmodels/glsqproc.h"
#include "gproc/gpcelsqIF.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gcoders/dvpteph405.h"
#include "gcoders/poleut1.h"
#include "gio/gfile.h"

using namespace std;
using namespace pugi;
using namespace great;
namespace gnut
{

	class  t_gcfg_pce : public t_gsetgen,
		public t_gsetinp,
		public t_gsetout,
		public t_gsetgnss,
		public t_gsetproc,
		public t_gsetpar,
		public t_gsetamb,
		public t_gsetrec,
		public t_gsetturboedit,
		public t_gsetflt
	{

	public:
		/** @brief default constructor. */
		t_gcfg_pce();
		/** @brief default constructor. */
		~t_gcfg_pce();

		/** @brief settings check. */
		void check();                                 // settings check
		/** @brief settings help. */
		void help();                                  // settings help

	protected:

	private:

	};

} // namespace

#endif