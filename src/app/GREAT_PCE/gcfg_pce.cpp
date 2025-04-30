#include <iomanip>
#include <sstream>

#include "gcfg_pce.h"

using namespace std;
using namespace pugi;

namespace gnut 
{

// Constructor
// ----------
	t_gcfg_pce::t_gcfg_pce()
 : t_gsetbase(),
   t_gsetgen(),
   t_gsetinp(),
   t_gsetout(),
   t_gsetgnss(),
   t_gsetproc(),
   t_gsetpar(),
   t_gsetrec(),
   t_gsetturboedit(),
   t_gsetflt()
{
  _IFMT_supported.insert(RINEXO_INP);
  _IFMT_supported.insert(RINEXC_INP);
  _IFMT_supported.insert(RINEXN_INP);
  _IFMT_supported.insert(ATX_INP);
  _IFMT_supported.insert(BLQ_INP);
  _IFMT_supported.insert(SP3_INP);
  _IFMT_supported.insert(BIABERN_INP);
  _IFMT_supported.insert(SINEX_INP);
  _IFMT_supported.insert(BIASINEX_INP);
  _IFMT_supported.insert(POLEUT1_INP);
  _IFMT_supported.insert(SATPARS_INP);
  _IFMT_supported.insert(LEAPSECOND_INP);
  _IFMT_supported.insert(DE_INP);
  _IFMT_supported.insert(EPODIR_INP);

  _OFMT_supported.insert(LOG_OUT);
  _OFMT_supported.insert(PPP_OUT);
  _OFMT_supported.insert(FLT_OUT);
  _OFMT_supported.insert(RECCLK_OUT);
  _OFMT_supported.insert(SATCLK_OUT);
  _OFMT_supported.insert(CLK_OUT);
  _OFMT_supported.insert(EPODIR_OUT);
}


// Destructor
// ----------
t_gcfg_pce::~t_gcfg_pce()
{}


// settings check
// ----------
void t_gcfg_pce::check()
{
  t_gsetgen::check();
  t_gsetinp::check();
  t_gsetout::check();
  t_gsetrec::check();
  t_gsetpar::check();
  t_gsetproc::check();
  t_gsetgnss::check();
  t_gsetamb::check();

}

// settings help
// ----------
void t_gcfg_pce::help()
{
  t_gsetbase::help_header();
  t_gsetgen::help();
  t_gsetinp::help();
  t_gsetout::help();
  t_gsetrec::help();
  t_gsetpar::help();
  t_gsetproc::help();
  t_gsetgnss::help();
  t_gsetamb::help();
  t_gsetbase::help_footer();
}

} // namespace
