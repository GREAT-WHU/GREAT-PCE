
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>

#include "gutils/gcommon.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {


// sampling synchronization filter for epochs (return true if the epoch fits sampling)
// ----------
bool time_sync(const t_gtime& epo, double smp, double scl, t_glog* glog)
{
  gtrace("time_sync");

  if( smp > 0 )
  {
    if( scl > 0 && smp <= 1 ){
      int smpl = (int)round( scl * smp);
      int iepo = (int)(round(epo.sod() *scl
			                     + epo.dsec()*scl)); // synced to .0 day i.e >=1Hz suggests .0 day epochs at least!

      if( smpl == 0 ) return false; // to be save if mixed high/low-rate sampling occurs
      int resi = iepo%smpl;

      if( resi != 0 ) return false;
    }
    else  if( int(round(epo.sod()+epo.dsec()))%int(smp) != 0 ){
      return false;
    }
  }
  return true;
}


} // namespace
