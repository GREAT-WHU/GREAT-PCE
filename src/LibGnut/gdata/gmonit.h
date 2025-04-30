
#ifndef GMONIT_H
#define GMONIT_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  Purpose: 
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <string>
#include <sstream>
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut {

class  LibGnut_LIBRARY_EXPORT  t_gmonit {

 public:
  t_gmonit( string id );
  virtual ~t_gmonit();

  virtual void  show( ostringstream& os, int verb );
//  virtual string show( int verb, int repeat, string& buff ) = 0;

 protected:
  string   _moni_id;

 private:

};

} // namespace

#endif
