/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <cstring>
#include <sstream>

#include "gio/gfile.h"
#include "gutils/gcommon.h"
#include "gutils/gfileconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gfile::t_gfile()
 : t_gio(),
   _irc(0),
   _gzip(false)
{
  gtrace("t_gfile::construct");
  
  _file  = 0;
  _size = FILEBUF_SIZE;
}


// destructor
// ----------
t_gfile::~t_gfile()
{
  gtrace("t_gfile::destruct");

  reset();
}


// set path
// ----------
int t_gfile::path( string str )
{
  gtrace("t_gfile::path(str)");

  if( str == "" ) return -1;

  size_t idx1 = 0, idx2 = 0;
  string prefix = GFILE_PREFIX;

  // check if file
  idx2 = str.find(prefix);
  if( idx2 == string::npos ){ str = prefix+str; idx2 = 0; }

  // file path
  idx1 = idx2+prefix.length();
  idx2 = str.length();
  if( idx2 != string::npos && idx2 > idx1  ){
      
    string name = str.substr(idx1,idx2-idx1);
    _set_gzip(name);

    if (_file == 0)
        _file = new t_giof(name.c_str());
    else
        return -1;

    if (_file)
        _file->mask(name);
    else
    {
        cerr << "cannot assign file to path!\n";
        _irc++;
        return -1;
    }

    ostringstream ltmp; 
    ltmp << "File: " << int(idx1) << ":" << int(idx2) << " = " << name;
    if( _log ) _log->comment(3, "gfile", ltmp.str());
    
  }else{
    if( _log ) _log->comment(0, "gfile", "warning: path does not contain file://dir/name  [check file://]");
    return -1; 
  }

  t_gio::path(str);
   
  return 1;
}


// get path
// ----------
string t_gfile::path(){
  return GFILE_PREFIX + mask();
}


// get name
// ----------
string t_gfile::name(){
   
  return mask();
}

// irc
// ----------
int t_gfile::irc()const
{
    if (_file) return (_file->irc() + _irc);

    return _irc;
}

// eof
// ----------
bool t_gfile::eof()
{
    if (_file)
        return _file->eof();

    return true;
}


// mask
// ----------
string t_gfile::mask()
{
    if (_file)
        return _file->mask();

    return "";
}

// get md5
// ----------
string t_gfile::md5sum()
{
    gtrace("t_gfile::md5sum");

    if (_file) return  _file->md5sum();

    return "";
}

   
// init read function (read header)
// ----------
int t_gfile::init_read()
{
  gtrace("t_gfile::init_read");
  
  if( ! _coder ){ return 0; }

  char* loc_buff = new char[FILEHDR_SIZE];

  int nbytes = 0;
  vector<string>  errmsg;
  while( ( nbytes = _gio_read( loc_buff, FILEHDR_SIZE ) ) > 0 && _stop != 1 ){
     
    if( _coder->decode_head( loc_buff, nbytes, errmsg ) < 0 ) break;
  }

  if( nbytes < 0 ){ // header not completed properly ?
    ++_irc; if( _coder ) _coder->mesg(GERROR, "Error: Incomplete header identified."); 
  }

  delete [] loc_buff;
  return 1; 
}


// init write function (write header)
// ----------
int t_gfile::init_write()
{ 
  gtrace("t_gfile::init_write");  

  if( ! _coder ){ return 0; }
   
  char* loc_buff = new char[FILEHDR_SIZE];

  vector<string> errmsg;
   
  int nbytes = 0;
  do{
    if( ( nbytes = _coder->encode_head( loc_buff, FILEHDR_SIZE, errmsg ) ) < 0 ) break;
  }while( (nbytes > 0) && (_gio_write( loc_buff, nbytes ) > 0) && _stop != 1  );

  delete [] loc_buff;
      
  return 1; 
}

// reset
// ----------
void t_gfile::reset()
{
    gtrace("t_gfile::reset");

    if (_file) { delete _file;   _file = 0; }
}

// read data
// ----------
int t_gfile::_gio_read( char* buff, int size )
{     
  gtrace("t_gfile::_gio_read");

  if( mask() == "" ) return -1;
  int nbytes =  this->_read(buff,size);

  if( nbytes == 0 && this->eof() ) return -1;
  return nbytes;

//  int nbytes = 0
//  nbytes = this->_read(buff,size);
//  if( this->eof() || nbytes == -1 ){
//    if( nbytes == 0 ) nbytes = -1;
//    }
//  }
//  return nbytes;
}


// send data
// ----------
int t_gfile::_gio_write( const char* buff, int size )
{
  gtrace("t_gfile::_gio_write");

  if( mask() == "" ) return 0;

  this->_write( buff, size );

  return size;
}


// common function for file close
// ----------
int t_gfile::_stop_common()
{ 
  gtrace("t_gfile::_stop_common");

  return t_gio::_stop_common();
}

   
// compressed or not ?
// ----------
void t_gfile::_set_gzip( string name )
{ 
  gtrace("t_gfile::_set_gzip");

  if( // name.compare(name.length()-2,name.length(),".Z")  == 0 || ///  NOT SUPPORTED
      name.compare(name.length()-3,name.length(),".gz") == 0 )

        _gzip = true;
   else _gzip = false;

}


// ----------
int t_gfile::_read( char* b, int s )
{
  gtrace("t_gfile::_read");

  if (_file)
      return _file->read(b, s);

  return -1;
}
   

// ----------
int t_gfile::_write( const char* b, int s )
{
  gtrace("t_gfile::_write");

  if(_file) return  _file->write( b, s );
  
  return -1;
}


} // namespace
