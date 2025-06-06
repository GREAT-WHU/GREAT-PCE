
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <sys/stat.h>
#include <iostream>
#include <vector>

#include "gutils/gcommon.h"
#include "gutils/gfileconv.h"

using namespace std;

namespace gnut {

// basename
// ----------
string base_name(const string& path)
{
//string delimiter = "/\\";
  string base_name =  path.substr(path.find_last_of(PATH_SEPARATOR) + 1);
  return base_name;
}


// dirname
// ----------
string dir_name(const string& path)
{
//string delimiter = "/\\";

  string dir_name =     path.substr(0,    path.find_last_of(PATH_SEPARATOR));

  if( dir_name == path || dir_name.empty() ) return "";

  return dir_name;
}


// filename
// ----------
string file_name(const string& path)
{
//string delimiter = "/\\";

  string file_name = path.substr(path.find_last_of(PATH_SEPARATOR)+1,path.length()
			       - path.find_last_of(PATH_SEPARATOR));

  return file_name;
}


// dirname
// ----------
bool dir_exists( const string& path)
{
  struct stat info;

  if(stat( path.c_str(), &info ) != 0) return false;
  else if(info.st_mode & S_IFDIR)      return true;

  return false;
}
   

// create single directory
// ----------
int make_dir(const string& path)
{
  string dir  = dir_name( path );
   
  if(      dir.empty()     ) return 0;
  else if( dir_exists(dir) ) return 0;

#if defined _WIN32 || defined _WIN64
  CreateDirectory(dir.c_str(), NULL);
#else
  mkdir(dir.c_str(), 0777);
#endif

//  cerr << "creating directory: " << dir << endl;

  return 1;
}


// create directory recursively
// ----------
int make_path(const string& path)
{
  string dir = path;
  vector<string> xxx;
   
  while( (dir = dir_name( dir )) != "" ){
     xxx.push_back(dir); 
  }

  int count = 0;
  for( auto it = xxx.rbegin(); it != xxx.rend(); ++it )
  {
    count += make_dir( *it+ PATH_SEPARATOR);
  }

  return count;
}

   
} // namespace
