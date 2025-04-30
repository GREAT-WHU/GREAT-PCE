/**
*
* @verbatim
    History
    2011-01-10  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gfile.h
* @brief      Purpose: implements file io based on gio class
*.
* @author     JD
* @version    1.0.0
* @date       2011-01-10
*
*/

#ifndef GFILE_H
#define GFILE_H


#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/gio.h"
#include "gio/giof.h"

// special buffer size for file reading
// --> must be bellow gcoder maximum limit !
#define FILEBUF_SIZE 20480
#define FILEHDR_SIZE    48

using namespace std;

namespace gnut {
    /** @brief class for t_gfile derive from t_gio. */
    class LibGnut_LIBRARY_EXPORT t_gfile : public t_gio {

     public:
         /** @brief default constructor. */
         t_gfile();
         virtual ~t_gfile();
         
         virtual int init_write();
         virtual int init_read();
         
         virtual int irc()const;                   // get irc status
         virtual bool eof();                       // integrate gzip/ascii
         virtual string md5sum();                  // get md5sum
         virtual string mask();                    // integrate gzip/ascii
         
         virtual bool compressed(){ return _gzip; }
         
         virtual int path( string str );           // set/get file://dir/name
         virtual string path();
         virtual string name();
         virtual void reset();                     // reset path

     protected:
         /** @brief write data.
        * @param[in]    buff    buffer of the data
        * @param[in]    size    buffer size of the data
        * @return
            @retval >0    number of bytes read
            @retval <=0    fail
        */
         virtual int _gio_write( const char* buff, int size );
         /** @brief read data.
        * @param[in]    buff    buffer of the data
        * @param[in]    size    buffer size of the data
        * @return
            @retval >0    number of bytes read
            @retval <=0    fail
        */
         virtual int _gio_read(  char* buff, int size );
         /**
        * @brief common function for file close.
        * @return        running status
        */
         virtual int _stop_common();
         /**
        * @brief compressed or not ?
        * @param[in]    name    file name
        */
         virtual void _set_gzip( string name );
         /**
         * @brief read.
         * @param[in]    b    buffer of the data
         * @param[in]    s    buffer size of the data
         * @return
             @retval <0    fail
         */
         virtual int _read( char* b, int s);
         /**
         * @brief write.
         * @param[in]    b    buffer of the data
         * @param[in]    s    buffer size of the data
         * @return
             @retval <0    fail
         */
         virtual int _write( const char* b, int s);

         int         _irc;
         bool        _gzip;                        // compressed
         t_giof*     _file;                        // ascii file
     private:
    };

} // namespace
 
#endif
