/*! \file pkdgrav3interface.cxx
 *  \brief this file contains stub routines that allow VELOCIraptor to interface with the pkdgrav3 N-body code.
 */

 #include "pkdgrav3interface.h"

 #ifdef PKDGRAV3INTERFACE
 
 Options* pkdgrav3_make_default_options() {
     // Safe for both aggregates and ctor-based types
     auto* opt = new Options{};      // value-initialize
     return opt;
 }
 
 void pkdgrav3_destroy_options(Options* opt) {
     delete opt; // legal here: we see the complete type
 }
 
 void pkdgrav3_load_options(const char* filename, Options &opt) {
     
     opt.pname = const_cast<char*>(filename); 
     fprintf(stdout, "Velociraptor Options->pname: %s\n", opt.pname);
 }
 
 #endif
 