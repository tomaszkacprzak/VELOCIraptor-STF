/*! \file pkdgrav3interface.cxx
 *  \brief this file contains stub routines that allow VELOCIraptor to interface with the pkdgrav3 N-body code.
 */

 #include "pkdgrav3interface.h"

 #ifdef PKDGRAV3INTERFACE
 
 Options* make_default_options() {
     // Safe for both aggregates and ctor-based types
     auto* o = new Options{};      // value-initialize
     return o;
 }
 
 void destroy_options(Options* p) {
     delete p; // legal here: we see the complete type
 }
 
 #endif
 