/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with the pkdgrav3 N-body code.
 */

 #ifndef PKDGRAV3INTERFACE_H
 #define PKDGRAV3INTERFACE_H
 
 #include "allvars.h"
 #include "proto.h"
 
 #ifdef PKDGRAV3INTERFACE
 
 Options* make_default_options();
 void destroy_options(Options* p);
 
 #endif
 
 #endif
 