/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with the pkdgrav3 N-body code.
 */

#ifndef PKDGRAV3INTERFACE_H
#define PKDGRAV3INTERFACE_H

#include "allvars.h"
#include "proto.h"

#ifdef PKDGRAV3INTERFACE

Options* pkdgrav3_make_default_options();
void pkdgrav3_destroy_options(Options* opt);
void pkdgrav3_load_options(const char* filename, Options &opt, const int numthreads);

#endif

#endif
