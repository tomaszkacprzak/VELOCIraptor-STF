/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with the pkdgrav3 N-body code.
 */

#ifndef PKDGRAV3INTERFACE_H
#define PKDGRAV3INTERFACE_H

#include "allvars.h"
#include "proto.h"

#ifdef PKDGRAV3INTERFACE

Options* Pkdgrav3MakeDefaultOptions();
void Pkdgrav3DestroyOptions(Options* opt);
int Pkdgrav3LoadOptions(const char* filename_options, const char* filename_output, Options &opt, const int numthreads, const double box_size, const int num_total_particles);
int Pkdgrav3InvokeVelociraptor(Options* opt);

#endif

#endif
