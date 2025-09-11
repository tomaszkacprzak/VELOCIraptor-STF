/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with the pkdgrav3 N-body code.
 */

#ifndef PKDGRAV3INTERFACE_H
#define PKDGRAV3INTERFACE_H

#include "allvars.h"
#include "proto.h"

#ifdef PKDGRAV3INTERFACE
/// Initialise VELOCIraptor when called from pkdgrav3.
extern "C" int InitVelociraptorPKD();

/// Invoke VELOCIraptor from pkdgrav3.
extern "C" int InvokeVelociraptorPKD();
#endif

#endif
