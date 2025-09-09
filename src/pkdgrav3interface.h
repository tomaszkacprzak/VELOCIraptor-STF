/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with pkdgrav3
 */

#ifndef PKDGRAV3INTERFACE_H
#define PKDGRAV3INTERFACE_H

#include <iostream>

#ifdef PKDGRAV3INTERFACE

extern "C" int InitVelociraptor(char *configname, const int numthreads);
extern "C" int InitVelociraptorExtra(const int iextra, char *configname, const int numthreads);
extern "C" int InvokeVelociraptor(const int snapnum, char *outputname, const int numthreads);
extern "C" int InvokeVelociraptorExtra(const int iextra, const int snapnum, char *outputname, const int numthreads);
void SetVelociraptorSimulationState();

#endif

#endif
