/*! \file pkdgrav3interface.h
 *  \brief this file contains definitions and routines for interfacing with pkdgrav3
 */

#ifndef PKDGRAV3INTERFACE_H
#define PKDGRAV3INTERFACE_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <getopt.h>
#include <sys/stat.h>

#ifdef USEMPI
#include <mpi.h>
#endif

// include the options structure via allvars
#include "allvars.h"
#include "proto.h"

///\name Include for NBodyFramework library.
//@{
#include <NBody.h>
#include <NBodyMath.h>
#include <KDTree.h>
//@}

using namespace std;
using namespace Math;
using namespace NBody;

#ifdef PKDGRAV3INTERFACE

///structure for interfacing with pkdgrav3 and extract cosmological information
namespace Pkdgrav3 {
    
    /// store all info from pkdgrav3
    struct allinfo {
      double atime;
      double littleh;
      double Omega_m;
      double Omega_r;
      double Omega_nu;
      double Omega_k;
      double Omega_b;
      double Omega_Lambda;
      double Omega_cdm;
      double w_de;
      double lengthtokpc;
      double velocitytokms;
      double masstosolarmass;
      double energyperunitmass;
      double gravity;
      double hubbleunit;
      int icosmologicalsim;
      int idarkmatter, igas, istar, ibh, iother;
      double mass_uniform_box;
    };
  
}

using namespace Pkdgrav3;

extern Options libvelociraptorOpt;

extern "C" int InitVelociraptor(const char* configname, Pkdgrav3::allinfo, const int numthreads);
extern "C" void InvokeVelociraptor(const int snapnum, char* outputname,
    Pkdgrav3::allinfo,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct pkdgrav3_vel_part *pkdgrav_parts, const int numthreads, const int ireturngroupinfoflag, 
    const int ireturnmostbound);
void SetVelociraptorSimulationState(Pkdgrav3::allinfo);
#endif

#endif
