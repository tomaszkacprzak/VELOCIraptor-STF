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
    /// store basic cosmological information
    struct cosmoinfo {
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
    };
    ///structure to store unit information
    struct unitinfo {
        double lengthtokpc,velocitytokms,masstosolarmass,energyperunitmass,gravity,hubbleunit;
    };

    struct cell_loc {
        double loc[3];
    };

    ///store simulation info
    struct siminfo {
        double period, zoomhigresolutionmass, interparticlespacing, spacedimension[3];
        int numcells;
        int numcellsperdim;
        cell_loc *cellloc;
        double cellwidth[3];
        double icellwidth[3];
        int *cellnodeids;
        int icosmologicalsim;
        int izoomsim;
        int idarkmatter, igas, istar, ibh, iother;
#ifdef NOMASS
        double mass_uniform_box;
#endif
    };

    struct groupinfo {
      int index;
      long long groupid;
    };

    struct vr_return_data {
      int num_gparts_in_groups;
      struct groupinfo *group_info;
      int num_most_bound;
      int *most_bound_index;
    };
}

struct pkdgrav3_vel_part;

using namespace Pkdgrav3;

extern Options libvelociraptorOpt;
int InitVelociraptor(Options &opt, const char* configname, unitinfo u, siminfo s, const int numthreads);
extern "C" int InitVelociraptor(const char* configname, unitinfo u, siminfo s, const int numthreads);
extern "C" vr_return_data InvokeVelociraptor(const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct pkdgrav3_vel_part *pkdgrav_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, const int ireturnmostbound);
void SetVelociraptorSimulationState(cosmoinfo c, siminfo s);
#endif

#endif
