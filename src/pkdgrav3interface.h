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

///include the options structure via allvars
#include "allvars.h"
#include "proto.h"

#ifdef PKDGRAV3INTERFACE

///structure for interfacing with pkdgrav3 and extract simulation information
namespace Pkdgrav3 {
    ///store simulation info
    struct siminfo {
        int dummy; ///< placeholder for simulation information
    };
}

#endif // PKDGRAV3INTERFACE

#endif // PKDGRAV3INTERFACE_H
