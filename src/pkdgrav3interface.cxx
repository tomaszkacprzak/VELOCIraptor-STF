/*! \file pkdgrav3interface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the pkdgrav3 N-body code.
 */

#include "ioutils.h"
#include "logging.h"
#include "pkdgrav3interface.h"
#include "timer.h"

#ifdef PKDGRAV3INTERFACE

Options libvelociraptorOpt;
Options libvelociraptorOptbackup;

/// \defgroup PKDGRAV3CONFIGERRORS Errors for pkdgrav3
//@{
#define PKDGRAV3CONFIGOPTMISSING 8
#define PKDGRAV3CONFIGOPTERROR 9
#define PKDGRAV3CONFIGOPTCONFLICT 10
//@}

///check configuration options with pkdgrav3
inline int ConfigCheckPkdgrav3(Options &opt, Pkdgrav3::siminfo &s)
{
    std::cout << "hello world from config ConfigCheckPkdgrav3" << std::endl;
    return 1;
}

#endif // PKDGRAV3INTERFACE
