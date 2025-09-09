/*! \file pkdgrav3interface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the pkdgrav3 N-body code.
 */

#include "pkdgrav3interface.h"

#ifdef PKDGRAV3INTERFACE

int InitVelociraptor(char* configname, unitinfo u, siminfo s, const int numthreads)
{
    std::cout << "hello world from velociraptor: InitVelociraptor" << std::endl;
    return 0;
}

vr_return_data InvokeVelociraptor(const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct pkdgrav3_vel_part *pkdgrav_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, const int ireturnmostbound)
{
    std::cout << "hello world from velociraptor: InvokeVelociraptor" << std::endl;
    vr_return_data data{};
    return data;
}

void SetVelociraptorSimulationState(cosmoinfo c, siminfo s)
{
    std::cout << "hello world from velociraptor: SetVelociraptorSimulationState" << std::endl;
}

#endif
