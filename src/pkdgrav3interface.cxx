/*! \file pkdgrav3interface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the pkdgrav3 N-body code.
 */

#include "ioutils.h"
#include "logging.h"
#include "pkdgrav3interface.h"
#include "timer.h"

#ifdef PKDGRAV3INTERFACE

Options libvelociraptorOpt;

/// \defgroup PKDGRAV3CONFIGERRORS Errors for pkdgrav3
//@{
#define PKDGRAV3CONFIGOPTMISSING 8
#define PKDGRAV3CONFIGOPTERROR 9
#define PKDGRAV3CONFIGOPTCONFLICT 10
//@}

///check configuration options with pkdgrav3
inline int ConfigCheckPkdgrav3(Options &opt, Pkdgrav3::siminfo &s)
{
    if (opt.iBaryonSearch && !(opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK)) {
        LOG_RANK0(error) << "Conflict in config file: both gas/star/etc particle type search AND the separate baryonic (gas,star,etc) search flag are on. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (opt.iBoundHalos && opt.iKeepFOF) {
        LOG_RANK0(error) << "Conflict in config file: Asking for Bound Field objects but also asking to keep the 3DFOF/then run6DFOF. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (opt.HaloMinSize==-1) opt.HaloMinSize=opt.MinSize;

    if (s.idarkmatter == false && opt.partsearchtype==PSTDARK) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but searching only dark matter. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (s.idarkmatter == false && opt.partsearchtype==PSTALL && opt.iBaryonSearch) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but using dark matter to define links when searching all particles . This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (s.igas == false && opt.partsearchtype==PSTGAS) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only gas. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if ((s.istar == false && opt.partsearchtype==PSTSTAR)) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only stars or using stars as basis for links. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }

#ifndef USEHDF
    if (opt.ibinaryout==OUTHDF){
        LOG_RANK0(error) << "Code not compiled with HDF output enabled. Recompile with this enabled or change Binary_output";
        return PKDGRAV3CONFIGOPTERROR;
    }
#endif

#ifndef USEADIOS
    if (opt.ibinaryout==OUTADIOS){
        LOG_RANK0(error) << "Code not compiled with ADIOS output enabled. Recompile with this enabled or change Binary_output";
        return PKDGRAV3CONFIGOPTERROR;
    }
#endif

    //call the general ConfigCheck now
    ConfigCheck(opt);
    return 1;
}

int InitVelociraptor(Options &opt, const char* configname, unitinfo u, siminfo s, const int numthreads)
{
    // if mpi invoked, init the velociraptor tasks and openmp threads
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    mpi_nlocal=new Int_t[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    mpi_nhalos=new Int_t[NProcs];
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    MinNumMPI=2;
#else
    int ThisTask = 0;
#endif
#ifdef USEOPENMP
    omp_set_num_threads(numthreads);
#endif

    // logging level is reset correctly later after reading the user configuration
    vr::init_logging(vr::LogLevel::trace);
    gsl_set_error_handler_off();

    int iconfigflag;
    ///read the parameter file
    opt.pname = const_cast<char*>(configname);
    LOG_RANK0(info) << "Initialising VELOCIraptor git revision " << velociraptor::git_sha1();
    LOG_RANK0(info) << "Reading VELOCIraptor config file...";
    GetParamFile(opt);
    //on the fly finding and using pkdgrav3's mesh mpi decomposition
    opt.iontheflyfinding = true;
    opt.impiusemesh = true;
    ///check configuration
    iconfigflag = ConfigCheckPkdgrav3(opt, s);
    if (iconfigflag != 1) return iconfigflag;

    LOG_RANK0(info) << "Setting cosmology, units, sim stuff";
    ///set units, here idea is to convert internal units so that have kpc, km/s, solar mass
    ///\todo switch this so run in reasonable pkdgrav3 units and store conversion
    opt.lengthtokpc=u.lengthtokpc;
    opt.velocitytokms=u.velocitytokms;
    opt.masstosolarmass=u.masstosolarmass;

    //run in pkdgrav3 internal units, don't convert units
    opt.lengthinputconversion=1.0;
    opt.massinputconversion=1.0;
    opt.velocityinputconversion=1.0;
    opt.energyinputconversion=1.0;
    opt.SFRinputconversion=1.0;
    opt.metallicityinputconversion=1.0;
    opt.istellaragescalefactor = 1;

    //set cosmological parameters that do not change
    ///these should be in units of kpc, km/s, and solar mass
    opt.G=u.gravity;
    opt.H=u.hubbleunit;

    //set if cosmological
    opt.icosmologicalin = s.icosmologicalsim;

    //store a general mass unit, useful if running uniform box with single mass
    //and saving memory by not storing mass per particle.
#ifdef NOMASS
    opt.MassValue = s.mass_uniform_box;
    NOMASSCheck(opt);
#endif

    //write velociraptor configuration info, appending .configuration to the input config file and writing every config option
    opt.outname = configname;

    //store list of names that
    WriteVELOCIraptorConfig(opt);

#ifdef USEMPI
    //initialize the mpi write communicator to comm world;
    MPIInitWriteComm();
#endif

    LOG_RANK0(info) << "Finished initialising VELOCIraptor";

    //return the configuration flag value
    return iconfigflag;

}

int InitVelociraptor(const char* configname, unitinfo u, siminfo s, const int numthreads)
{
    return InitVelociraptor(libvelociraptorOpt, configname, u, s, numthreads);
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
