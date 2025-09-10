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
inline int ConfigCheckPkdgrav3(Options &opt, Pkdgrav3::allinfo &pkdgrav3_info)
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

    if (pkdgrav3_info.idarkmatter == false && opt.partsearchtype==PSTDARK) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but searching only dark matter. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (pkdgrav3_info.idarkmatter == false && opt.partsearchtype==PSTALL && opt.iBaryonSearch) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but using dark matter to define links when searching all particles . This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (pkdgrav3_info.igas == false && opt.partsearchtype==PSTGAS) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only gas. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if ((pkdgrav3_info.istar == false && opt.partsearchtype==PSTSTAR)) {
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

int InitVelociraptor(const char* configname, allinfo pkdgrav3_info, const int numthreads)
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

    // Create options object
    Options opt;

    int iconfigflag;
    ///read the parameter file
    opt.pname = const_cast<char*>(configname);
    LOG_RANK0(info) << "Initialising VELOCIraptor git revision " << velociraptor::git_sha1();
    LOG_RANK0(info) << "Reading VELOCIraptor config file...       " << opt.pname;
    GetParamFile(opt);
    //on the fly finding and using pkdgrav3's mesh mpi decomposition
    opt.iontheflyfinding = true;
    opt.impiusemesh = true;
    ///check configuration
    iconfigflag = ConfigCheckPkdgrav3(opt, pkdgrav3_info);
    if (iconfigflag != 1) return iconfigflag;

    LOG_RANK0(info) << "Setting cosmology, units, sim stuff";
    ///set units, here idea is to convert internal units so that have kpc, km/s, solar mass
    ///\todo switch this so run in reasonable pkdgrav3 units and store conversion
    opt.lengthtokpc=pkdgrav3_info.lengthtokpc;
    opt.velocitytokms=pkdgrav3_info.velocitytokms;
    opt.masstosolarmass=pkdgrav3_info.masstosolarmass;

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
    opt.G=pkdgrav3_info.gravity;
    opt.H=pkdgrav3_info.hubbleunit;

    //set if cosmological
    opt.icosmologicalin = pkdgrav3_info.icosmologicalsim;

    //store a general mass unit, useful if running uniform box with single mass
    //and saving memory by not storing mass per particle.
#ifdef NOMASS
    opt.MassValue = pkdgrav3_info.mass_uniform_box;
    NOMASSCheck(opt);
#endif

    //write velociraptor configuration info, appending .configuration to the input config file and writing every config option
    opt.outname = const_cast<char*>(configname);

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


void InvokeVelociraptor(const int snapnum, char* outputname,
    allinfo pkdgrav3_info,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct pkdgrav3_vel_part *pkdgrav_parts,
    const int numthreads, const int ireturngroupinfoflag, const int ireturnmostbound)
{
    std::cout << "hello world from velociraptor: InvokeVelociraptor" << std::endl;
}

void SetVelociraptorSimulationState(allinfo pkdgrav3_info)
{
    std::cout << "hello world from velociraptor: SetVelociraptorSimulationState" << std::endl;
}

#endif
