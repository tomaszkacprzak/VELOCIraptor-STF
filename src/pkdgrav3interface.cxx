/*! \file pkdgrav3interface.cxx
 *  \brief this file contains stub routines that allow VELOCIraptor to interface with the pkdgrav3 N-body code.
 */

#include "pkdgrav3interface.h"
#include "ioutils.h"
#include "logging.h"
#include "timer.h"
#include <cstdint>

#ifdef PKDGRAV3INTERFACE

/// \defgroup PKDGRAV3CONFIGERRORS Errors for pkdgrav3
//@{
///
#define PKDGRAV3CONFIGOPTMISSING 8
#define PKDGRAV3CONFIGOPTERROR 9
#define PKDGRAV3CONFIGOPTCONFLICT 10
//@}

///check configuration options with swift
inline int ConfigCheckPkdgrav3(Options &opt, int idarkmatter, int igas, int istar)
{
    if (opt.iBaryonSearch && !(opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK)) {
        LOG_RANK0(error) << "Conflict in config file: both gas/star/etc particle type search AND the separate baryonic (gas,star,etc) search flag are on. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (opt.iBoundHalos && opt.iKeepFOF) {
        LOG_RANK0(error) << "Conflict in config file: Asking for Bound Field objects but also asking to keep the 3DFOF/then run 6DFOF. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (opt.HaloMinSize==-1) opt.HaloMinSize=opt.MinSize;

    if (idarkmatter == false && opt.partsearchtype==PSTDARK) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but searching only dark matter. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (idarkmatter == false && opt.partsearchtype==PSTALL && opt.iBaryonSearch) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but using dark matter to define links when searching all particles . This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if (igas == false && opt.partsearchtype==PSTGAS) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only gas. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }
    if ((istar == false && opt.partsearchtype==PSTSTAR)) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only stars or using stars as basis for links. This is incompatible. Check config";
        return PKDGRAV3CONFIGOPTCONFLICT;
    }


    //call the general ConfigCheck now
    ConfigCheck(opt);
    return 1;
}

Options* Pkdgrav3MakeDefaultOptions() {
    // Safe for both aggregates and ctor-based types
    auto* opt = new Options{};      // value-initialize
    return opt;
}

void Pkdgrav3DestroyOptions(Options* opt) {
    delete opt; // legal here: we see the complete type
}

int Pkdgrav3LoadOptions(const char* filename_options, const char* filename_output, Options &opt, const int numthreads, const double box_size, const int num_total_particles) {

    opt.pname = const_cast<char*>(filename_options); 
    opt.outname = const_cast<char*>(filename_output);
    fprintf(stdout, "Velociraptor Options->pname: %s  numthreads: %d\n", opt.pname, numthreads);

    // following swiftinterface.cxx

    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //mpi_domain=new MPI_Domain[NProcs];
    mpi_nlocal=new Int_t[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    mpi_nhalos=new Int_t[NProcs];
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    //store MinSize as when using mpi prior to stitching use min of 2;
    MinNumMPI=2;

    omp_set_num_threads(numthreads);


    // logging level is reset correctly later after reading the user configuration
    vr::init_logging(vr::LogLevel::trace);
    gsl_set_error_handler_off();

    int iconfigflag;

    GetParamFile(opt);

    fprintf(stdout, "Velociraptor read config file\n");
    fprintf(stdout, "Velociraptor opt.MinSize: %d\n", opt.MinSize);

    //on the fly finding and using swift's mesh mpi decomposition
    opt.iontheflyfinding = true;
    opt.impiusemesh = false;

    // for now only dark matter is supported
    int igas=0;
    int istar=0;
    int idarkmatter=1;
    iconfigflag = ConfigCheckPkdgrav3(opt, idarkmatter, igas, istar);
    if (iconfigflag != 1) return iconfigflag;

    LOG_RANK0(info) << "Setting cosmological parameters";

    opt.icosmologicalin = 1;
    opt.lengthtokpc=1000.0; // output in Gpc
    opt.velocitytokms=1.0;  // output in km/s
    opt.masstosolarmass=1.0e10; // output in 1e10 Msun
    opt.icomoveunit=1;

    // unit conversion
    // https://pkdgrav3.readthedocs.io/en/latest/units.html
    opt.lengthinputconversion=box_size;
    opt.massinputconversion=pow(box_size, 3) * 2.77536627208 * 1.0e11 / (double)num_total_particles; // in pkdgrav3, "mass" listed in the units docs is total mass of the box, so we need to divide by the number of particles
    opt.velocityinputconversion=box_size * sqrt(8.0/3.0*M_PI) / (box_size * 100.0); // needs to be divided by scale factor a later
    
    // cosmological parameters
    opt.G=1.0;
    opt.H=1.0;
    opt.h=1.0;
    opt.Omega_m=1.0;
    opt.Omega_Lambda=1.0;
    opt.Omega_de=1.0;
    opt.Omega_k=1.0;
    opt.Omega_cdm=1.0;
    opt.Omega_b=1.0;
    opt.Omega_r=1.0;
    opt.Omega_nu=1.0;
    opt.w_de=-1.0;
    opt.MassValue = 1.0;
    NOMASSCheck(opt);
    
    // WriteVELOCIraptorConfig(opt);

// #ifdef USEMPI
//     //initialize the mpi write communicator to comm world;
//     MPIInitWriteComm();
// #endif

//     LOG_RANK0(info) << "Finished initialising VELOCIraptor";

    return iconfigflag;


}

int Pkdgrav3InvokeVelociraptor(const char* filename_options, const char* filename_output, const int iStep, const double box_size, const long num_total_particles, std::vector<NBody::Particle> &vExportParticles, unsigned int &nExportParticleCount, const int numthreads) {

    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int sum_rank = 0;
    MPI_Allreduce(&rank, &sum_rank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   
    fprintf(stdout, "Velociraptor rank %d sum_rank %d invoked with options: options %s output %s step %d numthreads %d box_size %f nExportParticleCount %d num_total_particles %ld\n", rank, sum_rank, filename_options, filename_output, iStep, numthreads, box_size, nExportParticleCount, num_total_particles);

    /*
    Load Velociraptoroptions 
    */
    Options* opt;
    opt = Pkdgrav3MakeDefaultOptions();
    Pkdgrav3LoadOptions(filename_options, filename_output, *opt, numthreads, box_size, num_total_particles);

    /*
    Main Velociraptor interface code from swiftinterface.cxx InvokeVelociraptorHydro()
    */

    int NProcs=1;
    int ThisTask=0;
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    omp_set_num_threads(numthreads);

    Int_t Nlocal, Ntotal, Nmemlocal;
    Int_t *pfof = NULL, *pfofall = NULL, *pfofbaryons = NULL, *numingroup = NULL, **pglist = NULL;
    Int_t *nsub = NULL, *parentgid = NULL, *uparentgid =NULL, *stype = NULL;
    Int_t ndark, index;
    Int_t ngroup, nhalos;
    //to store information about the group
    PropData *pdata = NULL,*pdatahalos = NULL;
    Nlocal = Nmemlocal = nExportParticleCount; // to fit the name convention
    Ntotal = num_total_particles; // to fit the name convention

    // Add memory padding for MPI
    if(NProcs>1) {
        Nmemlocal*=(1+opt->mpipartfac);
        vExportParticles.resize(Nmemlocal);
    }
    
    fprintf(stdout, "Velociraptor rank %d NProcs %d ThisTask %d Nlocal %d Ntotal %d Nmemlocal %d mpipartfac %f\n", rank, NProcs, ThisTask, Nlocal, Ntotal, Nmemlocal, opt->mpipartfac);

    // Exchange particle counts
    MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgather(&Nlocal, 1, MPI_Int_t, mpi_nlocal, 1, MPI_Int_t, MPI_COMM_WORLD);



    /*
    Cleanup
    */
    Pkdgrav3DestroyOptions(opt);

    return 0;
}

#endif
