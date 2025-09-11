/*! \file pkdgrav3interface.cxx
 *  \brief this file contains stub routines that allow VELOCIraptor to interface with the pkdgrav3 N-body code.
 */

#include "pkdgrav3interface.h"
#include "ioutils.h"
#include "logging.h"
#include "timer.h"

#ifdef PKDGRAV3INTERFACE



Options* pkdgrav3_make_default_options() {
    // Safe for both aggregates and ctor-based types
    auto* opt = new Options{};      // value-initialize
    return opt;
}

void pkdgrav3_destroy_options(Options* opt) {
    delete opt; // legal here: we see the complete type
}

void pkdgrav3_load_options(const char* filename, Options &opt, const int numthreads) {

    opt.pname = const_cast<char*>(filename); 
    opt.outname = strdup("test.vr");
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

//     //on the fly finding and using swift's mesh mpi decomposition
//     opt.iontheflyfinding = true;
//     opt.impiusemesh = true;

//     LOG_RANK0(info) << "Setting cosmology, units, sim stuff";
//     ///set units, here idea is to convert internal units so that have kpc, km/s, solar mass
//     ///\todo switch this so run in reasonable swift units and store conversion
//     opt.lengthtokpc=u.lengthtokpc;
//     opt.velocitytokms=u.velocitytokms;
//     opt.masstosolarmass=u.masstosolarmass;

//     //run in swift internal units, don't convert units
//     opt.lengthinputconversion=1.0;
//     opt.massinputconversion=1.0;
//     opt.velocityinputconversion=1.0;
//     opt.energyinputconversion=1.0;
//     opt.SFRinputconversion=1.0;
//     opt.metallicityinputconversion=1.0;
//     opt.istellaragescalefactor = 1;

//     //set cosmological parameters that do not change
//     ///these should be in units of kpc, km/s, and solar mass
//     opt.G=u.gravity;
//     opt.H=u.hubbleunit;

//     //set if cosmological
//     opt.icosmologicalin = s.icosmologicalsim;

//     //store a general mass unit, useful if running uniform box with single mass
//     //and saving memory by not storing mass per particle.
// #ifdef NOMASS
//     opt.MassValue = s.mass_uniform_box;
//     NOMASSCheck(opt);
// #endif

//     //write velociraptor configuration info, appending .configuration to the input config file and writing every config option
//     opt.outname = configname;

//     //store list of names that
//     WriteVELOCIraptorConfig(opt);

// #ifdef USEMPI
//     //initialize the mpi write communicator to comm world;
//     MPIInitWriteComm();
// #endif

//     LOG_RANK0(info) << "Finished initialising VELOCIraptor";


}

#endif
