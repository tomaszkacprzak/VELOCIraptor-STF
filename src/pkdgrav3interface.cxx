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



}

#endif
