#include "pkdgrav3interface.h"

#ifdef PKDGRAV3INTERFACE

extern "C" int InitVelociraptor(char *configname, const int numthreads) {
    std::cout << "hello world from velociraptor: InitVelociraptor" << std::endl;
    return 0;
}

extern "C" int InitVelociraptorExtra(const int iextra, char *configname, const int numthreads) {
    std::cout << "hello world from velociraptor: InitVelociraptorExtra" << std::endl;
    return 0;
}

extern "C" int InvokeVelociraptor(const int snapnum, char *outputname, const int numthreads) {
    std::cout << "hello world from velociraptor: InvokeVelociraptor" << std::endl;
    return 0;
}

extern "C" int InvokeVelociraptorExtra(const int iextra, const int snapnum, char *outputname, const int numthreads) {
    std::cout << "hello world from velociraptor: InvokeVelociraptorExtra" << std::endl;
    return 0;
}

void SetVelociraptorSimulationState() {
    std::cout << "hello world from velociraptor: SetVelociraptorSimulationState" << std::endl;
}

#endif
