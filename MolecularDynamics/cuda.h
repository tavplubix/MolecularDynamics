#pragma once
#include "Molecule.h"
#include <vector>

typedef std::vector<Molecule> MoloculesList;

extern void calculateNewForces_GPU(MoloculesList &molecules1, MoloculesList &molecules2);

extern void allocateMemory(size_t size);
extern void freeMemory();


extern CUDASpace* copyToDevice(CUDASpace *h_cs/*, CUDASpace **d_cs*/);



