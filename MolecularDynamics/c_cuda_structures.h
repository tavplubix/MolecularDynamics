#pragma once
#include "cuda_runtime.h"





//==========================================================================
//						Structures
//==========================================================================


#define VECTOR_DIMENSION 3
struct CUDAVector
{
	double v[VECTOR_DIMENSION];
};

struct CUDAMolecule
{
	CUDAVector r, oldr;
	CUDAVector v;
	CUDAVector F, oldF, newF;
};

struct CUDAUnderspace
{
	CUDAMolecule *molecules;
	size_t numberOfMolecules;
};


#define LINEAR(p, x, y, z) z * p->Nx * p->Ny  +  y * p->Nx  +  x
#define BYTES(p) p->Nx * p->Ny * p->Nz * sizeof(CUDAUnderspace)
#define SIZE(p) p->Nx * p->Ny * p->Nz

struct CUDASpace
{
	//CUDAUnderspace *underspaces;
	CUDAUnderspace *underspaces;
	size_t Nx, Ny, Nz;
	double dt;
	int width, height;
};


//==========================================================================
//						CUDA device functions for Vector
//==========================================================================
__device__ extern void mov(const CUDAVector &a, CUDAVector& result);
__device__ extern void add(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void sub(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void mulv(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void mulc(const CUDAVector &a, const double b, CUDAVector& result);
__device__ extern void div(const CUDAVector &a, const double b, CUDAVector& result);
__device__ extern double square(const CUDAVector &a);

//==========================================================================
//						CUDA device functions for Molecule
//==========================================================================
//extern inline void Force_LennardJones(CUDAMolecule& m1, CUDAMolecule& m2);
__device__ inline void d_Force_LennardJones(const CUDAVector& r, double square, CUDAVector& F);		//r - distance between two molecules, square=r*r




//==========================================================================
//						CUDA device functions for Underspace
//==========================================================================
__device__ extern void d_recalculatePositions_Beeman(CUDAUnderspace *cus, double dt);
__device__ extern void d_recalculateSpeeds_Beeman(CUDAUnderspace *cus, double dt, int width, int height);
__device__ extern void d_calculateNewForcesForUnderspace(CUDASpace *cs, int nx, int ny, int nz);
__device__ extern void d_calculateNewForces(CUDAUnderspace *cus1, CUDAUnderspace *cus2);


//==========================================================================
//						CUDA global functions for Space
//==========================================================================
__global__ extern void cuda_recalculatePositions(CUDASpace *cs);
__global__ extern void cuda_recalculateSpeeds(CUDASpace *cs);
__global__ extern void cuda_recalculateForces(CUDASpace *cs);
__global__ extern void cuda_validate(CUDASpace *cs);


//==========================================================================
//						CUDA host functions for Space
//==========================================================================
extern void cuda_oneStep(CUDASpace *d_cs, int Nx, int Ny, int Nz);

extern CUDASpace* copyAndDeleteFromHost(CUDASpace *h_cs/*, CUDASpace **d_cs*/);
extern CUDASpace* copyAndDeleteFromDevice(CUDASpace *d_cs/*, CUDASpace *h_cs*/);

//extern void freeDeviceMem(CUDASpace *d_cs);
//extern void freeHostMem(CUDASpace *h_cs);
