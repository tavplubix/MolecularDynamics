#pragma once
#include "cuda_runtime.h"


#define ODINARY_PRECISION

#ifdef ODINARY_PRECISION
typedef float myfloat;
#else
typedef double myfloat;
#endif




#define LINEAR(p, x, y, z) z * p->Nx * p->Ny  +  y * p->Nx  +  x
#define BYTES(p) p->Nx * p->Ny * p->Nz * sizeof(CUDAUnderspace)
#define SIZE(p) p->Nx * p->Ny * p->Nz

#define GET_POINTER(type, ptr, shift) reinterpret_cast<type*>( reinterpret_cast<size_t>(ptr) + shift )
#define GET_SHIFT(ptr1, ptr2) (reinterpret_cast<size_t>(ptr2) - reinterpret_cast<size_t>(ptr1))
#define WHOLE_SIZE_OF_SPACE(space)  1*sizeof(CUDASpace) + space->Nx*space->Ny*space->Nz*sizeof(CUDAUnderspace) + space->numberOfAllMolecules*sizeof(CUDAMolecule);

//==========================================================================
//						Structures
//==========================================================================

typedef unsigned char byte;		//sizeof(byte) == 1
//5.3.3 - ... sizeof(char), sizeof(signed char) and sizeof(unsigned char) are 1 ...



#define VECTOR_DIMENSION 3
struct CUDAVector
{
	myfloat v[VECTOR_DIMENSION];
};

struct CUDAMolecule
{
	CUDAVector r, oldr;
	CUDAVector v;
	CUDAVector F, oldF, newF;
	short type;
	int id;
};

struct CUDAUnderspace
{
	//CUDAMolecule *molecules;
	size_t moleculesShift;
	size_t numberOfMolecules;
};



struct CUDASpace
{
	//CUDAUnderspace ***underspaces;
	//CUDAUnderspace *underspaces;
	size_t underspacesShift;
	size_t Nx, Ny, Nz;
	myfloat dt;
	int width, height, depth;
	size_t numberOfAllMolecules;
};




//==========================================================================
//						CUDA device functions for Vector
//==========================================================================
__device__ extern void mov(const CUDAVector &a, CUDAVector& result);
__device__ extern void add(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void sub(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void mul(const CUDAVector &a, const CUDAVector& b, CUDAVector& result);
__device__ extern void mul(const CUDAVector &a, const myfloat b, CUDAVector& result);
__device__ extern void div(const CUDAVector &a, const myfloat b, CUDAVector& result);
__device__ extern double square(const CUDAVector &a);

//==========================================================================
//						CUDA device functions for Molecule
//==========================================================================
//extern inline void Force_LennardJones(CUDAMolecule& m1, CUDAMolecule& m2);
__device__ inline void d_Force_LennardJones(const CUDAVector& r, myfloat square, CUDAVector& F);		//r - distance between two molecules, square=r*r




//==========================================================================
//						CUDA device functions for Underspace
//==========================================================================
__device__ extern void d_recalculatePositions_Beeman(CUDAUnderspace *cus, myfloat dt);
__device__ extern void d_recalculateSpeeds_Beeman(CUDAUnderspace *cus, myfloat dt, int width, int height, int depth);
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

extern CUDASpace* moveFromHost(CUDASpace *h_cs, size_t wholeSize = 0);
extern CUDASpace* moveFromDevice(CUDASpace *d_cs, size_t wholeSize = 0, byte *h_p = nullptr);

//extern void freeDeviceMem(CUDASpace *d_cs);
//extern void freeHostMem(CUDASpace *h_cs);
