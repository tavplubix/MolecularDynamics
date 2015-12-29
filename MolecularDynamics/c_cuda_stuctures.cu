#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <iostream>

#include "c_cuda_structures.h"


//==========================================================================
//					Operations with CUDAVector
//==========================================================================

__device__ void mov(const CUDAVector &a, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i];
}

__device__ void add(const CUDAVector &a, const CUDAVector& b, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i] + b.v[i];
}

__device__ void sub(const CUDAVector &a, const CUDAVector& b, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i] - b.v[i];
}

__device__ void mulv(const CUDAVector &a, const CUDAVector& b, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i] * b.v[i];
}

__device__ void mulc(const CUDAVector &a, const double b, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i] * b;
}

__device__ void div(const CUDAVector &a, const double b, CUDAVector& result)
{
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result.v[i] = a.v[i] / b;
}

__device__ double square(const CUDAVector &a)
{
	double result = 0;
	for (size_t i = 0; i < VECTOR_DIMENSION; ++i)
		result += a.v[i] * a.v[i];
	return result;
}


//==========================================================================
//					Operations with CUDAMolecule
//==========================================================================

#define Boltzmann  1.3806488e-23
#define Angstrom  10e-10
#define AtomicMassUnit  1.660538921e-27
#define sigma 2.74*10e-10	//m
#define epsilon 36.2*1.3806488e-23	//J
#define mass 20.1797 * 1.660538921e-27

#define maxDistSquare 4.0 * sigma * 4.0 * sigma


__device__ inline void d_Force_LennardJones(const CUDAVector& r, double square, CUDAVector& F)
{
	const double sigmaSquare = sigma * sigma;
	square = sigmaSquare / square;
	double U = 2.0*pow(square, 14 / 2) - pow(square, 8 / 2);
	const double c = -24.0 * epsilon / sigmaSquare;
	mulc(r, c * U, F);
}


//==========================================================================
//					Operations with CUDAUnderspace
//==========================================================================

__device__ void d_recalculatePositions_Beeman(CUDAUnderspace *cus, double dt)
{
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		mov(cus->molecules[i].r, cus->molecules[i].oldr);
		CUDAVector tmp;
		//i.r += i.v * dt;
		mulc(cus->molecules[i].v, dt, tmp);
		add(cus->molecules[i].r, tmp, cus->molecules[i].r);
		//i.r += 4.0 / 6.0 * (i.F / i.m) * (dt*dt);
		mulc(cus->molecules[i].F, (4.0 / 6.0) * (dt*dt) / mass, tmp);
		add(cus->molecules[i].r, tmp, cus->molecules[i].r);
		//i.r -= 1.0 / 6.0 * (i.oldF / i.m) * (dt*dt);
		mulc(cus->molecules[i].oldF, - (1.0/6.0) * (dt*dt) / mass, tmp);
		add(cus->molecules[i].r, tmp, cus->molecules[i].r);
	}
}

__device__ void d_recalculateSpeeds_Beeman(CUDAUnderspace *cus, double dt, int width, int height)
{
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		CUDAVector tmp;
		//i.v += 2.0 / 6.0 * (i.newF / i.m) * dt;
		mulc(cus->molecules[i].newF, (2.0 / 6.0) * dt / mass, tmp);
		add(cus->molecules[i].v, tmp, cus->molecules[i].v);
		//i.v += 5.0 / 6.0 * (i.F / i.m) * dt;
		mulc(cus->molecules[i].F, (5.0 / 6.0) * dt / mass, tmp);
		add(cus->molecules[i].v, tmp, cus->molecules[i].v);
		//i.v -= 1.0 / 6.0 * (i.oldF / i.m) * dt;
		mulc(cus->molecules[i].oldF, - (1.0 / 6.0) * dt / mass, tmp);
		add(cus->molecules[i].v, tmp, cus->molecules[i].v);

		if (cus->molecules[i].r.v[0] <= 0) {
			cus->molecules[i].v.v[0] = abs(cus->molecules[i].v.v[0]);
		}
		if (width * Angstrom <= cus->molecules[i].r.v[0]) {
			cus->molecules[i].v.v[0] = -abs(cus->molecules[i].v.v[0]);
		}

		if (cus->molecules[i].r.v[1] <= 0) {
			cus->molecules[i].v.v[1] = abs(cus->molecules[i].v.v[1]);
		}
		if (height * Angstrom <= cus->molecules[i].r.v[1]) {
			cus->molecules[i].v.v[1] = -abs(cus->molecules[i].v.v[1]);
		}
	}
}

__device__ void d_calculateNewForcesForUnderspace(CUDASpace *cs, int nx, int ny, int nz)
{
	CUDAUnderspace *centralSpace = &cs->underspaces[nx][ny][nz];
	int  closest = 1;
	for (int dx = -closest; dx <= closest; ++dx) {
		for (int dy = -closest; dy <= closest; ++dy) {
			for (int dz = -closest; dz <= closest; ++dz) {
				int x = nx + dx;
				int y = ny + dy;
				int z = nz + dz;
				if (x < 0 || y < 0 || z < 0) continue;
				if (x >= cs->Nx) continue;
				if (y >= cs->Ny) continue;
				if (z >= cs->Nz) continue;
				d_calculateNewForces(centralSpace, &cs->underspaces[x][y][z]);
			}
		}
	}
}

__device__ void d_calculateNewForces(CUDAUnderspace *cus1, CUDAUnderspace *cus2)
{
	for (size_t i = 0; i < cus1->numberOfMolecules; ++i) {
		for (size_t j = 0; j < cus2->numberOfMolecules; ++j) {
			if (i == j) continue;
			//Vector r = (*j).r - (*i).r;
			CUDAVector tmp;
			sub(cus2->molecules[j].r, cus1->molecules[i].r, tmp);
			double sq = square(tmp);
			if (maxDistSquare < sq) continue;
			//(*i).newF += d_Force_LennardJones(r, sq);
			d_Force_LennardJones(tmp, sq, tmp);
			add(cus1->molecules[i].newF, tmp, cus1->molecules[i].newF);
		}
	}
}


//==========================================================================
//					Operations with CUDASpace
//==========================================================================
__global__ void cuda_recalculatePositions(CUDASpace *cs)
{
	size_t nx = blockDim.x * blockIdx.x + threadIdx.x;
	if (nx >= cs->Nx) return;
	size_t ny = blockDim.y * blockIdx.y + threadIdx.y;
	if (ny >= cs->Ny) return;
	size_t nz = blockDim.z * blockIdx.z + threadIdx.z;
	if (ny >= cs->Ny) return;

	d_recalculatePositions_Beeman(&cs->underspaces[nx][ny][nz], cs->dt);
}

__global__ void cuda_recalculateSpeeds(CUDASpace *cs)
{
	size_t nx = blockDim.x * blockIdx.x + threadIdx.x;
	if (nx >= cs->Nx) return;
	size_t ny = blockDim.y * blockIdx.y + threadIdx.y;
	if (ny >= cs->Ny) return;
	size_t nz = blockDim.z * blockIdx.z + threadIdx.z;
	if (ny >= cs->Ny) return;

	d_recalculateSpeeds_Beeman(&cs->underspaces[nx][ny][nz], cs->dt, cs->width, cs->height);
}

__global__ void cuda_recalculateForces(CUDASpace *cs)
{
	size_t nx = blockDim.x * blockIdx.x + threadIdx.x;
	//if (nx >= cs->Nx) return;	
	size_t ny = blockDim.y * blockIdx.y + threadIdx.y;
	//if (ny >= cs->Ny) return;
	size_t nz = blockDim.z * blockIdx.z + threadIdx.z;
	//if (ny >= cs->Ny) return;

	d_calculateNewForcesForUnderspace(cs, nx, ny, nz);
}

__global__ void cuda_validate(CUDASpace *cs)
{
	size_t nx = blockDim.x * blockIdx.x + threadIdx.x;
	if (nx >= cs->Nx) return;
	size_t ny = blockDim.y * blockIdx.y + threadIdx.y;
	if (ny >= cs->Ny) return;
	size_t nz = blockDim.z * blockIdx.z + threadIdx.z;
	if (ny >= cs->Ny) return;

	CUDAUnderspace& cus = cs->underspaces[nx][ny][nz];
	for (size_t i = 0; i < cus.numberOfMolecules; ++i) {
		//t.oldF = t.F;
		mov(cus.molecules[i].oldF, cus.molecules[i].F);
		//t.F = t.newF;
		mov(cus.molecules[i].F, cus.molecules[i].newF);
	}
}

void cuda_oneStep(CUDASpace *d_cs,int Nx, int Ny, int Nz)
{
	int numberOfCores = 512;
	int coresPerDim = numberOfCores / 3;
	dim3 grid, blocks;
	grid = dim3(Nx / coresPerDim + 1, Ny / coresPerDim + 1, Nz / coresPerDim + 1);
	blocks = dim3(coresPerDim, coresPerDim, coresPerDim);

	cuda_recalculatePositions	<<<grid, blocks>>> (d_cs);
	cudaDeviceSynchronize();
	cuda_recalculateForces		<<<grid, blocks>>> (d_cs);
	cudaDeviceSynchronize();
	cuda_recalculateSpeeds		<<<grid, blocks>>> (d_cs);
	cudaDeviceSynchronize();
	cuda_validate				<<<grid, blocks>>> (d_cs);
	cudaDeviceSynchronize();
}



extern "C" CUDASpace* /*void*/ copyToDevice(CUDASpace *h_cs/*, CUDASpace **d_cs*/)
{
	CUDASpace *cs;
	cudaMalloc(&cs, sizeof(CUDASpace));
	cudaMalloc(&cs->underspaces, sizeof(CUDAUnderspace**) * h_cs->Nx);
	for (size_t i = 0; i < h_cs->Nx; ++i) {
		cudaMalloc(&cs->underspaces[i], sizeof(CUDAUnderspace*) * h_cs->Ny);
		for (size_t j = 0; j < h_cs->Ny; ++j) {
			cudaMalloc(&cs->underspaces[i][j], sizeof(CUDAUnderspace) * h_cs->Nz);
			for (size_t k = 0; k < h_cs->Nz; ++k) {
				cudaMalloc(&cs->underspaces[i][j][k].molecules, sizeof(CUDAMolecule) * h_cs->underspaces[i][j][k].numberOfMolecules);
				cudaMemcpy(cs->underspaces[i][j][k].molecules, h_cs->underspaces[i][j][k].molecules, h_cs->underspaces[i][j][k].numberOfMolecules, cudaMemcpyHostToDevice);
				cudaMemcpy(&cs->underspaces[i][j][k].numberOfMolecules, &h_cs->underspaces[i][j][k].numberOfMolecules, sizeof(h_cs->underspaces[i][j][k].numberOfMolecules), cudaMemcpyHostToDevice);
			}
		}
	}
	cudaMemcpy(&cs->dt, &h_cs->dt, sizeof(h_cs->dt), cudaMemcpyHostToDevice);
	cudaMemcpy(&cs->width, &h_cs->width, sizeof(h_cs->width), cudaMemcpyHostToDevice);
	cudaMemcpy(&cs->height, &h_cs->height, sizeof(h_cs->height), cudaMemcpyHostToDevice);
	cudaMemcpy(&cs->Nx, &h_cs->Nx, sizeof(h_cs->Nx), cudaMemcpyHostToDevice);
	cudaMemcpy(&cs->Ny, &h_cs->Ny, sizeof(h_cs->Ny), cudaMemcpyHostToDevice);
	cudaMemcpy(&cs->Nz, &h_cs->Nz, sizeof(h_cs->Nz), cudaMemcpyHostToDevice);
	//*d_cs = cs;
	return cs;
}

CUDASpace* copyFromDevice(CUDASpace *d_cs/*, CUDASpace *h_cs*/)
{
	CUDASpace *cs = new CUDASpace;
	cs->underspaces = new CUDAUnderspace**[d_cs->Nx];
	for (size_t i = 0; i < d_cs->Nx; ++i) {
		cs->underspaces[i] = new CUDAUnderspace*[d_cs->Ny];
		for (size_t j = 0; j < d_cs->Ny; ++j) {
			cs->underspaces[i][j] = new CUDAUnderspace[d_cs->Nz];
			for (size_t k = 0; k < d_cs->Nz; ++k) {
				cudaMalloc(&cs->underspaces[i][j][k].molecules, sizeof(CUDAMolecule) * d_cs->underspaces[i][j][k].numberOfMolecules);
				cudaMemcpy(cs->underspaces[i][j][k].molecules, d_cs->underspaces[i][j][k].molecules, d_cs->underspaces[i][j][k].numberOfMolecules, cudaMemcpyDeviceToHost);
				cudaMemcpy(&cs->underspaces[i][j][k].numberOfMolecules, &d_cs->underspaces[i][j][k].numberOfMolecules, sizeof(d_cs->underspaces[i][j][k].numberOfMolecules), cudaMemcpyDeviceToHost);
			}
		}
	}
	cudaMemcpy(&cs->dt, &d_cs->dt, sizeof(d_cs->dt), cudaMemcpyDeviceToHost);
	cudaMemcpy(&cs->width, &d_cs->width, sizeof(d_cs->width), cudaMemcpyDeviceToHost);
	cudaMemcpy(&cs->height, &d_cs->height, sizeof(d_cs->height), cudaMemcpyDeviceToHost);
	cudaMemcpy(&cs->Nx, &d_cs->Nx, sizeof(d_cs->Nx), cudaMemcpyDeviceToHost);
	cudaMemcpy(&cs->Ny, &d_cs->Ny, sizeof(d_cs->Ny), cudaMemcpyDeviceToHost);
	cudaMemcpy(&cs->Nz, &d_cs->Nz, sizeof(d_cs->Nz), cudaMemcpyDeviceToHost);

	return cs;
}

void freeDeviceMem(CUDASpace *d_cs)
{
	for (size_t i = 0; i < d_cs->Nx; ++i) {
		for (size_t j = 0; j < d_cs->Ny; ++j) {
			for (size_t k = 0; k < d_cs->Nz; ++k) {
				cudaFree(d_cs->underspaces[i][j][k].molecules);
			}
			cudaFree(d_cs->underspaces[i][j]);
		}
		cudaFree(d_cs->underspaces[i]);
	}
	cudaFree(d_cs->underspaces);
	cudaFree(d_cs);
}

void freeHostMem(CUDASpace *h_cs)
{
	for (size_t i = 0; i < h_cs->Nx; ++i) {
		for (size_t j = 0; j < h_cs->Ny; ++j) {
			for (size_t k = 0; k < h_cs->Nz; ++k) {
				delete h_cs->underspaces[i][j][k].molecules;
			}
			delete h_cs->underspaces[i][j];
		}
		delete h_cs->underspaces[i];
	}
	delete h_cs->underspaces;
	delete h_cs;
}


