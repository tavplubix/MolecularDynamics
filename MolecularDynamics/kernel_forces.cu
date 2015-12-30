
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <vector>

#include "Space.h"

#include "Vector.h"
#include "Molecule.h"
#include "cuda.h"

#define Boltzmann  1.3806488e-23d
#define Angstrom  10e-10d
#define AtomicMassUnit  1.660538921e-27d
#define sigma 2.74*10e-10	//m
#define epsilon 36.2*1.3806488e-23	//J

#define maxDistSquare 4.0 * sigma * 4.0 * sigma;




typedef std::vector<Molecule> MoloculesList;

const size_t maxSize = 1024 * 16;
double h_newF[3 * maxSize], h_r[3 * maxSize];
double *d_newF, *d_r, *d_squares;
double *space;


void allocateMemory(size_t size)
{
	if (size > maxSize)
		throw 0;		//WARNING
	cudaMalloc(&d_newF, sizeof(double) * 3 * size);
	cudaMalloc(&d_r, sizeof(double) * 3 * size);
	cudaMalloc(&d_squares, sizeof(double) * size);
}
void freeMemory()
{
	cudaFree(d_newF);
	cudaFree(d_r);
	cudaFree(d_squares);
}

__global__ void kernel_calculateSquares(double *squares, const double *r, int n)
{
	auto lambdaTest = [&]() -> int {
		return blockDim.x * blockIdx.x + threadIdx.x;
	};
	int molecule = lambdaTest();//blockDim.x * blockIdx.x + threadIdx.x;
	const double sigmaSquare = sigma * sigma;
	if (molecule >= n) return;
	int rIndex = molecule * 3;
	double sqX = r[rIndex + 0] * r[rIndex + 0];
	double sqY = r[rIndex + 1] * r[rIndex + 1];
	double sqZ = r[rIndex + 2] * r[rIndex + 2];
	double sqR = sqX + sqY + sqZ;		//sqR = r^2
	if (sqR == 0) {
		squares[molecule] = 0;
		return;
	}
	double sq = sigmaSquare / sqR;		//sq = (sigma/r)^2
	squares[molecule] = sq*sq*sq*sq;			// = (sigma/r)^8
}


__global__ void kernel_calculateNewForces_GPU(double *forces, const double *r, const double *squares, int n)
{
	int molecule = blockDim.x * blockIdx.x + threadIdx.x;
	int component = blockDim.y * blockIdx.y + threadIdx.y;
	if (molecule >= n) return;
	//if (component >= 3) return;
	const double c = -(24.0 * epsilon) / (sigma * sigma);
	int componetIndex = molecule * 3 + component;
	double sq = squares[molecule];
	forces[componetIndex] += c * ((2 * sq*sq) - (sq)) * r[componetIndex];
}





void calculateNewForces_GPU(MoloculesList &molecules1, MoloculesList &molecules2)
{
	auto end1 = molecules1.end();
	for (auto i = molecules1.begin(); i != end1; ++i) {
		auto end2 = molecules2.end();

		//copy data to C-array
		size_t index = 0;
		for (auto j = molecules2.begin(); j != end2; ++j) {
			h_r[index++] = (*j).r.x;
			h_r[index++] = (*j).r.y;
			h_r[index++] = (*j).r.z;
		}

		//allocate memory on GPU
		size_t size = molecules2.size();
		//if (size > maxSize)
		//	throw 0;		//WARNING
		//double *d_newF, *d_r, *d_squares;
		//cudaMalloc(&d_newF, sizeof(double) * 3 * size);
		//cudaMalloc(&d_r, sizeof(double) * 3 * size);
		//cudaMalloc(&d_squares, sizeof(double) * size);

		//copy array to GPU
		cudaMemcpy(d_r, h_r, 3 * sizeof(double)*size, cudaMemcpyHostToDevice);

		//calculations
		int numberOfCores = 512;
		dim3 grid, blocks;
		grid = dim3(size / numberOfCores + 1, 1, 1);
		blocks = dim3(numberOfCores, 1, 1);
		kernel_calculateSquares << <grid, blocks >> >(d_squares, d_r, size);
		cudaDeviceSynchronize();
		blocks = dim3(numberOfCores, 3, 1);
		kernel_calculateNewForces_GPU << <grid, blocks >> >(d_newF, d_r, d_squares, size);
		cudaDeviceSynchronize();

		//copy array from GPU;
		cudaMemcpy(h_r, d_r, 3 * sizeof(double)*size, cudaMemcpyDeviceToHost);

		//copy C-array to data
		index = 0;
		for (auto j = molecules2.begin(); j != end2; ++j) {
			(*j).r.x = h_r[index++];
			(*j).r.y = h_r[index++];
			(*j).r.z = h_r[index++];
		}

		//free memory on GPU
		//cudaFree(d_newF);
		//cudaFree(d_r);
		//cudaFree(d_squares);

	}
}


/*
__device__ inline Vector Force_LennardJones_GPU(Vector r, double square)
{
	register const double sigmaSquare = sigma * sigma;
	square = sigmaSquare / square;
	register double U = 2.0*pow(square, 14 / 2) - pow(square, 8 / 2);
	register const double c = -24.0 * epsilon / sigmaSquare;
	return c * U * r;
}

__device__ void calculateNewForces_GPU_rutine(MoloculesList &molecules1, MoloculesList &molecules2)
{
	auto end1 = molecules1.end();
	for (auto i = molecules1.begin(); i != end1; ++i) {
		auto end2 = molecules2.end();
		for (auto j = molecules2.begin(); j != end2; ++j) {
			if (i._Ptr == j._Ptr) continue;
			Vector r = (*j).r - (*i).r;
			register double square = r.square();
			//if (maxDistSquare < square) continue;
			(*i).newF += Force_LennardJones_GPU(r, square);
		}
	}
}

__global__ void calculateNewForcesForUnderspace_GPU(Space* space, int nx, int ny, int nz)
{
	//QThread::msleep(1);
	Underspace &centralSpace = space->underspaces[nx][ny][nz];
	auto closestSpaces = { -1, 0, 1 };
	for (const auto &dx : closestSpaces) {
		for (const auto &dy : closestSpaces) {
			for (const auto &dz : closestSpaces) {
				int x = nx + dx;
				int y = ny + dy;
				int z = nz + dz;
				if (x < 0 || y < 0 || z < 0) continue;
				if (x >= space->Nx) continue;
				if (y >= space->Ny) continue;
				if (z >= space->Nz) continue;
				calculateNewForces_GPU_rutine(centralSpace.molecules, space->underspaces[x][y][z].molecules);
			}
		}
	}
}*/



