
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <vector>



#include "Vector.h"
#include "Molecule.h"
#include "cuda.h"

#define Boltzmann  1.3806488e-23d
#define Angstrom  10e-10d
#define AtomicMassUnit  1.660538921e-27d
#define sigma 2.74*10e-10	//m
#define epsilon 36.2*1.3806488e-23	//J





typedef std::vector<Molecule> MoloculesList;

const size_t maxSize = 1024 * 16;
double h_newF[3 * maxSize], h_r[3 * maxSize];



__global__ void kernel_calculateSquares(double *squares, const double *r, int n)
{
	int molecule = blockDim.x * blockIdx.x + threadIdx.x;
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
		if (size > maxSize)
			throw 0;		//WARNING
		double *d_newF, *d_r, *d_squares;
		cudaMalloc(&d_newF, sizeof(double) * 3 * size);
		cudaMalloc(&d_r, sizeof(double) * 3 * size);
		cudaMalloc(&d_squares, sizeof(double) * size);

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
		cudaFree(d_newF);
		cudaFree(d_r);
		cudaFree(d_squares);

	}
}










