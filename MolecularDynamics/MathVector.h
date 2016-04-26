#pragma once
#include <math.h>
#include "c_cuda_structures.h"



class MathVector3D
{
public:
	friend MathVector3D operator*(const double I, const MathVector3D &V);

	double x, y, z;
	//double v;

	MathVector3D();
	MathVector3D(double endX, double endY, double endZ, double startX = 0, double startY = 0, double startZ = 0);
	MathVector3D(const MathVector3D &V);

	MathVector3D& operator=(const MathVector3D &V);
	MathVector3D& operator+=(const MathVector3D &V);
	MathVector3D& operator-=(const MathVector3D &V);
	MathVector3D& operator*=(const double I);
	MathVector3D operator/(const double I) const;
	MathVector3D operator+(const MathVector3D &V) const;
	MathVector3D operator-(const MathVector3D &V) const;
	MathVector3D operator*(const double I) const;
	double square() { return x*x + y*y + z*z; };
	operator double() { return sqrt(square()); }
	static double ScalarMultiply(const MathVector3D &V1, const MathVector3D &V2)
	{
		return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
	}

	void toCUDA(CUDAVector& cv) const;
	void fromCUDA(CUDAVector& cv);
};