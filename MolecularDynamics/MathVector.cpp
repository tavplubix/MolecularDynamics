#include "MathVector.h"

MathVector3D::MathVector3D()
{
	x = 0;
	y = 0;
	z = 0;
	//v = 0;
}
MathVector3D::MathVector3D(double endX, double endY, double endZ, double startX, double startY, double startZ)
{
	x = endX - startX;
	y = endY - startY;
	z = endZ - startZ;
	//v = sqrt(x*x + y*y + z*z);
}
MathVector3D::MathVector3D(const MathVector3D &V)
{
	x = V.x;
	y = V.y;
	z = V.z;
	//v = V.v;
}

MathVector3D& MathVector3D::operator=(const MathVector3D &V)
{
	x = V.x;
	y = V.y;
	z = V.z;
	//v = V.v;
	return *this;
}

MathVector3D MathVector3D::operator/(const double I) const
{
	return MathVector3D(x / I, y / I, z / I);
}

MathVector3D& MathVector3D::operator+=(const MathVector3D &V)
{
	*this = *this + V;
	return *this;
}
MathVector3D& MathVector3D::operator-=(const MathVector3D &V)
{
	*this = *this - V;
	return *this;
}
MathVector3D& MathVector3D::operator*=(const double I)
{
	*this = *this * I;
	return *this;
}

MathVector3D operator*(const double I, const MathVector3D &V)
{
	return V * I;
}

MathVector3D MathVector3D::operator+(const MathVector3D &V) const
{
	return MathVector3D(x + V.x, y + V.y, z + V.z);
}
MathVector3D MathVector3D::operator-(const MathVector3D &V) const
{
	return MathVector3D(x - V.x, y - V.y, z - V.z);
}
MathVector3D MathVector3D::operator*(const double I) const
{
	return  MathVector3D(x * I, y * I, z * I);
}


void MathVector3D::toCUDA(CUDAVector& cv) const
{
	cv.v[0] = x;
	cv.v[1] = y;
	cv.v[2] = z;
}

void MathVector3D::fromCUDA(CUDAVector& cv)
{
	x = cv.v[0];
	y = cv.v[1];
	z = cv.v[2];
}

