#pragma once
#include <math.h>

class Vector
{
public:
	double startX, startY, startZ;
	double endX, endY, endZ;
	double dx, dy, dz;
	double v;

	Vector();
	Vector(double startX, double startY, double startZ, double endX, double endY, double endZ);
	Vector(const Vector &V);

	/*Vector& operator=(Vector &v);*/

	double GetABS();
};