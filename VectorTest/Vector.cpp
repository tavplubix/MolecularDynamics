#include "Vector.h"

Vector::Vector(double startX, double startY, double startZ, double endX, double endY, double endZ)
	: startX(startX), startY(startY), startZ(startZ), endX(endX), endY(endY), endZ(endZ)
{
	dx = endX - startX;
	dy = endY - startY;
	dz = endZ - startZ;

	v = sqrt(dx*dx + dy*dy + dz*dz);
}

Vector::Vector(const Vector &V) {
	startX = V.startX;
	startY = V.startY;
	startZ = V.startZ;
	endX = V.endX;
	endY = V.endY;
	endZ = V.endZ;
	dx = V.dx;
	dy = V.dy;
	dz = V.dz;
	v = V.v;
}

double Vector::GetABS() {
	return v;
}
