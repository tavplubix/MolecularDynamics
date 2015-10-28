#include "Vector.h"

Point::Point()
{
	x = 0;
	y = 0;
	z = 0;
}
Point::Point(double x_, double y_, double z_)
	:x(x_), y(y_), z(z_)
{

}
Point::Point(const Point &P)
{
	x = P.x;
	y = P.y;
	z = P.z;
}

Point& Point::operator=(const Point &P)
{
	x = P.x;
	y = P.y;
	z = P.z;
	return *this;
}

Vector::Vector()
{
	x = 0;
	y = 0;
	z = 0;
	v = 0;
}
Vector::Vector(double endX, double endY, double endZ, double startX, double startY, double startZ)
{
	x = endX - startX;
	y = endY - startY;
	z = endZ - startZ;
	v = sqrt(x*x + y*y + z*z);
}
Vector::Vector(const Vector &V)
{
	x = V.x;
	y = V.y;
	z = V.z;
	v = V.v;
}

Vector& Vector::operator=(const Vector &V)
{
	x = V.x;
	y = V.y;
	z = V.z;
	v = V.v;
	return *this;
}

Vector Vector::operator/(const double I) const
{
	return (*this) * (1.0 / I);
}

Vector& Vector::operator+=(const Vector &V)
{
	*this = *this + V;
	return *this;
}
Vector& Vector::operator-=(const Vector &V)
{
	*this = *this - V;
	return *this;
}
Vector& Vector::operator*=(const double I)
{
	*this = *this * I;
	return *this;
}

Vector operator*(const double I, const Vector &V)
{
	return V * I;
}

Vector Vector::operator+(const Vector &V) const
{
	return Vector(x + V.x, y + V.y, z + V.z);
}
Vector Vector::operator-(const Vector &V) const
{
	return Vector(x - V.x, y - V.y, z - V.z);
}
Vector Vector::operator*(const double I) const
{
	return  Vector(x * I, y * I, z * I);
}


