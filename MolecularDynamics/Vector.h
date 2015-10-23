#pragma once
#include <math.h>

class Point;
class Vector;

class Point
{
private:
	double x, y, z;
public:
	Point();
	Point(double x_, double y_, double z_);
	Point(const Point &P);

	Point& operator=(const Point &P);

	double getX() {
		return x;
	}
	double getY() {
		return y;
	}
	double getZ() {
		return z;
	}

	void setX(double x_) {
		x = x_;
	}
	void setY(double y_) {
		y = y_;
	}
	void setZ(double z_) {
		z = z_;
	}
	void setXYZ(double x_, double y_, double z_)
	{
		setX(x_);
		setY(y_);
		setZ(z_);
	}

};

class Vector
{
public:
	double x, y, z;
	double v;

	Vector();
	Vector(Point S, Point E);
	Vector(double endX, double endY, double endZ, double startX = 0, double startY = 0, double startZ = 0);
	Vector(const Vector &V);

	Vector& operator=(const Vector &V);
	Vector& operator+=(const Vector &V);
	Vector& operator-=(const Vector &V);
	Vector& operator*=(const double &I);
	Vector operator+(const Vector &V);
	Vector operator-(const Vector &V);
	Vector operator*(const double &I);
	operator double() { return v; }

};