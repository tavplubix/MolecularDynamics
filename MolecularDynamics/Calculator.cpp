#include "Caculator.h"


void Calculator::oneStep()
{
	Space &space = *this->space;
	//forces

	for (int i = 0; i < space.molecules.size(); ++i) {
		space.molecules[i].oldF = space.molecules[i].F;
		space.molecules[i].F = Vector();
		for (int j = 0; j < space.molecules.size(); ++j) {
			if (i == j) continue;
			Vector dF = Force(space.molecules[i], space.molecules[j]);
			if (!std::isfinite(dF)) throw std::exception("infinite force");
			if (std::isnan(dF)) throw std::exception("dF == NaN");
			space.molecules[i].F += dF;
		}
	}

	for (auto &i : space.molecules) {
		//speeds
		Vector v = i.v + ((i.F + i.oldF) / (2.0*i.m)) * dt;
		//i.v += integrateWithTaylorAproximation(dt, i.F / i.m, (i.F - i.oldF)/dt);		//WARNING
		//cordinates
		i.r += (i.v + v) * 0.5 * dt;
		//i.r += integrateWithTaylorAproximation(dt, i.v, i.F / i.m);		
		i.v = v;
		if (i.r.x <= 0 || space.width * Angstrom <= i.r.x) {
			i.v.x = -i.v.x;
		}
		if (i.r.y <= 0 || space.height * Angstrom <= i.r.y) {
			i.v.y = -i.v.y;
		}
	}
	averageSpeed();
	//emit stateChanged();
}

Vector Calculator::Force(Molecule &m1, Molecule &m2)
{
	Vector r = m2.r - m1.r;
	double U = 2.0*pow(Molecule::sigma / r, 14) - pow(Molecule::sigma / r, 8);
	U *= 24.0 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return -U * r;
}

Calculator::Calculator(Space *space, QObject *parent /*= 0*/)
	:QObject(parent), space(space)
{
	calculationsRequired = false;
}

double Calculator::pow(double d, int i)
{
	double result = d;
	while (--i)  result *= d;
	return result;
}

void Calculator::start()
{
	calculationsRequired = true;
	//QMetaObject::invokeMethod(this, "modeling");
	modeling();
}

void Calculator::stop()
{
	calculationsRequired = false;
}

void Calculator::modeling()
{
	while (calculationsRequired) {
		space->mutex.lock();
		//space->mutex.lockForWrite();
		qDebug() << "	enter in modeling() cycle";
		for (int i = 0; i < 5; ++i) {
			oneStep();
		}
		qDebug() << "	return from modeling() cycle";
		space->mutex.unlock();	//WARNING мьютекс может не освободиться, если будет выкинуто исключение
	}
}

void Calculator::averageSpeed()
{
	space->averageV = 0;
	for (auto i : space->molecules)
		space->averageV += i.v;
	space->averageV /= space->molecules.size();
}


Vector Calculator::integrateWithTaylorAproximation(double h, const Vector &f, const Vector &d1f /*= Vector()*/, const Vector &d2f /*= Vector()*/)
{
	Vector result = f * h;
	result += d1f * (h*h) / 2.0;
	result += d2f * (h*h*h) / 6.0;
	return result;
}

