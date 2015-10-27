#include "Calculator.h"


void Calculator::oneStep()
{
	QThread::msleep(10);
	Space &space = *this->space;
	recalculatePositions_VelocityVerlet();
	recalculateForces_LennardJones();
	recalculateSpeeds_VelocityVerlet();
	averageSpeed();
	//emit stateChanged();
}

Vector Calculator::Force_LennardJones(Molecule &m1, Molecule &m2)
{
	Vector r = m2.r - m1.r;
	double U = 2.0*pow(Molecule::sigma / r, 14) - pow(Molecule::sigma / r, 8);
	U *= 24.0 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return -U * r;
}

void Calculator::recalculateForces_LennardJones()
{
	Space &space = *this->space;
	for (int i = 0; i < space.molecules.size(); ++i) {
		space.molecules[i].oldF = space.molecules[i].F;
		space.molecules[i].F = Vector();
		for (int j = 0; j < space.molecules.size(); ++j) {
			if (i == j) continue;
			Vector dF = Force_LennardJones(space.molecules[i], space.molecules[j]);
			if (!std::isfinite(dF)) throw std::exception("infinite force");
			if (std::isnan(dF)) throw std::exception("dF == NaN");
			space.molecules[i].F += dF;
		}
	}
}

void Calculator::recalculatePositions_VelocityVerlet()
{
	//Velocity Verlet
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.oldr = i.r;
		i.r += i.v * dt;
		i.r += (i.F / i.m) * (dt*dt) * 0.5;
	}
	//TODO использовать алгоритм Бимана для последующих итераций
}

void Calculator::recalculateSpeeds_VelocityVerlet()
{
	//Velocity Verlet
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.v = i.v + (i.F + i.oldF) / i.m * 0.5 * dt;
		//i.v = i.v + (i.F / i.m) * dt;
		if (i.r.x <= 0 || space.width * Angstrom <= i.r.x) {
			i.v.x = -i.v.x;
		}
		if (i.r.y <= 0 || space.height * Angstrom <= i.r.y) {
			i.v.y = -i.v.y;
		}
	}
}

void Calculator::recalculatePositions_Beeman()
{
	//Beeman's algorithm
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.oldr = i.r;
		i.r += i.v * dt;
		i.r += 4.0 / 6.0 * (i.F / i.m) * (dt*dt);
		i.r += 1.0 / 6.0 * (i.oldF / i.m) * (dt*dt);
	}
}

void Calculator::recalculateSpeed_Beeman()
{
	//Beeman's algorithm
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.v += 2.0 / 6.0 * (i.newF / i.m);
		i.v += 5.0 / 6.0 * (i.F / i.m);
		i.v -= 1.0 / 6.0 * (i.oldF / i.m);
	}
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
	recalculateForces_LennardJones();
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
	//space->averageV = space->molecules[0].v;
}


Vector Calculator::integrateWithTaylorAproximation(double h, const Vector &f, const Vector &d1f /*= Vector()*/, const Vector &d2f /*= Vector()*/)
{
	Vector result = f * h;
	result += d1f * (h*h) / 2.0;
	result += d2f * (h*h*h) / 6.0;
	return result;
}

Vector Calculator::VerletIntegration(const Vector &r, const Vector &oldr, const Vector &a)
{
	Vector newr = 2 * r;
	newr -= oldr;
	newr += a * (dt*dt);
	return newr;
}

