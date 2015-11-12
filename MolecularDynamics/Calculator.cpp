#include "Calculator.h"
#include "Space.h"
#include <QEventLoop>
#include <QEventLoopLocker>

void Calculator::oneStep()
{
	//QThread::msleep(1);
	recalculatePositions_Beeman();
	calculateNewForces();
	recalculateSpeeds_Beeman();
	for (auto &i : space->molecules) {
		i.oldF = i.F;
		i.F = i.newF;
	}

// 	recalculatePositions_VelocityVerlet();
// 	recalculateForces_LennardJones();
// 	recalculateSpeeds_VelocityVerlet();
 	averageSpeed();
	//emit stateChanged();
}

inline Vector Calculator::Force_LennardJones(Molecule &m1, Molecule &m2)
{
	Vector r = m2.r - m1.r;
	register double square = r.square();
	if (maxDistSquare < square) return Vector();
	square = (Molecule::sigma * Molecule::sigma) / square;
	double U = 2.0*pow(square, 14/2) - pow(square, 8/2);
	U *= 24.0 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return -U * r;
}

inline Vector Calculator::Force_LennardJones(Vector r, double square)
{
	square = (Molecule::sigma * Molecule::sigma) / square;
	double U = 2.0*pow(square, 14 / 2) - pow(square, 8 / 2);
	U *= 24.0 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return -U * r;
}



void Calculator::calculateNewForces()
{
	for (int i = 0; i < space->molecules.size(); ++i) {
		space->molecules[i].newF = Vector();
		int size = space->molecules.size();
		for (int j = 0; j < size; ++j) {
			if (i == j) continue;
			Vector r = space->molecules[j].r - space->molecules[i].r;
			register double square = r.square();
			if (maxDistSquare < square) continue;
			//Vector dF = Force_LennardJones(space->molecules[i], space->molecules[j]);
#ifdef DEBUG
			if (!std::isfinite(dF)) throw std::exception("infinite force");
			if (std::isnan(dF)) throw std::exception("dF == NaN");
#endif
			//space->molecules[i].newF += dF;
			space->molecules[i].newF += Force_LennardJones(r, square);
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
}

void Calculator::recalculateSpeeds_VelocityVerlet()
{
	//Velocity Verlet
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.v = i.v + (i.F + i.oldF) / i.m * 0.5 * dt;
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
		i.r -= 1.0 / 6.0 * (i.oldF / i.m) * (dt*dt);
	}
}

void Calculator::recalculateSpeeds_Beeman()
{
	//Beeman's algorithm
	Space &space = *this->space;
	for (auto &i : space.molecules) {
		i.v += 2.0 / 6.0 * (i.newF / i.m) * dt;
		i.v += 5.0 / 6.0 * (i.F / i.m) * dt;
		i.v -= 1.0 / 6.0 * (i.oldF / i.m) * dt;
		if (i.r.x <= 0 || space.width * Angstrom <= i.r.x) {
			i.v.x = -i.v.x;
		}
		if (i.r.y <= 0 || space.height * Angstrom <= i.r.y) {
			i.v.y = -i.v.y;
		}
	}
}


void Calculator::set_dt_precision(int precision)	//hot only
{
	bool cr = calculationsRequired;
	calculationsRequired = false;
	space->mutex.lock();
	dt = std::pow(10, -precision);
	space->mutex.unlock();
	start();
	calculationsRequired = cr;
}



Calculator::Calculator(Space *space, QObject *parent /*= 0*/)
	:QObject(parent), space(space)
{
	calculationsRequired = false;
	QMetaObject::invokeMethod(this, "modeling", Qt::QueuedConnection);		//modeling() does nothing if calculationsRequired == false
}

double Calculator::pow(double d, int i)
{
	double result = d;
	if (i > 0) while (--i)  result *= d;
	else while (++i <= 1) result /= d;
	return result;
}


double Calculator::pow(Vector v, int i)
{
	if (i && 1 == 0)
		return pow(v.x*v.x + v.y*v.y + v.z*v.z, i / 2);
	else
		return pow(double(v), i);
}

void Calculator::start()
{
	space->mutex.lock();
	//pre-init:
	calculateNewForces();
	for (auto &i : space->molecules) {
		i.F = i.newF;
	}
	//init: first iteration
	recalculatePositions_VelocityVerlet();
	calculateNewForces();
	for (auto &i : space->molecules) {
		i.F = i.newF;
	}
	recalculateSpeeds_VelocityVerlet();
	averageSpeed();
	//next iterations:
	calculationsRequired = true;
	space->mutex.unlock();
	//modeling();
}

void Calculator::pause()
{
	calculationsRequired = false;
}

void Calculator::modeling()
{
	while (calculationsRequired) {
		space->mutex.lock();
		//space->mutex.lockForWrite();
#ifdef DEBUG
		qDebug() << "	enter in modeling() cycle";
#endif
		for (int i = 0; i < 20; ++i) {
			oneStep();
		}
#ifdef DEBUG
		qDebug() << "	return from modeling() cycle";
#endif
		space->mutex.unlock();	//WARNING мьютекс может не освободиться, если будет выкинуто исключение
	}

	QMetaObject::invokeMethod(this, "modeling", Qt::QueuedConnection);
}

void Calculator::averageSpeed()
{
	space->averageV = 0;
	for (auto i : space->molecules)
		space->averageV += i.v;
	space->averageV /= space->molecules.size();
	if (space->averageV < space->minV)
		space->minV = space->averageV;
	if (space->averageV > space->maxV)
		space->maxV = space->averageV;
	space->deltaV = space->maxV - space->minV;
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

