#include "Calculator.h"
#include <QEventLoop>
#include <QEventLoopLocker>

//#define OLDCODE

#ifdef OLDCODE
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
#endif

void Calculator::_oneStep()
{
	
	forAllU(k, space->underspaces)
	{
		recalculatePositions_Beeman(k.molecules);
	}
	forAllM(t, space->underspaces)
	{
		t.newF = Vector();
	}
	forAllU(k, space->underspaces)
	{
		calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	}
	forAllU(k, space->underspaces)
	{
		recalculateSpeeds_Beeman(k.molecules);
	}
	forAllM(t, space->underspaces)
	{
		t.oldF = t.F;
		t.F = t.newF;
	}
	_averageSpeed();
}

#ifdef OLDCODE
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
#endif

inline Vector Calculator::Force_LennardJones(Vector r, double square)
{
	square = (Molecule::sigma * Molecule::sigma) / square;
	double U = 2.0*pow(square, 14 / 2) - pow(square, 8 / 2);
	U *= 24.0 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return -U * r;
}


//#ifdef OLDCODE
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
			if (!std::isfinite(Force_LennardJones(r, square))) throw std::exception("infinite force");
			if (std::isnan(Force_LennardJones(r, square))) throw std::exception("dF == NaN");
#endif
			//space->molecules[i].newF += dF;
			space->molecules[i].newF += Force_LennardJones(r, square);
		}
	}
}
// #endif

void Calculator::calculateNewForces(std::list<Molecule> &molecules1, std::list<Molecule> &molecules2)
{
	for (auto i = molecules1.begin(); i != molecules1.end(); ++i) {
		//(*i).newF = Vector();			//FIXME
		int size = space->molecules.size();
		for (auto j = molecules2.begin(); j != molecules2.end(); ++j) {
			if (i._Ptr == j._Ptr) continue;
			Vector r = (*j).r - (*i).r;
			register double square = r.square();
			if (maxDistSquare < square) continue;
			//Vector dF = Force_LennardJones(space->molecules[i], space->molecules[j]);
#ifdef DEBUG
			if (!std::isfinite(Force_LennardJones(r, square))) throw std::exception("infinite force");
			if (std::isnan(Force_LennardJones(r, square))) throw std::exception("dF == NaN");
#endif
			//space->molecules[i].newF += dF;
			(*i).newF += Force_LennardJones(r, square);
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

#ifdef OLDCODE
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
#endif

void Calculator::recalculatePositions_Beeman(std::list<Molecule> &molecules)
{
	for (auto &i : molecules) {
		i.oldr = i.r;
		i.r += i.v * dt;
		i.r += 4.0 / 6.0 * (i.F / i.m) * (dt*dt);
		i.r -= 1.0 / 6.0 * (i.oldF / i.m) * (dt*dt);
	}
}

void Calculator::recalculateSpeeds_Beeman(std::list<Molecule> &molecules)
{
	for (auto &i : molecules) {
		i.v += 2.0 / 6.0 * (i.newF / i.m) * dt;
		i.v += 5.0 / 6.0 * (i.F / i.m) * dt;
		i.v -= 1.0 / 6.0 * (i.oldF / i.m) * dt;
		if (i.r.x <= 0 || space->width * Angstrom <= i.r.x) {
			i.v.x = -i.v.x;
		}
		if (i.r.y <= 0 || space->height * Angstrom <= i.r.y) {
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
	//forAllU(k, space->underspaces) {
	//	calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	//}

	for (auto &t : space->molecules) {
		t.F = t.newF;
	}
	//forAllM(t, space->underspaces){
	//	t.F = t.newF;
	//}

	//init: first iteration
	recalculatePositions_VelocityVerlet();		//FIXME

	calculateNewForces();
	for (auto &i : space->molecules) {
		i.F = i.newF;
	}
	//forAllU(k, space->underspaces)
	//{
	//	calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	//}
	//forAllM(t, space->underspaces)
	//{
	//	t.F = t.newF;
	//}

	recalculateSpeeds_VelocityVerlet();		//FIXME

	space->toUnderspaces();
	//averageSpeed();
	//next iterations:
	calculationsRequired = true;
	space->mutex.unlock();
	//modeling();
}

void Calculator::pause()
{
	calculationsRequired = false;
}




void Calculator::calculateNewForcesForUnderspace(int nx, int ny, int nz)
{
	//QThread::msleep(1);
	Underspace &centralSpace = space->underspaces[nx][ny][nz];
	auto closestSpaces = { -1, 0, 1 };
	for (auto dx : closestSpaces) {
		for (auto dy : closestSpaces) {
			for (auto dz : closestSpaces) {
				int x = nx + dx; 
				int y = ny + dy;
				int z = nz + dz;
				if (x < 0 || y < 0 || z < 0) continue;
				if (x >= space->underspaces.size()) continue;
				if (y >= space->underspaces[x].size()) continue;
				if (z >= space->underspaces[x][y].size()) continue;
				Underspace &neighboringSpace = space->underspaces[x][y][z];
				calculateNewForces(centralSpace.molecules, neighboringSpace.molecules);
			}
		}
	}
}



void Calculator::validateUnderspace(Underspace &space)
{
	for (auto iter = space.molecules.begin(); iter != space.molecules.end(); ++iter) {
		Molecule molecule = *iter;
		int moveToX = molecule.r.x / space.size.x;
		int moveToY = molecule.r.y / space.size.y;
		int moveToZ = molecule.r.z / space.size.z;
		if (moveToX != space.nx || moveToY != space.ny || moveToZ != space.nz) {
			space.molecules.erase(iter);
			--iter;
			this->space->underspaces[moveToX][moveToY][moveToZ].molecules.push_back(molecule);
		}
	}
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
			_oneStep();
		}
#ifdef DEBUG
		qDebug() << "	return from modeling() cycle";
#endif
		space->mutex.unlock();	//WARNING мьютекс может не освободиться, если будет выкинуто исключение
	}

	QMetaObject::invokeMethod(this, "modeling", Qt::QueuedConnection);
}

#ifdef OLDCODE
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
}
#endif

void Calculator::_averageSpeed()
{
	space->averageV = 0;
	for (auto i : space->underspaces)
		for (auto j : i)
			for (auto k : j)
				for (auto t : k.molecules)
					space->averageV += t.v;
	space->averageV /= space->molecules.size();		//WARNING
	if (space->averageV < space->minV)
		space->minV = space->averageV;
	if (space->averageV > space->maxV)
		space->maxV = space->averageV;
	space->deltaV = space->maxV - space->minV;
}

#ifdef OLDCODE
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
#endif
