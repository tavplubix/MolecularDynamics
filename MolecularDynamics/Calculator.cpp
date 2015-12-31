#include "Calculator.h"
#include <QEventLoop>
#include <QEventLoopLocker>
#include <QtConcurrent>
#include <QFuture>

#include "cuda.h"

//#define OLDCODE


void Calculator::_oneStep()
{
	auto lambdaP = [&](std::vector<std::vector<Underspace>> & part) {
		for (auto &i : part)
			for (auto &k : i)
				recalculatePositions_Beeman(k.molecules);
	};
	auto futuresP = QtConcurrent::map(space->underspaces, lambdaP);
	forAllM(t, space->underspaces)
		t.newF = Vector();
	futuresP.waitForFinished();

#ifndef OLDCUDA
	auto lambdaF = [&](std::vector<std::vector<Underspace>> & part) {
		for (auto &i : part)
			for (auto &k : i)
				calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	};
	auto futuresF = QtConcurrent::map(space->underspaces, lambdaF);
	futuresF.waitForFinished();
#else
	forAllU(k, space->underspaces)
		calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
#endif

	auto lambdaS = [&](std::vector<std::vector<Underspace>> & part) {
		for (auto &i : part)
			for (auto &k : i)
				recalculateSpeeds_Beeman(k.molecules);
	};
	auto futuresS = QtConcurrent::map(space->underspaces, lambdaS);
	futuresS.waitForFinished();


	forAllM(t, space->underspaces)
	{
		t.oldF = t.F;
		t.F = t.newF;
	}

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
	register const double sigmaSquare = Molecule::sigma * Molecule::sigma;
	square = sigmaSquare / square;
	register double U = 2.0*pow(square, 14 / 2) - pow(square, 8 / 2);
	register const double c = -24.0 * Molecule::epsilon / sigmaSquare;
	return c * U * r;
}


void Calculator::calculateNewForces(MoloculesList &molecules1, MoloculesList &molecules2)
{
	auto end1 = molecules1.end();
	for (auto i = molecules1.begin(); i != end1; ++i) {
		auto end2 = molecules2.end();
		for (auto j = molecules2.begin(); j != end2; ++j) {
			if (i._Ptr == j._Ptr) continue;
			Vector r = (*j).r - (*i).r;
			register double square = r.square();
			if (maxDistSquare < square) continue;
			(*i).newF += Force_LennardJones(r, square);
		}
	}
}

void Calculator::recalculatePositions_VelocityVerlet(MoloculesList &molecules)
{
	//Velocity Verlet
	for (auto &i : molecules) {
		i.oldr = i.r;
		i.r += i.v * dt;
		i.r += (i.F / i.m) * (dt*dt) * 0.5;
	}
}

void Calculator::recalculateSpeeds_VelocityVerlet(MoloculesList &molecules)
{
	//Velocity Verlet
	for (auto &i : molecules) {
		i.v = i.v + (i.F + i.oldF) / i.m * 0.5 * dt;
		if (i.r.x <= 0 || space->width * Angstrom <= i.r.x) {
			i.v.x = -i.v.x;
		}
		if (i.r.y <= 0 || space->height * Angstrom <= i.r.y) {
			i.v.y = -i.v.y;
		}
	}
}


void Calculator::recalculatePositions_Beeman(MoloculesList &molecules)
{
//#pragma omp parallel 
	for (auto &i : molecules) {
		i.oldr = i.r;
		i.r += i.v * dt;
		i.r += 4.0 / 6.0 * (i.F / i.m) * (dt*dt);
		i.r -= 1.0 / 6.0 * (i.oldF / i.m) * (dt*dt);
	}
}

void Calculator::recalculateSpeeds_Beeman(MoloculesList &molecules)
{
//#pragma omp parallel 
	for (auto &i : molecules) {
		i.v += 2.0 / 6.0 * (i.newF / i.m) * dt;
		i.v += 5.0 / 6.0 * (i.F / i.m) * dt;
		i.v -= 1.0 / 6.0 * (i.oldF / i.m) * dt;
		if (i.r.x <= 0 ) {
			i.v.x = std::abs(i.v.x);
		}
		if (space->width * Angstrom <= i.r.x) {
			i.v.x = - std::abs(i.v.x);
		}
		if (i.r.y <= 0) {
			i.v.y = std::abs(i.v.y);
		}
		if (space->height * Angstrom <= i.r.y) {
			i.v.y = - std::abs(i.v.y);
		}
	}
}

void Calculator::set_dt_precision(int precision)	//hot only
{
	bool cr = calculationsRequired;
	calculationsRequired = false;
	//QThread::sleep(1);
	space->mutex.lock();
	dt = std::pow(10, -precision);
	space->mutex.unlock();
	//start();
	calculationsRequired = cr;
}



Calculator::Calculator(Space *space, QObject *parent /*= 0*/)
	:QObject(parent), space(space)
{
	allocateMemory(space->numberOfMolecules);
	calculationsRequired = false;
	QMetaObject::invokeMethod(this, "modeling", Qt::QueuedConnection);		//modeling() does nothing if calculationsRequired == false
}

Calculator::~Calculator()
{
	freeMemory();
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
	forAllU(k, space->underspaces) 
		calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	forAllM(t, space->underspaces)
		t.F = t.newF;

	//init: first iteration
	forAllU(k, space->underspaces)
		recalculatePositions_VelocityVerlet(k.molecules);
	//forAllU(k, space->underspaces)
	//	validateUnderspace(k);

	forAllU(k, space->underspaces)
		calculateNewForcesForUnderspace(k.nx, k.ny, k.nz);
	forAllM(t, space->underspaces)
		t.F = t.newF;
	forAllU(k, space->underspaces)
		recalculateSpeeds_VelocityVerlet(k.molecules);	

	_averageSpeed();
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
	for (const auto &dx : closestSpaces) {
		for (const auto &dy : closestSpaces) {
			for (const auto &dz : closestSpaces) {
				int x = nx + dx; 
				int y = ny + dy;
				int z = nz + dz;
				if (x < 0 || y < 0 || z < 0) continue;
				if (x >= space->Nx) continue;
				if (y >= space->Ny) continue;
				if (z >= space->Nz) continue;
#ifndef OLDCUDA
				calculateNewForces(centralSpace.molecules, space->underspaces[x][y][z].molecules);
#else
				calculateNewForces_GPU(centralSpace.molecules, space->underspaces[x][y][z].molecules);
#endif
			}
		}
	}
}



void Calculator::normalizeUnderspace(Underspace &space)
{
	auto iter = space.molecules.begin();
	//for (auto iter = space.molecules.begin(); iter != space.molecules.end(); ++iter) {
	while (iter != space.molecules.end()) {
		Molecule molecule = *iter;
		int moveToX = molecule.r.x / space.size.x;
		int moveToY = molecule.r.y / space.size.y;
		int moveToZ = molecule.r.z / space.size.z;
		if (moveToX != space.nx || moveToY != space.ny || moveToZ != space.nz) {
			iter = space.molecules.erase(iter);
			this->space->underspaces[moveToX][moveToY][moveToZ].molecules.push_back(molecule);
		}
		else {
			++iter;
		}
	}
}

void Calculator::normalizeUnderspaces_Vector()
{
	space->molecules.clear();
	forAllM(t, space->underspaces)
		space->molecules.push_back(t);

	space->toUnderspaces();
}


void Calculator::freeze(double vMul /*= 0.999*/)
{
	forAllM(m, space->underspaces)
		m.v *= vMul;
}

void Calculator::modeling()
{
	while (calculationsRequired) {
		space->mutex.lock();
		//space->mutex.lockForWrite();
#ifdef DEBUG
		qDebug() << "	enter in modeling() cycle";
#endif
#ifndef CUDA
		for (int i = 0; i < 30; ++i) {
			_oneStep();
			space->iterations++;
			space->time_s += dt;
		}
		//forAllU(k, space->underspaces)
		//	normalizeUnderspace(k);
#else
		CUDASpace *h_cs = space->toCUDA();
		h_cs->dt = dt;
		size_t wholeSize = WHOLE_SIZE_OF_SPACE(h_cs);
		//CUDASpace *d_cs = nullptr;
		CUDASpace *d_cs = moveFromHost(h_cs, wholeSize);
		//freeHostMem(h_cs);

		for (int i = 0; i < 30; ++i) {
			cuda_oneStep(d_cs, space->Nx, space->Ny, space->Nz);
			space->iterations++;
			space->time_s += dt;
		}

		h_cs = moveFromDevice(d_cs, wholeSize);
		//freeDeviceMem(d_cs);
		space->fromCuda(h_cs);

#endif
		freeze(heating);
		_averageSpeed();
		normalizeUnderspaces_Vector();
#ifdef DEBUG
		qDebug() << "	return from modeling() cycle";
#endif
		space->mutex.unlock();	//WARNING мьютекс может не освободиться, если будет выкинуто исключение
	}

	QMetaObject::invokeMethod(this, "modeling", Qt::QueuedConnection);
}


void Calculator::_averageSpeed()
{
	space->averageV = 0;
	forAllM (t, space->underspaces)
					space->averageV += t.v;
	space->averageV /= space->numberOfMolecules;		//WARNING
	if (space->averageV < space->minV)
		space->minV = space->averageV;
	if (space->averageV > space->maxV)
		space->maxV = space->averageV;
	space->deltaV = space->maxV - space->minV;
}

