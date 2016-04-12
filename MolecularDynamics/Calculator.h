#pragma once
#include <QObject>

#include "Molecule.h"
#include "Space.h"

#define CUDA


class Calculator : public QObject
{
	Q_OBJECT
private:
	Space* space;
	double dt = 1.0e-14;
	double maxDistSquare = pow(4 * Molecule::sigma, 2);

	void _averageSpeed();
	void _wholeEnergy();
	double calculatePotentialEnergyForUnderspace(int nx, int ny, int nz);
	double calculatePotentionalEnergy(MoloculesList &molecules1, MoloculesList &molecules2);

	inline Vector Force_LennardJones(Molecule &m1, Molecule &m2);
	inline Vector Force_LennardJones(Vector r, double square);		//r - distance between two molecules, square=r*r

	void recalculatePositions_VelocityVerlet(MoloculesList &molecules);
	void recalculateSpeeds_VelocityVerlet(MoloculesList &molecules);


	void recalculatePositions_Beeman(MoloculesList &molecules);
	void recalculateSpeeds_Beeman(MoloculesList &molecules);
	void calculateNewForces();		//for Beeman
	void calculateNewForces(MoloculesList &molecules1, MoloculesList &molecules2);

	void _oneStep();

	void calculateNewForcesForUnderspace(int nx, int ny, int nz);
	void normalizeUnderspace(Underspace &space);
	void normalizeUnderspaces_Vector();

	void freeze(double vMul = 0.999);

	Q_INVOKABLE void modeling();
public:
	volatile double heating = 1;
	Q_INVOKABLE void set_dt_precision(int precision);
	void setSpace(Space* s);
	double get_dt() { return dt; };
	volatile bool calculationsRequired;
	Calculator(Space *space, QObject *parent = 0);
	~Calculator();
	//static double LennardJonesPotential(Molecule &m1, Molecule &m2);
	static double pow(double d, int i);
	static double pow(Vector v, int i);
	public slots:
	void start();
	void pause();
signals:
	void stateChanged();
};




