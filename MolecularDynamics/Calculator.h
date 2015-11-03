#pragma once
#include "classes.h"



class Calculator : public QObject
{
	Q_OBJECT
private:
	Space* space;
	double dt = 10e-15;
	void averageSpeed();
	Vector Force_LennardJones(Molecule &m1, Molecule &m2);
	void recalculateForces_LennardJones();
	void recalculatePositions_VelocityVerlet();
	void recalculateSpeeds_VelocityVerlet();
	void recalculatePositions_Beeman();
	void recalculateSpeeds_Beeman();
	void calculateNewForces();		//for Beeman
	Vector VerletIntegration(const Vector &r, const Vector &oldR, const Vector &a);
	Vector integrateWithTaylorAproximation(double h, const Vector &f, const Vector &d1f = Vector(), const Vector &d2f = Vector());
	void oneStep();
	Q_INVOKABLE void modeling();
public:
	volatile bool calculationsRequired;
	Calculator(Space *space, QObject *parent = 0);
	//static double LennardJonesPotential(Molecule &m1, Molecule &m2);
	static double pow(double d, int i);
	public slots:
	void start();
	void pause();
signals:
	void stateChanged();
};




