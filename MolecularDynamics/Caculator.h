#pragma once
#include "classes.h"



class Calculator : public QObject
{
	Q_OBJECT
private:
	Space* space;
	double dt = 10e-14;
	bool calculationsRequired;
	void averageSpeed();
	static Vector Force(Molecule &m1, Molecule &m2);
	static Vector integrateWithTaylorAproximation(double h, const Vector &f, const Vector &d1f, const Vector &d2f = Vector());
	void oneStep();
public:
	void modeling();
	Calculator(Space *space, QObject *parent = 0);
	//static double LennardJonesPotential(Molecule &m1, Molecule &m2);
	static double pow(double d, int i);
	public slots:
	void start();
	void stop();
signals:
	void stateChanged();
};



