#pragma once
#include "classes.h"



class Calculator : public QObject
{
	Q_OBJECT
private:
	Space* space;
	double dt = 1.0e-15;
	double maxDistSquare = pow(4 * Molecule::sigma, 2);
	void averageSpeed();
	inline Vector Force_LennardJones(Molecule &m1, Molecule &m2);
	inline Vector Calculator::Force_LennardJones(Vector r, double square);		//r - distance between two molecules, square=r*r
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
	Q_INVOKABLE void set_dt_precision(int precision);
	void setSpace(Space* s);
	double get_dt() { return dt; };
	volatile bool calculationsRequired;
	Calculator(Space *space, QObject *parent = 0);
	//static double LennardJonesPotential(Molecule &m1, Molecule &m2);
	static double pow(double d, int i);
	static double pow(Vector v, int i);
	public slots:
	void start();
	void pause();
signals:
	void stateChanged();
};




