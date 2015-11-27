#pragma once
#include "Vector.h"


namespace constants {
	const double Boltzmann = 1.3806488e-23;
	const double Avogadro = 6.022140857e23;
	const double Angstrom = 10e-10;
	const double Electronvolt = 1.6021766208e-19;
	const double AtomicMassUnit = 1.660538921e-27;
	const double pi = 3.141592654;
}

using namespace constants;


class Molecule
{
public:
	//double x, y;
	//double vx, vy;
	//double Fx, Fy;
	Vector r, oldr, v, F, oldF, newF;
	const static double m, radius, sigma, epsilon;
};





