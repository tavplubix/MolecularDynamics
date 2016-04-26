#pragma once
#include "MathVector.h"
#include <vector>
#include <QString>


namespace constants {
	const double Boltzmann = 1.3806488e-23;
	const double Avogadro = 6.022140857e23;
	const double Angstrom = 1e-10;
	const double Electronvolt = 1.6021766208e-19;
	const double AtomicMassUnit = 1.660538921e-27;
	const double pi = 3.141592654;
}

using namespace constants;

class MoleculeType
{
	double m_sigma, m_epsilon;
	QString m_name;
public:
	MoleculeType(double sigma, double epsilon, const QString& name);
	double sigma();
	double sigma(const MoleculeType& mt);
	double epsilon();
	double epsilon(const MoleculeType& mt);
	QString name();
};


class Molecule
{
public:
	//double x, y;
	//double vx, vy;
	//double Fx, Fy;
	short type;
	int id;
	MathVector3D r, oldr, v, F, oldF, newF;
	const static double m, radius, sigma, epsilon;
	static std::vector<double> epsilonv, sigmav;
};

//==================================================================
//						Some types of molecules
//==================================================================

extern MoleculeType Ne;



