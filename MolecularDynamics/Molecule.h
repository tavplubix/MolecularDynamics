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
	double m_sigma, m_epsilon, m_mass;
	QString m_name;
public:
	MoleculeType(double sigma, double epsilon, double mass, const QString& name);
	double sigma();
	double sigma(const MoleculeType& mt);
	double epsilon();
	double epsilon(const MoleculeType& mt);
	double mass();
	QString name();
};

//==================================================================
//						Some types of molecules
//==================================================================

extern MoleculeType Ne;

class MoleculesGroup
{
public:
	MoleculeType type = Ne;
	double xshift = 0, yshift = 0, zshift = 0;
	int Nx = 0, Ny = 0, Nz = 0;
	double randomSpeedComponent = 10;
	double xspeed = 0, yspeed = 0, zspeed = 0;

	void setShift(double x, double y, double z) { xshift = x; yshift = y; zshift = z; }		//Shift components in m
	void setNumberOfMolecules(double x, double y, double z) { Nx = x; Ny = y; Nz = z; }
	void setSpeed(double x, double y, double z) { xspeed = x; yspeed = y; zspeed = z; }
	void setShiftInAngstroms(int x, int y, int z) { setShift(x*Angstrom, y*Angstrom, z*Angstrom); }

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





