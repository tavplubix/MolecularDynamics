#include "Molecule.h"
#include <math.h>
//Ne
const double Molecule::radius = Angstrom;	//m
const double Molecule::m = 20.1797 * AtomicMassUnit;	//kg
const double Molecule::sigma = 2.74 * Angstrom;		//m
const double Molecule::epsilon = 1.0 * 36.2 * Boltzmann;	//J

MoleculeType::MoleculeType(double sigma, double epsilon, const QString& name)
	: m_sigma(sigma), m_epsilon(epsilon), m_name(name)
{

}

double MoleculeType::sigma()
{
	return m_sigma;
}

double MoleculeType::sigma(const MoleculeType& mt)
{
	return 0.5*(m_sigma + mt.m_sigma);
}

double MoleculeType::epsilon()
{
	return m_epsilon;
}

double MoleculeType::epsilon(const MoleculeType& mt)
{
	return sqrt(m_epsilon * mt.m_epsilon);
}

QString MoleculeType::name()
{
	return m_name;
}

//==================================================================
//						Some types of molecules
//==================================================================

MoleculeType Ne(2.74 * Angstrom, 1.0 * 36.2 * Boltzmann, "Ne");


