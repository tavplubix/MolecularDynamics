#pragma once

#include <QMutex>
#include <QFile>
#include <list>
#include <vector>
#include <functional>
#include <iterator>


#include "Molecule.h"

#define forAllU(t, s) for (auto &i: s) for (auto &j: i) for (auto &t: j)
#define forAllM(t, s) for (auto &i: s) for (auto &j: i) for (auto &k: j) for (auto &t: k.molecules)


//class MoleculesContainer;


//typedef std::list<Molecule> MoloculesList;
typedef std::vector<Molecule> MoloculesList;

class Underspace
{
public:
	MoloculesList molecules;
	//const Vector minR, maxR;
	int nx, ny, nz;
	const static MathVector3D size;
	//Underspace(MoloculesList &&molecules);
	//Underspace(std::function<Molecule()> generator);
	//Underspace(std::function<Vector()> speedsGenerator, std::function<Vector(Vector)> positionsGenerator);
	//bool shouldContains(Molecule m);
	//bool actuallyContains(Molecule m);
	void toCUDA(CUDAUnderspace *cus, CUDAMolecule *placeForMolecules) const;
	void fromCUDA(CUDAUnderspace *cus);
};
template<typename T>
using  Matrix3D = std::vector < std::vector< std::vector<T> > >;

class DeprecatedSpace
{
private:
	void generateCoordinates();
	void generateSpeeds();
	void initializeUnderspaces();

	void generate2DWall();
	void generate2DBall();
	void generate2DRectangle(int xshift, int yshift, int xsize, int ysize, int type, int xspeed = 0, int yspeed = 0);
	void generate3DRectangle(int xshift, int yshift, int zshift, int xsize, int ysize, int zsize, int type, int xspeed = 0, int yspeed = 0, int zspeed = 0);
public:
	void toUnderspaces();
public:
	double time_s;
	unsigned long long iterations;
	//QReadWriteLock mutex;
	//TODO incapsulate this fields:
	mutable QMutex mutex;
	int width, height, depth;
	int Nx, Ny, Nz;
	double averageV;
	double K, U;
	double maxV, minV, deltaV;
	int numberOfMolecules;
	QFile trajektoryFile;
	long long int trajektoryTime = 103000;

	std::vector<std::vector<double>> sigma, epsilon;

	//===============================
	DeprecatedSpace(int width, int height, int n);
	DeprecatedSpace& operator=(const DeprecatedSpace&& s);
	std::vector<Molecule> molecules;
	Matrix3D<Underspace> underspaces;
	//public slots:
	//void saveCoordinates();
	void saveCoordinatesAndSpeeds(const QString& filename);
	void saveTrajektory();
	//void saveAll();
	//void loadStateC(const QString& filename);
	void loadStateCS(const QString& filename);	//load coordinates and speeds (CS)
	//void loadStateALL(const QString& filename);


	CUDASpace* toCUDA() const;
	void fromCuda(CUDASpace *cs);
};

/*
class SpaceIterator 
	: public std::iterator < std::forward_iterator_tag, Space >
{
	
public:

};
*/

/*
class MoleculesContainer
	: private std::vector<Molecule>
{
	static const double Reserve;
public:

};

const double MoleculesContainer::Reserve = 1.2;
*/
