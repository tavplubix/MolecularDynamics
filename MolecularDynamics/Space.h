#pragma once

#include <QObject>
#include <QMutex>
#include <QErrorMessage>
#include "classes.h"
#include <list>
#include <vector>
#include <functional>




class Underspace
{
public:
	std::list<Molecule> molecules;
	//const Vector minR, maxR;
	//int nx, ny, nz;
public:
	const static Vector size;
	//Underspace(std::list<Molecule> &&molecules);
	//Underspace(std::function<Molecule()> generator);
	//Underspace(std::function<Vector()> speedsGenerator, std::function<Vector(Vector)> positionsGenerator);
	//bool shouldContains(Molecule m);
	//bool actuallyContains(Molecule m);
};
template<typename T>
using  Matrix3D = std::vector < std::vector< std::vector<T> > >;

class Space
{
private:
	void generateCoordinates();
	void generateSpeeds();
	void initializeUnderspaces();
	void toUnderspaces();
public:
	//QReadWriteLock mutex;
	//TODO incapsulate this fields:
	mutable QMutex mutex;
	int width, height, depth;
	double averageV;
	double maxV, minV, deltaV;
	//===============================
	Space(int width, int height, int n);
	Space& operator=(const Space&& s);
	std::vector<Molecule> molecules;
	Matrix3D<Underspace> underspaces;
	public slots:
	//void saveCoordinates();
	void saveCoordinatesAndSpeeds(const QString& filename);
	//void saveAll();
	//void loadStateC(const QString& filename);
	void loadStateCS(const QString& filename);	//load coordinates and speeds (CS)
	//void loadStateALL(const QString& filename);
};