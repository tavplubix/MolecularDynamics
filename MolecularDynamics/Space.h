#pragma once

#include <QObject>
#include <QMutex>
#include "classes.h"
#include <vector>


class Space
{
private:
	void generateCoordinates();
	void generateSpeeds();
public:
	//QReadWriteLock mutex;
	//TODO incapsulate this fields:
	mutable QMutex mutex;
	int width, height;
	double averageV;
	double maxV, minV, deltaV;
	//===============================
	Space(int width, int height, int n);
	Space& operator=(const Space&& s);
	std::vector<Molecule> molecules;
	public slots:
	//void saveCoordinates();
	void saveCoordinatesAndSpeeds(const QString& filename);
	//void saveAll();
	//void loadStateC(const QString& filename);
	void loadStateCS(const QString& filename);	//load coordinates and speeds (CS)
	//void loadStateALL(const QString& filename);
};