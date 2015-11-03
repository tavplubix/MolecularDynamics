#pragma once

#include <QObject>
#include <QMutex>
#include "classes.h"
#include <vector>


class Space
{
public:
	//QReadWriteLock mutex;
	QMutex mutex;
	int width, height;
	double averageV;
	Space(int width, int height, int n);
	std::vector<Molecule> molecules;
	public slots:
	//void saveCoordinates();
	void saveCoordinatesAndSpeeds(const QString& filename);
	//void saveAll();
	//void loadStateC(const QString& filename);
	void loadStateCS(const QString& filename);	//load coordinates and speeds (CS)
	//void loadStateALL(const QString& filename);
};