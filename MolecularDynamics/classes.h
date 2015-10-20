#pragma once
#include <exception>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <QWidget>
#include <QTimer>
#include <QPainter>
#include <QObject>
#include <QDebug>
#include <QThread>
#include <QReadWriteLock>
#include <QMutex>


class PaintWidget;
class Molecule;
class Calculator;
class Space;

namespace constants {
	const double Boltzmann = 1.3806488e-23;
	const double Avogadro = 6.022140857e23;
	const double Angstrom = 10e-10;
	const double Electronvolt = 1.6021766208e-19;
	const double AtomicMassUnit = 1.660538921e-27;
	const double pi = 3.141592654;
}

using namespace constants;

class PaintWidget : public QWidget
{
	Q_OBJECT
private:
	Space *space;
	const int zoom = 1;
	const int hIndent = 35;
public:
	PaintWidget(Space *space, QWidget *parent = 0);
	void paintEvent(QPaintEvent *);
};

class Molecule
{
public:
	double x, y;
	double vx, vy;
	double Fx, Fy;
	const static double m, r, sigma, epsilon;
};

class Space
{
public:
	//QReadWriteLock mutex;
	QMutex mutex;
	int width, height;
	double averageV;
	Space(int width, int height, int n);
	std::vector<Molecule> molecules;
};

class Calculator : public QObject
{
	Q_OBJECT
private:
	Space* space;
	double dt = 10e-16;
	bool calculationsRequired;
	void averageSpeed();
	static double Force(Molecule &m1, Molecule &m2);
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



