#pragma once
#include <QObject>
#include "Space.h"

class AbstractCalculator : QObject
{
	Q_OBJECT
private:
	Space* space;
	virtual void oneStep() = 0;
	Q_INVOKABLE virtual void modeling();
	double timeStep; //dt
public:
	AbstractCalculator(Space* space, double timeStep, QObject* parent = 0);
	~AbstractCalculator();
	virtual void calculateAverageMoleculesSpeeds() = 0; //Root mean square
	virtual void calculateWholeSystemEnergy() = 0;
public slots:
	virtual void start();
	virtual void stop();
	virtual void pause();
};