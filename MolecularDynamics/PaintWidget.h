#pragma once

#include <QWidget>
#include <QTimer>
#include <QPainter>
#include <QObject>
#include <QDebug>
#include <QThread>
#include <QReadWriteLock>
#include <QMutex>

#include "DeprecatedSpace.h"



class PaintWidget : public QWidget
{
	Q_OBJECT
private:
	DeprecatedSpace *space;
	const int hIndent = 25;

	double maxE = -1000, minE = 1000;		//CRUTCH
public:
	PaintWidget(DeprecatedSpace *space, QWidget *parent = 0);
	void paintEvent(QPaintEvent *);
	double zoom = 1;
};

