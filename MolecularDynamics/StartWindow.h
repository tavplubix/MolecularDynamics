#pragma once
#include <QMainWindow>
#include <functional>
#include "ui_StartWindow.h"

class StartWindow : public QMainWindow
{
	Q_OBJECT

	Ui_StartWindow uisw;

public:
	StartWindow();

signals:
	void createSpaceRequest(int width, int height, int depth = 0);
	void addMoleculesToSpaceRequest(/*type: (gas, liquid, solid); constants: (sigma, epsilon)*/);
	void startModelingRequest(/*algorithm*/);
	void pauseContinueRequest(std::function<void()> callItBackWhenPausedOrContinued);
	void saveRequest(const QString& filename);
	void loadRequest(const QString& filename);
	void stopRequest();
};
