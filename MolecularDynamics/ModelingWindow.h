#pragma once
#include "ui_ModelingWindow.h"
#include "PaintWidget.h"
#include <QWindow>
#include <functional>

class ModelingWindow : public QWindow
{
	Q_OBJECT

	Ui_ModelingWindow uimw;

public:
	ModelingWindow();

};
