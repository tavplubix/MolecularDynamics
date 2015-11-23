#ifndef MAINGUI_H
#define MAINGUI_H

#include <QtWidgets/QMainWindow>
#include "ui_maingui.h"
#include "classes.h"

class MainGui : public QMainWindow
{
	Q_OBJECT

public:
	MainGui(QWidget *parent = 0);
	~MainGui();

private:
	unsigned long long lastIter, lastTime;
	inline void setButtons();
	QThread *calculationsThread;
	QTimer *timer;
	Calculator *calculator;
	Space *space;
	PaintWidget *pw;
	Ui::MainGuiClass ui;
};

#endif // MAINGUI_H
