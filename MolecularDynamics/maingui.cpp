#include "maingui.h"
#include "classes.h"
#include "Calculator.h"
#include "Space.h"
#include <QFileDialog>

MainGui::MainGui(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	space = new Space(300, 160, 50);
	calculator = new Calculator(space);

	//Create PaintWidget for space
	pw = new PaintWidget(space, this);
	ui.mainLayout->addWidget(pw, ui.mainLayout->rowCount(), 0, 1, -1);
	ui.centralWidget->setLayout(ui.mainLayout);

	//Set actions for buttons
	ui.saveButton->setHidden(true);
	connect(ui.restartButton, &QPushButton::clicked, [&](){ 
		if (calculator->calculationsRequired == false) ui.pauseContinueButton->click();		//CRUTCH
		calculator->calculationsRequired = false; 
		//TODO recreate space
		QMetaObject::invokeMethod(calculator, "start"); 
	});
	connect(ui.pauseContinueButton, &QPushButton::clicked, [&](){
		if (calculator->calculationsRequired) {		//pause
			calculator->calculationsRequired = false;
			ui.saveButton->setHidden(false);
			ui.pauseContinueButton->setText("Continue");
		}
		else {		//continue
			calculator->calculationsRequired = true;
			ui.saveButton->setHidden(true);
			ui.pauseContinueButton->setText("Pause");
		}
	});
	connect(ui.saveButton, &QPushButton::clicked, [&]() {
		QString filename = QFileDialog::getSaveFileName(this, "Save state", "./../state.mdcs", "MD (*.mdcs)");
		space->saveCoordinatesAndSpeeds(filename);
	});
	connect(ui.loadButton, &QPushButton::clicked, [&](){
		QString filename = QFileDialog::getOpenFileName(this, "Load state", "./../state.mdcs", "MD (*.mdcs)");
		if (calculator->calculationsRequired == false) ui.pauseContinueButton->click();		//CRUTCH
		calculator->calculationsRequired = false;
		space->loadStateCS(filename);
		QMetaObject::invokeMethod(calculator, "start");
	});
	


	//Use timer to update gui every 100 ms
	timer = new QTimer(this);
	timer->setInterval(100);
	connect(timer, &QTimer::timeout, pw, static_cast<void(PaintWidget::*)(void)>(&PaintWidget::update));
	timer->start();

	//Start calculations in another thread
	calculationsThread = new QThread(this);
	calculator->moveToThread(calculationsThread);
	calculationsThread->start();
	QMetaObject::invokeMethod(calculator, "start");
}

MainGui::~MainGui()
{
	calculator->calculationsRequired = false;
	QEventLoop loop;
	connect(calculator, &QObject::destroyed, &loop, &QEventLoop::quit);
	QMetaObject::invokeMethod(calculator, "deleteLater");
	loop.exec();
	delete space;
}
