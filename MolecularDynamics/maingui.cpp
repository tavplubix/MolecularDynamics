#include "maingui.h"
#include "classes.h"
#include "Calculator.h"

MainGui::MainGui(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	space = new Space(300, 160, 50);
	pw = new PaintWidget(space, this);
	ui.verticalLayout->addWidget(pw);
	ui.centralWidget->setLayout(ui.verticalLayout);
	calculator = new Calculator(space);
	calculationsThread = new QThread(this);
	calculator->moveToThread(calculationsThread);

	//connect(calculator, &Calculator::stateChanged, pw, static_cast<void(PaintWidget::*)(void)>(&PaintWidget::update));
	timer = new QTimer(this);
	timer->setInterval(100);
	connect(timer, &QTimer::timeout, pw, static_cast<void(PaintWidget::*)(void)>(&PaintWidget::update));
	//connect(timer, &QTimer::timeout, calculator, [=](){ calculator->oneStep(*space); });
	timer->start();
	calculationsThread->start();
	QTimer::singleShot(0, calculator, &Calculator::start);
	//calculator->start();
}

MainGui::~MainGui()
{
	delete space;
}
