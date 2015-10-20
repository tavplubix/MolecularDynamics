#include "maingui.h"
#include "classes.h"

MainGui::MainGui(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	space = new Space(300, 50, 2);
	pw = new PaintWidget(space, this);
	ui.verticalLayout->addWidget(pw);
	ui.centralWidget->setLayout(ui.verticalLayout);
	calculator = new Calculator();
	connect(calculator, &Calculator::stateChanged, pw, static_cast<void(PaintWidget::*)(void)>(&PaintWidget::update));
	timer = new QTimer(this);
	timer->setInterval(1);
	connect(timer, &QTimer::timeout, calculator, [=](){ calculator->oneStep(*space); });
	timer->start();
}

MainGui::~MainGui()
{
	delete space;
}
