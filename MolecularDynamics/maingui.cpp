#include "maingui.h"
#include <functional>
#include <QFileDialog>
#include <ctime>


MainGui::MainGui(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	this->resize(1850, 990);
	ui.numberOfMoleculesSpinBox->setValue(3000);
	ui.widthSpinBox->setValue(300);
	ui.heightSpinBox->setValue(430);
	ui.precisionSpinBox->setValue(14);

	space = new Space(ui.widthSpinBox->value(), ui.heightSpinBox->value(), ui.numberOfMoleculesSpinBox->value());
	calculator = new Calculator(space);

	//Create PaintWidget for space
	pw = new PaintWidget(space, this);
	ui.mainLayout->addWidget(pw, ui.mainLayout->rowCount(), 0, 1, -1);
	ui.centralWidget->setLayout(ui.mainLayout);

	setButtons();

	//Precision
	connect(ui.precisionSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&](int precision){
		calculator->set_dt_precision(precision);
	});

	//Zoom
	connect(ui.zoomSpinBox, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [&](double zoom) {
		pw->zoom = zoom;
	});

	//Heating/Freezing
	connect(ui.heatingSpinBox, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [&](double heating) {
		calculator->heating = heating;
	});

	//Borders
	connect(ui.heightSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&](int height){
		space->height = height;		//FIXME resize Space::underspaces
	});
	connect(ui.widthSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&](int width){
		space->width = width;		//FIXME resize Space::underspaces
	});


	//Use timer to update GUI every 100 ms
	timer = new QTimer(this);
	timer->setInterval(100);
	connect(timer, &QTimer::timeout, pw, static_cast<void(PaintWidget::*)(void)>(&PaintWidget::update));
	lastIter = 0;
	lastTime = std::time(nullptr);
	connect(timer, &QTimer::timeout, this, [&](){
		unsigned long long iter = space->iterations, t = time(nullptr);
		this->setWindowTitle(QString::number(calculator->get_dt()));
		ui.timeLabel->setText("Time: " + QString::number(space->time_s));
		ui.iterationsLabel->setText("Iterations: " + QString::number(iter));
		ui.iterPerSecLabel->setText("Iterations per second: " + QString::number((iter - lastIter) / double(t - lastTime)));
		//int dt = t - lastTime;
		//if (dt > 30) {
		//	lastIter = iter;
		//	lastTime = t;
		//}
	});
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

void MainGui::setButtons()
{
	//Set actions for buttons
	//Disable some widgets
	ui.saveButton->setHidden(true);
	ui.numberOfMoleculesSpinBox->setEnabled(false);

	//Stop/Restart Button
	connect(ui.stopStartButton, &QPushButton::clicked, [&](){
		//Stop Button
		if (ui.stopStartButton->text() == "Stop") {
			calculator->calculationsRequired = false;	//stop calculations
			//hide other buttons
			ui.pauseContinueButton->setHidden(true);
			ui.loadButton->setHidden(true);
			ui.saveButton->setHidden(true);
			//remove space
			Space s(0, 0, 0);
			*space = std::move(s);
			
			ui.stopStartButton->setText("Restart");
			ui.numberOfMoleculesSpinBox->setEnabled(true);
		}
		//Restart Button
		else if (ui.stopStartButton->text() == "Restart") {
			//create new space
			int height = ui.heightSpinBox->value();
			int width = ui.widthSpinBox->value();
			int numberOfMolecules = ui.numberOfMoleculesSpinBox->value();
			Space s(width, height, numberOfMolecules);
			*space = std::move(s);
			//disable some widgets
			ui.numberOfMoleculesSpinBox->setEnabled(false);
			//restore other widgets
			ui.pauseContinueButton->setHidden(false);
			ui.loadButton->setHidden(false);
			ui.stopStartButton->setText("Stop");
			ui.pauseContinueButton->setText("Pause");
			calculator->calculationsRequired = true;
		}
		else throw std::exception("something strange has happend in slot connected to stopStartButton");
		lastIter = space->iterations;
		lastTime = std::time(nullptr);
	});

	//Pause/Continue Button
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
		lastIter = space->iterations;
		lastTime = std::time(nullptr);
	});

	//Save Button
	connect(ui.saveButton, &QPushButton::clicked, [&]() {
		QString filename = QFileDialog::getSaveFileName(this, "Save state", "./../state.mdcs", "MD (*.mdcs)");
		space->saveCoordinatesAndSpeeds(filename);
	});

	//Load Button
	connect(ui.loadButton, &QPushButton::clicked, [&](){
		QString filename = QFileDialog::getOpenFileName(this, "Load state", "./../state.mdcs", "MD (*.mdcs)");
		if (calculator->calculationsRequired == false) ui.pauseContinueButton->click();		//CRUTCH
		calculator->calculationsRequired = false;
		space->loadStateCS(filename);
		QMetaObject::invokeMethod(calculator, "start");
	});
}
