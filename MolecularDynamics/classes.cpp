#include "classes.h"


//Ne
const double Molecule::radius = Angstrom;	//m
const double Molecule::m = 20.1797 * AtomicMassUnit;	//kg
const double Molecule::sigma = 2.74 * Angstrom;		//m
const double Molecule::epsilon = 1.0 * 36.2 * Boltzmann;	//J



Space::Space(int width, int height, int n)
	:width(width), height(height)
{
	std::srand(std::time(nullptr));
	molecules.resize(n);

	for (auto &i : molecules) {
		//Molecule m;
		i.r.x = std::rand() % width;
		i.r.x *= Angstrom;
		i.r.y = std::rand() % height;
		i.r.y *= Angstrom;
		double v = 300 + (std::rand() % 200 - 100);
		double alphaDeg = std::rand() % 360;
		double alpha = alphaDeg / 360.0 * 2 * pi;
		i.v.x = v * std::cos(alpha);
		i.v.y = v * std::sin(alpha);
	}
}

void Calculator::oneStep()
{
	Space &space = *this->space;
	//forces
	
	for (int i = 0; i < space.molecules.size(); ++i) {
		space.molecules[i].F = Vector();
		for (int j = 0; j < space.molecules.size(); ++j) {
			if (i == j) continue;
			Vector dF = Force(space.molecules[i], space.molecules[j]);
			//if (!std::isfinite(dF)) throw std::exception("infinite force");
			//if (std::isnan(dF)) throw std::exception("dF == NaN");
			space.molecules[i].F += dF;
		}
	}
	
	for (auto &i : space.molecules) {
		//speeds
 		Vector v = i.v + (i.F *(1.0 / i.m)) * dt;
		if (std::abs(v - i.v) > 1000.0)
			//throw std::exception("strange speed");
			qDebug() << "strange speed, F = " + QString::number(i.F);
		//cordinates
		i.r += (i.v + v) * 0.5 * dt;
		i.v = v;
		if (i.r.x <= 0 || space.width * Angstrom <= i.r.x) {
			i.v.x = -v.x;
		}
		if (i.r.y <= 0 || space.height * Angstrom <= i.r.y) {
			i.v.y = -v.y;
		}
	}
	averageSpeed();
	//emit stateChanged();
}

Vector Calculator::Force(Molecule &m1, Molecule &m2)
{
	Vector dr = m2.r - m1.r;
	//double U = std::pow(Molecule::sigma / r, 12) - std::pow(Molecule::sigma / r, 6);
	double U = pow(Molecule::sigma / dr, 12) - pow(Molecule::sigma / dr, 6);
	U *= 4 * Molecule::epsilon;
	return dr * (-U);
}

Calculator::Calculator(Space *space, QObject *parent /*= 0*/)
	:QObject(parent), space(space)
{
	calculationsRequired = false;
}

double Calculator::pow(double d, int i)
{
	double result = d;
	while (--i)  result *= d;
	return result;
}

void Calculator::start()
{
	calculationsRequired = true;
	//QMetaObject::invokeMethod(this, "modeling");
	modeling();
}

void Calculator::stop()
{
	calculationsRequired = false;
}

void Calculator::modeling()
{
	while (calculationsRequired) {
		space->mutex.lock();
		//space->mutex.lockForWrite();
		qDebug() << "	enter in modeling() cycle";
		for (int i = 0; i < 5; ++i) {
			oneStep();
		}
		qDebug() << "	return from modeling() cycle";
		space->mutex.unlock();	//WARNING мьютекс может не освободиться, если будет выкинуто исключение
	}
}

void Calculator::averageSpeed()
{
	space->averageV = 0;
	for (auto i : space->molecules)
		space->averageV += i.v;
	space->averageV /= space->molecules.size();
}

PaintWidget::PaintWidget(Space *space, QWidget *parent)
	:QWidget(parent), space(space)
{

}

void PaintWidget::paintEvent(QPaintEvent *)
{
	space->mutex.lock();
	//space->mutex.lockForRead();
	qDebug() << "\n	enter in paintEvent()\n";
	static std::vector<int> oldx, oldy;
	QPainter painter(this);
	painter.setPen(Qt::SolidLine);
	painter.setPen(Qt::red);
	painter.drawRect(0, 0, space->width * zoom, space->height * zoom);
	painter.setPen(Qt::green);
	for (auto i : space->molecules) {
		int x = i.r.x / Angstrom;
		int y = i.r.y / Angstrom;
		int r = 1 + 6 * i.radius / Angstrom;
		painter.drawEllipse(x * zoom, y * zoom, r * zoom, r * zoom);
		oldx.push_back(x);
		oldy.push_back(y);
	}
	painter.setPen(Qt::blue);
	for (int i = 0; i < oldx.size(); ++i) {
		painter.drawPoint(oldx[i], oldy[i]);
	}
	painter.setPen(Qt::black);
	painter.drawText(5, space->height * zoom + hIndent, QString("Average speed: ") + QString::number(space->averageV, 'f', 3));
	qDebug() << "\n	return from paintEvent()\n";
	space->mutex.unlock();
}
