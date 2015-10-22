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
		double x = std::rand() % width;
		x *= Angstrom;
		double y = std::rand() % height;
		y *= Angstrom;
		i.r = Vector(x, y, 0);
		double v = 300 + (std::rand() % 200 - 100);
		double alphaDeg = std::rand() % 360;
		double alpha = alphaDeg / 360.0 * 2 * pi;
		double vy = v * std::sin(alpha);
		double vx = v * std::cos(alpha);
		i.v = Vector(vx, vy, 0);
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
			if (!std::isfinite(dF)) throw std::exception("infinite force");
			if (std::isnan(dF)) throw std::exception("dF == NaN");
			space.molecules[i].F += dF;
		}
	}
	
	for (auto &i : space.molecules) {
		//speeds
		Vector v = i.v + (i.F / i.m) * dt;
		if (std::abs(v - i.v) > 1000.0 )
			qDebug() << "strange speed, Fx = " + QString::number(i.F.x) + ", Fy = " + QString::number(i.F.y);
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
	Vector r = m2.r - m1.r;
	double U = pow(Molecule::sigma / r, 14) - pow(Molecule::sigma / r, 8);
	U *= 4 * Molecule::epsilon / pow(Molecule::sigma, 2);
	return - U * r;
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
	//static std::vector<int> oldx, oldy;
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
		//oldx.push_back(x);
		//oldy.push_back(y);
	}
	painter.setPen(Qt::blue);
	//for (int i = 0; i < oldx.size(); ++i) {
	//	painter.drawPoint(oldx[i], oldy[i]);
	//}
	painter.setPen(Qt::black);
	painter.drawText(5, space->height * zoom + hIndent, QString("Average speed: ") + QString::number(space->averageV, 'f', 3));
	qDebug() << "\n	return from paintEvent()\n";
	space->mutex.unlock();
}

Molecule::Molecule(const Vector &_r, const Vector &_v)
	:r(_r), v(_v)
{

}
