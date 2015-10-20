#include "classes.h"


//Ne
const double Molecule::r = Angstrom;	//m
const double Molecule::m = 20.1797 * AtomicMassUnit;	//kg
const double Molecule::sigma = 2.74 * Angstrom;		//m
const double Molecule::epsilon = 1.0 * 36.2 * Boltzmann;	//J



Space::Space(int width, int height, int n)
	:width(width), height(height)
{
	std::srand(std::time(nullptr));
	molecules.resize(n);

	for (auto &i : molecules) {
		Molecule m;
		m.x = std::rand() % width;
		m.x *= Angstrom;
		m.y = std::rand() % height;
		m.y *= Angstrom;
		double v = 300 + (std::rand() % 200 - 100);
		double alphaDeg = std::rand() % 360;
		double alpha = alphaDeg / 360.0 * 2 * pi;
		m.vx = v * std::cos(alpha);
		m.vy = v * std::sin(alpha);
		i = m;
	}
}



void Calculator::oneStep(Space &space)
{
	//forces
	
	for (int i = 0; i < space.molecules.size(); ++i) {
		space.molecules[i].Fx = 0;
		space.molecules[i].Fy = 0;
		for (int j = 0; j < space.molecules.size(); ++j) {
			if (i == j) continue;
			double dF = Force(space.molecules[i], space.molecules[j]);
			if (!std::isfinite(dF)) throw std::exception("infinite force");
			if (std::isnan(dF)) throw std::exception("dF == NaN");
			double dx = space.molecules[j].x - space.molecules[i].x;
			double dy = space.molecules[j].y - space.molecules[i].y;
			double hip = std::sqrt(dx*dx + dy*dy);
			double dFx = dF * dx / hip;
			double dFy = dF * dy / hip;
			space.molecules[i].Fx += dFx;
			space.molecules[i].Fy += dFy;
		}
	}
	
	for (auto &i : space.molecules) {
		//speeds
		double vx = i.vx + (i.Fx / i.m) * dt;
		double vy = i.vy + (i.Fy / i.m) * dt;
		if (std::abs(vx - i.vx) > 1000.0 || std::abs(vy - i.vy) > 1000.0)
			//throw std::exception("strange speed");
			qDebug() << "strange speed, Fx = " + QString::number(i.Fx) + ", Fy = " + QString::number(i.Fy);
		//cordinates
		i.x += (i.vx + vx) * 0.5 * dt;
		i.y += (i.vy + vy) * 0.5 * dt;
		i.vx = vx;
		i.vy = vy;
		if (i.x <= 0 || space.width * Angstrom <= i.x) {
			i.vx = -vx;
		}
		if (i.y <= 0 || space.height * Angstrom <= i.y) {
			i.vy = -vy;
		}
	}
	average(space);
	emit stateChanged();
}

double Calculator::Force(Molecule &m1, Molecule &m2)
{
	double dx = m2.x - m1.x;
	double dy = m2.y - m1.y;
	double r = std::sqrt(dx*dx + dy*dy);
	if (std::isnan(r)) throw std::exception("r == NaN");
	//double U = std::pow(Molecule::sigma / r, 12) - std::pow(Molecule::sigma / r, 6);
	double U = pow(Molecule::sigma / r, 12) - pow(Molecule::sigma / r, 6);
	U *= 4 * Molecule::epsilon;
	return - U / r;
}

double Calculator::pow(double d, int i)
{
	double result = d;
	while (--i)  result *= d;
	return result;
}

void Calculator::average(Space &space)
{
	space.averageV = 0;
	for (auto i : space.molecules) 
		space.averageV += std::sqrt(i.vx*i.vx + i.vy*i.vy);
	space.averageV /= space.molecules.size();
}

PaintWidget::PaintWidget(Space *space, QWidget *parent)
	:QWidget(parent), space(space)
{

}

void PaintWidget::paintEvent(QPaintEvent *)
{
	static std::vector<int> oldx, oldy;
	QPainter painter(this);
	painter.setPen(Qt::SolidLine);
	painter.setPen(Qt::red);
	painter.drawRect(0, 0, space->width * zoom, space->height * zoom);
	painter.setPen(Qt::green);
	for (auto i : space->molecules) {
		int x = i.x / Angstrom;
		int y = i.y / Angstrom;
		int r = 1 + 6 * i.r / Angstrom;
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
}
