#include "classes.h"


//Ne
const double Molecule::r = Angstrom;	//m
const double Molecule::m = 20.1797 * AtomicMassUnit;	//kg
const double Molecule::sigma = 2.74 * Angstrom;		//m
const double Molecule::epsilon = 10.0 * 36.2 * Boltzmann;	//J




Space::Space(int width, int height, int n)
	:width(width), height(height)
{
	srand(time(nullptr));
	molecules.resize(n);
	for (auto &i : molecules) {
		Molecule m;
		m.x = rand() % width;
		m.x *= Angstrom;
		m.y = rand() % height;
		m.y *= Angstrom;
		double v = 300 + (rand() % 200 - 100);
		double alphaDeg = rand() % 360;
		double alpha = alphaDeg / 360.0 * 2 * pi;
		m.vx = v * cos(alpha);
		m.vy = v * sin(alpha);
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
			double dx = space.molecules[i].x - space.molecules[j].x;
			double dy = space.molecules[i].y - space.molecules[j].y;
			double hip = sqrt(dx*dx + dy*dy);
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
	emit stateChanged();
}

double Calculator::Force(Molecule &m1, Molecule &m2)
{
	double r = sqrt(m1.x*m1.x + m1.y*m2.y);
	double U = pow(Molecule::sigma / r, 12) - pow(Molecule::sigma / r, 6);
	U *= 4 * Molecule::epsilon;
	return - U / r;
}

PaintWidget::PaintWidget(Space *space, QWidget *parent /*= 0*/)
	:QWidget(parent), space(space)
{

}

void PaintWidget::paintEvent(QPaintEvent *)
{
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
	}
}
