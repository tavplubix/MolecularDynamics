#include "classes.h"
#include "Space.h"


//Ne
const double Molecule::radius = Angstrom;	//m
const double Molecule::m = 20.1797 * AtomicMassUnit;	//kg
const double Molecule::sigma = 2.74 * Angstrom;		//m
const double Molecule::epsilon = 1.0 * 36.2 * Boltzmann;	//J





PaintWidget::PaintWidget(Space *space, QWidget *parent)
	:QWidget(parent), space(space)
{

}

void PaintWidget::paintEvent(QPaintEvent *)
{
	space->mutex.lock();
		//space->mutex.lockForRead();
		qDebug() << "\n	enter in paintEvent()\n";
		std::vector<Molecule> copy = space->molecules;
		double averageV = space->averageV;
		double deltaV = space->deltaV;
		double height = space->height;
		double width = space->width;
	space->mutex.unlock();

//	static std::vector<int> oldx, oldy;
	QPainter painter(this);
	painter.setPen(Qt::SolidLine);
	painter.setPen(Qt::red);
	painter.drawRect(0, 0, width * zoom, height * zoom);
	painter.setPen(Qt::green);
	for (auto i : copy) {
		int x = i.r.x / Angstrom;
		int y = i.r.y / Angstrom;
		int r = 1 + 6 * i.radius / Angstrom;
		painter.drawEllipse(x * zoom, y * zoom, r * zoom, r * zoom);
// 		oldx.push_back(x);
// 		oldy.push_back(y);
	}
	painter.setPen(Qt::blue);
// 	for (int i = 0; i < oldx.size(); ++i) {
// 		painter.drawPoint(oldx[i], oldy[i]);
// 	}
	painter.setPen(Qt::black);
	painter.drawText(5, height * zoom + hIndent, QString("Average speed: ") + QString::number(averageV, 'f', 3));
	painter.drawText(5, height * zoom + 2*hIndent, QString("Delta: ") + QString::number(deltaV, 'f', 3));
	qDebug() << "\n	return from paintEvent()\n";
}

Molecule::Molecule(const Vector &_r, const Vector &_v)
	:r(_r), v(_v)
{

}
