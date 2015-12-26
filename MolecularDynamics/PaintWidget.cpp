#include "PaintWidget.h"


PaintWidget::PaintWidget(Space *space, QWidget *parent)
	:QWidget(parent), space(space)
{

}

void PaintWidget::paintEvent(QPaintEvent *)
{
	space->mutex.lock();
	//space->mutex.lockForRead();
#ifdef DEBUG
	qDebug() << "\n	enter in paintEvent()\n";
#endif
	std::vector<Molecule> copy = space->molecules;
	double averageV = space->averageV;
	double deltaV = space->deltaV;
	double height = space->height;
	double width = space->width;
#ifdef DEBUG
	qDebug() << "\n	return from paintEvent()\n";
#endif
	space->mutex.unlock();

	//static std::vector<int> oldx, oldy;
	QPainter painter(this);
	painter.setPen(Qt::SolidLine);
	painter.setPen(Qt::red);
	painter.drawRect(0, 0, width * zoom, height * zoom);
	painter.setPen(Qt::blue);
	painter.setBrush(Qt::green);

	//for (auto i : copy) {
	//	int x = i.r.x / Angstrom;
	//	int y = i.r.y / Angstrom;
	//	int r = 1 + 6 * i.radius / Angstrom;
	//	painter.drawEllipse(x * zoom, y * zoom, r * zoom, r * zoom);
	//		//oldx.push_back(x);
	//		//oldy.push_back(y);
	//}
	//space->toUnderspaces();

	//forAllM(t, space->underspaces) {
	for (auto &t : copy) {
		int x = t.r.x / Angstrom;
		int y = t.r.y / Angstrom;
		if (space->numberOfMolecules <= 0/*5000*/) {
			int r = 1 + 6 * t.radius / Angstrom;
			painter.drawEllipse((x - r / 2) * zoom, (y - r / 2) * zoom, r * zoom, r * zoom);
		}
		else {
			painter.drawEllipse((x - 1)*zoom, (y - 1)*zoom, 3, 3);
			//painter.drawPoint(x * zoom, y * zoom);
		}
		//oldx.push_back(x);
		//oldy.push_back(y);
	}
	//space->mutex.unlock();


	//painter.setPen(Qt::blue);
	//for (int i = 0; i < oldx.size(); ++i) {
	//	painter.drawPoint(oldx[i], oldy[i]);
	//}
	painter.setPen(Qt::black);
	painter.drawText(5, height * zoom + hIndent, QString("Average speed: ") + QString::number(averageV, 'f', 3));
	painter.drawText(5, height * zoom + 2 * hIndent, QString("Delta: ") + QString::number(deltaV, 'f', 3));


}



