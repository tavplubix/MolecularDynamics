#include "distributiontest.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	DistributionTest w;
	w.show();
	return a.exec();
}
