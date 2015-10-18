#include "maingui.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MainGui w;
	w.show();
	return a.exec();
}
