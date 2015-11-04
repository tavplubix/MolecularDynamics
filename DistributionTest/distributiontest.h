#ifndef DISTRIBUTIONTEST_H
#define DISTRIBUTIONTEST_H

#include <QtWidgets/QMainWindow>
#include "ui_distributiontest.h"

class DistributionTest : public QMainWindow
{
	Q_OBJECT

public:
	DistributionTest(QWidget *parent = 0);
	~DistributionTest();

private:
	Ui::DistributionTestClass ui;
};








#endif // DISTRIBUTIONTEST_H
