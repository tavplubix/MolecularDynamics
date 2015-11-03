#include "Space.h"
#include <QFile>




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
		i.F = Vector();
		//i.oldr = i.r;
	}
}

void Space::saveCoordinatesAndSpeeds(const QString& filename)
{
	mutex.lock();
	QFile file(filename);
	file.open(QIODevice::WriteOnly);
	file.write(QByteArray::number((unsigned long long)molecules.size()));
	file.write("\n");
	for (auto molecule : molecules) {
		file.write((char*)&molecule.r.x, sizeof(molecule.r.x));
		file.write((char*)&molecule.r.y, sizeof(molecule.r.y));
		file.write((char*)&molecule.r.z, sizeof(molecule.r.z));

		file.write((char*)&molecule.v.x, sizeof(molecule.v.x));
		file.write((char*)&molecule.v.y, sizeof(molecule.v.y));
		file.write((char*)&molecule.v.z, sizeof(molecule.v.z));
	}
	mutex.unlock();
}

void Space::loadStateCS(const QString& filename)
{
	mutex.lock();
	qDebug() << "			enter in load()";
	QFile file(filename);
	file.open(QIODevice::ReadOnly);
	molecules.clear();
	QByteArray arrSize = file.readLine();
	arrSize.chop(1);
	molecules.resize(arrSize.toULongLong());
	for (auto &molecule : molecules) {
		const size_t size = 6 * sizeof(double);		//WARNING
		char buf[size];
		size_t actuallySize = file.read(buf, size);
		if (actuallySize != size && actuallySize != 0) throw std::exception("invalid input file");
		double *t = reinterpret_cast<double*>(buf);
		molecule.r = Vector(t[0], t[1], t[2]);
		molecule.v = Vector(t[3], t[4], t[5]);
	}
	qDebug() << "			return from load()";
	mutex.unlock();
}
