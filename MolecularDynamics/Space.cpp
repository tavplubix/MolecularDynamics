#include "Space.h"
#include <QFile>

#include <map>
#include <set>
#include <random>




void Space::generateCoordinates()
{
	size_t size = molecules.size();
	const int minDistance = 2 + int(Molecule::sigma / Angstrom);
	std::map<int, std::set<int>> used;		//CRUTCH
	for (auto &i : molecules) {

		int tmpX = (std::rand() * minDistance) % width;
		int tmpY;
		do {
			tmpY = (std::rand() * minDistance) % height;
		} while (used[tmpX].find(tmpY) != used[tmpX].end());		//WARNING
		used[tmpX].insert(tmpY);
		i.r.x = double(tmpX) * Angstrom;
		i.r.y = double(tmpY) * Angstrom;
	}
}

void Space::generateSpeeds()
{
	std::random_device rd;
	std::default_random_engine generator(rd());
	double averageSpeed = 400.0;
	double sigma = averageSpeed / std::sqrt(3.0);
	std::normal_distribution<double> normal(0, sigma);
	for (auto &i : molecules) {
		i.v.x = normal(generator);
		i.v.y = normal(generator);
		//i.v.z = normal(generator);
	}
}

Space::Space(int width, int height, int n)
	:width(width), height(height)
{
	maxV = 0;
	minV = std::numeric_limits<double>::infinity();
	deltaV = 0;
	std::srand(std::time(nullptr));
	molecules.resize(n);

	generateCoordinates();
	generateSpeeds();

#ifdef DEBUG
	saveCoordinatesAndSpeeds("./../last.log.mdcs");
#endif
}

Space& Space::operator=(const Space&& s)
{
	if (this == &s) return *this;
	s.mutex.lock();
	mutex.lock();
	molecules = std::move(s.molecules);
	width = s.width;
	height = s.height;
	averageV = s.averageV;
	minV = s.minV;
	maxV = s.maxV;
	deltaV = s.deltaV;
	mutex.unlock();
	s.mutex.unlock();
	return *this;
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
	mutex.unlock();
}
