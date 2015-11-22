#include "Space.h"
#include <QFile>

#include <map>
#include <set>
#include <random>


const Vector Underspace::size = Vector(300 * Angstrom, 300 * Angstrom, 300 * Angstrom);

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
		i.r.z = 0;
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
		i.v.z = 0;
	}
}

void Space::initializeUnderspaces()
{

	int Nx = 1 + width * Angstrom / Underspace::size.x;
	int Ny = 1 + height * Angstrom / Underspace::size.y;
	int Nz = 1 + depth * Angstrom / Underspace::size.z;
	underspaces.resize(Nx);
	for (auto &i : underspaces)
		i.resize(Ny);
	for (auto &i : underspaces) {
		for (auto &j : i) {
			j.resize(Nz);
		}
	}
	for (size_t x = 0; x < Nx; ++x) {
		for (size_t y = 0; y < Ny; ++y) {
			for (size_t z = 0; z < Nz; ++z) {
				underspaces[x][y][z].nx = x;
				underspaces[x][y][z].ny = y;
				underspaces[x][y][z].nz = z;
			}
		}
	}

}

void Space::toUnderspaces()
{
	for (auto &i : underspaces)
		for (auto &j : i)
			for (auto &k : j)
				k.molecules.clear();

	for (auto i : molecules) {
		int nx = i.r.x / Underspace::size.x;
		int ny = i.r.y / Underspace::size.y;
		int nz = i.r.z / Underspace::size.z;
		try {
			underspaces[nx][ny][nz].molecules.push_back(i);
		} catch (...) {
			qDebug() << "Space::toUnderspaces(): Vector: out of range";
			QErrorMessage err;
			err.showMessage("Space::toUnderspaces(): Vector: out of range");
			err.exec();
		}
	}

}

Space::Space(int width, int height, int n)
	:width(width), height(height), depth(1)
{
	maxV = 0;
	minV = std::numeric_limits<double>::infinity();
	deltaV = 0;
	std::srand(std::time(nullptr));
	molecules.resize(n);

	generateCoordinates();
	generateSpeeds();

	initializeUnderspaces();
	toUnderspaces();

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
