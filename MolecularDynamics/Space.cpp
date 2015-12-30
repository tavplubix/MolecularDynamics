#include "Space.h"
#include <QFile>
#include <QDebug>
#include <QErrorMessage>
#include <QObject>
#include <QMutex>

#include <map>
#include <set>
#include <random>
#include <ctime>


const Vector Underspace::size = Vector(3 * Molecule::sigma, 3 * Molecule::sigma, 3 * Molecule::sigma);
//const Vector Underspace::size = Vector(100 * Angstrom, 100 * Angstrom, 100 * Angstrom);

void Space::generateCoordinates()
{
	size_t size = molecules.size();
	const double minDistance = 1.0 * Molecule::sigma;
	double distance = sqrt((width * Angstrom * height * Angstrom) / double(numberOfMolecules));
	if (distance < minDistance) {
		qDebug() << "WARNING: in Space::generateCoordinates(): distance < minDistance";
		throw 0;
	}
	int NX = 1 + width*Angstrom / distance;
	int NY = 1 + height*Angstrom / distance;
	//std::map<int, std::set<int>> used;		//CRUTCH
	int nx = 0, ny = 0;
	for (auto &m : molecules) {
		m.r.x = nx * distance;
		m.r.y = ny * distance;
		m.r.z = 0;

		++nx;
		if (nx >= NX) {
			nx = 0;
			++ny;
		}
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

	Nx = 1 + width * Angstrom / Underspace::size.x;
	Ny = 1 + height * Angstrom / Underspace::size.y;
	Nz = 1 + depth * Angstrom / Underspace::size.z;
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
		int nx = int(i.r.x / Underspace::size.x);
		int ny = int(i.r.y / Underspace::size.y);
		int nz = int(i.r.z / Underspace::size.z);
		if (nx < 0) nx = 0; if (Nx <= nx) nx = Nx - 1;
		if (ny < 0) ny = 0; if (Ny <= ny) ny = Ny - 1;
		if (nz < 0) nz = 0; if (Nz <= nz) nz = Nz - 1;
		//try {
			underspaces[nx][ny][nz].molecules.push_back(i);
		//} catch (...) {
			//qDebug() << "Space::toUnderspaces(): Vector: out of range";
			//QErrorMessage err;
			//err.showMessage("Space::toUnderspaces(): Vector: out of range");
			//err.exec();
		//}
	}

}

Space::Space(int width, int height, int n)
	:width(width), height(height), depth(1), numberOfMolecules(n)
{
	maxV = 0;
	minV = std::numeric_limits<double>::infinity();
	deltaV = 0;
	time_s = 0;
	iterations = 0;
	std::srand(std::time(nullptr));
	molecules.resize(n);

	generateCoordinates();
	generateSpeeds();

	initializeUnderspaces();
	toUnderspaces();
	molecules.clear();
	//molecules.shrink_to_fit();

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
	underspaces = std::move(s.underspaces);
	Nx = s.Nx;
	Ny = s.Ny;
	Nz = s.Nz;
	time_s = s.time_s;
	iterations = s.iterations;
	width = s.width;
	height = s.height;
	averageV = s.averageV;
	minV = s.minV;
	maxV = s.maxV;
	deltaV = s.deltaV;
	numberOfMolecules = s.numberOfMolecules;
	mutex.unlock();
	s.mutex.unlock();
	return *this;
}

void Space::saveCoordinatesAndSpeeds(const QString& filename)
{
	mutex.lock();
	QFile file(filename);
	file.open(QIODevice::WriteOnly);
	file.write(QByteArray::number((unsigned long long)numberOfMolecules));
	file.write("\n");
	forAllM (molecule, underspaces) {
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
	numberOfMolecules = arrSize.toULongLong();
	molecules.resize(numberOfMolecules);
	for (auto &molecule : molecules) {
		const size_t size = 6 * sizeof(double);		//WARNING
		char buf[size];
		size_t actuallySize = file.read(buf, size);
		if (actuallySize != size && actuallySize != 0) throw std::exception("invalid input file");
		double *t = reinterpret_cast<double*>(buf);
		molecule.r = Vector(t[0], t[1], t[2]);
		molecule.v = Vector(t[3], t[4], t[5]);
	}
	toUnderspaces();
	mutex.unlock();
}


CUDASpace* Space::toCUDA() const
{
	CUDASpace *cs = new CUDASpace;
	cs->width = width;
	cs->height = height;
	cs->Nx = Nx;
	cs->Ny = Ny;
	cs->Nz = Nz;
	cs->underspaces = new CUDAUnderspace[Nx*Ny*Nz];
	for (size_t i = 0; i < cs->Nx; ++i) {
		for (size_t j = 0; j < cs->Ny; ++j) {
			for (size_t k = 0; k < cs->Nz; ++k) {
				underspaces[i][j][k].toCUDA(cs->underspaces + LINEAR(cs, i, j, k));
			}
		}
	}
	return cs;
}


void Space::fromCuda(CUDASpace *cs)
{
	for (size_t i = 0; i < cs->Nx; ++i) {
		for (size_t j = 0; j < cs->Ny; ++j) {
			for (size_t k = 0; k < cs->Nz; ++k) {
				underspaces[i][j][k].fromCUDA(cs->underspaces + LINEAR(cs, i, j, k));
			}
		}
	}
}

void Underspace::toCUDA(CUDAUnderspace *cus) const
{
	cus->numberOfMolecules = molecules.size();
	cus->molecules = new CUDAMolecule[cus->numberOfMolecules];
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		molecules[i].oldr	.toCUDA(cus->molecules[i].oldr);
		molecules[i].r		.toCUDA(cus->molecules[i].r);
		molecules[i].oldF	.toCUDA(cus->molecules[i].oldF);
		molecules[i].F		.toCUDA(cus->molecules[i].F);
		molecules[i].newF	.toCUDA(cus->molecules[i].newF);
		molecules[i].v		.toCUDA(cus->molecules[i].v);
	}
}

void Underspace::fromCUDA(CUDAUnderspace *cus)
{
	molecules.resize(cus->numberOfMolecules);
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		molecules[i].oldr	.fromCUDA(cus->molecules[i].oldr);
		molecules[i].r		.fromCUDA(cus->molecules[i].r);
		molecules[i].oldF	.fromCUDA(cus->molecules[i].oldF);
		molecules[i].F		.fromCUDA(cus->molecules[i].F);
		molecules[i].newF	.fromCUDA(cus->molecules[i].newF);
		molecules[i].v		.fromCUDA(cus->molecules[i].v);
	}
}


