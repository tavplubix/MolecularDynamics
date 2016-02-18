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

#define WALL


const Vector Underspace::size = Vector(3 * Molecule::sigma, 3 * Molecule::sigma, 3 * Molecule::sigma);
//const Vector Underspace::size = Vector(100 * Angstrom, 100 * Angstrom, 100 * Angstrom);

void Space::generateCoordinates()
{
	size_t size = molecules.size();
	const double minDistance = 1.0 * Molecule::sigma;
	double distance = pow((width * height *depth) / double(numberOfMolecules), 1.0/3.0) * Angstrom;
	if (distance < minDistance) {
		qDebug() << "WARNING: in Space::generateCoordinates(): distance < minDistance";
		throw 0;
	}
	int NX = 1 + width*Angstrom / distance;
	int NY = 1 + height*Angstrom / distance;
	int NZ = 1 + depth*Angstrom / distance;
	//std::map<int, std::set<int>> used;		//CRUTCH
	int nx = 0, ny = 0, nz = 0;
	int n = 1;
	for (auto &m : molecules) {
		m.r.x = nx * distance;
		m.r.y = ny * distance;
		m.r.z = nz * distance;
		m.type = 1;
		m.id = n++;

		++nx;
		if (nx >= NX) {
			nx = 0;
			++ny;
			if (ny >= NY) {
				ny = 0;
				++nz;
				if (nz >= NZ)
					throw 0;
			}
		}
	}
}

void Space::generateSpeeds()
{
	std::random_device rd;
	std::default_random_engine generator(rd());
	double averageSpeed = 300.0;
	double sigma = averageSpeed / std::sqrt(3.0);
	std::normal_distribution<double> normal(0, sigma);
	for (auto &i : molecules) {
		i.v.x = normal(generator);
		i.v.y = normal(generator);
		i.v.z = normal(generator);
	}
}


void Space::generate2DWall()
{
	//hardcoded settings
	int xshift = 100;
	int yshift = 40;
	int zshift = 10;
	int NX = 5;
	int NY = 120;
	int NZ = 1;

	double distance = 1.11 * Molecule::sigma;
	double xStep = distance * sqrt(3.0) * 0.5;
	double yStep = distance;

	//initialize random generator with normal distribution
	std::random_device rd;
	std::default_random_engine generator(rd());
	double averageSpeed = 30.0;
	double sigma = averageSpeed / std::sqrt(3.0);
	std::normal_distribution<double> normal(0, sigma);

	for (int nx = 0; nx < NX; ++nx) {
		for (int ny = 0; ny < NY; ++ny) {
			for (int nz = 0; nz < NZ; ++nz) {
				Molecule m;
				//set positions
				m.r.x = xshift*Angstrom + nx*xStep;
				m.r.y = yshift*Angstrom + ny*yStep;
				if (nx % 2 == 1)
					m.r.y += 0.5 * distance;
				m.r.z = 0;// zshift*Angstrom + nz*distance;		//WARNING

				//set coordinates
				m.v.x = normal(generator);
				m.v.y = normal(generator);
				m.v.z = 0;// normal(generator);

				m.type = 1;
				m.id = molecules.size() + 1;

				molecules.push_back(m);
			}
		}
	}

}

void Space::generate2DBall()
{
	//hardcoded settings
	const int xshift = 0;
	const int yshift = 190;
	const int zshift = 0;
	const int NX = 5;
	const int NY = 5;
	const int NZ = 1;
	const int xSpeed = 2000;

	double distance = 1.11 * Molecule::sigma;
	double xStep = distance * sqrt(3.0) * 0.5;
	double yStep = distance;

	//initialize random generator with normal distribution
	std::random_device rd;
	std::default_random_engine generator(rd());
	double averageSpeed = 30.0;
	double sigma = averageSpeed / std::sqrt(3.0);
	std::normal_distribution<double> normal(0, sigma);

	for (int nx = 0; nx < NX; ++nx) {
		for (int ny = 0; ny < NY; ++ny) {
			for (int nz = 0; nz < NZ; ++nz) {
				Molecule m;
				//set positions
				m.r.x = xshift*Angstrom + nx*xStep;
				m.r.y = yshift*Angstrom + ny*yStep;
				if (nx % 2 == 1)
					m.r.y += 0.5 * distance;
				m.r.z = 0;// zshift*Angstrom + nz*distance;		//WARNING

				//set coordinates
				m.v.x = normal(generator) + xSpeed;
				m.v.y = normal(generator);
				m.v.z = 0;// normal(generator);

				m.type = 2;
				m.id = molecules.size() + 1;

				molecules.push_back(m);
			}
		}
	}

	
}







void Space::generate2DRectangle(int xshift, int yshift, int xsize, int ysize, int type, int xspeed /*= 0*/, int yspeed /*= 0*/)
{
	//hardcoded settings
	const int NX = xsize;
	const int NY = ysize;
	const int NZ = 1;

	double distance = 1.11 * Molecule::sigma;
	double xStep = distance * sqrt(3.0) * 0.5;
	double yStep = distance;

	//initialize random generator with normal distribution
	std::random_device rd;
	std::default_random_engine generator(rd());
	double averageSpeed = 30.0;
	double sigma = averageSpeed / std::sqrt(3.0);
	std::normal_distribution<double> normal(0, sigma);

	for (int nx = 0; nx < NX; ++nx) {
		for (int ny = 0; ny < NY; ++ny) {
			for (int nz = 0; nz < NZ; ++nz) {
				Molecule m;
				//set positions
				m.r.x = xshift*Angstrom + nx*xStep;
				m.r.y = yshift*Angstrom + ny*yStep;
				if (nx % 2 == 1)
					m.r.y += 0.5 * distance;
				m.r.z = 0;// zshift*Angstrom + nz*distance;		//WARNING

						  //set coordinates
				m.v.x = normal(generator) + xspeed;
				m.v.y = normal(generator) + yspeed;
				m.v.z = 0;// normal(generator);

				m.type = type;
				m.id = molecules.size() + 1;

				molecules.push_back(m);
			}
		}
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
#ifndef WALL
	molecules.resize(n);

	generateCoordinates();
	generateSpeeds();
#else
	//generate2DWall();
	//generate2DBall();
	generate2DRectangle(300, 120, 190, 190, 1, 0, 0);
	generate2DRectangle(10, 370, 100, 11, 2, 1000, 0);
	numberOfMolecules = molecules.size();
#endif

	initializeUnderspaces();
	toUnderspaces();
	molecules.clear();
	//molecules.shrink_to_fit();

	//TODO: full VMD format support
	if (n != 0)
	{
		char buff[40];
		time_t now = time(NULL);
		strftime(buff, 40, "%d_%m_%Y__%H_%M_%S.lammpstrj", localtime(&now));
		trajektoryFile.setFileName(buff);
		trajektoryFile.open(QIODevice::Append);
	}

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

void Space::saveTrajektory()
{
	/*TODO other output formats, now only VMD trajektory*/
	int _i = 1;
	QTextStream out(&trajektoryFile);
	out.setCodec("UTF-8");
	//out << "abc";
	
	out << "ITEM: TIMESTEP\n" << trajektoryTime+10 << "\nITEM: NUMBER OF ATOMS\n" << numberOfMolecules << "\n" 
		<< "ITEM: BOX BOUNDS pp pp pp\n0 50.25\n0 50.25\n0 50.25\nITEM: ATOMS type id xs ys zs \n";
	//ITEM : ATOMS type id xs ys zs
	//1 1 0.306349 0.483313 0.084611
	//TODO delete kostil
	forAllM(molecule, underspaces) {
		out << molecule.type << " " << molecule.id << " " << molecule.r.x /** std::pow(10,8)*/ << " " << molecule.r.y /** std::pow(10, 8)*/ << " " << molecule.r.z /** std::pow(10, 8)*/ << "\n";
		_i++;
	}
	out.flush();
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
	//allocate memory for all data
	size_t wholeSize = 1 * sizeof(CUDASpace) + Nx*Ny*Nz*sizeof(CUDAUnderspace) + numberOfMolecules*sizeof(CUDAMolecule);
	byte *p = new byte[wholeSize];		
	auto cs =  reinterpret_cast<CUDASpace*>(p);
	cs->width = width;
	cs->height = height;
	cs->depth = depth;
	cs->Nx = Nx;
	cs->Ny = Ny;
	cs->Nz = Nz;
	cs->numberOfAllMolecules = numberOfMolecules;
	p += 1*sizeof(CUDASpace);
	auto cus = reinterpret_cast<CUDAUnderspace*>(p);
	cs->underspacesShift = GET_SHIFT(cs, p);
	p += Nx*Ny*Nz*sizeof(CUDAUnderspace);
	for (size_t i = 0; i < cs->Nx; ++i) {
		for (size_t j = 0; j < cs->Ny; ++j) {
			for (size_t k = 0; k < cs->Nz; ++k) {
				auto placeForMolecules = reinterpret_cast<CUDAMolecule*>(p);
				auto currentUnderspace = cus + LINEAR(cs, i, j, k);
				underspaces[i][j][k].toCUDA(currentUnderspace, placeForMolecules);
				p += underspaces[i][j][k].molecules.size() * sizeof(CUDAMolecule);
			}
		}
	}
	return cs;
}


void Space::fromCuda(CUDASpace *cs)
{
	auto cus = GET_POINTER(CUDAUnderspace, cs, cs->underspacesShift);
	for (size_t i = 0; i < cs->Nx; ++i) {
		for (size_t j = 0; j < cs->Ny; ++j) {
			for (size_t k = 0; k < cs->Nz; ++k) {
				underspaces[i][j][k].fromCUDA(cus + LINEAR(cs, i, j, k));
			}
		}
	}
}

void Underspace::toCUDA(CUDAUnderspace *cus, CUDAMolecule *placeForMolecules) const
{
	cus->numberOfMolecules = molecules.size();
	cus->moleculesShift = GET_SHIFT(cus, placeForMolecules);
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		molecules[i].oldr	.toCUDA(placeForMolecules[i].oldr);
		molecules[i].r		.toCUDA(placeForMolecules[i].r);
		molecules[i].oldF	.toCUDA(placeForMolecules[i].oldF);
		molecules[i].F		.toCUDA(placeForMolecules[i].F);
		molecules[i].newF	.toCUDA(placeForMolecules[i].newF);
		molecules[i].v		.toCUDA(placeForMolecules[i].v);
		placeForMolecules[i].type = molecules[i].type;
		placeForMolecules[i].id = molecules[i].id;
	}
}




void Underspace::fromCUDA(CUDAUnderspace *cus)
{
	auto cm = GET_POINTER(CUDAMolecule, cus, cus->moleculesShift);
	molecules.resize(cus->numberOfMolecules);
	for (size_t i = 0; i < cus->numberOfMolecules; ++i) {
		molecules[i].oldr	.fromCUDA(cm[i].oldr);
		molecules[i].r		.fromCUDA(cm[i].r);
		molecules[i].oldF	.fromCUDA(cm[i].oldF);
		molecules[i].F		.fromCUDA(cm[i].F);
		molecules[i].newF	.fromCUDA(cm[i].newF);
		molecules[i].v		.fromCUDA(cm[i].v);
		molecules[i].type = cm[i].type;
		molecules[i].id = cm[i].id;
	}
}


