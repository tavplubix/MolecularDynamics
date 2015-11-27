#pragma once
#include "Molecule.h"
#include <vector>

typedef std::vector<Molecule> MoloculesList;

extern void calculateNewForces_GPU(MoloculesList &molecules1, MoloculesList &molecules2);




