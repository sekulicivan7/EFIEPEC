#pragma once
#ifndef _EFIEMATRIXSYS_HEADER_
#define _EFIEMATRIXSYS_HEADER_
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "Trianinfo.h"
#include <complex>

#define COMPLEX complex<double>



namespace EFIE{

void assemble_system_matrix(vector<COMPLEX> &Alocal, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, unsigned int Nt, unsigned int maxele, int numprocs, int rank, COMPLEX k);

}

#endif
