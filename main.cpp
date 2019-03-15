#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include <mpi.h>
#include "Points.h"
#include "math.h"
#include "Trianinfo.h"
#include "Products.h"
#include "EFIEmatrixsys.h"
#include "excvec.h"
#include <complex>
#include <time.h>
#include <Eigen/Dense>



#define COMPLEX complex<double>

#define PI           double(3.14159265358979323846)  /* pi */
#define I             COMPLEX (0,1)
#define eta          double(119.9169832*PI)
#define eps         double(2.2204e-16)
#define lam         double(1) // lambda equals to 1m
#define k           COMPLEX ((2 * PI) / lam)

using namespace std;

using namespace Eigen;

typedef Matrix<COMPLEX, Dynamic, Dynamic> MatrixXCPL;
typedef Matrix<COMPLEX,Dynamic, 1>VectorXCPL;


void send_data(vector<COMPLEX> &local_data, int n, int numprocs, int my_rank) {

	int SIZE = n*n;

	MPI_Send(&local_data[0], SIZE, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
}

void receive_data(vector<COMPLEX> &A, int n, int numprocs) {

	int SIZE = n*n;

	vector<COMPLEX> temp(SIZE);


for (int rank = 1; rank < numprocs; ++rank) {


	MPI_Recv(&temp[0], SIZE, MPI_DOUBLE_COMPLEX, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


 for (int j1 = 0; j1 < n; ++j1){


	for (int j2 = 0; j2 != n; ++j2) {


		A[j1*n + j2] += temp[j1*n + j2];


	}

}

	}
}



int main(int args, char *argv[]) {
	int my_rank, numprocs;

	// MPI initialization
	MPI_Init(&args, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


	vector<double>coord;
	vector<int> topol;
	vector<int> trian;

	vector<int>sizes(3);
	double num1 = 0.0;
	int num2 = 0;
	int num3 = 0;

       clock_t start, end;
       double cpu_time_used;



	if(my_rank==0){

		ifstream file1("coord.txt");
		ifstream file2("topol.txt");
		ifstream file3("trian.txt");

		while (file1 >> num1) {
			coord.emplace_back(num1);
		}

		while (file2 >> num2) {
			topol.emplace_back(num2);
		}

		while (file3 >> num3) {
			trian.emplace_back(num3);
		}

		int sizeC = coord.size();
		int sizeTO = topol.size();
		int sizeTR = trian.size();

      sizes = {sizeC,sizeTO,sizeTR};

	}

  MPI_Bcast(&sizes[0], sizes.size(), MPI_INT, 0, MPI_COMM_WORLD);

	   coord.resize(sizes[0]);
	   topol.resize(sizes[1]);
	   trian.resize(sizes[2]);



  MPI_Bcast(&coord[0], coord.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&topol[0], topol.size(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&trian[0], trian.size(), MPI_INT, 0, MPI_COMM_WORLD);



	for (int i = 0; i < topol.size(); ++i) {

		topol[i] = topol[i] - 1;
	}

	Mesh mesh(coord, topol, trian);

	unsigned int Nt = (topol.size()) / 3;  //number of triangles



	vector<Trianinfo> Triangles;
	int maxele = 1;

	for (int i = 0; i < trian.size(); ++i) {

		if (abs(trian[i]) >= maxele)
			maxele = abs(trian[i]);
	}


	for (int i = 0; i < Nt; ++i) {

		int* vertind = mesh.getNOvertex(i);

		double* vertcoord1 = mesh.getCoord(vertind[0]);
		double* vertcoord2 = mesh.getCoord(vertind[1]);
		double* vertcoord3 = mesh.getCoord(vertind[2]);

		Trianinfo trian(vertcoord1, vertcoord2, vertcoord3);
		Triangles.push_back(trian);

	}



	Points points;

	int SIZE = maxele*maxele;
	vector<COMPLEX> Aglobal(SIZE);
	vector<COMPLEX> Alocal(SIZE);

	fill(Aglobal.begin(), Aglobal.end(), COMPLEX(0));
	fill(Alocal.begin(), Alocal.end(), COMPLEX(0));

	if (my_rank != 0) {

		EFIE::assemble_system_matrix(Alocal, mesh, Triangles, points, Nt, maxele, numprocs, my_rank, k);

		send_data(Alocal, maxele, numprocs, my_rank);
	}
	else {

	start = clock();

	receive_data(Aglobal, maxele, numprocs);

	end = clock();
	
	vector<COMPLEX> E(maxele);
	fill(E.begin(), E.end(), COMPLEX(0));
	
	EFIE::excvecE::assemble_exic_vector(E, mesh, Triangles, points, Nt, k);
	
	MatrixXCPL A(maxele, maxele);
	VectorXCPL C(maxele);
	VectorXCPL B(maxele);

   COMPLEX zbroj = COMPLEX(0);

   for (int i = 0; i < maxele; ++i) {
   
  			 C(i)= E[i];
   
		for (int j = 0; j < maxele; ++j) {

			A(i,j) = Aglobal[i*maxele + j];
		}
	}
	
	
	B = A.colPivHouseholderQr().solve(C); //solve system matrix

	cout << B.sum() << endl;

	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        cout << cpu_time_used << endl;

	}

	MPI_Finalize();

	return 0;
}
