// Main.cpp - точка входа в программу

#define _CRT_SECURE_NO_WARNINGS

// дл€ PARDISO
#include <omp.h>
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_spblas.h>


// дл€ остального
#include <time.h>
//#include <math.h>
//#include <cstdio>
#include <cstdlib>
#include <cstdint>

// файлы  Ёћ
//#include "Mesh.h"
#include "fem.h"
//#include "geometry.h"

ofstream logfile;
string output_folder = "test";

int Read_Long_From_Txt_File(const char* fname, int* number)
{
	FILE* fp;
	int temp;
	int retcode;

	if ((fp = fopen(fname, "r")) == 0)
	{
		char str[100];
		sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
		return 1;
	}

	retcode = fscanf(fp, "%ld", &temp);
	if (retcode != 1)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\".\n", fname);
		fclose(fp);
	}
	*number = temp;

	fclose(fp);
	return 0;
}

int Read_Bin_File_Of_Double(const char* fname, double* massiv, int n_of_records, int len_of_record)
{
	int temp;
	FILE* fp;

	// открываем файл на чтение
	if ((fp = fopen(fname, "r+b")) == 0)
	{
		char str[100];
		sprintf(str, "Error: Cannot open file \"%s\" for reading.\n", fname);
		return 1;
	}

	temp = fread(massiv, sizeof(double) * len_of_record, n_of_records, fp);
	if (temp != n_of_records)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
		fclose(fp);
		return 1;
	}

	fclose(fp);

	return 0;
}

int Read_Bin_File_Of_Long(const char* fname, int* massiv, int n_of_records, int len_of_record)
{
	int temp;
	FILE* fp;

	// открываем файл на чтение
	if ((fp = fopen(fname, "r+b")) == 0)
	{
		char str[100];
		sprintf(str, "Cannot open file %s.\n", fname);
		return 1;
	}

	// чтение
	temp = fread(massiv, sizeof(int) * len_of_record, n_of_records, fp);
	if (temp != n_of_records)
	{
		char str[100];
		sprintf(str, "Error reading file \"%s\". %ld of %ld records was read.\n", fname, temp, n_of_records);
		fclose(fp);
		return 1;
	}


	fclose(fp);
	return 0;
}

void FromRSFToCSR_Real_1_Sym(int nb, int* ig, int* sz_ia, int* sz_ja)
{
	*sz_ia = nb + 1;
	*sz_ja = ig[nb] + nb;
}

void FromRSFToCSR_Real_2_Sym(int nb, int* ig, int* jg, double* di, double* gg,
	MKL_INT64* ia, MKL_INT64* ja, double* a)
{
	int i, j, k;
	vector<MKL_INT64> adr;

	// подсчитываем число элементов в каждой строчке
	adr.resize(nb, 0);

	for (i = 0; i < nb; i++)
	{
		adr[i] += 1; // диагональ

		// верхний треугольник
		for (j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			k = jg[j];
			adr[k]++;
		}
	}

	// ia
	ia[0] = 0;
	for (i = 0; i < nb; i++)
		ia[i + 1] = ia[i] + adr[i];

	// ja,  a
	for (i = 0; i <= ig[nb] + nb; i++)
		a[i] = 0;

	for (i = 0; i < nb; i++)
		adr[i] = ia[i]; // в какую позицию заносить значение

	// диагональ
	for (i = 0; i < nb; i++)
	{
		ja[adr[i]] = i;
		a[adr[i]] = di[i];
		adr[i]++;
	}

	// верхний треугольник
	for (i = 0; i < nb; i++)
	{
		for (j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			k = jg[j];
			ja[adr[k]] = i;

			a[adr[k]] = gg[j];

			adr[k]++;
		}
	}
}



int WriteResultTxt(const char* fname, double* result_vector, int n_of_records, int len_of_record) {
	FILE* fp;
	if ((fp = fopen(fname, "w")) == 0) {
		printf("ERROR: Cannot open file \"@%s\" for writing.\n", fname);
		return 1;
	}

	printf("Writing %s...", fname);
	for (int i = 0; i < n_of_records; i++) {
		for (int j = 0; j < len_of_record; j++)
			fprintf(fp, "%25.13e\t", result_vector[i * len_of_record + j]);
		fprintf(fp, "\n");
	}
	printf("Done!\n");
	fclose(fp);
	return 0;
}

void solve_pardiso_symm(MKL_INT64 n, MKL_INT64* ia, MKL_INT64* ja, double* a, double* b, double* x) {
	logfile.open("pardiso64.log");
	if (!logfile) {
		cerr << "Cannot open pardiso64.log" << endl;
		exit(1);
	}
	MKL_INT64 mtype = 2;	// 2 - Real symmetric positive defined matrix
							// -2 - Real symmetric indefinited matrix
	MKL_INT64 nrhs = 1;	// Number of right hand sides.
	// Internal solver memory pointer pt,
	// 32-bit: int pt[64]; 64-bit: long int pt[64]
	// or void *pt[64] should be OK on both architectures
	void* pt[64];
	//pt[0] = 0;
	for (int i = 0; i < 64; i++)
		pt[i] = 0;

	// PARDISO control parameters
	MKL_INT64 iparm[64];
	iparm[0] = 0;
	for (int i = 1; i < 64; i++)
		iparm[i] = 0;
	//double dparm[64];

	// Numbers of processors, value of OMP_NUM_THREADS
	
	MKL_INT64 maxfct = 1; // Maximum number of numerical factorizations.
	MKL_INT64 mnum = 1;   // Which factorization to use.
	MKL_INT64 msglvl = 1; // Print statistical information
	MKL_INT64 error = 0;  // Initialize error flag
	MKL_INT64 phase = 13;
	MKL_INT64 info = -100;



	MKL_INT64* perm = new MKL_INT64[n];

	cout << "Pardiso start..." << endl << flush;
	PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, perm, &nrhs,
		iparm, &msglvl, b, x, &info);

	if (info != 0) {
		printf("\nERROR during solution: %d\n", info);
		system("pause");
		exit(3);
	}
	WriteResultTxt((output_folder + "/u.txt").c_str(), x, n, 1);

	logfile.close();
	logfile.clear();
}


int main() {
	bool is_hole = false;
	cout << "Define the presence of hole (0/1):" << endl;
	cin >> is_hole;

	Mesh mesh(is_hole);

	string	filename_nodes = output_folder + "\\nodes.txt",
			filename_elements = output_folder + "\\elements.txt";
	CreateMesh(mesh, filename_nodes, filename_elements);


	FEM fem;
	fem.AssembleGlobalStiffnessMatrix(mesh);
	fem.AssembleGlobalLoadVector(mesh);
	
	//ofstream out;


	//out.open(output_folder + "\\out_stiffness_matrix.txt");
	//cout << "Printing results...\n";

	//for (int i = 0; i < fem.GlobalStiffnessMatrix.size(); i++) {
	//	for (int j = 0; j < fem.GlobalStiffnessMatrix.size(); j++) {
	//		out << fem.GlobalStiffnessMatrix[i][j] << "\t";
	//	}
	//	out << "\n";
	//}


	//out.close();
	//out.open(output_folder + "\\out_loads_vec.txt");

	//for (int i = 0; i < fem.GlobalLoadVector.size(); i++) {
	//	out << fem.GlobalLoadVector[i] << "\n";
	//}
	//cout << "Results has been printed\n";
	//out.close();


	int n = fem.GlobalStiffnessMatrix.size();
	Portrait portrait(n);
	int ig_n_1 = 0;
	fem.GeneratePortrait(portrait, fem.GlobalStiffnessMatrix, ig_n_1);

	double* b = new double[n];
	for (int i = 0; i < n; i++) {
		b[i] = fem.GlobalLoadVector[i];
	}

	int sz_ia = 0;
	int sz_ja = 0;


	double* x = new double[n];		// вектор неизвестных

	FromRSFToCSR_Real_1_Sym(n, portrait.ig, &sz_ia, &sz_ja);
	MKL_INT64* ia = new MKL_INT64[sz_ia];	// массив адресов начала ненулевых строк
	MKL_INT64* ja = new MKL_INT64[sz_ja];	// массив адресов столбцов ненулевых элементов
	double* ggl = new double[ig_n_1];
	double* a = new double[sz_ja];


	FromRSFToCSR_Real_2_Sym(n, portrait.ig, portrait.jg, portrait.di, portrait.ggl, ia, ja, a);
	for (int i = 0; i < sz_ia; i++) {
		ia[i]++;
	}
	for (int i = 0; i < sz_ja; i++) {
		ja[i]++;
	}



	solve_pardiso_symm(n, ia, ja, a, b, x);


	cout << "\n Solution is done!\n";

	fem.Get_X_Stresses(x, mesh, output_folder);

	delete[] b;
	delete[] x;
	delete[] ia;
	delete[] ja;


	return 0;
}



