#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <omp.h>
#include <chrono>
#include <random>
#include "CRS.h"

int main(int argc, char * argv[])
{
	string name;
	if (argc == 1)
	{
		name = argv[1];
	}
	else
	{
		cout << "Wrong count of arguments" << endl;
		return 1;
	}
	CRSmtx A, B, Bt, C;
	FILE* fileName = fopen((name ).c_str(), "rb");
	int n, m, nnz;
	fread(&n, sizeof(int), 1, fileName);
	fread(&m, sizeof(int), 1, fileName);
	fread(&nnz, sizeof(int), 1, fileName);
	InitMtx(n, m, nnz, A);
	ReadFromBinaryFile(fileName,A);
	fread(&n, sizeof(int), 1, fileName);
	fread(&m, sizeof(int), 1, fileName);
	fread(&nnz, sizeof(int), 1, fileName);
	InitMtx(n, m, nnz, B);
	ReadFromBinaryFile(fileName,B);

	TransposeMtx(B, Bt);
	GetMtxTemplate(A, Bt, C);
	double start, end;
	start = omp_get_wtime();
	CRSMatrMult(A, Bt, C);
	end = omp_get_wtime();
	double time = end - start;
	FILE* sol = fopen((name + ".sol").c_str(), "wb");
	WriteInBinaryFile(sol,C);
	fwrite(&time, sizeof(time), 1, sol);
	return 0;
}