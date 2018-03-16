#ifndef GVOZDEVA_VE_TEST_VERSION_CRSMULTIPLY_H
#define GVOZDEVA_VE_TEST_VERSION_CRSMULTIPLY_H
#include <iostream>
#include <vector>
using namespace std;
struct CRSmtx
{
	int n, m;
	int nnz;
	double* val;
	int* col_ind;
	int* row_ptr;
};

	void InitMtx(int n, int m, int nnz, CRSmtx &Mtx);
	void DeleteMtx(CRSmtx &Mtx);
	int CompareMtx(CRSmtx A, CRSmtx B, double& diff);
	void TransposeMtx(CRSmtx &A, CRSmtx &At); 
	void PrintMtx(CRSmtx Mtx);
	void GetMtxTemplate(CRSmtx A, CRSmtx B, CRSmtx &C); 
	void CRSMatrMult(CRSmtx A, CRSmtx Bt, CRSmtx &C);
	void WriteInBinaryFile(FILE* fileName, CRSmtx A);
	void ReadFromBinaryFile(FILE* fileName, CRSmtx A);
#endif //GVOZDEVA_VE_TEST_VERSION_CRSMULTIPLY_H