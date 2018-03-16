#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <random>
#include <ctime>
#include <chrono>
#include "CRS.h"
#include <omp.h>
using namespace std;
//n1,m1 - dimension1; m2,n2 - dimension2
default_random_engine gen(chrono::system_clock::now().time_since_epoch().count());
uniform_real_distribution <double> distr_val(1, 1e2);
uniform_int_distribution <int> col_ind(0,10);
void GenerateRegularCRS(int seed, int n, int m, int cntInRow, CRSmtx& mtx)
{
	int i, j, k, f, tmp, notNull, c;
	notNull = cntInRow * n;
	InitMtx(n, m, n*cntInRow, mtx);
	col_ind = uniform_int_distribution<int>(0, m - 1);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < cntInRow; j++)
		{
			do
			{
				//mtx.col_ind[i * cntInRow + j] = col_ind(gen) % (m-1);
				mtx.col_ind[i * cntInRow + j] = col_ind(gen);
				f = 0;
				for (k = 0; k < j; k++)
					if (mtx.col_ind[i * cntInRow + j] ==
						mtx.col_ind[i * cntInRow + k])
						f = 1;
			} while (f == 1);
		}
		for (j = 0; j < cntInRow - 1; j++)
			for (k = 0; k < cntInRow - 1; k++)
				if (mtx.col_ind[i * cntInRow + k] >
					mtx.col_ind[i * cntInRow + k + 1])
				{
					tmp = mtx.col_ind[i * cntInRow + k];
					mtx.col_ind[i * cntInRow + k] =
						mtx.col_ind[i * cntInRow + k + 1];
					mtx.col_ind[i * cntInRow + k + 1] = tmp;
				}
	}
	for (i = 0; i < cntInRow * n; i++)
		mtx.val[i] = distr_val(gen);
	c = 0;
	for (i = 0; i <= n; i++)
	{
		mtx.row_ptr[i] = c;
		c += cntInRow;
	}
}

int main(int argc, char * argv[])
{
	CRSmtx A, B, Bt, C;
	int i, j, k, f, tmp, notNull, c;
	string fileName, answerName, name;

	int n1, n2, m1, m2, per1, per2, nnz_per_row1, nnz_per_row2, nnz1, nnz2;
	if (argc > 2)
	{
		//n1 = n_tests[atoi(argv[1])];
		//m1 = m_tests[atoi(argv[2])];
		//n2 = n_tests[atoi(argv[3])];
		//m2 = m_tests[atoi(argv[4])];
		//nnz_per_row1 = nnz_per_row[atoi(argv[5])];
		//nnz_per_row2 = nnz_per_row[atoi(argv[6])];
		n1 = atoi(argv[1]);
		nnz_per_row1 = atoi(argv[2]);
		nnz_per_row2= atoi(argv[3]);
		name = argv[4]; 
		fileName = (name + ".test").c_str();
		answerName = (name + ".ans").c_str();
		if ((n1 < nnz_per_row1) || (n1 < nnz_per_row2))
		{
			cout << "Wrong arguments" << endl;
			return 1;
		}
	}
	else
	{
		cout << "Wrong number of arguments" << endl;
		return 1;
	}
	n2 = n1;
	m1 = n1;
	m2 = n1;
	nnz1 = nnz_per_row1*n1;
	nnz2 = nnz_per_row2*n2;

	GenerateRegularCRS(10000, n1, m1, nnz_per_row1, A);
	GenerateRegularCRS(10, n2, m2, nnz_per_row2, B);
	FILE* test = fopen((name).c_str(), "wb");

	WriteInBinaryFile(test,A); 
	WriteInBinaryFile(test,B); 
	TransposeMtx(B, Bt);
	GetMtxTemplate(A, Bt, C);
	double start, end;
	start = omp_get_wtime();
	CRSMatrMult(A, Bt, C); 
	end = omp_get_wtime();

	double time = end - start;
	FILE* perfect = fopen((name+".ans").c_str(), "wb");
	WriteInBinaryFile(perfect,C); 
	fwrite(&time, sizeof(time), 1, perfect);
	return 0;
}