#include "omp_mult.h"
void CRSMatrMultParallel(CRSmtx A, CRSmtx B, CRSmtx &C)
{

#pragma omp parallel for
	for (int i = 0; i < C.n; ++i)
	{
		int *x = new int[A.m];
		memset(x, -1, sizeof(int)*A.m);

		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
		{
			x[A.col_ind[j]] = j;
		}

		for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1]; ++j)
		{
			double sum = 0.0;

			int col = C.col_ind[j];
			for (int q = B.row_ptr[col]; q < B.row_ptr[col + 1]; ++q)
			{
				int ind = B.col_ind[q];
				if (x[ind] != -1)
					sum += A.val[x[ind]] * B.val[q];
			}
			C.val[j] = sum;
		}
		delete[] x;

	}
}