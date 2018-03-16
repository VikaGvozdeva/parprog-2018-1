#include "CRS.h"

void InitMtx(int n, int m, int nnz, CRSmtx &Mtx)
{
	Mtx.n = n;
	Mtx.m = m;
	Mtx.nnz = nnz;
	Mtx.val = new double[nnz];
	Mtx.col_ind = new int[nnz];
	Mtx.row_ptr = new int[n + 1];
	for (int i = 0; i < nnz; ++i) 
	{
		Mtx.col_ind[i] = 0;
		Mtx.val[i] = 0.0;
 	}
	for (int i = 0; i < n + 1; ++i)
	{
		Mtx.row_ptr[i] = 0;
	}
}
void DeleteMtx(CRSmtx &Mtx)
{
	delete[] Mtx.val;
	delete[] Mtx.col_ind;
	delete[] Mtx.row_ptr;
}
int CompareMtx(CRSmtx A, CRSmtx B, double& diff)
{
	if (A.n != B.n || A.m != B.m)
		return 1;
	int n = A.n;
	int m = A.m;

	vector<vector<double>> p(n, vector<double>(m));


	for (int i = 0; i < n; ++i)
	{
		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
		{
			p[i][A.col_ind[j]] = A.val[j];
		}
	}


	double max = -1;
	double tmp, tmp2;
	for (int i = 0; i < n; ++i)
	{
		for (int j = B.row_ptr[i]; j < B.row_ptr[i + 1]; ++j)
		{
			tmp = abs(p[i][B.col_ind[j]] - B.val[j]);
			if (tmp > max) max = tmp;
		}
	}

	diff = max;
	return 0;
}
void TransposeMtx(CRSmtx &A, CRSmtx &At)
{
	int n = A.n;
	int m = A.m;
	int nnz = A.nnz;
	InitMtx(m, n, nnz, At);

	for (int i = 0; i < nnz; ++i) 
	{
		At.row_ptr[A.col_ind[i] + 1]++;
	}

	int ind = 0;
	int tmp;
	for (int i = 1; i < m + 1; ++i) 
	{
		tmp = At.row_ptr[i];
		At.row_ptr[i] = ind;
		ind += tmp;
	}
	int col, row;
	double val;
	for (int i = 0; i < n; ++i) 
	{
		col = i;
		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) 
		{
			val = A.val[j];
			row = A.col_ind[j];
			ind = At.row_ptr[row + 1];
			At.val[ind] = val;
			At.col_ind[ind] = col;
			At.row_ptr[row + 1]++;
		}
	}
}
void PrintMtx(CRSmtx Mtx)
{
	int n = Mtx.n;
	int m = Mtx.m;
	vector<vector<double>> p(Mtx.n, vector<double>(Mtx.m));


	for (int i = 0; i < n; ++i)
	{
		for (int j = Mtx.row_ptr[i]; j < Mtx.row_ptr[i + 1]; ++j)
		{
			p[i][Mtx.col_ind[j]] = Mtx.val[j];
		}
	}
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < m; ++j) 
		{
			cout << p[i][j] << " ";
			
		}
		cout << endl;
	}
}
void GetMtxTemplate(CRSmtx A, CRSmtx B, CRSmtx &C)
{
	int **templ = new int*[A.n];
	int nz2 = 0;
	for (int i = 0; i < A.n; ++i)
	{
		int *x = new int[A.m];
		templ[i] = new int[A.m];

		memset(x, 0, sizeof(int)*A.m);
		memset(templ[i], sizeof(int) * 0, A.m);

		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
		{
			x[A.col_ind[j]] = 1;
		}

		for (int j = 0; j < B.n; ++j)
		{
			int nz = 0;

			for (int q = B.row_ptr[j]; q < B.row_ptr[j + 1]; ++q)
			{
				if (x[B.col_ind[q]] == 1)
				{
					nz++;
					break;
				}
			}
			if (nz > 0)
			{
				templ[i][j] = -1;
				nz2++;
			}

		}
		delete[] x;
	}
	InitMtx(A.n, A.m, nz2, C);

	int cnt = 0;

	for (int i = 0; i < A.n; ++i)
	{
		for (int j = 0; j < A.m; ++j)
		{
			//if (templ[i][j]>0)
			if (templ[i][j] ==  -1)
			{
				C.col_ind[cnt] = j;
				C.val[cnt] = 1;
				cnt++;
			}

		}
		C.row_ptr[i + 1] = cnt;
		//delete[] templ[i];
	}

	if (cnt != nz2) cout << "fail!" << endl;

	for (int i = 0; i < A.n; ++i)
		delete templ[i];

	delete[] templ;
}
void GetMtxTemplate1(CRSmtx A, CRSmtx Bt, CRSmtx &C)
{
	int **templ = new int*[A.n];
	int nz2 = 0;
	for (int i = 0; i < A.n; ++i)
	{
		int *x = new int[A.m];
		templ[i] = new int[A.m];

		memset(x, 0, sizeof(int)*A.m);
		memset(templ[i], sizeof(int) * 0, A.m);

		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
		{
			x[A.col_ind[j]] = 1;
		}

		for (int j = 0; j < Bt.n; ++j)
		{
			int nz = 0;

			for (int q = Bt.row_ptr[j]; q < Bt.row_ptr[j + 1]; ++q)
			{
				if (x[Bt.col_ind[q]] == 1)
				{
					nz++;
					break;
				}
			}
			if (nz > 0)
			{
				templ[i][j] = 1;
				nz2++;
			}

		}
		delete[] x;
	}
	InitMtx(A.n, Bt.n, nz2, C);
	//cout << nz;
	int cnt = 0;


	for (int i = 0; i < A.n; ++i)
	{
		for (int j = 0; j < Bt.n; ++j)
		{
			if (templ[i][j] >  0)
			{
				C.col_ind[cnt] = j;
				C.val[cnt] = 1;
				cnt++;
			}

		}
		C.row_ptr[i + 1] = cnt;
		//delete[] templ[i];
	}

	if (cnt != nz2) cout << "fail!" << endl;

} 
void CRSMatrMult(CRSmtx A, CRSmtx B, CRSmtx &C)
{
	for (int i = 0; i < C.n; ++i)
	{
		int *x = new int[A.m];
		memset(x, -1, sizeof(int)*A.m);

		for (int j = A.row_ptr[i]; j <A.row_ptr[i + 1]; ++j)
		{
			x[A.col_ind[j]] = j;
		}

		for (int j = C.row_ptr[i]; j < C.row_ptr[i + 1]; ++j)
		{
			double sum = 0.0;

			int col = C.col_ind[j];
			for (int q =B.row_ptr[col]; q < B.row_ptr[col + 1]; ++q)
			{
				int ind = B.col_ind[q];
				if (x[ind] != -1)
					sum += A.val[x[ind]] * B.val[q];
			}
			C.val[j] = sum;
		}
		//delete[] x;

	}
}

void ReadFromBinaryFile(FILE* fileName, CRSmtx A)
{
	fread(A.val, sizeof(double), A.nnz, fileName);
	fread(A.col_ind, sizeof(int), A.nnz, fileName);
	fread(A.row_ptr, sizeof(int), A.n+1, fileName);
}

void WriteInBinaryFile(FILE* fileName, CRSmtx A)
{
	fwrite(&A.n, sizeof(int), 1, fileName);
	fwrite(&A.m, sizeof(int), 1, fileName);
	fwrite(&A.nnz, sizeof(int), 1, fileName);
	fwrite(A.val, sizeof(double), A.nnz, fileName);
	fwrite(A.col_ind, sizeof(int), A.nnz, fileName);
	fwrite(A.row_ptr, sizeof(int), A.n+1, fileName);

}
