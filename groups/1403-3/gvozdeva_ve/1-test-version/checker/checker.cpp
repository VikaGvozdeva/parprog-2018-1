#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>
#include <string>
#include <omp.h>
#include "CRS.h"
using namespace std;
// ������������ ��� �������������� � ����������� ��������
//////////////////////////////////////////////////////////////////////////////////////////
/*
 Checker ����� ������������� ��� ��� ��� ��������:
AC = Accepted = ������� ����� ���������� ��������� �� ������ �����
WA = Wrong Answer = ������� ����� ������������ ��������� �� ������ �����
PE = Presentation Error = ������ ������� �������� ������
 ��������� �������� checker �� ����� �������������
NO = No verdict = ������� �����������
CE = Compilation Error = ������ ����������
ML = Memory Limit Exceeded = ��������� ����������� �� ������
TL = Time Limit Exceeded = ��������� ����������� �� ������� ������
RE = Runtime Error = ������ ������� ���������� ���������
IL = Idle Limit Exceeded = ��������� ����� ������� (�����������) ���������
DE = Deadly Error = ������ ����������� �������
*/
enum verdict { NO = 1, AC, WA, CE, ML, TL, RE, IL, PE, DE };
class result
{
private:
	FILE * bur;
public:
	enum ext_cls { NO = 1, VERDICT, MESSAGE, TIME, MEMORY };
	result(bool read = false)
	{
		if (read) bur = fopen("result.txt", "r"); 
		else bur = fopen("result.txt", "w");
	}
	~result() { fclose(bur); }
	void write_type(ext_cls t) 
	{
		fwrite(&t, sizeof(t), 1, bur); 
	}
	// �������� ����������� �������, ��� ������� �������� ���� �� ��������� verdict
	void write_verdict(verdict v)
	{
		write_type(ext_cls::VERDICT); 
		fwrite(&v, sizeof(v), 1, bur);
	}
	// �������� ��������� �� checker'a ������������.
	 //��������, ��� ������� ������, ��� ��������.
	 //������������ ������ ��������� ����� � ����� ����������
	void write_message(string str)
	{
		write_type(ext_cls::MESSAGE); 
		int l = str.size(); 
		fwrite(&l, sizeof(l), 1, bur);
		fwrite(&str[0], sizeof(str[0]), l, bur);
	}
	/* �������� ����������� ������� ����� ������ ��������� ���������,
	 ����������� � ������� before_code
	 x ����� ����������� 100 �� = 10 ^ (-7) ���*/
	void write_time(long long x)
	{
		write_type(ext_cls::TIME); 
		fwrite(&x, sizeof(x), 1, bur);
	}
	 //�������� ����������� �������, ������ ������������� ���������� ���������
	void write_memory(unsigned long long x)
	{
		write_type(ext_cls::MEMORY); 
		fwrite(&x, sizeof(x), 1, bur);
	}
} checker_result;
//////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
	string name;
	if (argc > 1)
	{
		name = argv[1];
	}
	else
	{
		cout << "Wrong count of arguments" << endl;
		return 1;
	}
	int n, m, nnz;
	int n1, m1, nnz1;
	CRSmtx C, C_;
	double res_time, ans_time;
	double diff, diff_time;
	FILE* ans = fopen((name + ".ans").c_str(), "rb");
	FILE* sol = fopen((name + ".sol").c_str(), "rb");
	fread(&n, sizeof(int), 1, ans);
	fread(&m, sizeof(int), 1, ans);
	fread(&nnz, sizeof(int), 1, ans);
	fread(&n1, sizeof(int), 1, sol);
	fread(&m1, sizeof(int), 1, sol);
	fread(&nnz1, sizeof(int), 1, sol);
	bool result = true;
	if ((n != n1) || (m != m1) || (nnz != nnz1))
	{
		result = false;
	}
	if (result)
	{
		InitMtx(n, m, nnz, C);
		InitMtx(n, m, nnz, C_);
		ReadFromBinaryFile(ans, C);
		ReadFromBinaryFile(sol, C_);
		CompareMtx(C, C_, diff);
	}

	fread(&ans_time, sizeof(double), 1, ans);
	fread(&res_time, sizeof(double), 1, sol);
	fclose(sol); fclose(ans);
	if (diff < 1e-6)
	{
		checker_result.write_message("AC. Numbers are equal.");
		checker_result.write_verdict(verdict::AC);
	}
	else
	{
		checker_result.write_message("WA. Output is not correct.");
		checker_result.write_verdict(verdict::WA);
	}
//	 ���������� ����� � ���������� ����������� 
	checker_result.write_time(ans_time*1e7);
	return 0;
}