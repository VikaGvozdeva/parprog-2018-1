#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>
#include <string>
#include <omp.h>
#include "CRS.h"
using namespace std;
// Используется для взаимодействия с тестирующей системой
//////////////////////////////////////////////////////////////////////////////////////////
/*
 Checker может устанавливать вот эти три вердикта:
AC = Accepted = Решение выдаёт корректный результат на данном тесте
WA = Wrong Answer = Решение выдаёт некорректный результат на данном тесте
PE = Presentation Error = Ошибка формата выходных данных
 Остальные вердикты checker не может устанавливать
NO = No verdict = Вердикт отсутствует
CE = Compilation Error = Ошибка компиляции
ML = Memory Limit Exceeded = Превышено ограничение по памяти
TL = Time Limit Exceeded = Превышено ограничение по времени работы
RE = Runtime Error = Ошибка времени исполнения программы
IL = Idle Limit Exceeded = Превышено время простоя (бездействия) программы
DE = Deadly Error = Ошибка тестирующей системы
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
	// Сообщить тестирующей системе, что решение получило один из вердиктов verdict
	void write_verdict(verdict v)
	{
		write_type(ext_cls::VERDICT); 
		fwrite(&v, sizeof(v), 1, bur);
	}
	// Написать сообщение от checker'a пользователю.
	 //Например, что решение верное, или неверное.
	 //Использовать только латинские буквы и знаки препинания
	void write_message(string str)
	{
		write_type(ext_cls::MESSAGE); 
		int l = str.size(); 
		fwrite(&l, sizeof(l), 1, bur);
		fwrite(&str[0], sizeof(str[0]), l, bur);
	}
	/* Сообщить тестирующей системе время работы программы участника,
	 вычисленное с помощью before_code
	 x имеет размерность 100 нс = 10 ^ (-7) сек*/
	void write_time(long long x)
	{
		write_type(ext_cls::TIME); 
		fwrite(&x, sizeof(x), 1, bur);
	}
	 //Сообщить тестирующей системе, память затребованную программой участника
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
//	 Записываем время в правильной размерности 
	checker_result.write_time(ans_time*1e7);
	return 0;
}