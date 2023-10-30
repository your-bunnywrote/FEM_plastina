#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

// ����� ��� �������� ������� � ����������� �������: ������ ����� ������ ����� ������������ � PARDISO.
// ������������� ������ ��� CSR ����� ������ �� ��������� � ������ Rectangle �������
class Portrait {
public:
	double* di;
	double* ggl;
	double* ggu;
	int* ig;
	int* jg;
	int N;
	Portrait();
	Portrait(int N);
	~Portrait();
	void PushBack(double*& arr, int &size, const double val);	// ��� ��������� ������� �������
	void PushBack(int*& arr, int& size, const int val);
};



#endif // !MATRIX_H


