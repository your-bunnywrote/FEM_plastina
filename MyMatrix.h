#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

// класс для хранения матрицы в разреженном формате: данные этого класса будут отправляться в PARDISO.
// формироваться данные для CSR будут исходя из собранной в классе Rectangle матрицы
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
	void PushBack(double*& arr, int &size, const double val);	// для изменения размера массива
	void PushBack(int*& arr, int& size, const int val);
};



#endif // !MATRIX_H


