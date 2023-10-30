#include "MyMatrix.h"


Portrait::Portrait() {};

Portrait::Portrait(int N) {
	this->N = N;
	di = new double[N];
	ggl = new double;
	ggu = new double;
	ig = new int[N+1];
	jg = new int;
}

Portrait::~Portrait() {
	delete[] di;
	delete[] ggl;
	delete[] ggu;
	delete[] ig;
	delete[] jg;
	std::cout << "Arrays was deleted" << std::endl;
}

void Portrait::PushBack(double*& arr, int &size, const double val) {
	double* newArray = new double[size + 1];
	for (int i = 0; i < size; i++) {
		newArray[i] = arr[i];
	}
	newArray[size++] = val;	// сначала добавляем новое значение в текущий элемент массива, затем смещаем индекс
	delete[] arr;
	arr = newArray;		// меняем указатель со старого удаленного массива на новый созданный
}

void Portrait::PushBack(int*& arr, int& size, const int val) {
	int* newArray = new int[size + 1];
	for (int i = 0; i < size; i++) {
		newArray[i] = arr[i];
	}
	newArray[size++] = val;	// сначала добавляем новое значение в текущий элемент массива, затем смещаем индекс
	delete[] arr;
	arr = newArray;		// меняем указатель со старого удаленного массива на новый созданный
}