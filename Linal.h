// Linal.h
// здесь представлены методы для упрощенной работы с матрицами
// например, объявление матрицы размером 4х4 и заполнение ее нулями выглядит как A(4,4)
// обращение к элементу матрицы

#ifndef LINAL_H
#define LINAL_H

#include <iostream>
#include <fstream>
#include <vector>

class MyArray;
class Matrix;

// Вектор
class MyArray {
public:
	MyArray();
	MyArray(int array_size);	// инициализация вектора
	MyArray(const MyArray &a);
	~MyArray();
	void Show();	// вывод содержимого вектора на экран
	void zap();		// заполнение вектора случайными числами от 0 до 9
	void Set(int index, double value);	// присваивание значения элементу вектора
	int get_size();		// получение размер вектора
	double & operator [](int index);
	MyArray operator =(const MyArray &a);
	double norma();

	void grad(double (*f)(MyArray x), MyArray &res, double eps);
	void Hessian(double (*f)(MyArray x), Matrix &matr, double eps);

	void WriteToFile();
    void ReadFromFile();
private:
	int array_size;
	double *p;
};

class Matrix{
public:
	Matrix();
	Matrix(int n);
	Matrix(int row, int col);
	Matrix(const Matrix &a);
	Matrix(int row, int col, bool isDiag);
	~Matrix();
	void Show();	// вывод матрицы на экран
	double & operator ()(int i,int j);
	void Set(int index1,int index2,double value);	// присваивание значения элементу матрицы
	int get_row();	// получение размера строки
	int get_col();	// получение размера столбца
	void zap();		// заполниение матрицы случайными числами в диапазоне от 0 до 9
	Matrix operator =(const Matrix &a);
	void LU_decomposition(Matrix &L, Matrix &U, int n);
	void Solve_Gauss_reverse(Matrix matr, MyArray B, MyArray &res, int n, bool isDown);
	void LU_solve(Matrix A, MyArray B, MyArray &result, int n);
	void CGM_solve(Matrix A, MyArray B, MyArray &result, int n);
	void scale(double value);	// умножение матрицы на число
    Matrix Sum(Matrix &a);	// матричное сложение
	Matrix Difference(Matrix &a);	// матричное вычитание
	Matrix Product(Matrix &a);	//	матричное умножение
	Matrix & transpose();	// транспонирование
	Matrix & transpose2();
	Matrix & Gauss();
	double det_gauss();
	double det(int n);
	void Get_matrix(int n, Matrix &temp_matr, int indRow, int indCol);	
	void inverse(Matrix &matr, int n, bool isDiag);

	void WriteToFile();
	void ReadFromFile();

private:
	int row;
	int col;
	bool isDiag;
	bool isTrans;
    double *m;
};

struct Triplet {
	Triplet(int x_value, int y_value, double value) {
		this->x_value = x_value;
		this->y_value = y_value;
		this->value = value;
	}
    int get_x() { return x_value; }
    int get_y() { return y_value; }
    double get_value() { return value; }

	void Show() {
		std::cout << x_value << " " << y_value << " " << value << "\n";
	}
	int x_value;
	int y_value;
	double value;
};

class SparseMatrix {
public:
	SparseMatrix(int sparse_size);
	~SparseMatrix();
    void ConvertTripletToSparse(std::vector<Triplet> t);
    void ConvertToMatrix(Matrix& M);
	void SortIt();
    void SparseLU();
    int get_size();
	void Show();
private:
	std::vector<Triplet> v;
	int nonzero;
	int sparse_size;
	int *x;
	int *y;
	double *data;
};
int CountNonZero(std::vector<Triplet> t);

#endif
