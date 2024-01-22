// Linal.h
// ����� ������������ ������ ��� ���������� ������ � ���������
// ��������, ���������� ������� �������� 4�4 � ���������� �� ������ �������� ��� A(4,4)
// ��������� � �������� �������

#ifndef LINAL_H
#define LINAL_H

#include <iostream>
#include <fstream>
#include <vector>

class MyArray;
class Matrix;

// ������
class MyArray {
public:
	MyArray();
	MyArray(int array_size);	// ������������� �������
	MyArray(const MyArray &a);
	~MyArray();
	void Show();	// ����� ����������� ������� �� �����
	void zap();		// ���������� ������� ���������� ������� �� 0 �� 9
	void Set(int index, double value);	// ������������ �������� �������� �������
	int get_size();		// ��������� ������ �������
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
	void Show();	// ����� ������� �� �����
	double & operator ()(int i,int j);
	void Set(int index1,int index2,double value);	// ������������ �������� �������� �������
	int get_row();	// ��������� ������� ������
	int get_col();	// ��������� ������� �������
	void zap();		// ����������� ������� ���������� ������� � ��������� �� 0 �� 9
	Matrix operator =(const Matrix &a);
	void LU_decomposition(Matrix &L, Matrix &U, int n);
	void Solve_Gauss_reverse(Matrix matr, MyArray B, MyArray &res, int n, bool isDown);
	void LU_solve(Matrix A, MyArray B, MyArray &result, int n);
	void CGM_solve(Matrix A, MyArray B, MyArray &result, int n);
	void scale(double value);	// ��������� ������� �� �����
    Matrix Sum(Matrix &a);	// ��������� ��������
	Matrix Difference(Matrix &a);	// ��������� ���������
	Matrix Product(Matrix &a);	//	��������� ���������
	Matrix & transpose();	// ����������������
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
