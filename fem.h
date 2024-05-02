// fem.h - ����� ����������� ����������� �� �������� ����� �������������, ������� ��� �����������
// ���� � ������ ��� ������ � �������-���������� �������
#pragma once
#ifndef FEM_H
#define FEM_H

#include "Mesh.h"
#include "MyMatrix.h"

class block2x2;
class block1x2;


class Load;
class Rectangle : public Element {
public:
	static void init();
	Rectangle();
	// ���������� �������� ������� �����
	// �������� �� ��������� x0, x1, x, ������� ����� ����������� ��� ���������� ������� �� Y, ��� ��� ����� ���������� ���
	double bfunc1D(size_t fucn_num, double x0, double x1, double x) override;

	// ����������� ������� ����
	double dbfunc1D(size_t func_num, double x0, double x1, double x) override;

	double param_func_x(double t, Point& from, Point& to);		// ��������������� ������������� ������� ���������
	double param_func_y(double t, Point& from, Point& to);		
	// ��������� ���������� ������� �����
	// � �������� ����������� ���������� x ���� ����� ���������� ����� ������, �� ������� ����� �������������
	// �� ���� ��� ������� ������ ���� � ����� �� ������ ������
	double phi(size_t func_num, Point& from, Point& to, double x, double y);

	// ����������� ��������� ������� ����
	//double dphi(size_t var, size_t num1, size_t num2, Point& from, Point& to, double x, double y);
	double dphi(size_t var, size_t num1, size_t num2, Point from, Point to, double x, double y);

	double gauss_points_local[3];
	double gauss_weights[3];
	Point integrate_points[9];
	//static double gauss_points_local[3];
	//static Point integrate_points[9];
	//static double gauss_weights[3];
	// �������������� ������� ������ �� ��������
	double Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var);
	// �������������� ������� ������ �� �����
	double Edge_IntegrateGauss3(Point& from, Point& to, size_t num);
	// ���������� ��������� ������� ��������� ��� ��������
	void Calculate_LocalStiffnessMatrix(Element& element);
	vector<int> fixed_nodes;
	vector<int> loaded_nodes;
	// ���� ����� ������� �����, ��� ��� ������� ��������� ����� ������� ���������
	vector<vector<block2x2>> LocalStiffnessMatrix_block;
	vector<vector<block2x2>> GlobalStiffnessMatrix_block;
	// ������������ ������� ������� ���������
	vector<vector<double>> LocalStiffnessMatrix;
	vector<vector<double>> GlobalStiffnessMatrix;
	// ������ ���������� ������� ���������
	void Assemble_GlobalStiffnessMatrix(Mesh& mesh);
	
	// ���������� ���������� ������� ��������
	void Calculate_LocalLoadVector(Element& element, Load P);
	vector<block1x2> LocalLoadVector_block;
	vector<block1x2> GlobalLoadVector_block;

	// ������������ ������� �������� ��������
	vector<double> LocalLoadVector;
	vector<double> GlobalLoadVector;

	void Assemble_GlobalLoadVector(Mesh& mesh);

	// ��������������� ������� � ����������� ������ ��� ������������� � PARDISO
	
	void GeneratePortrait(Portrait & Matrix, vector<vector<double>> GSM, int &ja_sz);
};

// ����� ��������
class Load {
public:
	Load();

	double Px, Py;	// ���������� ������������� ���
};

// ���� ���� ����� ��������� ��������� ������� ���������
// ���������� ����
class block2x2 {
public:
	block2x2();
	block2x2(double val11, double val12, double val21, double val22);
	double val11, val12,
		   val21, val22;
	// �������� ������ (����, � �������� ���������� - �����, � ������������ ���� - ������)
	block2x2 operator + (block2x2 block2);
	block2x2 operator += (block2x2 block2);
	block2x2 operator = (double val);


};

// ���� ���� ����� ��������� ������� �������� ��� ����������, ��� � �����������

class block1x2 {
public:
	block1x2();
	block1x2(double val1, double val2);
	double val1, val2;
	block1x2 operator +(block1x2 block2);
	block1x2 operator += (block1x2 block2);
	block1x2 operator = (double val);
};


#endif // !FEM_H

