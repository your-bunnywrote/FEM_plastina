#ifndef FEM_H
#define FEM_H

#include "Mesh.h"

class Rectangle : public Element {
public:
	void init();
	// ����� ��� ��������� ������ ���������� ������� ����� �� ������ ���������
	void TwoByOne(size_t two, size_t& one1, size_t& one2);
	// ���������� ������� �����
	double bfunc1D(size_t func_n, double ksi);
	// ����������� ������� ����
	double dbfunc1D(size_t func_n, double ksi);
	// ���������� ������� �����
	double bfunc2D_local(size_t num, const Point& p);
	double bfunc2D(size_t num, const Point& p);
	Point GaussPoints[12];
	double GaussWeights[12];


};











#endif // !FEM_H

