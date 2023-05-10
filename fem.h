#ifndef FEM_H
#define FEM_H

#include "Mesh.h"

class Rectangle : public Element {
public:
	void init();
	// Метод для получения номера одномерной функции формы по номеру двумерной
	void TwoByOne(size_t two, size_t& one1, size_t& one2);
	// Одномерные функции формы
	double bfunc1D(size_t func_n, double ksi);
	// Производные функций форм
	double dbfunc1D(size_t func_n, double ksi);
	// Двухмерные функции формы
	double bfunc2D_local(size_t num, const Point& p);
	double bfunc2D(size_t num, const Point& p);
	Point GaussPoints[12];
	double GaussWeights[12];


};











#endif // !FEM_H

