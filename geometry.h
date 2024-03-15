
#ifndef GEOMETRY_H
#define GEOMETRY_H

#define _USE_MATH_DEFINES


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <math.h>

using namespace std;

// Точка/узел
class Point {
public:
	double x, y;
	size_t num;
	Point();
	Point(const double x, const double y, const int num);
	Point(const double x, const double y) :
		x(x),
		y(y),
		num()
	{}
	Point operator + (Point p);
	Point operator + (double val);
	Point operator / (double val);
	Point operator - (Point p);
	Point operator * (double t);
};

enum curve_type { line, arc };
// Кривая/линия
class Curve {
protected:

	curve_type type;
public:
	int num;
	Point begin, end;
	Point center;
	double radius;

	Curve();
	// по умолчанию (без указания типа) кривая - линия, у которой всегда есть начало и конец, поэтому конструктор по умолчанию на входе получает две точки - начало и конец, а также номер
	Curve(Point begin_, Point end_, int num_);
	// если кривая - дуга, то в конструктор подается тип (дуга), центр, радиус, начало и конец, а также номер
	Curve(curve_type type, Point center_, double radius_, Point begin_, Point end_, int num_);
};



// Ребро
//class Edge :public Curve{};



// Расчетная область
class comp_domain {
public:
	int Nx, Ny;			// количество линий, ограничивающих подобласти
	vector<Point> coords;	// вектор координат линий, ограничивающих подобласти
	double length;
	double width;
	double hole_radius;
	vector<pair<Point, Point >> rect_domains;	// содержит пару поинтов - противоположные углы прямоугольника (x1,y1) и (x2,y2), фактически описывающие прямоугольник
	void read_noholegeom_info();
	void create_holegeom_info();
	void nonsymmetric_hole_geom();
	bool is_contain(const Point& node);	// проверяет, попадает ли узел в истинную подобласть
	comp_domain();
};



#endif // !GEOMETRY_H

