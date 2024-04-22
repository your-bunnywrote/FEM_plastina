
#pragma once
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

// папка дл€ входных данных
extern string input_folder;

// “очка/узел
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
//  рива€/лини€
class Curve {

public:
	curve_type type;
	int num;
	Point begin, end;
	Point center;
	double radius;

	Curve();
	// по умолчанию (без указани€ типа) крива€ - лини€, у которой всегда есть начало и конец, поэтому конструктор на входе получает две точки - начало и конец, а также номер
	Curve(Point begin_, Point end_, int num_);
	// если крива€ - дуга, то в конструктор подаетс€ тип (дуга), центр, радиус, начало и конец, а также номер
	Curve(curve_type type, Point center_, double radius_, Point begin_, Point end_, int num_);
};



// –ебро
//class Edge :public Curve{};



// –асчетна€ область
class comp_domain {
public:
	int Nx, Ny;			// количество линий, ограничивающих подобласти
	bool is_hole;
	vector<Point> coords;	// вектор координат линий, ограничивающих подобласти (дл€ пр€моугольной геометрии)
	vector<vector<Curve>> vertical_curves;
	vector<vector<Curve>> horizontal_curves;
	double length;
	double width;
	double hole_radius;
	Point hole_center;
	vector<pair<Point, Point >> rect_domains;	// содержит пару поинтов - противоположные углы пр€моугольника (x1,y1) и (x2,y2), фактически описывающие пр€моугольник
	void read_noholegeom_info();
	void create_holegeom_info();
	void nonsymmetric_hole_geom();
	bool is_contain(const Point& node);	// провер€ет, попадает ли узел в истинную подобласть
	comp_domain();
	comp_domain(bool is_hole);
};



#endif // !GEOMETRY_H

