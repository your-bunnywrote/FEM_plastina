
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

// �����/����
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
// ������/�����
class Curve {
protected:

	curve_type type;
public:
	int num;
	Point begin, end;
	Point center;
	double radius;

	Curve();
	// �� ��������� (��� �������� ����) ������ - �����, � ������� ������ ���� ������ � �����, ������� ����������� �� ��������� �� ����� �������� ��� ����� - ������ � �����, � ����� �����
	Curve(Point begin_, Point end_, int num_);
	// ���� ������ - ����, �� � ����������� �������� ��� (����), �����, ������, ������ � �����, � ����� �����
	Curve(curve_type type, Point center_, double radius_, Point begin_, Point end_, int num_);
};



// �����
//class Edge :public Curve{};



// ��������� �������
class comp_domain {
public:
	int Nx, Ny;			// ���������� �����, �������������� ����������
	vector<Point> coords;	// ������ ��������� �����, �������������� ����������
	double length;
	double width;
	double hole_radius;
	vector<pair<Point, Point >> rect_domains;	// �������� ���� ������� - ��������������� ���� �������������� (x1,y1) � (x2,y2), ���������� ����������� �������������
	void read_noholegeom_info();
	void create_holegeom_info();
	void nonsymmetric_hole_geom();
	bool is_contain(const Point& node);	// ���������, �������� �� ���� � �������� ����������
	comp_domain();
};



#endif // !GEOMETRY_H
