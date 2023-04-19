#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>

using namespace std;

// Класс узла
class Point {
public:
	double x, y;
	size_t num;
	Point();
	Point(const double x, const double y, const int num);
};

// Класс "Расчетная область"
class comp_domain {
public:
	int Nx, Ny;			// количество линий, ограничивающих подобласти
	vector<Point> coords;	// вектор координат линий, ограничивающих подобласти
	vector<pair<Point, Point >> domains;	// содержит пару поинтов - противоположные углы прямоугольника (x1,y1) и (x2,y2)
	void readfile_domains();
	bool is_contain(const Point& node);	// проверяет, попадает ли узел в истинную подобласть
};

class Phys_area {
public:
	double E;
	double mu;
};

class Element {
public:
	Point loc_nodes[4];
	int material;
	size_t num;
};

class Mesh {
public:
	comp_domain subdomain;
	vector<Point> nodes;
	vector<Element> elements;
	vector<int> nx,		// количество узлов на подынтервалах
				ny;
	vector<double> kx,	// коэффициенты разрядки
				   ky;
	void readfile_partition_info();	// чтение информации о разбиениии подобластей на элементы
	// считаем координаты узлов с учетом коэффициента разрядки
	vector<double> calc_nodes_x_coords(const vector<double>& k, const vector<int>& n, vector<Point>& coordXw);
	vector<double> calc_nodes_y_coords(const vector<double>& k, const vector<int>& n, vector<Point>& coordYw);
	void fill_nodes();
	void fill_elements();
	vector<Point> edit_nodes(vector<Point>& nodes);
	vector<Element> edit_elements(vector<Element>& elements);
	void output(const string& filename);
	void output(const string& filename, const vector<Point>& nodes);
	void output(const string& filename, const vector<Element>& elements);
};

// для разделения считываемых параметров во время чтения
vector<string> split(string& s, char delimeter);




#endif // !MESH_H

