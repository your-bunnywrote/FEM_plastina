// Mesh.h - здесь объявляются классы узла, области, материала и сетки,
// а также метода для ее построения
#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>

using namespace std;

// Класс узла
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

// Класс "Расчетная область"
class comp_domain {
public:
	int Nx, Ny;			// количество линий, ограничивающих подобласти
	vector<Point> coords;	// вектор координат линий, ограничивающих подобласти
	vector<pair<Point, Point >> domains;	// содержит пару поинтов - противоположные углы прямоугольника (x1,y1) и (x2,y2)
	void readfile_domains();
	bool is_contain(const Point& node);	// проверяет, попадает ли узел в истинную подобласть
	bool is_match(const Point& node);	// проверяет, лежит ли узел на прямой
};

struct Material {
public:
	Material();
	size_t num;
	double E;
	double mu;
	double c;
	double thickness;

};

class Element {
public:
	vector<Point> loc_nodes;
	size_t num;
	Material mat;
	vector<vector<double>> D{ {mat.c, mat.c * mat.mu, 0},
							  {mat.c * mat.mu, mat.c * 1, 0},
							  {0, 0, mat.c * (1 - mat.mu)/2}};	// матрица упругости, для изотропного материала имеющая вид
	//												    E   |1 mu	 0	  |
	//											D  = ------	|mu 1	 0	  |
	//												 1-mu^2 |0  0 (1-mu)/2|
	Element();
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
	vector<int> num_nodes_in_new_mesh;	// содержит новую нумерацию узлов, индекс этого узла (с нуля) - номер узла до удаления

};

// для разделения считываемых параметров во время чтения
vector<string> split(string& s, char delimeter);
// генерация сетки (чтение данных о подобластях, вычисление коэффициента разрядки, координат узлов и формирование списков узлов и элементов)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

