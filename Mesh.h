
#pragma once
#ifndef MESH_H
#define MESH_H

#include "geometry.h"

using namespace std;

struct Material {
public:
	Material();
	size_t num;
	double E;
	double mu;
	double c;
	double thickness;

};


//  ласс €чейки, эквивалентный элементу, но содержащий в себе только информацию о геометрии
enum etype {RECT,QUAD};
class Cell {
public:
	vector<Point> loc_nodes;
	size_t num;
	Material mat;
	Cell();
};

class Mesh {
public:
	Mesh();
	Mesh(bool is_hole);

	comp_domain subdomain;
	vector<Point> nodes;
	vector<Cell> elements;
	vector<int> nx,		// количество узлов на подынтервалах
				ny;
	vector<double> kx,	// коэффициенты разр€дки
				   ky;
	vector<int> num_nodes_in_new_mesh;	// содержит новую нумерацию узлов, индекс этого узла (с нул€) - номер узла до удалени€


	// вычисление координат узлов сетки с круговым отверстием
	// так как дл€ вычислени€ координат нам нужны координаты интервалов горизонтальных и вертикальных линий, эта функци€ должна быть членом класса comp_domain или Mesh
	void calculate_coords(vector<double>& x, vector<double>& y);
	etype check_element_type(Cell& element);
};

// дл€ разделени€ считываемых параметров во врем€ чтени€
vector<string> split(string& s, char delimeter);

// генераци€ сетки (чтение данных о подобласт€х, вычисление коэффициента разр€дки, координат узлов и формирование списков узлов и элементов)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

