// Mesh.h - здесь объ€вл€ютс€ классы узла, области, материала и сетки,
// а также метода дл€ ее построени€
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


class Element {
public:
	enum etype { RECTANGLE, QUADRILATERAL };
	Point* loc_nodes;
	size_t num;
	Material mat;
	// –≈ƒј “»–ќ¬ј“№: можно заменить на обычный двумерный массив
	vector<vector<double>> D = {{mat.c, mat.c * mat.mu, 0},
							   {mat.c * mat.mu, mat.c * 1, 0},
							   {0, 0, mat.c * (1 - mat.mu)/2}};	// матрица упругости, дл€ изотропного материала имеюща€ вид		
	//												    E   |1 mu	 0	  |
	//											D  = ------	|mu 1	 0	  |
	//												 1-mu^2 |0  0 (1-mu)/2|
	Element();
	~Element();
};



class Mesh {
public:
	comp_domain subdomain;
	vector<Point> nodes;
	vector<Element> elements;
	vector<int> nx,		// количество узлов на подынтервалах
				ny;
	vector<double> kx,	// коэффициенты разр€дки
				   ky;
	vector<int> num_nodes_in_new_mesh;	// содержит новую нумерацию узлов, индекс этого узла (с нул€) - номер узла до удалени€

};

// дл€ разделени€ считываемых параметров во врем€ чтени€
vector<string> split(string& s, char delimeter);
// генераци€ сетки (чтение данных о подобласт€х, вычисление коэффициента разр€дки, координат узлов и формирование списков узлов и элементов)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

