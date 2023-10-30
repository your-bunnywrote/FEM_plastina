// Mesh.h - ����� ����������� ������ ����, �������, ��������� � �����,
// � ����� ������ ��� �� ����������
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

// ����� ����
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

// ����� "��������� �������"
class comp_domain {
public:
	int Nx, Ny;			// ���������� �����, �������������� ����������
	vector<Point> coords;	// ������ ��������� �����, �������������� ����������
	vector<pair<Point, Point >> domains;	// �������� ���� ������� - ��������������� ���� �������������� (x1,y1) � (x2,y2)
	void readfile_domains();
	bool is_contain(const Point& node);	// ���������, �������� �� ���� � �������� ����������
	bool is_match(const Point& node);	// ���������, ����� �� ���� �� ������
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
							  {0, 0, mat.c * (1 - mat.mu)/2}};	// ������� ���������, ��� ����������� ��������� ������� ���
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
	vector<int> nx,		// ���������� ����� �� �������������
				ny;
	vector<double> kx,	// ������������ ��������
				   ky;
	vector<int> num_nodes_in_new_mesh;	// �������� ����� ��������� �����, ������ ����� ���� (� ����) - ����� ���� �� ��������

};

// ��� ���������� ����������� ���������� �� ����� ������
vector<string> split(string& s, char delimeter);
// ��������� ����� (������ ������ � �����������, ���������� ������������ ��������, ��������� ����� � ������������ ������� ����� � ���������)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

