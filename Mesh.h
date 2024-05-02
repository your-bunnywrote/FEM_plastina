
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

enum etype { RECTANGLE, QUADR };
class Element {
public:
	etype type;
	vector<Point> loc_nodes;
	size_t num;
	Material mat;
	// �������������: ����� �������� �� ������� ��������� ������
	vector<vector<double>> D = {{mat.c, mat.c * mat.mu, 0},
							   {mat.c * mat.mu, mat.c * 1, 0},
							   {0, 0, mat.c * (1 - mat.mu)/2}};	// ������� ���������, ��� ����������� ��������� ������� ���		
	//												    E   |1 mu	 0	  |
	//											D  = ------	|mu 1	 0	  |
	//												 1-mu^2 |0  0 (1-mu)/2|
	Element();
protected:
	void twoD_to_oneD(size_t i, size_t& mu, size_t& nu);
	virtual double bfunc1D(size_t func_num, double x0, double x1, double x) {
		return 1;
	};
	virtual double dbfunc1D(size_t func_num, double x0, double x1, double x) {
		return 1;
	};


	


};



class Mesh {
public:
	Mesh();
	Mesh(bool is_hole);

	comp_domain subdomain;
	vector<Point> nodes;
	vector<Element> elements;
	vector<int> nx,		// ���������� ����� �� �������������
				ny;
	vector<double> kx,	// ������������ ��������
				   ky;
	vector<int> num_nodes_in_new_mesh;	// �������� ����� ��������� �����, ������ ����� ���� (� ����) - ����� ���� �� ��������


	// ���������� ��������� ����� ����� � �������� ����������
	// ��� ��� ��� ���������� ��������� ��� ����� ���������� ���������� �������������� � ������������ �����, ��� ������� ������ ���� ������ ������ comp_domain ��� Mesh
	void calculate_coords(vector<double>& x, vector<double>& y);

};

// ��� ���������� ����������� ���������� �� ����� ������
vector<string> split(string& s, char delimeter);

// ��������� ����� (������ ������ � �����������, ���������� ������������ ��������, ��������� ����� � ������������ ������� ����� � ���������)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

