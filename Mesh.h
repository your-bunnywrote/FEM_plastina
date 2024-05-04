
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


// ����� ������, ������������� ��������, �� ���������� � ���� ������ ���������� � ���������
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
	vector<int> nx,		// ���������� ����� �� �������������
				ny;
	vector<double> kx,	// ������������ ��������
				   ky;
	vector<int> num_nodes_in_new_mesh;	// �������� ����� ��������� �����, ������ ����� ���� (� ����) - ����� ���� �� ��������


	// ���������� ��������� ����� ����� � �������� ����������
	// ��� ��� ��� ���������� ��������� ��� ����� ���������� ���������� �������������� � ������������ �����, ��� ������� ������ ���� ������ ������ comp_domain ��� Mesh
	void calculate_coords(vector<double>& x, vector<double>& y);
	etype check_element_type(Cell& element);
};

// ��� ���������� ����������� ���������� �� ����� ������
vector<string> split(string& s, char delimeter);

// ��������� ����� (������ ������ � �����������, ���������� ������������ ��������, ��������� ����� � ������������ ������� ����� � ���������)
void CreateMesh(Mesh& mesh, string& filename_nodes,string& filename_elements);




#endif // !MESH_H

