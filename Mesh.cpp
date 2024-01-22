// Mesh.cpp - определение всех функций и операторов, модуль для построения сетки

#include "Mesh.h"

vector<string> split(string& s, char delimeter) {
		stringstream ss(s);
		string item;
		vector<string> tokens;
		while (getline(ss, item, delimeter))
		{
			tokens.push_back(item);
		}
		return tokens;
}


string input_folder = "input";

// =======================================================================

Point::Point() {
	x = 0.;
	y = 0.;
	num = 0;
}

Point::Point(double x, double y, int num) {
	this->x = x;
	this->y = y;
	this->num = num;
}

Point Point::operator*(double t) {
	double x_res, y_res;
	x_res = x * t;
	y_res = y * t;
	return Point(x_res, y_res);
}

Point Point::operator+(Point p) {
	double x_res, y_res;
	x_res = x + p.x;
	y_res = y + p.y;
	return Point(x_res, y_res);
}

Point Point::operator-(Point p) {
	double x_res, y_res;
	x_res = x - p.x;
	y_res = y - p.y;
	return Point(x_res, y_res);
}

Point Point::operator+(double val) {
	double res_x, res_y;
	res_x = x + val;
	res_y = y + val;
	return Point(res_x, res_y);
}

Point Point::operator/(double val) {
	double res_x, res_y;
	res_x = x / val;
	res_y = y / val;
	return Point(res_x, res_y);
}

Element::Element() {
	loc_nodes.resize(4);
}

//========================================================================

void comp_domain::readfile_domains() {
	ifstream W;
	W.open( input_folder + "\\subdomains.txt");

	W >> Nx;
	double x;
	for (size_t i = 0; i < Nx; i++) {
		W >> x;
		coords.push_back(Point(x, -1., 0));
	}
	W >> Ny;
	for (size_t i = 0; i < Ny; i++) {
		double y;
		W >> y;
		if (i < Nx) {
			coords[i].y = y;
		}
		else {
			coords.push_back(Point(-1., y, 0));
		}
	}
	int W_count;
	W >> W_count;
	domains.resize(W_count);
	double x1, x2, y1, y2;
	int mat;
	for (size_t i = 0; i < W_count; i++) {
		W >> mat >> x1 >> x2 >> y1 >> y2;
		domains[i].first.x = coords[x1].x;
		domains[i].first.y = coords[y1].y;
		domains[i].second.x = coords[x2].x;
		domains[i].second.y = coords[y2].y;
	}
	W.close();
}

bool comp_domain::is_contain(const Point& node) {
	for (size_t i = 0; i < domains.size(); i++) {
		if ((node.x >= domains[i].first.x) && (node.x <= domains[i].second.x + 1) && (node.y >= domains[i].first.y) && (node.y <= domains[i].second.y)) 
			return true;
	}
	return false;
}



// =======================================================================

void CreateMesh(Mesh& mesh, string& filename_nodes, string& filename_elements) {
	// информация о разбиении областей
	mesh.subdomain.readfile_domains();
	ifstream Xmsh, Ymsh;
	// x
	Xmsh.open( input_folder + "\\partition_info_x.txt");
	string line;
	while (getline(Xmsh, line)) {
		vector<string> tokens = split(line, '\t');
		mesh.nx.push_back(stoi(tokens.at(0)));
		mesh.kx.push_back(stod(tokens.at(1)));		// k>0 - каждый следующий интервал больше предыдущего в k раз
													// k<0 - каждый следующий интервал меньше предыдущего в k раз (т. е. больше предыдущего в 1/k раз)
		if (mesh.kx[mesh.kx.size() - 1] < 0) {
			mesh.kx[mesh.kx.size() - 1] = 1 / abs(mesh.kx[mesh.kx.size() - 1]);
		}
	}
	Xmsh.close();
	line.clear();
	// y
	Ymsh.open( input_folder + "\\partition_info_y.txt");
	while (getline(Ymsh, line)) {
		vector<string> tokens = split(line, '\t');
		mesh.ny.push_back(stoi(tokens.at(0)));
		mesh.ky.push_back(stod(tokens.at(1)));
		if (mesh.ky[mesh.ky.size() - 1] < 0) {
			mesh.ky[mesh.ky.size() - 1] = 1 / abs(mesh.ky[mesh.ky.size() - 1]);
		}
	}
	Ymsh.close();
	line.clear();
	// считаем координаты узлов
	// по x
	vector<double> x, y;
	double hx, hy, sum;
	x.resize(accumulate(mesh.nx.begin(), mesh.nx.end(), 0) + 1);
	for (size_t j = 0; j < mesh.nx.size(); j++) {
		sum = 0;
		int g_index = accumulate(mesh.nx.begin(), mesh.nx.begin() + j, 0);
		for (int i = 0; i < mesh.nx[j]; i++) {
			sum += pow(mesh.kx[j], i);
		}
		hx = (mesh.subdomain.coords[j + 1].x - mesh.subdomain.coords[j].x) / sum;
		x[g_index] = (mesh.subdomain.coords[j].x);
		for (int i = 1; i <= mesh.nx[j]; i++) {
			x[i + g_index] = x[i - 1 + g_index] + hx;
			hx *= mesh.kx[j];
		}
		x[g_index + mesh.nx[j]] = mesh.subdomain.coords[j+1].x;	// округляем посчитанную координату узла на границе подобластей до точной координаты границы (чтобы избежать значения в формате .000000001)
	}
	// по y
	y.resize(accumulate(mesh.ny.begin(), mesh.ny.end(), 0) + 1);
	for (size_t j = 0; j < mesh.ny.size(); j++) {
		sum = 0;
		int g_index = accumulate(mesh.ny.begin(), mesh.ny.begin() + j, 0);
		for (int i = 0; i < mesh.ny[j]; i++) {
			sum += pow(mesh.ky[j], i);
		}
		hy = (mesh.subdomain.coords[j + 1].y - mesh.subdomain.coords[j].y) / sum;
		y[g_index] = (mesh.subdomain.coords[j].y);
		for (int i = 1; i <= mesh.ny[j]; i++) {
			y[i + g_index] = y[i - 1 + g_index] + hy;
			hy *= mesh.ky[j];
		}
		y[g_index + mesh.ny[j]] = mesh.subdomain.coords[j + 1].y;
	}

	// заполняем вектор узлов
	mesh.nodes.resize(x.size() * y.size());
	size_t nodenum = 0;
	for (size_t j = 0; j < y.size(); j++) {
		for (size_t i = 0; i < x.size(); i++) {
			mesh.nodes[i + j * x.size()].x = x[i];
			mesh.nodes[i + j * x.size()].y = y[j];
			mesh.nodes[i + j * x.size()].num = nodenum + 1;
			nodenum++;
		}
	}



	// заполняем вектор элементов
	Element el;
	size_t elemnum = 0;
	mesh.elements.resize((x.size() - 1) * (y.size() - 1));
	for (size_t j = 0; j < y.size() - 1; j++) {
		for (size_t i = 0; i < x.size() - 1; i++) {
			mesh.elements[i + j * (x.size() - 1)].loc_nodes[0] = Point(x[i], y[j], mesh.nodes[i + j * x.size()].num);
			mesh.elements[i + j * (x.size() - 1)].loc_nodes[1] = Point(x[i + 1], y[j], mesh.nodes[i + 1 + j * x.size()].num);
			mesh.elements[i + j * (x.size() - 1)].loc_nodes[2] = Point(x[i + 1], y[j + 1], mesh.nodes[i + 1 + (j + 1) * x.size()].num);
			mesh.elements[i + j * (x.size() - 1)].loc_nodes[3] = Point(x[i], y[j + 1], mesh.nodes[i + (j + 1) * x.size()].num);
			mesh.elements[i + j * (x.size() - 1)].mat.num = 1;
			mesh.elements[i + j * (x.size() - 1)].num = elemnum + 1;
			elemnum++;
		}
	}
	// учет пустот в геометрии: удаление ненужных узлов и элементов, их перенумерация
	// удаление узлов и перенумерация оставшихся
	Mesh NewMesh;
	bool is_remove_node = false;
	uint32_t removed_nodes = 0;
	mesh.num_nodes_in_new_mesh.resize(mesh.nodes.size());
	fill(mesh.num_nodes_in_new_mesh.begin(), mesh.num_nodes_in_new_mesh.end(), 0);
	// помечаем узлы на удаление
	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		is_remove_node = !mesh.subdomain.is_contain(mesh.nodes[i]);
		if (is_remove_node) {
			mesh.nodes[i].num = 0;
		}
	}
	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i].num == 0) {
			removed_nodes++;
			continue;
		}
		else {
			NewMesh.nodes.push_back(mesh.nodes[i]);
			NewMesh.nodes[NewMesh.nodes.size() - 1].num -= removed_nodes;
			mesh.num_nodes_in_new_mesh[mesh.nodes[i].num - 1] = mesh.nodes[i].num - removed_nodes;
		}
	}
	mesh.nodes = NewMesh.nodes;

	// удаление лишних элементов, перенумерация оставшихся и перенумерация локальных узлов с учетом удаленных
	removed_nodes = 0;
	uint32_t old_num_node;
	uint32_t new_num_node;
	uint32_t removed_elements = 0;
	bool is_remove_elem;
	// помечаем элементы на удаление
	for (size_t i = 0; i < mesh.elements.size(); i++) {
		is_remove_elem = false;
		for (size_t j = 0; j < 4; j++) {
			is_remove_elem = !mesh.subdomain.is_contain(mesh.elements[i].loc_nodes[j]);
			if (is_remove_elem) {
				is_remove_elem = true;
				mesh.elements[i].mat.num = 0;
				break;
			}
		}
	}

	// перенумерация локальных узлов элементов
	for (size_t i = 0; i < mesh.elements.size(); i++) {
		if (mesh.elements[i].mat.num == 0) {
			for (size_t j = 2; j < 4; j++) {
				if (!mesh.subdomain.is_contain(mesh.elements[i].loc_nodes[j])) {
					removed_nodes++;
					break;
				}
			}
			removed_elements++;
		}
		else {
			NewMesh.elements.push_back(mesh.elements[i]);
			NewMesh.elements[NewMesh.elements.size() - 1].num -= removed_elements;
			for (size_t j = 0; j < 4; j++) {
				old_num_node = NewMesh.elements[NewMesh.elements.size() - 1].loc_nodes[j].num;
				new_num_node = mesh.num_nodes_in_new_mesh[old_num_node - 1];
				NewMesh.elements[NewMesh.elements.size() - 1].loc_nodes[j].num = new_num_node;
			}
		}
	}
	mesh.elements = NewMesh.elements;

	// выгружаем данные сетки в файлы
	ofstream out;
	out.open(filename_nodes);
	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		out << mesh.nodes[i].num << "\t" << mesh.nodes[i].x << "\t" << mesh.nodes[i].y << "\n";
	}
	out.close();
	out.clear();
	out.open(filename_elements);
	for (size_t i = 0; i < mesh.elements.size(); i++) {
		out << mesh.elements[i].num;
		for (size_t j = 0; j < 4; j++) {
			out << "\t" << mesh.elements[i].loc_nodes[j].num;
		}
		out << "\n";
	}
	out.close();
	out.clear();
}





