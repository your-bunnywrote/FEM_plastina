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

//========================================================================

void comp_domain::readfile_domains() {
	ifstream W;
	W.open(".\\text_files\\subdomains.txt");

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
		if ((node.x >= domains[i].first.x) && (node.x <= domains[i].second.x) && (node.y >= domains[i].first.y) && (node.y <= domains[i].second.y+1))
			return true;
	}
	return false;
}

// =======================================================================

void Mesh::readfile_partition_info() {
	subdomain.readfile_domains();
	ifstream Xmsh, Ymsh;
	Xmsh.open(".\\text_files\\partition_info_x.txt");
	string line;
	while (getline(Xmsh, line)) {
		vector<string> tokens = split(line, '\t');
		nx.push_back(stoi(tokens.at(0)));
		kx.push_back(stod(tokens.at(1)));
	}
	Xmsh.close();
	line.clear();
	Ymsh.open(".\\text_files\\partition_info_y.txt");
	while (getline(Ymsh, line)) {
		vector<string> tokens = split(line, '\t');
		ny.push_back(stoi(tokens.at(0)));
		ky.push_back(stod(tokens.at(1)));
	}
	Ymsh.close();
	line.clear();
}

vector<double> Mesh::calc_nodes_x_coords(const vector<double>& k, const vector<int>& n, vector<Point>& coordXw) {
	vector<double> x;
	double h, sum;
	x.resize(accumulate(n.begin(), n.end(), 0) + 1);
	for (size_t j = 0; j < n.size(); j++) {
		sum = 0;
		int g_index = accumulate(n.begin(), n.begin() + j, 0);
		for (int i = 0; i < n[j]; i++) {
			sum += pow(k[j], i);
		}
		h = (coordXw[j + 1].x - coordXw[j].x) / sum;
		x[g_index] = coordXw[j].x;
		for (int i = 1; i <= n[j]; i++) {
			x[i + g_index] = x[i - 1 + g_index] + h;
			h *= k[j];
		}
	}
	return x;
}

vector<double> Mesh::calc_nodes_y_coords(const vector<double>& k, const vector<int>& n, vector<Point>& coordYw) {
	vector<double> y;
	double h, sum;

	y.resize(accumulate(n.begin(), n.end(), 0) + 1);
	for (int j = 0; j < n.size(); j++) {
		sum = 0;
		int g_index = accumulate(n.begin(), n.begin() + j, 0);
		for (int i = 0; i < n[j]; i++) {
			sum += pow(k[j], i);
		}
		h = (coordYw[j + 1].y - coordYw[j].y) / sum;
		y[g_index] = coordYw[j].y;
		for (int i = 1; i <= n[j]; i++) {
			y[i + g_index] = y[i - 1 + g_index] + h;
			h *= k[j];
		}
	}
	return y;
}

vector<Point> Mesh::fill_nodes() {
	readfile_partition_info();
	vector<double> x = calc_nodes_x_coords(kx, nx, subdomain.coords);
	vector<double> y = calc_nodes_y_coords(ky, ny, subdomain.coords);
	nodes.resize(x.size() * y.size());
	size_t nodenum = 0;
	for (size_t j = 0; j < y.size(); j++) {
		for (size_t i = 0; i < x.size(); i++) {
			nodes[i + j * x.size()].x = x[i];
			nodes[i + j * x.size()].y = y[j];
			nodes[i + j * x.size()].num = nodenum + 1;
			nodenum++;
		}
	}
	return nodes;
}

vector<Element> Mesh::fill_elements() {
	//readfile_partition_info();
	vector<double> x = calc_nodes_x_coords(kx, nx, subdomain.coords);
	vector<double> y = calc_nodes_y_coords(ky, ny, subdomain.coords);
	Element el;
	size_t elemnum = 0;
	elements.resize((x.size() - 1) * (y.size() - 1));
	for (size_t j = 0; j < y.size() - 1; j++) {
		for (size_t i = 0; i < x.size() - 1; i++) {
			elements[i + j * (x.size() - 1)].loc_nodes[0] = Point(x[i], y[j], nodes[i + j * x.size()].num);
			elements[i + j * (x.size() - 1)].loc_nodes[1] = Point(x[i + 1], y[j], nodes[i + 1 + j * x.size()].num);
			elements[i + j * (x.size() - 1)].loc_nodes[2] = Point(x[i + 1], y[j + 1], nodes[i + 1 + (j + 1) * x.size()].num);
			elements[i + j * (x.size() - 1)].loc_nodes[3] = Point(x[i], y[j + 1], nodes[i + (j + 1) * x.size()].num);
			elements[i + j * (x.size() - 1)].material = 1;
			elements[i + j * (x.size() - 1)].num = elemnum + 1;
			elemnum++;
		}
	}
	return elements;
}

void Mesh::output(const string& filename, const vector<Point>& nodes) {
	ofstream out;
	out.open(filename);
	for (size_t i = 0; i < nodes.size(); i++) {
		out << nodes[i].num << "\t" << nodes[i].x << "\t" << nodes[i].y << "\n";
	}
	out.close();
}

void Mesh::output(const string& filename, const vector<Element>& elements) {
	ofstream out;
	out.open(filename);
	for (size_t i = 0; i < elements.size(); i++) {
		out << elements[i].num;
		for (size_t j = 0; j < 4; j++) {
			out << "\t" << elements[i].loc_nodes[j].num;
		}
		out << "\n";
	}
	out.close();
}

vector<Point> Mesh::edit_nodes(vector<Point>& nodes) {
	Mesh NewMesh;
	bool is_remove_node = false;
	int removed_nodes = 0;
	num_nodes_in_new_mesh.resize(nodes.size());
	fill(num_nodes_in_new_mesh.begin(), num_nodes_in_new_mesh.end(), 0);
	// помечаем узлы на удаление
	for (size_t i = 0; i < nodes.size(); i++) {
		is_remove_node = !subdomain.is_contain(nodes[i]);
		if (is_remove_node) {
			nodes[i].num = 0;
		}
	}
	for (size_t i = 0; i < nodes.size(); i++) {
		if (nodes[i].num == 0) {
			removed_nodes++;
			continue;
		}
		else {
			NewMesh.nodes.push_back(nodes[i]);
			NewMesh.nodes[NewMesh.nodes.size() - 1].num -= removed_nodes;
			num_nodes_in_new_mesh[nodes[i].num - 1] = nodes[i].num - removed_nodes;
		}
	}
	return NewMesh.nodes;
}

vector<Element> Mesh::edit_elements(vector<Element>& elements) {
	Mesh NewMesh;
	uint32_t old_num_node;
	uint32_t new_num_node;
	uint32_t removed_nodes = 0;;
	uint32_t removed_elements = 0;
	bool is_remove_elem;
	// помечаем элементы на удаление
	for (size_t i = 0; i < elements.size(); i++) {
		is_remove_elem = false;
		for (size_t j = 0; j < 4; j++) {
			is_remove_elem = !subdomain.is_contain(elements[i].loc_nodes[j]);
			if (is_remove_elem) {
				is_remove_elem = true;
				elements[i].material = 0;
				break;
			}
		}
	}
	for (size_t i = 0; i < elements.size(); i++) {
		if (elements[i].material == 0) {
			for (size_t j = 2; j < 4; j++) {
				if (!subdomain.is_contain(elements[i].loc_nodes[j])) {
					removed_nodes++;
					break;
				}
			}
			removed_elements++;
		}
		else {
			NewMesh.elements.push_back(elements[i]);
			NewMesh.elements[NewMesh.elements.size() - 1].num -= removed_elements;
			for (size_t j = 0; j < 4; j++) {
				old_num_node = NewMesh.elements[NewMesh.elements.size() - 1].loc_nodes[j].num;
				new_num_node = num_nodes_in_new_mesh[old_num_node - 1];
				NewMesh.elements[NewMesh.elements.size() - 1].loc_nodes[j].num = new_num_node;
			}
		}
	}
	return NewMesh.elements;

	
}



