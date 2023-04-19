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
		domains[i].second.x = coords[x1].x;
		domains[i].second.y = coords[y1].y;
	}
	W.close();
}

bool comp_domain::is_contain(const Point& node) {
	for (size_t i = 0; i < domains.size(); i++) {
		if ((node.x >= domains[i].first.x) && (node.x >= domains[i].second.x) && (node.y >= domains[i].first.y) && (node.y >= domains[i].second.y))
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