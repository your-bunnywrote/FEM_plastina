// Mesh.cpp - определение всех функций и операторов, модуль для построения сетки
// Основная литература - "Метод конечных элементов для решения скалярных и векторных задач" [1]


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
	loc_nodes = new Point[4];
}
Element::~Element() {
	delete[] loc_nodes;
}

//========================================================================

void Mesh::calculate_coords(vector<double>& x, vector<double>& y) {
	// вычисление координат по x
	int x_size = (accumulate(nx.begin(), nx.end(), 0) + 1) + (ny[1] - 1) * nx[0];	// сумма кол-ва горизонтальнх подынтервалов + 1, и еще (ny[1] - 1) * nx[0] координат, необходимых из-за смещения узлов дугой
	x.resize(x_size);
	double hx=0,
		hy=0,
		sum = 0;
	double phi = 0;
	// необходимо учитывать, какой интервал бьем: основной или дочерний (образованный параметрами разбиения y - линий), скорее всего, координаты дочерних интервалов придется считать после основных
	// бьем основные горизонтальные интервалы
	for (int k = 0; k < subdomain.horizontal_curves.size(); k++) {
		for (int j = 0; j < subdomain.horizontal_curves[k].size(); j++) {
			Curve current_hor_interval = subdomain.horizontal_curves[k][j];
			sum = 0;
			int g_index = accumulate(nx.begin(), nx.begin() + j, 0);
			for (int i = 0; i < nx[j]; i++) {
				sum += pow(kx[j], i);
			}
			if (current_hor_interval.type == arc) {
				double r = subdomain.hole_radius;
				phi = M_PI_4 / nx[j];
				hx = ((current_hor_interval.end.x - current_hor_interval.begin.x) / sum) * cos(phi);
			}
			else
				hx = (current_hor_interval.end.x - current_hor_interval.begin.x) / sum;
			x[g_index] = current_hor_interval.begin.x;
			for (int i = 1; i <= nx[j]; i++) {
				x[i + g_index] = x[i - 1 + g_index] + hx;
				hx *= kx[j];
			}
			x[g_index + nx[j]] = current_hor_interval.end.x;	// округляем посчитанную координату узла на границе подобластей до точной координаты границы (чтобы избежать значения в формате .000000001)

		}
	}

	int y_size = (accumulate(ny.begin(), ny.end(), 0) + 1) + (nx[1] - 1) * ny[0] * 2;	// сумма кол-ва вертикальных подынтервалов + 1, и еще (nx[1] - 1) * ny[0] * 2 координат, необходимых из-за симметричного смещения узлов дугами
	y.resize(y_size);
	// бьем основные вертикальные интервалы, необходимо учесть тип кривой
	for (int k = 0; k < subdomain.vertical_curves.size(); k++) {
		for (int j = 0; j < subdomain.vertical_curves[k].size(); j++) {
			sum = 0;
			int g_index = accumulate(ny.begin(), ny.begin() + j, 0);
			for (int i = 0; i < ny[j]; i++) {
				sum += pow(ky[j], i);
			}
			hx = (subdomain.vertical_curves[k][j].end.y - subdomain.vertical_curves[k][j].begin.y) / sum;
			y[g_index] = subdomain.vertical_curves[k][j].begin.y;
			for (int i = 1; i <= ny[j]; i++) {
				y[i + g_index] = y[i - 1 + g_index] + hx;
				hx *= ky[j];
			}
			y[g_index + ny[j]] = subdomain.vertical_curves[k][j].end.y;	// округляем посчитанную координату узла на границе подобластей до точной координаты границы (чтобы избежать значения в формате .000000001)

		}
	}
	
	// вычисляем координаты узлов, образованных дочерними ортогональным линиями

}


// =======================================================================

void CreateMesh(Mesh& mesh, string& filename_nodes, string& filename_elements) {
	// информация о разбиении областей
	//mesh.subdomain;
	ifstream Xmsh, Ymsh;
	string line;
	vector<double> x, y;
	double hx, hy, sum;
	// чтение параметров разбиения и вычисление координат узлов с круглым отверстием
	if (mesh.subdomain.is_hole) {
		// чтение файлов с параметрами разбиения
		Xmsh.open(input_folder + "/partition_info_x_holed.txt");
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
		Ymsh.open(input_folder + "/partition_info_y_holed.txt");
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

		// вычисление координат узлов
		mesh.calculate_coords(x, y);





	}
	// чтение параметров разбиения и вычисление координат узлов для панели с прямоугольным отверстием
	else {
		// чтение файлов с разбиением
		// x
		Xmsh.open( input_folder + "/partition_info_x.txt");
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
		Ymsh.open( input_folder + "/partition_info_y.txt");
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
			x[g_index + mesh.nx[j]] = mesh.subdomain.coords[j + 1].x;	// округляем посчитанную координату узла на границе подобластей до точной координаты границы (чтобы избежать значения в формате .000000001)
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
	Mesh NewMesh;		// входит в конструктор comp_domain и снова запрашивает наличие отверстия - ИСПРАВИТЬ!
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

	cout << "Meshing complete!" << endl;
}





