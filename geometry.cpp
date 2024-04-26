#include "geometry.h"

string input_folder = "test";

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


Curve::Curve() {};

Curve::Curve(Point begin_, Point end_, int num_) {
	type = line;
	begin = begin_;
	end = end_;
	num = num_;
	center = Point(0, 0);
	radius = 0;
}

Curve::Curve(curve_type type_, Point center_, double radius_, Point begin_, Point end_, int num_) {
	type = type_;
	center = center_;
	radius = radius_;
	begin = begin_;
	end = end_;
	num = num_;
}


comp_domain::comp_domain(){}
comp_domain::comp_domain(bool is_hole) {

	// ��� ��� ������� ��������� ��������� ���������� � ������ ������� �� ���������, ����� ����������� ������� ��������� ����� �������
	// ���� ��������� ���, �� ��������� ������� ������ ����� �����������, ������� ������������ � ����� "no_hole_geom.txt" ("subdomains.txt")
	// ���� ��������� ����, �� ��������� ������� ������ ����� ����� ������ ��� � ������������� ����� ����������

	// ��������� ������ ��� ��������� � ���������� ����� �������� �����:
	// ������ �������� ����� - ���� ��������������, ������������ ��� ����� � ������ - �������� � �����
	// ������ ���������� - �������� � �����
	// (�� ������ � ���������) ����� ���������� - ����������� ����� ������� ��������������
	// ����� ��������� ��������� ����������� ���������� ��������� �������� ����� (�������� ���, �� �������� �� �������������� � ������������ �����) � �������� � ���������
	// 


	if (is_hole)
		create_holegeom_info();
	else
		read_noholegeom_info();


}

// ������ �����, ������������ ������������� ���������
void comp_domain::read_noholegeom_info() {
	ifstream W;
	W.open(input_folder + "/subdomains.txt");

	W >> Nx;
	double x;
	for (size_t i = 0; i < Nx; i++) {
		W >> x;
		coords.push_back(Point(x, -1., 0));
	}

	length = (coords[Nx - 1].x - coords[0].x);

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

	width = (coords[Ny - 1].y - coords[0].y);

	int W_count;
	W >> W_count;
	rect_domains.resize(W_count);
	double x1, x2, y1, y2;
	int mat;
	for (size_t i = 0; i < W_count; i++) {
		W >> mat >> x1 >> x2 >> y1 >> y2;
		rect_domains[i].first.x = coords[x1].x;
		rect_domains[i].first.y = coords[y1].y;
		rect_domains[i].second.x = coords[x2].x;
		rect_domains[i].second.y = coords[y2].y;
	}
	W.close();
}




// �������� ���������� � ��������� �������� � ����������
void comp_domain::create_holegeom_info() {

	fstream plategeomfile;
	plategeomfile.open(input_folder + "/plate_geom.txt");
	rect_domains.resize(1);		// ��� ��� ���������� � ��� ���� ������������� �������
	double x1, y1, x2, y2;		// ���������� ������ ������� � ������� �������� ����� ��������������
	plategeomfile >> x1 >> y1 >> x2 >> y2;
	plategeomfile >> hole_radius;
	plategeomfile.close();

	// �������� ��������� ��������� ������, ����� ������� ��� ����������� ���������� ��������������� � ��������� � �� ����� � �����
	bool is_symmetric = false;
	cout << "Define task symmetry (0/1):" << endl;
	cin >> is_symmetric;

	if (is_symmetric) {
		hole_center = Point(x2, (y2 - y1) / 2, 0);
		//plategeomfile.open(input_folder + "/curves_and_domains.txt");
		// �������� ����������, ��� ��� ��� ����� ��������� ���������� ������ ���� ������� ��������� ������������ �������
		int horizontal_keypoints_count = 3,		// ���������� ����� ����� ��� x (������ ����. �����)
			vertical_keypoints_count = 4;		// ���������� ����� ����� ��� y (������ ���. �����)

		vector<double> kp_x(horizontal_keypoints_count * vertical_keypoints_count, 0.);
		vector<double> kp_y(horizontal_keypoints_count * vertical_keypoints_count, 0.);


		for (int i = 0; i < vertical_keypoints_count * horizontal_keypoints_count; i++) {

			int k = i / horizontal_keypoints_count;

			kp_x[3 * k + 0] = 0.0;
			kp_x[3 * k + 1] = hole_center.x - hole_radius * cos(M_PI_4);
			kp_x[3 * k + 2] = x2;


			int l = i / vertical_keypoints_count;
			double phi = (1 + (l + 1) / 3) * M_PI_4;	// pi/4 ��� ����� ������ ���� "������������ �����", ����� ����� ��������� �� �������� �������� ���� ����������. pi/2 ��� ��������� �����, ����� ����� ����� �� ��� ���������
			kp_y[4 * l + 0] = 0.0;
			kp_y[4 * l + 1] = hole_center.y - hole_radius * sin(phi);
			kp_y[4 * l + 1] = hole_center.y - hole_radius * sin(phi);
			kp_y[4 * l + 2] = hole_center.y + hole_radius * sin(phi);
			kp_y[4 * l + 3] = y2;
		}

		vector<Point> keypoints(horizontal_keypoints_count * vertical_keypoints_count);

		for (int j = 0; j < vertical_keypoints_count; j++) {
			for (int i = 0; i < horizontal_keypoints_count; i++) {
				keypoints[i + j * horizontal_keypoints_count].x = kp_x[i];
				keypoints[i + j * horizontal_keypoints_count].y = kp_y[4 * i + j];
			}
		}



		// test
		curve_type type = arc;

		horizontal_curves.resize(vertical_keypoints_count);
		for (int i = 0; i < horizontal_curves.size(); i++) {
			horizontal_curves[i].resize(horizontal_keypoints_count - 1);
		}

		int hor_interval_num = 0;
		for (int i = 0; i < horizontal_curves.size(); i++) {
			int j = 0;
			int kp_index = i * horizontal_keypoints_count;
			horizontal_curves[i][j] = Curve(keypoints[kp_index], keypoints[kp_index + 1], hor_interval_num++);
			for (j = 1; j < horizontal_curves[i].size(); j++) {
				if ((i == 1 || i == 2) && j == 1)
					horizontal_curves[i][j] = Curve(arc, hole_center, hole_radius, horizontal_curves[i][j - 1].end, keypoints[kp_index + 1 + j], hor_interval_num++);
				else
					horizontal_curves[i][j] = Curve(horizontal_curves[i][j - 1].end, keypoints[kp_index + 1 + j], hor_interval_num++);	// �������� ������� �������������� �������� �����


			}
		}

		vertical_curves.resize(horizontal_keypoints_count);
		for (int i = 0; i < vertical_curves.size(); i++)
			vertical_curves[i].resize(vertical_keypoints_count - 1);

		int ver_interval_num = 0;
		for (int i = 0; i < vertical_curves.size(); i++) {
			int j = 0;
			vertical_curves[i][j] = Curve(keypoints[i], keypoints[i+horizontal_keypoints_count], ver_interval_num++);
			for (j = 1; j < vertical_curves[i].size(); j++) {
				if (i == 1 && j == 1)
					vertical_curves[i][j] = Curve(arc, hole_center, hole_radius, keypoints[i + j * horizontal_keypoints_count], keypoints[i + (j + 1) * horizontal_keypoints_count], ver_interval_num++);
				else
				vertical_curves[i][j] = Curve(keypoints[i + j * horizontal_keypoints_count], keypoints[i + (j + 1) * horizontal_keypoints_count], ver_interval_num++);
			}

		}


		// �������� ������������ �����������. ��� ��������������� ��������� ������ ���������:
		// [1,���. 461]
		// ���������� �����������
		// ������ ����� - ����� ���������, ������ ����� - ����� ������������ ������, �������������� ���������� �����, ������ - ����� ������������ ������, �������������� ���������� ������
		// ��������� - ����� �������������� ������, �������������� ���������� �����, ����� - ����� �������������� ������, �������������� ���������� ������ 
		int W_count = 3;
		int mat = 1;
		vector<pair<int, int>> vertical_curves_indeces{ {0,1},{1,2},{1,2} };
		vector<pair<int, int>> horizontal_curves_indeces{ {0,3},{0,1},{2,3} };

		ofstream ofile(input_folder + "/curves_and_domains.txt");
		ofile << horizontal_keypoints_count << " " << vertical_keypoints_count << endl;
		for (int j = 0; j < vertical_keypoints_count; j++) {
			for (int i = 0; i < horizontal_keypoints_count; i++) {
				ofile << keypoints[i + horizontal_keypoints_count * j].x << " " << keypoints[i + horizontal_keypoints_count * j].y << ((i == horizontal_keypoints_count - 1) ? ("\n") : (" "));
			}
		}
		ofile << W_count << endl;
		int k = 0;
		for (int i = 0; i < W_count; i++) {
			ofile << mat << " " << vertical_curves_indeces[k].first << " " << vertical_curves_indeces[k].second << " " << horizontal_curves_indeces[k].first << " " << horizontal_curves_indeces[k].second << endl;
			k++;
		}




		ofile.close();

	}
	else {

		nonsymmetric_hole_geom();
	}

}

void comp_domain::nonsymmetric_hole_geom() {


}


bool comp_domain::is_contain(const Point& node) {
	if (!is_hole) {
		for (size_t i = 0; i < rect_domains.size(); i++) {
			if ((node.x >= rect_domains[i].first.x) && (node.x <= rect_domains[i].second.x + 1) && (node.y >= rect_domains[i].first.y) && (node.y <= rect_domains[i].second.y))
				return true;
		}
		return false;
	}
	else {
		// ����������� �� ��������� ����������: (x - center.x)^2 + (y-center.y)^2 = R^2
		// ���� ����� ����� ������ ������, ���� �������� � ����������, �������������, �� �������� � ���������� � �������� ��������
		// ���� ����� ����� ������ ���� ����� ������, ���� �� �������� � ���� � �������� � ������

		// ��������� ����� �����, ��� ��� ���������� ��� ����� ���� double
		int number_of_dicimal = 1e4;
		double rounded = round( ((node.x - hole_center.x) * (node.x - hole_center.x) + (node.y - hole_center.y) * (node.y - hole_center.y)) * number_of_dicimal)/number_of_dicimal;
		if (rounded < hole_radius * hole_radius)
			return false;
		return true;
	}

}


