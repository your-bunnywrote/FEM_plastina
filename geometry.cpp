#include "geometry.h"

string input_folder = "test";

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



comp_domain::comp_domain() {

	// так как процесс обработки геометрии нелинейный и сильно зависит от отверстия, будем запрашивать наличие отверстия через консоль
	// если отверстия нет, то структура входных данных будет стандартная, заранее определенная в файле "no_hole_geom.txt" ("subdomains.txt")
	// если отверстие есть, то структура входных данных будет иметь другой вид и формироваться будет программой

	// структура данных для геометрии с отверстием будет примерно такая:
	// четыре ключевые точки - углы прямоугольника, определяющие его длину и ширину - задаются в файле
	// радиус окружности - задается в файле
	// (не входит в структуру) центр окружности - вычисляется через размеры прямоугольника
	// через имеющиеся параметры вычисляются координаты остальных ключевых точек (середины дуг, их проекции на горизонтальные и вертикальные линии) и вносятся в структуру
	// 

	bool is_hole = false;
	cout << "Define the presence of hole (0/1):" << endl;
	cin >> is_hole;
	if (is_hole)
		create_holegeom_info();
	else
		read_noholegeom_info();


}

// чтение файла, описывающего прямоугольную геометрию
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




// создание информации о геометрии пластины с отверстием
void comp_domain::create_holegeom_info() {

	fstream plategeomfile;
	plategeomfile.open(input_folder + "/plate_geom.txt");
	rect_domains.resize(1);		// так как изначально у нас одна прямоугольная область
	double x1, y1, x2, y2;		// координаты левого нижнего и правого верхнего углов прямоугольника
	plategeomfile >> x1 >> y1 >> x2 >> y2;
	plategeomfile >> hole_radius;
	plategeomfile.close();

	// включить обработку симметрии задачи, чтобы считать все необходимые координаты непосредственно в программе и не лезть в файлы
	bool is_symmetric = false;
	cout << "Define task symmetry (0/1):" << endl;
	cin >> is_symmetric;

	if (is_symmetric) {
		Point hole_center(x2, (y2 - y1) / 2, 0);
		//plategeomfile.open(input_folder + "/curves_and_domains.txt");
		// величины постоянные, так как для такой геометрии существует только один вариант обработки используемым методом
		int horizontal_keypoints_count = 3,		// количество точек вдоль оси x (задают верт. линии)
			vertical_keypoints_count = 4;		// количество точек вдоль оси y (задают гор. линии)

		vector<double> kp_x(horizontal_keypoints_count * vertical_keypoints_count, 0.);
		vector<double> kp_y(horizontal_keypoints_count * vertical_keypoints_count, 0.);


		for (int i = 0; i < vertical_keypoints_count * horizontal_keypoints_count; i++) {

			int k = i / horizontal_keypoints_count;

			kp_x[3 * k + 0] = 0.0;
			kp_x[3 * k + 1] = hole_center.x - hole_radius * cos(M_PI_4);
			kp_x[3 * k + 2] = x2;


			int l = i / vertical_keypoints_count;
			double phi = (1 + (l + 1) / 3) * M_PI_4;	// pi/4 для точек первых двух "вертикальных линий", когда точка находится на середине четверти дуги окружности. pi/2 для последней линии, когда точка лежит на оси симметрии
			kp_y[4 * l + 0] = 0.0;
			kp_y[4 * l + 1] = hole_center.y - hole_radius * sin(phi);	// придумать зависимость угла от l такую, чтобы в подсчете последней четверки чисел угол был pi/2
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

		vector<vector<Point>> horizontal_curves(vertical_keypoints_count);
		for (int i = 0; i < horizontal_curves.size(); i++)
			horizontal_curves[i].resize(horizontal_keypoints_count);


		vector<vector<Point>> vertical_curves(horizontal_keypoints_count);
		for (int i = 0; i < vertical_curves.size(); i++) {
			vertical_curves[i].resize(vertical_keypoints_count);
		}

		// строим по ключевым точкам горизонтальные кривые
		for (int i = 0; i < horizontal_curves.size(); i++) {
			for (int j = 0; j < horizontal_curves[i].size(); j++) {
				horizontal_curves[i][j] = keypoints[j + i * horizontal_keypoints_count];
			}
		}

		// строим по клчючевым точкам вертикальные кривые
		for (int i = 0; i < vertical_curves.size(); i++) {
			for (int j = 0; j < vertical_curves[i].size(); j++) {
				vertical_curves[i][j] = keypoints[i + j * horizontal_keypoints_count];
			}
		}


		// test
		curve_type type = arc;

		// забабахать цикл
		vector<Curve> hor_intervals(vertical_keypoints_count * (horizontal_keypoints_count - 1));
		hor_intervals[0] = Curve(keypoints[0], keypoints[1], 0);
		hor_intervals[1] = Curve(keypoints[1], keypoints[2], 1);
		hor_intervals[2] = Curve(keypoints[3], keypoints[4], 2);
		hor_intervals[3] = Curve(arc, hole_center, hole_radius, keypoints[4], keypoints[5], 3);
		hor_intervals[4] = Curve(keypoints[6], keypoints[7], 4);
		hor_intervals[5] = Curve(arc, hole_center, hole_radius, keypoints[7], keypoints[8], 5);
		hor_intervals[6] = Curve(keypoints[9], keypoints[10], 6);
		hor_intervals[7] = Curve(keypoints[10], keypoints[11], 7);

		//забабахать цикл
		vector<vector<Curve>> horizontal_lines(vertical_keypoints_count);
		for (int i = 0; i < horizontal_curves.size(); i++) {
			horizontal_curves[i].resize(horizontal_keypoints_count - 1);
		}

		int num = 0;
		for (int i = 0; i < horizontal_lines.size(); i++) {
			int j = 0;
			horizontal_lines[i][j] = Curve(keypoints[i * horizontal_keypoints_count], keypoints[i * horizontal_keypoints_count + 1], num);
			for (j = 1; j < horizontal_curves[i].size(); j++) {
				if ((i == 1 || i == 2) && j == 1)
					horizontal_lines[i][j] = Curve(arc, hole_center, hole_radius, horizontal_lines[i][j - 1].end, keypoints[j + 1], num++);
				else
					horizontal_lines[i][j] = Curve(horizontal_lines[i][j - 1].end, keypoints[j + 1], num++);	// изменить формулу индексирования конечной точки
			}
		}
		horizontal_lines[0][0] = Curve(keypoints[0], keypoints[1], 0);					// номер интервала - сумма индексов
		horizontal_lines[0][1] = Curve(keypoints[1], keypoints[2], 1);
		horizontal_lines[1][0] = Curve(keypoints[3], keypoints[4], 2);
		horizontal_lines[1][1] = Curve(arc, hole_center, hole_radius, keypoints[4], keypoints[5], 3);
		horizontal_lines[2][0] = Curve(keypoints[6], keypoints[7], 4);
		horizontal_lines[2][1] = Curve(arc, hole_center, hole_radius, keypoints[7], keypoints[8], 5);
		horizontal_lines[3][0] = Curve(keypoints[9], keypoints[10], 6);
		horizontal_lines[3][1] = Curve(keypoints[10], keypoints[11], 7);



		vector<Curve> vert_intervals(horizontal_keypoints_count*(vertical_keypoints_count - 1));
		vert_intervals[0] = Curve(keypoints[0], keypoints[0 + horizontal_keypoints_count], 0);
		vert_intervals[1] = Curve(vert_intervals[0].end, keypoints[2 * horizontal_keypoints_count], 1);
		vert_intervals[2] = Curve(vert_intervals[1].end, keypoints[3 * horizontal_keypoints_count], 2);
		vert_intervals[3] = Curve(keypoints[1], keypoints[1 + horizontal_keypoints_count], 3);
		vert_intervals[4] = Curve(arc, hole_center, hole_radius, keypoints[1 + horizontal_keypoints_count], keypoints[1 + 2 * horizontal_keypoints_count], 4);
		vert_intervals[5] = Curve(vert_intervals[4].end, keypoints[1 + 3 * horizontal_keypoints_count], 5);
		vert_intervals[6] = Curve(vert_intervals[5].end, keypoints[2 + horizontal_keypoints_count], 6);
		vert_intervals[7] = Curve(vert_intervals[6].end, keypoints[2 + 2 * horizontal_keypoints_count], 7);
		vert_intervals[8] = Curve(vert_intervals[7].end, keypoints[2 + 3 * horizontal_keypoints_count], 8);

		// Описание элементарных подобластей. Для рассматриваемой геометрии данные постоянны:
		// [1,стр. 461]
		// количество подобластей
		// первое число - номер материала, второе число - номер вертикальной кривой, ограничивающей подобласть слева, третье - номер вертикальной кривой, ограничивающей подобласть справа
		// четвертое - номер горизонтальной кривой, ограничивающей подобласть снизу, пятое - номер горизонтальной кривой, ограничивающей подобласть сверху 
		int W_count = 3;
		int mat = 1;
		vector<pair<int, int>> vertical_curves_indeces{ {0,1},{1,2},{1,2} };
		vector<pair<int, int>> horizontal_lines_indeces{ {0,3},{0,1},{2,3} };

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
			ofile << mat << " " << vertical_curves_indeces[k].first << " " << vertical_curves_indeces[k].second << " " << horizontal_lines_indeces[k].first << " " << horizontal_lines_indeces[k].second << endl;
			k++;
		}




		ofile.close();

	}
	else {

		nonsymmetric_hole_geom();
	}

	// придумать, как можно в этот массив внести координаты вертикальных и горизонтальных границ подобластей, но прежде нужно обработать информацию об окружности
	this->coords;

}

void comp_domain::nonsymmetric_hole_geom() {


}


bool comp_domain::is_contain(const Point& node) {
	for (size_t i = 0; i < rect_domains.size(); i++) {
		if ((node.x >= rect_domains[i].first.x) && (node.x <= rect_domains[i].second.x + 1) && (node.y >= rect_domains[i].first.y) && (node.y <= rect_domains[i].second.y))
			return true;
	}
	return false;
}


