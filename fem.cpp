// fem.cpp - ����������� ���� ����� ����������� �������

#include "fem.h"

double Rectangle::gauss_points_local[3] = { -0.77459666924148337703585307995648, 0, 0.77459666924148337703585307995648 };
double Rectangle::gauss_weights[3] = { 5. / 9., 8. / 9., 5. / 9. };
Point Rectangle::integrate_points[9] = {};


Material::Material() {
	E = 2.1e5;		// [���] ��� ��� ��� ��������� �������� � ����������, ������ ���� �������� � ���
	mu = 0.3;
	c = E / (1 - mu * mu);
	thickness = 3; // [��]
}

Load::Load() {
	LineLength = 0.0;
	Px = -1000;							// ��������, ������� ������������ � �����, [�/��]
	Py = 0.0;
}
void Load::GetLineLength(Mesh mesh) {
	double x1, x2, y1, y2;				// ���������� �����, ����������� �����, ����� ������� ����� ���������
	x1 = mesh.subdomain.coords.begin()->x;
	y1 = mesh.subdomain.coords.begin()->y;
	x2 = (mesh.subdomain.coords.end() - 1)->x;
	y2 = (mesh.subdomain.coords.end() - 1)->y;
	LineLength = sqrt((x1 - x1) * (x1 - x1) + (y2 - y1) * (y2 - y1));		// ����� ����� ���� ���������� �� �����
}

block2x2::block2x2() {
	this->val11 = 0.0;
	this->val12 = 0.0;
	this->val21 = 0.0;
	this->val22 = 0.0;
}

block2x2::block2x2(double val11, double val12, double val21, double val22) {
	this->val11 = val11;
	this->val12 = val12;
	this->val21 = val21;
	this->val22 = val22;
}

block2x2 block2x2::operator+(block2x2 block2) {
	double val11_res, val12_res,
		   val21_res, val22_res;
	val11_res = val11 + block2.val11;
	val12_res = val12 + block2.val12;
	val21_res = val21 + block2.val21;
	val22_res = val22 + block2.val22;
	return block2x2(val11_res, val12_res, val21_res, val22_res);
}

block2x2 block2x2::operator += (block2x2 block2) {
	val11 += block2.val11;
	val12 += block2.val12;
	val21 += block2.val21;
	val22 += block2.val22;
	return block2x2(val11,val12,val21,val22);
}

block2x2 block2x2::operator = (double val) {
	val11 = val;
	val12 = val;
	val21 = val;
	val22 = val;
	return block2x2(val11, val12, val21, val22);
}

block1x2::block1x2() {
	this->val1 = 0.0;
	this->val2 = 0.0;
}

block1x2::block1x2(double val1, double val2) {
	this->val1 = val1;
	this->val2 = val2;
}

block1x2 block1x2::operator+(block1x2 block2) {
	double val1_res, val2_res;
	val1_res = val1 + block2.val1;
	val2_res = val2 + block2.val2;
	return block1x2(val1_res, val2_res);
}

block1x2 block1x2::operator += (block1x2 block2) {
	val1 += block2.val1;
	val2 += block2.val2;
	return block1x2(val1, val2);
}

block1x2 block1x2::operator = (double val) {
	val1 = val;
	val2 = val;
	return block1x2(val1, val2);
}

Rectangle::Rectangle() {
	
	LocalStiffnessMatrix_block.resize(4);
	for (int i = 0; i < 4; i++) {
		LocalStiffnessMatrix_block[i].resize(4);
	}
	LocalStiffnessMatrix.resize(8);
	for (int i = 0; i < 8; i++) {
		LocalStiffnessMatrix[i].resize(8);
	}
	LocalLoadVector_block.resize(4);
	LocalLoadVector.resize(8);
}

ostream& operator<<(ostream& os, const block2x2& block) {
	os << block.val11 << "\t" << block.val12 << "\n"
		<< block.val21 << "\t" << block.val22;
	return os;
};

ostream& operator<<(ostream& os, const block1x2& block) {
	os << block.val1 << "\n"
		<< block.val2 << "\n\n";
	return os;
};


void Rectangle::init()
{
	int p_id = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			integrate_points[p_id] = Point(gauss_points_local[j], gauss_points_local[i]);
			p_id++;
		}
	}


}

double Rectangle::bfunc1D_1(double x0, double x1, double x) {
	return (x1 - x) / (x1 - x0);
}

double Rectangle::bfunc1D_2(double x0, double x1, double x) {
	return (x - x0) / (x1 - x0);
}

double Rectangle::dbfunc1D_1(double x0, double x1, double x) {
	return -1.0 / (x1 - x0);
}
double Rectangle::dbfunc1D_2(double x0, double x1, double x) {
	return 1.0 / (x1 - x0);
}

double Rectangle::param_func_x(double t, Point& from, Point& to) {
	return from.x + t* (to.x - from.x);
}
double Rectangle::param_func_y(double t, Point& from, Point& to) {
	return from.y + t * (to.y - from.y);
}



double Rectangle::phi(size_t func_num, Point& from, Point& to, double x, double y) {

	switch (func_num) {
	case 0:
		return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y);	// phi1(x,y) = X1(x)*Y1(y)
	case 1:
		return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y);	// phi2(x,y) = X2(x)*Y1(y)
	case 2:
		return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y);	// phi3(x,y) = X2(x)*Y2(y)
	case 3:
		return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y);	// phi4(x,y) = X1(x)*Y2(y)
	}
}

double Rectangle::parametric_phi(size_t func_num, double t) {
	switch (func_num) {
	case 0: 
		return (1 - t) * (1 - t);
	case 1:
		return t * (1 - t);
	case 2:
		return (1 - t) * t;
	case 3:
		return t * t;
	}
}


 // �������� ��� ������������ Kvar[num1][num2], ��� var - ����������, �� ������� �������������� ������� phi(x,y)
 // var = 0 - x, var = 1 - y, var = 2 - xy, var = 3 - yx;
double Rectangle::dphi(size_t var, size_t num1, size_t num2, Point& from, Point& to, double x, double y) {
	switch (var) {
	// Kx[num1][num2] = dphi[num1]/dx * dphi[num2]/dx
	case 0:
		switch (num1) {
		case 0:
			switch (num2) {
			case 0: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi1/dx = dX1/dx * dX1/dx * Y1 * Y1
			case 1: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi2/dx = dX1/dx * dX2/dx * Y1 * Y1
			case 2: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi1/dx * dphi3/dx = dX1/dx * dX2/dx * Y1 * Y2
			case 3: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dx = dX1/dx * dX1/dx * Y1 * Y2
			}
		case 1:
			switch (num2) {
			case 0: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi1/dx = dX2/dx * dX1/dx * Y1 * Y1
			case 1: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi2/dx = dX2/dx * dX2/dx * Y1 * Y1
			case 2: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi3/dx = dX2/dx * dX2/dx * Y1 * Y2
			case 3: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dx = dX2/dx * dX2/dx * Y2 * Y2
			}
		case 2:
			switch (num2) {
			case 0: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi1/dx = dX2/dx * dX1/dx * Y2 * Y1
			case 1: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi2/dx = dx2/dx * dX2/dx * Y2 * Y1
			case 2: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi3/dx = dX2/dx * dX2/dx * Y2 * Y2
			case 3: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi4/dx = dX2/dx * dX1/dx * Y2 * Y2
			}
		case 3:
			switch (num2) {
			case 0: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi1/dx = dX1/dx * dX1/dx * Y2 * Y1
			case 1: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi2/dx = dX1/dx * dX2/dx * Y2 * Y1
			case 2: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi3/dx = dX1/dx * dX2/dx * Y2 * Y2
			case 3: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi4/dx * dX1/dx * dX1/dx * Y2 * Y2
			}
		}

	// Ky[num1][num2] = dphi[num1]/dy * dphi[num2]dy
	case 1:
		switch (num1) {
		case 0:
			switch (num2) {
			case 0: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi1/dy * dphi1/dy = X1 * X1 * dY1 * dY1
			case 1: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi1/dy * dphi2/dy = X1 * X2 * dY1 * dY1
			case 2: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi3/dy = X1 * X2 * dY1 * dY2
			case 3: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi4/dy = X1 * X1 * dY1 * dY2
			}
		case 1:
			switch (num2) {
			case 0: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi1/dy = X2 * X1 * dY1 * dY1
			case 1: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi2/dy = X2 * X2 * dY1 * dY1
			case 2: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi3/dy = X2 * X2 * dY1 * dY2
			case 3: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi4/dy = X2 * X1 * dY1 * dY2
			}
		case 2:
			switch (num2) {
			case 0: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi3/dy * dphi1/dy = X2 * X1 * dY2 * dY1
			case 1: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi3/dy * dphi2/dy = X2 * X2 * dY2 * dY1
			case 2: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi3/dy = X2 * X2 * dY2 * dY2
			case 3: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi4/dy = X2 * X1 * dY2 * dY2
			}
		case 3:
			switch (num2) {
			case 0: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi4/dy * dphi1/dy = X1 * X1 * dY2 * dY1
			case 1: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi2/dy = X1 * X2 * dY2 * dY1
			case 2: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi3/dy = X1 * X2 * dY2 * dY2
			case 3: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi4/dy = X1 * X1 * dY2 * dY2
			}
		}
	// Kxy[num1][num2] = dphi[num1]/dx * dphi[num2]/dy
	case 2:
		switch (num1) {
		case 0:
			switch (num2) {
			case 0: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi1/dx * dphi1/dy = dX1 * X1 * Y1 * dY1
			case 1: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi2/dx = dX1 * X2 * Y1 * dY1
			case 2: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dx * dphi3/dy = dX1 * X2 * Y1 * dY2
			case 3: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y); // dphi1/dx * dphi4/dy = dX1 * X1 * Y2 * dY2
			}
		case 1:
			switch (num2) {
			case 0: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi1/dy = dX2 * X1 * Y1 * dY1
			case 1: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi2/dy = dX2 * X2 * Y1 * dY1
			case 2: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi3/dy = dX2 * X2 * Y1 * dY2
			case 3: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dy = dX1 * X1 * Y1 * dY2
			}
		case 2:
			switch (num2) {
			case 0: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi3/dx * dphi1/dy = dX2 * X1 * Y2 * dY1
			case 1: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi2/dy = dX2 * X2 * Y1 * dY2
			case 2: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi3/dy = dX2 * X2 * Y2 * dY2
			case 3: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y); // dphi3/dx * dphi4/dy = dX2 * X1 * Y1 * dY2
			}
		case 3:
			switch (num2) {
			case 0: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi1/dy = dX1 * X1 * Y2 * dY1
			case 1: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi2/dy = dX1 * X2 * Y2 * dY1
			case 2: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi3/dy = dX1 * X2 * Y2 * dY2
			case 3: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi4/dy = dX1 * X1 * Y2 * dY2
			}
		}
	// Kyx[num1][num2] = dphi[num1]/dy * dphi[num2]/dx
	case 3:
		switch (num1) {
		case 0:
			switch (num2) {
			case 0: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y); // dphi1/dy * dphi1/dx = X1 * dX1 * dY1 * Y1
			case 1: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dy * dphi2/dx = X1 * dX2 * dY1 * Y1
			case 2: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi3/dx = X1 * dX2 * dY1 * Y2
			case 3: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y); // dphi1/dy * dphi4/dx = X1 * dX1 * dY2 * Y2
			}
		case 1:
			switch (num2) {
			case 0: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi1/dx = X2 * dX1 * dY1 * Y1
			case 1: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi2/dx = X2 * dX2 * dY1 * Y1
			case 2: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi3/dx = X2 * dX2 * dY1 * Y2
			case 3: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi4/dx = X1 * dX1 * dY1 * Y2
			}
		case 2:
			switch (num2) {
			case 0: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y); // dphi3/dy * dphi1/dx = X2 * dX1 * dY2 * Y1
			case 1: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dy * dphi2/dx = X2 * dX2 * dY1 * Y2
			case 2: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi3/dx = X2 * dX2 * dY2 * Y2
			case 3: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y); // dphi3/dy * dphi4/dx = X2 * dX1 * dY1 * Y2
			}
		case 3:
			switch (num2) {
			case 0: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi1/dx = X1 * dX1 * dY2 * Y1
			case 1: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi2/dx = X1 * dX2 * dY2 * Y1
			case 2: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi3/dx = X1 * dX2 * dY2 * Y2
			case 3: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi4/dx = X1 * dX1 * dY2 * Y2
			}
		}
	}
}

double Rectangle::Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) {
	this->init();
	int p_id = 0;
	double res = 0.0;
	double x_centre = (to.x + from.x) / 2., 
		   y_centre = (to.y + from.y) / 2.;
	double hx = to.x - from.x,
		   hy = to.y - from.y;
	double ksi, eta;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			// ��������� � ���������� ���������� ��������
			ksi = x_centre + integrate_points[p_id].x * hx / 2.;
			eta = y_centre + integrate_points[p_id].y * hy / 2.;
			// � ����������� � ���������� �����������
			res += gauss_weights[i] * gauss_weights[j] * dphi(var, num1, num2, from, to, ksi, eta);
			p_id++;
		}
	}
	return res * hx * hy / 4.0;	// hx*hy/4 - ������� �������� �� ��������� ��������� �������������� � ���������� � ��������� ������
}



double Rectangle::Edge_IntegrateGauss3(Point& from, Point& to, size_t num) {
	double res = 0.0;
	//double hx = to.x - from.x;
	double hy = to.y - from.y;
	//double x_centre = (to.x+from.x) / 2;
	double y_centre = (to.y+from.y) / 2.;
	double ksi,eta;
	for (int g = 0; g < 3; g++) {
		//ksi = x_centre + gauss_points_local[g] * hx/2;
		eta = y_centre + gauss_points_local[g] * hy / 2.;
		switch (num) {
		case 0: res += hy / 2. * gauss_weights[g] * bfunc1D_1(from.y, to.y, eta); break;	// ���� �������� ��������� � ����� ��� ������ ������ (����������� �) �� ����������� ���������� �.�. Y1
		//case 1: res += hx/2. * gauss_weights[g] * bfunc1D_1(from.x, to.x, ksi); break;	/// ���� �������� ��������� � ������� ��� ������ ������ (����������� �), ����������� ���������� �.�. X1
		//case 2: res += hx/2. * gauss_weights[g] * bfunc1D_2(from.x, to.x, ksi); break;	/// � X2
		case 3: res += hy / 2. * gauss_weights[g] * bfunc1D_2(from.y, to.y, eta); break;	// Y2
		}
	}
	return res;
}


// ���� ����� ��������� ������� �� ����� ���������� ������ ���������
// �������������� ��� ����������� ������������� K ��� ������� ���� ����� ����������� ��������������� � ���� �������
void Rectangle::Calculate_LocalStiffnessMatrix(Element& element) {
	block2x2 block;
	double Kx,	// ��� ��� ����� ������� �������� �� �������� ������� �� ������ ���� �����
		Ky,
		Kxy,
		Kyx;
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			Kx = Element_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i, j, 0);
			Ky = Element_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i, j, 1);
			Kxy = Element_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i, j, 2);
			Kyx = Element_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i, j, 3);
			// ����� ����������� ���:
			block.val11 = mat.thickness * (D[0][0] * Kx + D[2][2] * Ky);
			block.val12 = mat.thickness * (D[0][1] * Kxy + D[2][2] * Kyx);
			block.val21 = mat.thickness * (D[1][0] * Kyx + D[2][2] * Kxy);
			block.val22 = mat.thickness * (D[2][2] * Kx + D[1][1] * Ky);
			// � �������� � ��������� ��, ��� ��� ��� ����� ������� ��������� (�� �������� - �����)
			LocalStiffnessMatrix_block[i][j] = block;

		}
	}


	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			// ���� �� ������ ��������� � �������� ������ �������, �� ��������� ������ ���������� ������ ����� ������
			if ((i + 1) % 2 != 0) {
				if ((j + 1) % 2 != 0)
					LocalStiffnessMatrix[i][j] = LocalStiffnessMatrix_block[i / 2][j / 2].val11;
				else 
					LocalStiffnessMatrix[i][j] = LocalStiffnessMatrix_block[i / 2][j / 2].val12;
			}
			// ���� � ������, �� ��������� ������ ���������� ������ ����� ������
			else {
				if ((j + 1) % 2 != 0)
					LocalStiffnessMatrix[i][j] = LocalStiffnessMatrix_block[i / 2][j / 2].val21;
				else
					LocalStiffnessMatrix[i][j] = LocalStiffnessMatrix_block[i / 2][j / 2].val22;
			}
		}
	}

}


void Rectangle::Assemble_GlobalStiffnessMatrix(Mesh& mesh) {

	cout << "Assembling Global Stiffness Matrix...\n";

	int nodes_count = mesh.nodes.size();
	block2x2 block;
	// ���������� ������� ��������� ����� ����������� NxN, ��� N = ����������_��������_�������_����(2) * ����������_�����
	// ���� � ������� �������, �� K*K, ��� K = ����������_�����, ��� ��� � ������ ����� ����� ��� 4 �������� (2 �������, 2 ������)
	// ��� ��� ��������� � ���������� ������� ��������� ������� �� ������, ������ �������� ��� ������
	GlobalStiffnessMatrix_block.resize(nodes_count);
	for (int k = 0; k < nodes_count; k++) {
		GlobalStiffnessMatrix_block[k].resize(nodes_count);
	}

	// ��������
	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
		Element current_element = mesh.elements[i_elem];
		Calculate_LocalStiffnessMatrix(current_element);
		for (int i = 0; i < current_element.loc_nodes.size(); i++) {
			for (int j = 0; j < current_element.loc_nodes.size(); j++) {
				GlobalStiffnessMatrix_block[current_element.loc_nodes[i].num - 1][current_element.loc_nodes[j].num - 1] += LocalStiffnessMatrix_block[i][j];
			}
		}
	}


	// �������� ������ ��� �������
	GlobalStiffnessMatrix.resize(nodes_count * 2);
	for (int i = 0; i < nodes_count * 2; i++) {
		GlobalStiffnessMatrix[i].resize(nodes_count * 2, 0);
	}


	// ������������������ ���������� ������� � ������������ ���
	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
		Element current_element = mesh.elements[i_elem];
		Calculate_LocalStiffnessMatrix(current_element);
		for (int i = 0; i < current_element.loc_nodes.size(); i++) {
			for (int j = 0; j < current_element.loc_nodes.size(); j++) {
				int i_1 = current_element.loc_nodes[i].num - 1,
					j_1 = current_element.loc_nodes[j].num - 1;
				GlobalStiffnessMatrix[2 * i_1][2 * j_1] += LocalStiffnessMatrix_block[i][j].val11;
				GlobalStiffnessMatrix[2 * i_1][2 * j_1 + 1] += LocalStiffnessMatrix_block[i][j].val12;
				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1] += LocalStiffnessMatrix_block[i][j].val21;
				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1 + 1] += LocalStiffnessMatrix_block[i][j].val22;
			}
		}
	}


	
	//	(�������) ����������� ����� (�������� ������-������� ������������� ���� � ������� ���������, �� ��������� ������ �������), ����� ����� ����� ������� ��� ����� �������������� ���������� ��������


	// =========== �������� ������ ������������ ����� ==============
	// ������ � ������������� ������ ����� ������������, ����������� �� ���������� ��������� ����� � ������������ ������������ ������ (������ ����� �� ���������� � ��������� �������)

	comp_domain domain = mesh.subdomain;
	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i].x == mesh.subdomain.coords[2].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� �����������:
			fixed_nodes.push_back(mesh.nodes[i].num);					// coords[0].x - ���������� ����� ������ ��������
																		// coords[1].x - ���������� ��������� ��������
																		// coords[2].x - ���������� ��� ��������� �������� � ������ ��������������� ������
	}

	// =========== �������� ������ ����������� ����� ==============

	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		if (mesh.nodes[i].x == mesh.subdomain.coords[0].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� ��������:
			loaded_nodes.push_back(mesh.nodes[i].num);
	}



	for(int k = 0; k < fixed_nodes.size();k++) {		// !!!!! k �� 0 �� ������_�������_������������_�����
		for (int i = 0; i < GlobalStiffnessMatrix_block.size(); i++) {
			if (i == (fixed_nodes[k] - 1)) {	// ���� ����� ������ ��������� � ������� ������������� ����
				for (int j = 0; j < GlobalStiffnessMatrix_block.size(); j++) {
					GlobalStiffnessMatrix_block[i][j] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix_block[j][i] = 0;		// �������� ��������������� �������
					GlobalStiffnessMatrix_block[i][i].val11 = 1;		// ������ �� ��������� 1
					GlobalStiffnessMatrix_block[i][i].val22 = 1;
				}
			}
		}
	}

	// ���������� ������� ����������� � ����������� ��������� �������
	for (int k = 0; k < fixed_nodes.size(); k++) {
		for (int i = 0; i < GlobalStiffnessMatrix.size()/2; i++) {
			if (i == (fixed_nodes[k] - 1)) {	// ���� ����� ������ ��������� � ������� ������������� ����
				for (int j = 0; j < GlobalStiffnessMatrix.size()/2; j++) {
					GlobalStiffnessMatrix[2 * i][2 * j] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i][2 * j +1] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i + 1][2 * j] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i + 1][2 * j + 1] = 0;		// �������� ��������������� ������


					GlobalStiffnessMatrix[2 * j][2 * i] = 0;		// �������� ��������������� �������
					GlobalStiffnessMatrix[2 * j][2 * i +1] = 0;		// �������� ��������������� �������
					GlobalStiffnessMatrix[2 * j + 1][2 * i] = 0;		// �������� ��������������� �������
					GlobalStiffnessMatrix[2 * j + 1][2 * i + 1] = 0;		// �������� ��������������� �������

					GlobalStiffnessMatrix[2 * i][2 * i] = 1;		// ������ �� ��������� 1
					//GlobalStiffnessMatrix[2 * i][2 * i +1] = 1;		// ������ �� ��������� 1
					//GlobalStiffnessMatrix[2 * i + 1][2 * i] = 1;		// ������ �� ��������� 1
					GlobalStiffnessMatrix[2 * i + 1][2 * i + 1] = 1;		// ������ �� ��������� 1
				}
			}
		}
	}






	cout << "Assembling Global Stiffness Matrix complete\n";
}

void Rectangle::Calculate_LocalLoadVector(Element& element, Load P) {
	block1x2 block;

	for (int i = 0; i < 4; i++) {
		block.val1 = P.Px * Edge_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i);
		block.val2 = P.Py * Edge_IntegrateGauss3(element.loc_nodes[0], element.loc_nodes[2], i);

		LocalLoadVector_block[i] = block;
	}
		// ������������������ ���������� ������� �������� � ������������ ������
	for (int i = 0; i < LocalLoadVector.size(); i++) {
		if ((i + 1) % 2 != 0)
			LocalLoadVector[i] = LocalLoadVector_block[i / 2].val1;
		else
			LocalLoadVector[i] = LocalLoadVector_block[i / 2].val2;
	}

}

void Rectangle::Assemble_GlobalLoadVector(Mesh& mesh) {
	cout << "\nAssembling Global Load Vector...\n";
	int nodes_count = mesh.nodes.size();
	Load P;
	P.GetLineLength(mesh);
	block1x2 block;
	// ������ ������� �������� - 1�N, ��� N = ����������_��������_�������_����(2) * ����������_�����
	// � ������� ������� - 1xK, ��� K = ����������_�����
	GlobalLoadVector_block.resize(nodes_count);

	// ��������
	// ������ �������� ������ ��������� ������ ��� ��� ���������, ����� ������� ���������� ��������
	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
		Element current_element = mesh.elements[i_elem];
		for (int i_node = 0; i_node < loaded_nodes.size(); i_node++) {
			// ���� ����� ����� �������� ��������� (�. �. ��� ������ ��� ��������� ���� ����� �� ����������� ������), �� �� ����� ����� � ������ �������� � ������� ��� ���� ��������� ������
			if (current_element.loc_nodes[0].num == loaded_nodes[i_node] || current_element.loc_nodes[3].num == loaded_nodes[i_node]) {
				Calculate_LocalLoadVector(current_element, P);
				if (current_element.loc_nodes[0].num == loaded_nodes[i_node])
					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += LocalLoadVector_block[0];
				if (current_element.loc_nodes[3].num == loaded_nodes[i_node])
					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += LocalLoadVector_block[3];
			}
		}
	}

	GlobalLoadVector.resize(nodes_count * 2,0);

	// ������ ����������� ������� �������� � ����������� �������
	for (int i = 0; i < GlobalLoadVector_block.size(); i++) {
		GlobalLoadVector[2 * i] = GlobalLoadVector_block[i].val1;
		GlobalLoadVector[2 * i + 1] = GlobalLoadVector_block[i].val2;
	}



	// ��������� ������� ����������� � ����� � ��������� ������� ��������
	for (int k = 0; k <fixed_nodes.size(); k++) {		
		for (int i = 0; i < GlobalStiffnessMatrix_block.size(); i++) {
			if (i == (fixed_nodes[k] - 1)) {
				GlobalLoadVector_block[i] = 0;		// �������� ��������������� ������
			}
		}
	}

	// ���������� ������� ����������� � ����������� ������� �������� � ������������ �������
	for (int k = 0; k < fixed_nodes.size(); k++) {
		for (int i = 0; i < GlobalLoadVector.size()/2; i++) {
			if (i == (fixed_nodes[k] - 1)) {
				GlobalLoadVector[2 * i] = 0;		// �������� ��������������� ������
				GlobalLoadVector[2 * i + 1] = 0;
			}
		}
	}

	cout << "Assembling Global load Vector complete\n";
}


void Rectangle::GeneratePortrait(Portrait& portrait,vector<vector<double>> GSM, int &ja_sz) {
	int ggl_size = 0;
	int temp = 0;
	portrait.ig[0] = 0;
	for (int i = 0; i < GSM.size(); i++) {
		int nonzero_in_row = 0;
		portrait.di[i] = GSM[i][i];
		for (int j = 0; j < i; j++) {
			if (GSM[i][j] != 0) {
				portrait.PushBack(portrait.ggl, ggl_size, GSM[i][j]);		// ������� ��������� ������� � ������ ggl, �������� ���
				temp = ggl_size - 1;						// ��� ��� ������ ggl ���������� �� 1, ��������� ���
				portrait.PushBack(portrait.jg, temp, j);					// � ������� �� ����� ������� ����� �������
				nonzero_in_row++;							// ����������� ����� ��������� ��������� � ������
			}

		}
		// ig[i+1] = ig[i] + ���-�� ��������� ��������� (������, ���. 496)
		
			portrait.ig[i+1] = portrait.ig[i] + nonzero_in_row;

	}


	ja_sz = portrait.ig[GSM.size()];

	//ofstream outfile;
	//outfile.open("CSRmatrix\\jg.txt");
	//for (int i = 0; i < ggl_size; i++) {
	//	outfile << portrait.jg[i] << endl;
	//}
	//outfile.close();

	//outfile.open("CSRmatrix\\ig.txt");
	//for (int i = 0; i < GSM.size() + 1; i++)
	//	outfile << portrait.ig[i] << endl;
	//outfile.close();

	//outfile.open("CSRmatrix\\ggl.txt");
	//for (int i = 0; i < ggl_size; i++)
	//	outfile << portrait.ggl[i] << endl;
	//outfile.close();
}
