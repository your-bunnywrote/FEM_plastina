// fem.cpp - ����������� ���� ����� ����������� �������

#include "fem.h"




Material::Material() {
	E = 2.1e5;		// [���] ��� ��� ��� ��������� �������� � ����������, ������ ���� �������� � ���
	mu = 0.3;
	c = E / (1 - mu * mu);
	thickness = 3; // [��]
}

Load::Load() {
	Px = -1000;							// ��������, ������� ������������ � �����, [�/��]
	Py = 0.0;
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

void Element::twoD_to_oneD(size_t i, size_t& mu, size_t& nu) {
	mu = ((i + 1) / 2) % 2;		// mu(i): mu(0) = 0, mu(1) = 1, mu(2) = 1, mu(3) = 0
	nu = i / 2;					// nu(i): nu(0) = 0, nu(1) = 0, nu(2) = 1, nu(3) = 1
}



Element::Element() {
	loc_nodes.resize(4);

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

	// ����� ������ �� ��������� ������-�������� [-1,1]x[-1,1]
	gauss_points_local[0] = -0.77459666924148337703585307995648;
	gauss_points_local[1] = 0.;
	gauss_points_local[2] = 0.77459666924148337703585307995648;

	gauss_weights[0] = 5. / 9.;
	gauss_weights[1] = 8. / 9.;
	gauss_weights[2] = 5. / 9.;
}


Rectangle::Rectangle() {
	loc_nodes = Element::loc_nodes;
	int p_id = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			integrate_points[p_id] = Point(gauss_points_local[j], gauss_points_local[i]);
			p_id++;
		}
	}

}



void Element::init() {

}




double Rectangle::bfunc1D(size_t func_num, double x0, double x1, double x) {
	switch (func_num) {
	case 0:
		return (x1 - x) / (x1 - x0);	// = 1-ksi - � ��������� (������������) �����������, ��� ksi = (x-x0)/(x1-x0)
	case 1:
		return (x - x0) / (x1 - x0);	// = ksi
	}
}



double Rectangle::dbfunc1D(size_t func_num, double x0, double x1, double x) {
	switch (func_num) {
	case 0:
		return -1.0 / (x1 - x0);
	case 1:
		return 1.0 / (x1 - x0);

	}
}




double Rectangle::phi(size_t func_num, Point& from, Point& to, double x, double y) {

	size_t mu = 0, nu = 0;
	twoD_to_oneD(func_num, mu, nu);
	return bfunc1D(mu, from.x, to.x, x) * bfunc1D(nu, from.y, to.y, y);

	//switch (func_num) {
	//case 0:
	//	return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y);	// phi1(x,y) = X1(x)*Y1(y)
	//case 1:
	//	return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y);	// phi2(x,y) = X2(x)*Y1(y)
	//case 2:
	//	return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y);	// phi3(x,y) = X2(x)*Y2(y)
	//case 3:
	//	return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y);	// phi4(x,y) = X1(x)*Y2(y)
	//}
}

// dphi
double Rectangle::dphi(size_t var, size_t i, size_t j, Point from, Point to, double x, double y) {
	size_t mu1 = 0, nu1 = 0;
	size_t mu2 = 0, nu2 = 0;
	twoD_to_oneD(i, mu1, nu1);
	twoD_to_oneD(j, mu2, nu2);
	switch (var) {
	// Kx[i][j] = dphi[i]/dx * dphi[j]/dx
	case 0: return dbfunc1D(mu1, from.x, to.x, x) * dbfunc1D(mu2, from.x, to.x, x) * bfunc1D(nu1, from.y, to.y, y) * bfunc1D(nu2, from.y, to.y, y);

	// Ky[i][j] = dphi[i]/dy * dphi[j]/dy
	case 1: return bfunc1D(mu1, from.x, to.x, x) * bfunc1D(mu2, from.x, to.x, x) * dbfunc1D(nu1, from.y, to.y, y) * dbfunc1D(nu2, from.y, to.y, y);

	// Kxy[i][j] = dphi[i]dx * dphi[j]/dy
	case 2: return dbfunc1D(mu1, from.x, to.x, x) * bfunc1D(mu2, from.x, to.x, x) * bfunc1D(nu1, from.y, to.y, y) * dbfunc1D(nu2, from.y, to.y, y);

	// Kyx[i][j] = dphi[i]dy * dphi[j]/dy
	case 3: return bfunc1D(mu1, from.x, to.x, x) * dbfunc1D(mu2, from.x, to.x, x) * dbfunc1D(nu1, from.y, to.y, y) * bfunc1D(nu2, from.y, to.y, y);
	}

}



Quadrilateral::Quadrilateral() {
	loc_nodes = Element::loc_nodes;
	alpha0 = (loc_nodes[1].x - loc_nodes[0].x) * (loc_nodes[3].y - loc_nodes[1].y) - (loc_nodes[1].y - loc_nodes[0].y) * (loc_nodes[3].x - loc_nodes[0].x);
	alpha1 = (loc_nodes[1].x - loc_nodes[0].x) * (loc_nodes[3].y - loc_nodes[2].y) - (loc_nodes[1].y - loc_nodes[0].y) * (loc_nodes[2].x - loc_nodes[3].x);
	alpha2 = (loc_nodes[3].y - loc_nodes[0].y) * (loc_nodes[2].x - loc_nodes[1].x) - (loc_nodes[3].x - loc_nodes[0].x) * (loc_nodes[2].y - loc_nodes[1].y);
	beta1 = loc_nodes[3].x - loc_nodes[0].x;
	beta2 = loc_nodes[1].x - loc_nodes[0].x;
	beta3 = loc_nodes[3].y - loc_nodes[0].y;
	beta4 = loc_nodes[1].y - loc_nodes[0].y;
	beta5 = loc_nodes[0].x - loc_nodes[1].x - loc_nodes[3].x + loc_nodes[2].x;
	beta6 = loc_nodes[0].y - loc_nodes[1].y - loc_nodes[3].y + loc_nodes[2].y;

	// ������� � �������� ������-�������� [-1,1]x[-1,1] �� ��������� [0,1]x[0,1]
	int p_id = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			integrate_points[i] = Point((gauss_points_local[i] + 1.0) / 2, (gauss_points_local[j] + 1.0) / 2);
		}
	}

}

double Quadrilateral::bfunc1D(size_t func_num, double ksi) {
	switch (func_num) {
	case 0:
		return 1 - ksi;
	case 1:
		return ksi;
	}
}

double Quadrilateral::dbfunc1D(size_t func_num, double ksi) {
	switch (func_num) {
	case 0:
		return -1;
	case 1:
		return 1;
	}
}

Point Quadrilateral::dphi(size_t num, double ksi, double eta) {
	size_t mu = 0, nu = 0;
	twoD_to_oneD(num, mu, nu);

	double W_mu = bfunc1D(mu, ksi);
	double dW_mu = dbfunc1D(mu, ksi);
	double W_nu = bfunc1D(nu, eta);
	double dW_nu = dbfunc1D(nu, eta);

	return Point(dW_mu * W_nu, W_mu * dW_nu);
}


Point Quadrilateral::to_global(Point& p) {
	double x = (1 - p.x) * (1 - p.y) * loc_nodes[0].x + p.x * (1 - p.y) * loc_nodes[1].x + p.x * p.y * loc_nodes[2].x + (1 - p.x) * p.y * loc_nodes[3].x;
	double y = (1 - p.x) * (1 - p.y) * loc_nodes[0].y + p.x * (1 - p.y) * loc_nodes[1].y + p.x * p.y * loc_nodes[2].y + (1 - p.x) * p.y * loc_nodes[3].y;
	return Point(x, y);
}

double Quadrilateral::det_J(double ksi, double eta) {
	double J = alpha0 + alpha1 * ksi + alpha2 * eta;
	return J;
}


double Quadrilateral::grad_bfunc_2d(size_t num,size_t var, double ksi, double eta) {
	//size_t mu = 0, nu = 0;
	//twoD_to_oneD(num, mu, nu);

	double dphi_dksi = dphi(num, ksi, eta).x;
	double dphi_deta = dphi(num, ksi, eta).y;

	switch (var) {
	// ��� Kx
	case 0:
		return dphi_dksi * (beta6 * ksi + beta3) - dphi_deta * (beta6 * eta + beta4);
	// ��� Ky
	case 1:
		return -dphi_dksi * (beta5 * ksi + beta1) + dphi_deta * (beta5 * eta + beta2);

	}

}


 // �������� ��� ������������ Kvar[num1][num2], ��� var - ����������, �� ������� �������������� ������� phi(x,y)
 // var = 0 - x, var = 1 - y, var = 2 - xy, var = 3 - yx;

//double Rectangle::dphi(size_t var, size_t num1, size_t num2, Point& from, Point& to, double x, double y) {
//	switch (var) {
//	// Kx[num1][num2] = dphi[num1]/dx * dphi[num2]/dx
//	case 0:
//		switch (num1) {
//		case 0:
//			switch (num2) {
//			case 0: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi1/dx = dX1/dx * dX1/dx * Y1 * Y1
//			case 1: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi2/dx = dX1/dx * dX2/dx * Y1 * Y1
//			case 2: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi1/dx * dphi3/dx = dX1/dx * dX2/dx * Y1 * Y2
//			case 3: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dx = dX1/dx * dX1/dx * Y1 * Y2
//			}
//		case 1:
//			switch (num2) {
//			case 0: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi1/dx = dX2/dx * dX1/dx * Y1 * Y1
//			case 1: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi2/dx = dX2/dx * dX2/dx * Y1 * Y1
//			case 2: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi3/dx = dX2/dx * dX2/dx * Y1 * Y2
//			case 3: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dx = dX2/dx * dX2/dx * Y2 * Y2
//			}
//		case 2:
//			switch (num2) {
//			case 0: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi1/dx = dX2/dx * dX1/dx * Y2 * Y1
//			case 1: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi2/dx = dx2/dx * dX2/dx * Y2 * Y1
//			case 2: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi3/dx = dX2/dx * dX2/dx * Y2 * Y2
//			case 3: return dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi4/dx = dX2/dx * dX1/dx * Y2 * Y2
//			}
//		case 3:
//			switch (num2) {
//			case 0: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi1/dx = dX1/dx * dX1/dx * Y2 * Y1
//			case 1: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi2/dx = dX1/dx * dX2/dx * Y2 * Y1
//			case 2: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi3/dx = dX1/dx * dX2/dx * Y2 * Y2
//			case 3: return dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi4/dx * dX1/dx * dX1/dx * Y2 * Y2
//			}
//		}
//
//	// Ky[num1][num2] = dphi[num1]/dy * dphi[num2]dy
//	case 1:
//		switch (num1) {
//		case 0:
//			switch (num2) {
//			case 0: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi1/dy * dphi1/dy = X1 * X1 * dY1 * dY1
//			case 1: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi1/dy * dphi2/dy = X1 * X2 * dY1 * dY1
//			case 2: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi3/dy = X1 * X2 * dY1 * dY2
//			case 3: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi4/dy = X1 * X1 * dY1 * dY2
//			}
//		case 1:
//			switch (num2) {
//			case 0: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi1/dy = X2 * X1 * dY1 * dY1
//			case 1: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi2/dy = X2 * X2 * dY1 * dY1
//			case 2: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi3/dy = X2 * X2 * dY1 * dY2
//			case 3: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi4/dy = X2 * X1 * dY1 * dY2
//			}
//		case 2:
//			switch (num2) {
//			case 0: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi3/dy * dphi1/dy = X2 * X1 * dY2 * dY1
//			case 1: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi3/dy * dphi2/dy = X2 * X2 * dY2 * dY1
//			case 2: return bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi3/dy = X2 * X2 * dY2 * dY2
//			case 3: return bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi4/dy = X2 * X1 * dY2 * dY2
//			}
//		case 3:
//			switch (num2) {
//			case 0: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi4/dy * dphi1/dy = X1 * X1 * dY2 * dY1
//			case 1: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi2/dy = X1 * X2 * dY2 * dY1
//			case 2: return bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi3/dy = X1 * X2 * dY2 * dY2
//			case 3: return bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi4/dy = X1 * X1 * dY2 * dY2
//			}
//		}
//	// Kxy[num1][num2] = dphi[num1]/dx * dphi[num2]/dy
//	case 2:
//		switch (num1) {
//		case 0:
//			switch (num2) {
//			case 0: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi1/dx * dphi1/dy = dX1 * X1 * Y1 * dY1
//			case 1: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi1/dx * dphi2/dx = dX1 * X2 * Y1 * dY1
//			case 2: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi1/dx * dphi3/dy = dX1 * X2 * Y1 * dY2
//			case 3: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y); // dphi1/dx * dphi4/dy = dX1 * X1 * Y2 * dY2
//			}
//		case 1:
//			switch (num2) {
//			case 0: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi1/dy = dX2 * X1 * Y1 * dY1
//			case 1: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi2/dx * dphi2/dy = dX2 * X2 * Y1 * dY1
//			case 2: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi3/dy = dX2 * X2 * Y1 * dY2
//			case 3: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi2/dx * dphi4/dy = dX1 * X1 * Y1 * dY2
//			}
//		case 2:
//			switch (num2) {
//			case 0: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y); // dphi3/dx * dphi1/dy = dX2 * X1 * Y2 * dY1
//			case 1: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi3/dx * dphi2/dy = dX2 * X2 * Y1 * dY2
//			case 2: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi3/dx * dphi3/dy = dX2 * X2 * Y2 * dY2
//			case 3: return dbfunc1D_2(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y); // dphi3/dx * dphi4/dy = dX2 * X1 * Y1 * dY2
//			}
//		case 3:
//			switch (num2) {
//			case 0: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi1/dy = dX1 * X1 * Y2 * dY1
//			case 1: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_1(from.y, to.y, y);	// dphi4/dx * dphi2/dy = dX1 * X2 * Y2 * dY1
//			case 2: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi3/dy = dX1 * X2 * Y2 * dY2
//			case 3: return dbfunc1D_1(from.x, to.x, x) * bfunc1D_1(from.x, to.x, x) * bfunc1D_2(from.y, to.y, y) * dbfunc1D_2(from.y, to.y, y);	// dphi4/dx * dphi4/dy = dX1 * X1 * Y2 * dY2
//			}
//		}
//	// Kyx[num1][num2] = dphi[num1]/dy * dphi[num2]/dx
//	case 3:
//		switch (num1) {
//		case 0:
//			switch (num2) {
//			case 0: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y); // dphi1/dy * dphi1/dx = X1 * dX1 * dY1 * Y1
//			case 1: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi1/dy * dphi2/dx = X1 * dX2 * dY1 * Y1
//			case 2: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi1/dy * dphi3/dx = X1 * dX2 * dY1 * Y2
//			case 3: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y); // dphi1/dy * dphi4/dx = X1 * dX1 * dY2 * Y2
//			}
//		case 1:
//			switch (num2) {
//			case 0: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi1/dx = X2 * dX1 * dY1 * Y1
//			case 1: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi2/dy * dphi2/dx = X2 * dX2 * dY1 * Y1
//			case 2: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi3/dx = X2 * dX2 * dY1 * Y2
//			case 3: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi2/dy * dphi4/dx = X1 * dX1 * dY1 * Y2
//			}
//		case 2:
//			switch (num2) {
//			case 0: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y); // dphi3/dy * dphi1/dx = X2 * dX1 * dY2 * Y1
//			case 1: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi3/dy * dphi2/dx = X2 * dX2 * dY1 * Y2
//			case 2: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi3/dy * dphi3/dx = X2 * dX2 * dY2 * Y2
//			case 3: return bfunc1D_2(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y); // dphi3/dy * dphi4/dx = X2 * dX1 * dY1 * Y2
//			}
//		case 3:
//			switch (num2) {
//			case 0: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi1/dx = X1 * dX1 * dY2 * Y1
//			case 1: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_1(from.y, to.y, y);	// dphi4/dy * dphi2/dx = X1 * dX2 * dY2 * Y1
//			case 2: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi3/dx = X1 * dX2 * dY2 * Y2
//			case 3: return bfunc1D_1(from.x, to.x, x) * dbfunc1D_1(from.x, to.x, x) * dbfunc1D_2(from.y, to.y, y) * bfunc1D_2(from.y, to.y, y);	// dphi4/dy * dphi4/dx = X1 * dX1 * dY2 * Y2
//			}
//		}
//	}
//}

double Rectangle::Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) {
	//this->init();
	int p_id = 0;
	double res = 0.0;
	double x_centre = (to.x + from.x) / 2., 
		   y_centre = (to.y + from.y) / 2.;
	double hx = to.x - from.x,
		   hy = to.y - from.y;
	double x, y;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			// ��������� � ���������� ���������� ��������
			x = x_centre + integrate_points[p_id].x * hx / 2.;
			y = y_centre + integrate_points[p_id].y * hy / 2.;
			// � ����������� � ���������� �����������
			res += gauss_weights[i] * gauss_weights[j] * dphi(var, num1, num2, from, to, x, y);
			p_id++;
		}
	}
	return res * hx * hy / 4.0;	// hx*hy/4 - ������� �������� �� ��������� ��������� �������������� � ���������� � ��������� ������
}



double Quadrilateral::Element_IntegrateGauss3(size_t i, size_t j, size_t dvar1, size_t dvar2) {
	
	int p_id = 0;
	double res = 0.0;
	double func = 0.0;
	double ksi = 0.0, eta = 0.0;
	double J = 0.0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			ksi = integrate_points[p_id].x;
			eta = integrate_points[p_id].y;
			func = grad_bfunc_2d(i, dvar1, ksi, eta) * grad_bfunc_2d(j, dvar2, ksi, eta);
			J = det_J(ksi, eta);
			func /= J;
			res += gauss_weights[i] * gauss_weights[j] * func;


		}
	}
	return res / 4.0;
}



double Rectangle::Edge_IntegrateGauss3(Point& from, Point& to, size_t num) {
	double res = 0.0;
	//double hx = to.x - from.x;
	double hy = to.y - from.y;
	//double x_centre = (to.x+from.x) / 2;
	double y_centre = (to.y+from.y) / 2.;
	double ksi,eta;
	size_t mu = 0, nu = 0;
	twoD_to_oneD(num, mu, nu);
	for (int g = 0; g < 3; g++) {
		//ksi = x_centre + gauss_points_local[g] * hx/2;
		eta = y_centre + gauss_points_local[g] * hy / 2.;
		res += hy / 2. * gauss_weights[g] * bfunc1D(mu,from.y, to.y, eta);
		//switch (num) {
		//case 0: res += hy / 2. * gauss_weights[g] * bfunc1D_1(from.y, to.y, eta); break;	// ���� �������� ��������� � ����� ��� ������ ������ (����������� �) �� ����������� ���������� �.�. Y1
		//case 1: res += hx/2. * gauss_weights[g] * bfunc1D_1(from.x, to.x, ksi); break;	/// ���� �������� ��������� � ������� ��� ������ ������ (����������� �), ����������� ���������� �.�. X1
		//case 2: res += hx/2. * gauss_weights[g] * bfunc1D_2(from.x, to.x, ksi); break;	/// � X2
		//case 3: res += hy / 2. * gauss_weights[g] * bfunc1D_2(from.y, to.y, eta); break;	// Y2
		//}
	}
	return res;
}

double Quadrilateral::Edge_IntegrateGauss3(Point& from, Point& to, size_t num) {
	return 1;
}


// ���� ����� ��������� ������� �� ����� ���������� ������ ���������
// �������������� ��� ����������� ������������� K ��� ������� ���� ����� ����������� ��������������� � ���� �������
void Quadrilateral::CalculateLocalStiffnessMatrix() {
	cout << "Calculate local stiffness matrix for Quad\n";
	block2x2 block;
	double Kx,	// ��� ��� ����� ������� �������� �� �������� ������� �� ������ ���� �����
		Ky,
		Kxy,
		Kyx;
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			Kx = Element_IntegrateGauss3(i, j, 0, 0);
			Ky = Element_IntegrateGauss3(i, j, 1, 1);
			Kxy = Element_IntegrateGauss3(i, j, 0, 1);
			Kyx = Element_IntegrateGauss3(i, j, 1, 0);
			// ����� ����������� ���:
			block.val11 = mat.thickness * (D[0][0] * Kx + D[2][2] * Ky);
			block.val12 = mat.thickness * (D[0][1] * Kxy + D[2][2] * Kyx);
			block.val21 = mat.thickness * (D[1][0] * Kyx + D[2][2] * Kxy);
			block.val22 = mat.thickness * (D[2][2] * Kx + D[1][1] * Ky);
			// � �������� � ��������� ��, ��� ��� ��� ����� ������� ��������� (�� �������� - �����)
			LocalStiffnessMatrix_block[i][j] = block;

		}
	}

	// ������������������ ������� ��������� ������� � ������������ ������
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

void Rectangle::CalculateLocalStiffnessMatrix() {
	cout << "Calculate local stiffness matrix for Rectangle\n";
	block2x2 block;
	double Kx,	// ��� ��� ����� ������� �������� �� �������� ������� �� ������ ���� �����
		Ky,
		Kxy,
		Kyx;
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			Kx = Element_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i, j, 0);
			Ky = Element_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i, j, 1);
			Kxy = Element_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i, j, 2);
			Kyx = Element_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i, j, 3);
			// ����� ����������� ���:
			block.val11 = mat.thickness * (D[0][0] * Kx + D[2][2] * Ky);
			block.val12 = mat.thickness * (D[0][1] * Kxy + D[2][2] * Kyx);
			block.val21 = mat.thickness * (D[1][0] * Kyx + D[2][2] * Kxy);
			block.val22 = mat.thickness * (D[2][2] * Kx + D[1][1] * Ky);
			// � �������� � ��������� ��, ��� ��� ��� ����� ������� ��������� (�� �������� - �����)
			LocalStiffnessMatrix_block[i][j] = block;

		}
	}

	// ������������������ ������� ��������� ������� � ������������ ������
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






//void Rectangle::Assemble_GlobalStiffnessMatrix(Mesh& mesh) {
//
//	cout << "Assembling Global Stiffness Matrix...\n";
//
//	int nodes_count = mesh.nodes.size();
//	block2x2 block;
//	// ���������� ������� ��������� ����� ����������� NxN, ��� N = ����������_��������_�������_����(2) * ����������_�����
//	// ���� � ������� �������, �� K*K, ��� K = ����������_�����, ��� ��� � ������ ����� ����� ��� 4 �������� (2 �������, 2 ������)
//	// ��� ��� ��������� � ���������� ������� ��������� ������� �� ������, ������ �������� ��� ������
//	GlobalStiffnessMatrix_block.resize(nodes_count);
//	for (int k = 0; k < nodes_count; k++) {
//		GlobalStiffnessMatrix_block[k].resize(nodes_count);
//	}
//
//	// ��������
//	//for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
//	//	Element current_element = mesh.elements[i_elem];
//	//	//if ((current_element.loc_nodes[0].y == current_element.loc_nodes[1].y)
//	//	//	&& (current_element.loc_nodes[0].x == current_element.loc_nodes[3].x)
//	//	//	&& (current_element.loc_nodes[1].x == current_element.loc_nodes[2].x)
//	//	//	&& (current_element.loc_nodes[2].y == current_element.loc_nodes[3].y))
//	//	//	current_element.type = RECTANGLE;
//	//	//else
//	//	//	current_element.type = QUADR;
//	//	Calculate_LocalStiffnessMatrix(current_element);
//	//	for (int i = 0; i < 4; i++) {
//	//		for (int j = 0; j < 4; j++) {
//	//			GlobalStiffnessMatrix_block[current_element.loc_nodes[i].num - 1][current_element.loc_nodes[j].num - 1] += LocalStiffnessMatrix_block[i][j];
//	//		}
//	//	}
//	//}
//
//
//	// �������� ������ ��� �������
//	GlobalStiffnessMatrix.resize(nodes_count * 2);
//	for (int i = 0; i < nodes_count * 2; i++) {
//		GlobalStiffnessMatrix[i].resize(nodes_count * 2, 0);
//	}
//
//
//	// ������������������ ���������� ������� � ������������ ���
//	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
//
//
//
//		Element current_element = mesh.elements[i_elem];
//		Calculate_LocalStiffnessMatrix(current_element);
//		for (int i = 0; i < 4; i++) {
//			for (int j = 0; j < 4; j++) {
//				int i_1 = current_element.loc_nodes[i].num - 1,
//					j_1 = current_element.loc_nodes[j].num - 1;
//				GlobalStiffnessMatrix[2 * i_1][2 * j_1] += LocalStiffnessMatrix_block[i][j].val11;
//				GlobalStiffnessMatrix[2 * i_1][2 * j_1 + 1] += LocalStiffnessMatrix_block[i][j].val12;
//				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1] += LocalStiffnessMatrix_block[i][j].val21;
//				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1 + 1] += LocalStiffnessMatrix_block[i][j].val22;
//			}
//		}
//	}
//
//
//	
//
//
//	// =========== �������� ������ ������������ ����� ==============
//	// ������ � ������������� ������ ����� ������������, ����������� �� ���������� ��������� ����� � ������������ ������������ ������ (������ ����� �� ���������� � ��������� �������)
//
//	comp_domain domain = mesh.subdomain;
//	for (size_t i = 0; i < mesh.nodes.size(); i++) {
//		//if (mesh.nodes[i].x == mesh.subdomain.coords[1].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� �����������:
//		if(mesh.nodes[i].x == mesh.subdomain.vertical_curves[2][0].begin.x)		// ����� ���������� ���������� ���� � ������
//			fixed_nodes.push_back(mesh.nodes[i].num);					// coords[0].x - ���������� ����� ������ ��������
//																		// coords[1].x - ���������� ��������� ��������
//																		// coords[2].x - ���������� ��� ��������� �������� � ������ ��������������� ������
//	}
//
//	// =========== �������� ������ ����������� ����� ==============
//
//	for (size_t i = 0; i < mesh.nodes.size(); i++) {
//		//if (mesh.nodes[i].x == mesh.subdomain.coords[0].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� ��������:
//		if(mesh.nodes[i].x == mesh.subdomain.vertical_curves[0][0].begin.x)	// ���������� ���������� ���� � ������
//			loaded_nodes.push_back(mesh.nodes[i].num);
//	}
//
//
//
//	for(int k = 0; k < fixed_nodes.size();k++) {		// !!!!! k �� 0 �� ������_�������_������������_�����
//		for (int i = 0; i < GlobalStiffnessMatrix_block.size(); i++) {
//			if (i == (fixed_nodes[k] - 1)) {	// ���� ����� ������ ��������� � ������� ������������� ����
//				for (int j = 0; j < GlobalStiffnessMatrix_block.size(); j++) {
//					GlobalStiffnessMatrix_block[i][j] = 0;		// �������� ��������������� ������
//					GlobalStiffnessMatrix_block[j][i] = 0;		// �������� ��������������� �������
//					GlobalStiffnessMatrix_block[i][i].val11 = 1;		// ������ �� ��������� 1
//					GlobalStiffnessMatrix_block[i][i].val22 = 1;
//				}
//			}
//		}
//	}
//
//	// ���������� ������� ����������� � ����������� ��������� �������
//	for (int k = 0; k < fixed_nodes.size(); k++) {
//		for (int i = 0; i < GlobalStiffnessMatrix.size()/2; i++) {
//			if (i == (fixed_nodes[k] - 1)) {	// ���� ����� ������ ��������� � ������� ������������� ����
//				for (int j = 0; j < GlobalStiffnessMatrix.size()/2; j++) {
//					GlobalStiffnessMatrix[2 * i][2 * j] = 0;		// �������� ��������������� ������
//					GlobalStiffnessMatrix[2 * i][2 * j +1] = 0;		// �������� ��������������� ������
//					GlobalStiffnessMatrix[2 * i + 1][2 * j] = 0;		// �������� ��������������� ������
//					GlobalStiffnessMatrix[2 * i + 1][2 * j + 1] = 0;		// �������� ��������������� ������
//
//
//					GlobalStiffnessMatrix[2 * j][2 * i] = 0;		// �������� ��������������� �������
//					GlobalStiffnessMatrix[2 * j][2 * i +1] = 0;		// �������� ��������������� �������
//					GlobalStiffnessMatrix[2 * j + 1][2 * i] = 0;		// �������� ��������������� �������
//					GlobalStiffnessMatrix[2 * j + 1][2 * i + 1] = 0;		// �������� ��������������� �������
//
//					GlobalStiffnessMatrix[2 * i][2 * i] = 1;		// ������ �� ��������� 1
//					//GlobalStiffnessMatrix[2 * i][2 * i +1] = 1;		// ������ �� ��������� 1
//					//GlobalStiffnessMatrix[2 * i + 1][2 * i] = 1;		// ������ �� ��������� 1
//					GlobalStiffnessMatrix[2 * i + 1][2 * i + 1] = 1;		// ������ �� ��������� 1
//				}
//			}
//		}
//	}
//
//
//
//
//
//
//	cout << "Assembling Global Stiffness Matrix complete\n";
//}

void FEM::AssembleGlobalStiffnessMatrix(Mesh& mesh) {

	cout << "Assembling Global Stiffness Matrix...\n";

	int nodes_count = mesh.nodes.size();
	//block2x2 block;
	// ���������� ������� ��������� ����� ����������� NxN, ��� N = ����������_��������_�������_����(2) * ����������_�����
	// ���� � ������� �������, �� K*K, ��� K = ����������_�����, ��� ��� � ������ ����� ����� ��� 4 �������� (2 �������, 2 ������)
	// ��� ��� ��������� � ���������� ������� ��������� ������� �� ������, ������ �������� ��� ������


	// �������� ������ ��� �������
	GlobalStiffnessMatrix.resize(nodes_count * 2);
	for (int i = 0; i < nodes_count * 2; i++) {
		GlobalStiffnessMatrix[i].resize(nodes_count * 2, 0);
	}

	Element element;
	element.init();
	Element* el;
	//Rectangle rect;
	//Quadrilateral quad;
	// ������������������ ���������� ������� � ������������ ���
	element.loc_nodes[0].x = 1.;
	Rectangle rect1;
	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
		element.loc_nodes = mesh.elements[i_elem].loc_nodes;
		element.type = mesh.check_element_type(mesh.elements[i_elem]);
		Rectangle rect;
		Quadrilateral quad;
		if (element.type == RECT)
			el = &rect;
		else
			el = &quad;
		el->CalculateLocalStiffnessMatrix();	// ���������, ��� �������� �������

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				int i_1 = el->loc_nodes[i].num - 1,
					j_1 = el->loc_nodes[j].num - 1;
				GlobalStiffnessMatrix[2 * i_1][2 * j_1] += el->LocalStiffnessMatrix_block[i][j].val11;
				GlobalStiffnessMatrix[2 * i_1][2 * j_1 + 1] += el->LocalStiffnessMatrix_block[i][j].val12;
				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1] += el->LocalStiffnessMatrix_block[i][j].val21;
				GlobalStiffnessMatrix[2 * i_1 + 1][2 * j_1 + 1] += el->LocalStiffnessMatrix_block[i][j].val22;
			}
		}
	}





	// =========== �������� ������ ������������ ����� ==============
	// ������ � ������������� ������ ����� ������������, ����������� �� ���������� ��������� ����� � ������������ ������������ ������ (������ ����� �� ���������� � ��������� �������)

	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		//if (mesh.nodes[i].x == mesh.subdomain.coords[1].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� �����������:
		if (mesh.nodes[i].x == mesh.subdomain.vertical_curves[2][0].begin.x)		// ����� ���������� ���������� ���� � ������
			fixed_nodes.push_back(mesh.nodes[i].num);					// coords[0].x - ���������� ����� ������ ��������
		// coords[1].x - ���������� ��������� ��������
		// coords[2].x - ���������� ��� ��������� �������� � ������ ��������������� ������
	}

	// =========== �������� ������ ����������� ����� ==============

	for (size_t i = 0; i < mesh.nodes.size(); i++) {
		//if (mesh.nodes[i].x == mesh.subdomain.coords[0].x)	// ���� ���������� x ���� ��������� � ����������� ������������ ������, �� �������� ��� �� ��������:
		if (mesh.nodes[i].x == mesh.subdomain.vertical_curves[0][0].begin.x)	// ���������� ���������� ���� � ������
			loaded_nodes.push_back(mesh.nodes[i].num);
	}



	// ���������� ������� ����������� � ����������� ��������� �������
	for (int k = 0; k < fixed_nodes.size(); k++) {
		for (int i = 0; i < GlobalStiffnessMatrix.size() / 2; i++) {
			if (i == (fixed_nodes[k] - 1)) {	// ���� ����� ������ ��������� � ������� ������������� ����
				for (int j = 0; j < GlobalStiffnessMatrix.size() / 2; j++) {
					GlobalStiffnessMatrix[2 * i][2 * j] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i][2 * j + 1] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i + 1][2 * j] = 0;		// �������� ��������������� ������
					GlobalStiffnessMatrix[2 * i + 1][2 * j + 1] = 0;		// �������� ��������������� ������


					GlobalStiffnessMatrix[2 * j][2 * i] = 0;		// �������� ��������������� �������
					GlobalStiffnessMatrix[2 * j][2 * i + 1] = 0;		// �������� ��������������� �������
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


void Quadrilateral::CalculateLocalLoadVector(Load& P) {
	block1x2 block;

	for (int i = 0; i < 4; i++) {
		block.val1 = P.Px * Edge_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i);
		block.val2 = P.Py * Edge_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i);

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

void Rectangle::CalculateLocalLoadVector(Load& P) {
	block1x2 block;

	for (int i = 0; i < 4; i++) {
		block.val1 = P.Px * Edge_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i);
		block.val2 = P.Py * Edge_IntegrateGauss3(loc_nodes[0], loc_nodes[2], i);

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


//void Rectangle::Assemble_GlobalLoadVector(Mesh& mesh) {
//	cout << "\nAssembling Global Load Vector...\n";
//	int nodes_count = mesh.nodes.size();
//	Load P;
//	//P.GetLineLength(mesh);
//	block1x2 block;
//	// ������ ������� �������� - 1�N, ��� N = ����������_��������_�������_����(2) * ����������_�����
//	// � ������� ������� - 1xK, ��� K = ����������_�����
//	GlobalLoadVector_block.resize(nodes_count);
//
//	// ��������
//	// ������ �������� ������ ��������� ������ ��� ��� ���������, ����� ������� ���������� ��������
//	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
//		Element current_element = mesh.elements[i_elem];
//		for (int i_node = 0; i_node < loaded_nodes.size(); i_node++) {
//			// ���� ����� ����� �������� ��������� (�. �. ��� ������ ��� ��������� ���� ����� �� ����������� ������), �� �� ����� ����� � ������ �������� � ������� ��� ���� ��������� ������
//			if (current_element.loc_nodes[0].num == loaded_nodes[i_node] || current_element.loc_nodes[3].num == loaded_nodes[i_node]) {
//				Calculate_LocalLoadVector(current_element, P);
//				if (current_element.loc_nodes[0].num == loaded_nodes[i_node])
//					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += LocalLoadVector_block[0];
//				if (current_element.loc_nodes[3].num == loaded_nodes[i_node])
//					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += LocalLoadVector_block[3];
//			}
//		}
//	}
//
//	GlobalLoadVector.resize(nodes_count * 2,0);
//
//	// ������ ����������� ������� �������� � ����������� �������
//	for (int i = 0; i < GlobalLoadVector_block.size(); i++) {
//		GlobalLoadVector[2 * i] = GlobalLoadVector_block[i].val1;
//		GlobalLoadVector[2 * i + 1] = GlobalLoadVector_block[i].val2;
//	}
//
//
//
//	// ��������� ������� ����������� � ����� � ��������� ������� ��������
//	for (int k = 0; k <fixed_nodes.size(); k++) {		
//		for (int i = 0; i < GlobalStiffnessMatrix_block.size(); i++) {
//			if (i == (fixed_nodes[k] - 1)) {
//				GlobalLoadVector_block[i] = 0;		// �������� ��������������� ������
//			}
//		}
//	}
//
//	// ���������� ������� ����������� � ����������� ������� �������� � ������������ �������
//	for (int k = 0; k < fixed_nodes.size(); k++) {
//		for (int i = 0; i < GlobalLoadVector.size()/2; i++) {
//			if (i == (fixed_nodes[k] - 1)) {
//				GlobalLoadVector[2 * i] = 0;		// �������� ��������������� ������
//				GlobalLoadVector[2 * i + 1] = 0;
//			}
//		}
//	}
//
//	cout << "Assembling Global load Vector complete\n";
//}

void FEM::AssembleGlobalLoadVector(Mesh& mesh) {
	cout << "\nAssembling Global Load Vector...\n";
	int nodes_count = mesh.nodes.size();
	Load P;
	//P.GetLineLength(mesh);
	block1x2 block;
	// ������ ������� �������� - 1�N, ��� N = ����������_��������_�������_����(2) * ����������_�����
	// � ������� ������� - 1xK, ��� K = ����������_�����
	GlobalLoadVector_block.resize(nodes_count);
	GlobalLoadVector.resize(2 * nodes_count);

	// ��������
	// ������ �������� ������ ��������� ������ ��� ��� ���������, ����� ������� ���������� ��������
	Element element;
	Element* el;
	Rectangle rect;
	Quadrilateral quad;
	for (int i_elem = 0; i_elem < mesh.elements.size(); i_elem++) {
		Cell current_element = mesh.elements[i_elem];
		element.type = mesh.check_element_type(current_element);
		for (int i_node = 0; i_node < loaded_nodes.size(); i_node++) {
			// ���� ����� ����� �������� ��������� (�. �. ��� ������ ��� ��������� ���� ����� �� ����������� ������), �� �� ����� ����� � ������ �������� � ������� ��� ���� ��������� ������
			if (current_element.loc_nodes[0].num == loaded_nodes[i_node] || current_element.loc_nodes[3].num == loaded_nodes[i_node]) {
				if (element.type == RECT)
					el = &rect;
				else
					el = &quad;
				el->CalculateLocalLoadVector(P);


				if (current_element.loc_nodes[0].num == loaded_nodes[i_node])
					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += el->LocalLoadVector_block[0];
				if (current_element.loc_nodes[3].num == loaded_nodes[i_node])
					GlobalLoadVector_block[loaded_nodes[i_node] - 1] += el->LocalLoadVector_block[3];
			}
		}
	}

	GlobalLoadVector.resize(nodes_count * 2, 0);

	// ������ ����������� ������� �������� � ����������� �������
	for (int i = 0; i < GlobalLoadVector_block.size(); i++) {
		GlobalLoadVector[2 * i] = GlobalLoadVector_block[i].val1;
		GlobalLoadVector[2 * i + 1] = GlobalLoadVector_block[i].val2;
	}



	// ��������� ������� ����������� � ����� � ��������� ������� ��������
	for (int k = 0; k < fixed_nodes.size(); k++) {
		for (int i = 0; i < GlobalStiffnessMatrix_block.size(); i++) {
			if (i == (fixed_nodes[k] - 1)) {
				GlobalLoadVector_block[i] = 0;		// �������� ��������������� ������
			}
		}
	}

	// ���������� ������� ����������� � ����������� ������� �������� � ������������ �������
	for (int k = 0; k < fixed_nodes.size(); k++) {
		for (int i = 0; i < GlobalLoadVector.size() / 2; i++) {
			if (i == (fixed_nodes[k] - 1)) {
				GlobalLoadVector[2 * i] = 0;		// �������� ��������������� ������
				GlobalLoadVector[2 * i + 1] = 0;
			}
		}
	}

	cout << "Assembling Global load Vector complete\n";

}



void FEM::GeneratePortrait(Portrait& portrait,vector<vector<double>> GSM, int &ja_sz) {
	cout << "Generating Portrait...\n";
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

	cout << "Portrait generated.\n";
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
