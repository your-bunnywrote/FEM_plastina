// fem.h - здесь объявляется производный от Элемента класс Прямоугольник, имеющий все необходимые
// поля и методы для работы с конечно-элементной моделью
#pragma once
#ifndef FEM_H
#define FEM_H

#include "Mesh.h"
#include "MyMatrix.h"

class block2x2;
class block1x2;
class Load;


// если экземпляры класса нигде не создаются, то он абстрактный
class Element {
public:
	etype type;
	vector<Point> loc_nodes;
	size_t num;
	Material mat;
	// РЕДАКТИРОВАТЬ: можно заменить на обычный двумерный массив
	vector<vector<double>> D = { {mat.c, mat.c * mat.mu, 0},
							   {mat.c * mat.mu, mat.c * 1, 0},
							   {0, 0, mat.c * (1 - mat.mu) / 2} };	// матрица упругости, для изотропного материала имеющая вид		
	//												    E   |1 mu	 0	  |
	//											D  = ------	|mu 1	 0	  |
	//												 1-mu^2 |0  0 (1-mu)/2|
	Element();

	vector<vector<double>> LocalStiffnessMatrix;
	vector<vector<block2x2>> LocalStiffnessMatrix_block;
	virtual void CalculateLocalStiffnessMatrix() {};
	vector<double> LocalLoadVector;
	vector<block1x2> LocalLoadVector_block;
	virtual void CalculateLocalLoadVector(Load& P) {};

	virtual double Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) {
		return 1;
	};
	virtual double Edge_IntegrateGauss3(Point& from, Point& to, size_t num) {
		return 1;
	};

	// вычисление деформаций по элементу
	virtual Point GetElementCentroid() = 0;
	virtual double GetElementalStrainX(double* u, Point& p) = 0;
	virtual double GetElementalStrainY(double* u, Point& p) = 0;
	virtual double GetElementalStrainXY(double* u, Point& p) = 0;
	

	// узловые деформации
	//virtual double GetNodalStrainX(double* u, Point& p) = 0;
	//virtual double GetNodalStrainY(double* u, Point& p) = 0;
	//virtual double GetNodalStrainXY(double* u, Point& p) = 0;

	virtual double GetElementalStressX(double x_strain, double y_strain) = 0;
	virtual double GetElementalStressY(double x_strain, double y_strain) = 0;
	virtual double GetElementalStressXY(double xy_strain) = 0;

	// определение содержания точки внутри элемента
	bool is_inside(Point& p);

	//Element* GetElement(Point& p);

protected:
	// получение номера одномерной бф через номер двухмерной
	void twoD_to_oneD(size_t i, size_t& mu, size_t& nu);

	// базисные функции и их производные для прямоугольного элемента в глобальных координатах
	virtual double bfunc1D(size_t func_num, double x0, double x1, double x) {
		return 1;
	};
	virtual double dbfunc1D(size_t func_num, double x0, double x1, double x) {
		return 1;
	};

	// базисные функции и их производные для четырехугольного элемента в шаблонных координатах
	virtual double bfunc1D(size_t func_num, double ksi) {
		return 1;
	}
	virtual double dbfunc1D(size_t func_num, double ksi) {
		return 1;
	}
	double gauss_points_local[3];
	double gauss_weights[3];




};




class Rectangle : public Element {
public:
	//etype type = RECT;
	//static void init();
	//Rectangle();
	Rectangle(vector<Point>& local_nodes);
	// Одномерные линейные функции формы
	// Несмотря на аргументы x0, x1, x, функции также применяются для вычислений функций по Y, так как имеют одинаковый вид
	double bfunc1D(size_t fucn_num, double x0, double x1, double x) override;

	// Производные функций форм
	double dbfunc1D(size_t func_num, double x0, double x1, double x) override;

	// Двумерные билинейные функции формы
	// в качестве независимой переменной x сюда будут посылаться точки Гаусса, по которым будем интегрировать
	// то есть эта функция должна быть в цикле по точкам Гаусса
	double phi(size_t func_num, Point& from, Point& to, double x, double y);

	// производные двумерных функции форм

	double dphi(size_t var, size_t i, size_t j, Point& from, Point& to, double x, double y);

	Point dphi(size_t num, Point& from, Point& to, Point& p);


	// интегрирование методом Гаусса по элементу
	double Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) override;
	// интегрирование методом Гаусса по ребру
	double Edge_IntegrateGauss3(Point& from, Point& to, size_t num) override;
	// вычисление локальной матрицы жесткости для элемента
	void CalculateLocalStiffnessMatrix() override;

	
	// вычисление локального вектора нагрузок
	void CalculateLocalLoadVector( Load& P) override;


	Point GetElementCentroid() override;
	double GetElementalStrainX(double* u, Point& p) override;
	double GetElementalStrainY(double* u, Point& p) override;
	double GetElementalStrainXY(double* u, Point& p) override;


	double GetElementalStressX(double x_strain, double y_strain) override;
	double GetElementalStressY(double x_strain, double y_strain) override;
	double GetElementalStressXY(double xy_strain) override;


	Point integrate_points[9];
};

class Quadrilateral : public Element {
public:
	Quadrilateral(vector<Point>& local_nodes);
	Point integrate_points[9];
	//etype type = QUAD;
	// одномерные базисные функции в шаблонных координатах
	double bfunc1D(size_t fucn_num, double ksi) override;
	// производные одномерных бф в шаблонных координатах
	double dbfunc1D(size_t func_num, double ksi) override;

	// двумрная бф в шаблонных координатах
	double phi(size_t func_num, double ksi, double eta);

	// производные двумерных бф в шаблонных координатах
	Point dphi(size_t func_num, double ksi, double eta);

	// производные двумерных бф в глобальных координатах
	Point dphi_global(size_t func_num, Point& p);

	double grad_bfunc_2d(size_t var, size_t num, double ksi, double eta);

	// Якобиан преобразования координат шаблонного элемента [ksi,eta] в координаты глобального [x,y]
	double det_J(double ksi, double eta);

	double Element_IntegrateGauss3(size_t i, size_t j, size_t var1, size_t var2);
	double Edge_IntegrateGauss3(size_t num);

	void CalculateLocalStiffnessMatrix() override;
	void CalculateLocalLoadVector(Load& P) override;

	// перевод шаблонных координат [0,1] в глобальные
	Point to_global(Point& p);
	// перевод глобальных координат в шаблонные [0,1]
	Point to_local(Point& p);

	Point GetElementCentroid() override;
	double GetElementalStrainX(double* u, Point& p) override;
	double GetElementalStrainY(double* u, Point& p) override;
	double GetElementalStrainXY(double* u, Point& p) override;

	double GetElementalStressX(double x_strain, double y_strain) override;
	double GetElementalStressY(double x_strain, double y_strain) override;
	double GetElementalStressXY(double xy_strain) override;


private:
	double alpha0, alpha1, alpha2;
	double beta1, beta2, beta3, beta4, beta5, beta6;


};


// класс нагрузки
class Load {
public:
	Load();

	double Px, Py;	// компоненты поверхностных сил
};

// этот блок будет элементом локальной матрицы жесткости
// глобальной тоже
class block2x2 {
public:
	block2x2();
	block2x2(double val11, double val12, double val21, double val22);
	double val11, val12,
		   val21, val22;
	// сложение блоков (блок, к которому прибавляем - слева, а прибавляемый блок - справа)
	block2x2 operator + (block2x2 block2);
	block2x2 operator += (block2x2 block2);
	block2x2 operator = (double val);


};

// этот блок будет элементом вектора нагрузок как локального, так и глобального

class block1x2 {
public:
	block1x2();
	block1x2(double val1, double val2);
	double val1, val2;
	block1x2 operator +(block1x2 block2);
	block1x2 operator += (block1x2 block2);
	block1x2 operator = (double val);
};

class FEM {
public:

	void AssembleGlobalStiffnessMatrix(Mesh& mesh);
	vector<vector<double>> GlobalStiffnessMatrix;
	vector<vector<block2x2>> GlobalStiffnessMatrix_block;

	void AssembleGlobalLoadVector(Mesh& mesh);
	vector<double> GlobalLoadVector;
	vector<block1x2> GlobalLoadVector_block;

	vector<int> loaded_nodes;
	vector<int> fixed_nodes;

	void GeneratePortrait(Portrait& portrait, vector<vector<double>> GSM, int &ig_n_1);


	vector<double> elemental_stressX;
	vector<double> elemental_stressY;
	vector<double> elemental_stressXY;
	vector<double> Elemental_VonMises_Stress;

	Element* GetElement(Point& p, Mesh& mesh);

	void Get_X_Stresses(double* u, Mesh& mesh, string output_folder);




private:

};

#endif // !FEM_H

