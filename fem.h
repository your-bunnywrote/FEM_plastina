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
	void init();
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

protected:
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
	Point integrate_points[9];



};




class Rectangle : public Element {
public:
	etype type = RECT;
	//static void init();
	Rectangle();
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
	//double dphi(size_t var, size_t num1, size_t num2, Point& from, Point& to, double x, double y);
	double dphi(size_t var, size_t num1, size_t num2, Point from, Point to, double x, double y);


	//static double gauss_points_local[3];
	//static Point integrate_points[9];
	//static double gauss_weights[3];
	// интегрирование методом Гаусса по элементу
	double Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) override;
	// интегрирование методом Гаусса по ребру
	double Edge_IntegrateGauss3(Point& from, Point& to, size_t num) override;
	// вычисление локальной матрицы жесткости для элемента
	void CalculateLocalStiffnessMatrix() override;
	// vector<int> fixed_nodes;
	// vector<int> loaded_nodes;
	// сюда будем вносить блоки, так как матрица жесткости имеет блочную структуру
	//vector<vector<block2x2>> LocalStiffnessMatrix_block;
	//vector<vector<block2x2>> GlobalStiffnessMatrix_block;
	// поэлементные форматы матрицы жесткости
	//vector<vector<double>> LocalStiffnessMatrix;
	//vector<vector<double>> GlobalStiffnessMatrix;
	// сборка глобальной матрицы жесткости
	// void Assemble_GlobalStiffnessMatrix(Mesh& mesh);
	
	// вычисление локального вектора нагрузок
	void CalculateLocalLoadVector( Load& P) override;
	// vector<block1x2> LocalLoadVector_block;
	// vector<block1x2> GlobalLoadVector_block;

	// поэлементные форматы векторов нагрузок
	//vector<double> LocalLoadVector;
	//vector<double> GlobalLoadVector;

	//void Assemble_GlobalLoadVector(Mesh& mesh);

	// конвертирование матрицы в разреженный формат для использования в PARDISO
	
	//void GeneratePortrait(Portrait & Matrix, vector<vector<double>> GSM, int &ja_sz);
};

class Quadrilateral : public Element {
public:
	etype type = QUAD;
	// одномерные базисные функции в шаблонных координатах
	double bfunc1D(size_t fucn_num, double ksi) override;
	// производные одномерных бф в шаблонных координатах
	double dbfunc1D(size_t func_num, double ksi) override;

	double phi(size_t func_num, Point& ksi);

	double dphi(size_t func_num, Point& ksi);

	// Якобиан преобразования координат шаблонного элемента [ksi,eta] в координаты глобального [x,y]
	double det_J(Point& p);

	double Element_IntegrateGauss3(Point& from, Point& to, size_t num1, size_t num2, size_t var) override;
	double Edge_IntegrateGauss3(Point& from, Point& to, size_t num) override;

	void CalculateLocalStiffnessMatrix() override;
	void CalculateLocalLoadVector(Load& P) override;

	Point to_global(Point& p);
	Quadrilateral();
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
	Rectangle rect;
	Quadrilateral quad;

	void AssembleGlobalStiffnessMatrix(Mesh& mesh);
	vector<vector<double>> GlobalStiffnessMatrix;
	vector<vector<block2x2>> GlobalStiffnessMatrix_block;

	void AssembleGlobalLoadVector(Mesh& mesh);
	vector<double> GlobalLoadVector;
	vector<block1x2> GlobalLoadVector_block;

	vector<int> loaded_nodes;
	vector<int> fixed_nodes;

	void GeneratePortrait(Portrait& portrait, vector<vector<double>> GSM, int &ig_n_1);

private:

};

#endif // !FEM_H

