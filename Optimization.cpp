#include "Optimization.h"
#include "Matrix.h"
#include <stdio.h>
#define _CRT_SECURE_NO_WARNINGS

const int maxStep = 1000;
const double h = 1e-5;

double f3_sphere(const vector<double>& x) {
	double result = 0;
	for (int i = 0; i < x.size(); i++)
		result += pow(x[i], 2);
	return result;
}

double f3_hyper_ellipsoid(const vector<double>& x) {
	return pow(x[0], 2) + pow(x[0] + x[1], 2);
}

double f3_test(const vector<double>& x) {
	return pow(x[0], 2) - x[0] * x[1] + pow(x[1], 2) + 9 * x[0] - 6 * x[1] + 20;
}

double f_diff(const vector<double>& x, double(*func)(const vector<double>&), int i) {
	vector<double> xplus = x;
	xplus[i] = x[i] + h;
	return (func(xplus) - func(x)) / h;
}

vector<double> gradient(const vector<double>& x, double (*func)(const vector<double>&)) {
	vector<double> grad(x.size(), 0);
	for (int i = 0; i < x.size(); i++)
		grad[i] = f_diff(x, func, i);
	return grad;
}

Matrix<double> hessian(const vector<double>& x, double (*func)(const vector<double>&), vector<double>& grad) {
	Matrix<double> hessian(x.size(), x.size());
	vector<double> xplus = x;
	vector<double> gplus = grad;
	for (int i = 0; i < x.size(); i++) {
		xplus[i] += h;
		gplus = gradient(xplus, func);
		for (int j = 0; j < x.size(); j++)
			hessian[i][j] = (gplus[j] - grad[j]) / h;
		xplus[i] = x[i];
	}
	
	return hessian.LUInverse();
}

vector<double> gradientDescent(const vector<double>& x0, double (*func)(const vector<double>&)) {
	vector<double> x = x0;
	vector<double> grad(x0.size(), 0);
	FILE* hfile;
	fopen_s(&hfile, "trajectory.dat", "w");
	const double alpha = 0.1;
	int counter = 0;
	for (int i = 0; i < maxStep; i++) {
		counter++;

		for (int j = 0; j < x.size(); j++) {
			fprintf(hfile, "%f ", x[j]);
		}

		fprintf(hfile, "\n");
		grad = gradient(x, func);
		if (vectorLength(grad) < precision) break;
		x = x - alpha * grad;
	}
	fclose(hfile);
	cout << "The Gradient Descent method find the solution in " << counter << " iterations.\n";
	return x;
}

vector<double> newton(const vector<double>& x0, double (*func)(const vector<double>&)) {
	vector<double> x = x0;
	vector<double> grad(x0.size(), 0);

	Matrix<double> hess(x.size(), x.size());

	const double alpha = 1;
	int counter = 0;
	for (int i = 0; i < maxStep; i++) {
		counter++;
		grad = gradient(x, func);
		if (vectorLength(grad) < precision) break;
		hess = hessian(x, func, grad);
		x = x - alpha * (hess * grad);
	}
	cout << "The Newton method find the solution in " << counter << " iterations.\n";
	return x;
}

vector<double> BFGS(const vector<double>& x0, double (*func)(const vector<double>&)) {
	// Начальное приближение квазиньютоновской матрицы
	Matrix<double> H(x0.size(), x0.size());

	// Единичная матрица для вычислений квазиньютоновской
	Matrix<double> I(x0.size(), x0.size());

	// k - результат скалярного умножения векторов
	// double k = 0;

	// alpha - специальный коэффициент, определяемый из условий Вульфе, но в нашем случае положенный просто единицей
	double alpha = 1;

	// Специальные матрица с вектором для промежуточных результатов из-за особенности реализованной мною перегрузки
	Matrix<double> interM1(x0.size(), x0.size());

	// Вектора для изменения градиента итерации
	vector<double> gnow = x0;
	vector<double> gprev = x0;
	vector<double> y = x0;

	// Вектора для шага алгоритма итерации
	vector<double> xnow = x0;
	vector<double> xprev = x0;
	vector<double> s = x0;

	// Точка, в направлении которой будет производить поиск
	vector<double> p = x0;
	int counter = 0;
	for (int i = 0; i < maxStep; i++) {
		counter++;
		// Считаем градиент в точке x
		gnow = gradient(xnow, func);
		if (vectorLength(gnow) < precision) break;
		interM1 = -H;
		// Определяем точку, в направлении которой будем производить поиск
		p = interM1 * gnow;
		// Изменяем значения x
		xnow = xprev + alpha * p;
		// Определяем вектор s
		s = xnow - xprev;
		// После этого можно запомнить нынешнее значение точки как прошлое
		xprev = xnow;
		// Запоминаем предыдущее значение градиента в предыдущей точке
		gprev = gnow;
		// Рассчитываем значение градиента в текущей (измененной в строке 130) точке
		gnow = gradient(xnow, func);
		// Определяем вектор y
		y = gnow - gprev;
		// Считаем k как скалярное произведение
		// k = 1.0 / (transpose(y) * s);
		// Здесь очень замудренный расчет новой квазиньютоновской матрицы
		// из-за особенностей реализованной мною перегрузки операторов для класса Matrix
		// interM1 = s * y;
		// interM2 = y * s;
		// interM3 = s * s;
		H = ((I - (1.0 / (transpose(y) * s)) * s * transpose(y)) *
			H * (I - (1.0 / (transpose(y) * s)) * y * transpose(s))) +
			((1.0 / (transpose(y) * s)) * s * transpose(s));
	}
	cout << "The BFGS method find the solution in " << counter << " iterations.\n";
	// возвращаем результат
	return xnow;
}