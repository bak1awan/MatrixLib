#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;



class Matrix {
public:
	// конструктор
	int rows;
	int cols;
	Matrix(int nRow, int nCol);

	// транспонирование матрицы
	Matrix transpose();

	Matrix operator- ();

	// перегрузка оператора +
	Matrix operator+ (const Matrix&);

	// перегрузка оператора -
	Matrix operator- (const Matrix&);

	// перегрузка оператора * с матрицей
	Matrix operator* (const Matrix&);

	// перегрузка оператора * с числом справа
	friend Matrix operator* (Matrix&, double);

	// перегрузка оператора * с числом слева
	friend Matrix operator* (double, Matrix&);

	// перегрузка оператора * с вектором справа
	friend vector<double> operator* (Matrix&, const vector<double>&);

	// перегрузка оператора * с вектором слева
	friend vector<double> operator* (const vector<double>&, Matrix&);

	// перегрузка оператора + с числом справа
	friend Matrix operator+ (Matrix&, double);

	// перегрузка оператора + с числом слева
	friend Matrix operator+ (double, Matrix&);

	// перегрузка <<
	friend std::ostream& operator << (std::ostream&, const Matrix&);

	// перегрузка >>
	friend std::istream& operator >> (std::istream&, Matrix&);

	// QR-разложение
	void QRDecomp(Matrix&, Matrix&);

	// LU-разложение
	void LUDecomp(Matrix&, Matrix&);

	// Определитель через QR-разложение
	double LUDeterminant();

	// Обратная матрица через LU-разложение
	Matrix LUInverse();

	// Обратная матрица через QR-разложение
	Matrix QRInverse();

	// Решение СЛАУ методом Гаусса – Зейделя
	void GaussSeidelSolution(vector<double>&, const vector<double>&, double = 0.001);

	// Норма матрицы
	double norm();

	// Проверка матрицы на диагональное преобладание
	bool diagonal();

	// Обратная матрица методом Шульца (Newton-Schulz-Hotelling)
	Matrix ShultzInverse(double = 0.001);

	// Решение СЛАУ через LU-разложение
	void LUSolution(vector<double>&, const vector<double>&);

	// Решение СЛАУ через QR-разложение
	void QRSolution(vector<double>&, const vector<double>&);

	void cholesky(Matrix&);

	// Собственные числа через QR-разложение
	void QREigen(vector<double>&);

	// Собственное число через соотношение Рэлея
	double RayleighEigen(double = 0.01);

	// Перегрузка оператора [] для матриц
	vector<double>& operator[] (const int i);

	const vector<double>& operator[] (const int i) const;

	// поле для инициализации
	vector<vector<double>> arr;


};

// Умножение вектора на число
vector<double> operator* (double, const vector<double>&);
vector<double> operator* (const vector<double>&, double);

// Сложение векторов
vector<double> operator+ (const vector<double>&, const vector<double>&);

// Вычитание векторов
vector<double> operator- (const vector<double>&, const vector<double>&);

// Вычисление длины вектора
double vectorLength(const vector<double>&);

// Тензорное произведение векторов
Matrix operator*(const vector<double>&, const vector<double>&);

// Скалярное произведение векторов
double scalarOp(const vector<double>&, const vector<double>&);

// Генерация единичной матрицы
Matrix generateE(int);

#endif