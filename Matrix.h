#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;

class VectorT;

class Matrix {
public:
	// конструктор
	int rows;
	int cols;
	Matrix(int nRow = 1, int nCol = 1);

	// транспонирование матрицы
	Matrix transpose();

	// Унарный минус
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

class VectorT {
public:
	VectorT(int size = 1) : v(size) {}

	// Конструктор копирования
	VectorT(vector<double>& v1) : v(v1) {}

	// Конструктор перемещения
	VectorT(vector<double>&& v1) : v(v1) {}

	// Конструктор присваивания копированием
	VectorT operator= (vector<double>& v1) {
		this->v = v1;
	}

	// Конструктор присваивания перемещением
	VectorT operator= (vector<double>&& v1) {
		this->v = v1;
	}

	int size();

	// Умножение транспонированного вектора на обычный - число
	double operator* (vector<double>&);

	// Сумма транспонированных векторов
	VectorT operator+ (VectorT&);

	// Разность транспонированных векторов
	VectorT operator- (VectorT&);

	// Оператор []
	double& operator[](int);
	const double& operator[](int i) const;

	// Сложение с числом
	VectorT operator+ (double);
	friend VectorT operator+ (double, VectorT&);

	// Вычитание числа
	VectorT operator- (double);

	// Умножение на число
	VectorT operator* (double);
	friend VectorT operator* (double, VectorT&);

	// Деление на число
	VectorT operator/ (double);

	// Вывод транспонированного вектора
	friend ostream& operator<< (ostream&, VectorT&);

	// Считывание транспонированного вектора
	friend istream& operator>> (istream&, VectorT&);

private:
	// Поле, содержащее веткор
	vector<double> v;
};

// Умножение вектора на число
vector<double> operator* (double, const vector<double>&);
vector<double> operator* (const vector<double>&, double);

// Сложение вектора с числом
vector<double> operator+ (double, const vector<double>&);
vector<double> operator+ (const vector<double>&, double);

vector<double> operator- (const vector<double>&, double);

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

VectorT transpose(vector<double>&);

// Генерация единичной матрицы
Matrix generateE(int);

// Умножение обычного вектора на транспонированный - матрица
Matrix operator* (const vector<double>&, const VectorT&);

#endif