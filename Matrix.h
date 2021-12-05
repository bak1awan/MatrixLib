#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;

const double precision = 1e-5;
const int maxIteration = 40;

// Предварительное объявление класса VectorT
template <class variableType>
class VectorT;



// Далее следует функционал класса матриц Matrix



template <class variableType>
class Matrix {
public:
	// конструктор, по умолчанию генерирует единичную матрицу
	Matrix(int nRow = 1, int nCol = 1) : rows(nRow), cols(nCol) {
		arr = vector<vector<variableType>>(nRow, vector<variableType>(nCol, 0));
	};

	// транспонирование матрицы
	Matrix<variableType> transpose() const;

	// Унарный минус
	Matrix<variableType> operator- ();

	// перегрузка оператора + с матрицей
	Matrix<variableType> operator+ (const Matrix<variableType>&);

	// перегрузка оператора - с матрицей
	Matrix<variableType> operator- (const Matrix<variableType>&);

	// перегрузка оператора * с матрицей
	Matrix<variableType> operator* (const Matrix<variableType>&);

	// Умножение матрицы на число справа
	Matrix<variableType> operator* (variableType);

	// Умножение матрицы на число слева
	template <class T> friend Matrix<T> operator* (T, Matrix<T>&);

	// Умножение матрицы на вектор справа
	vector<variableType> operator* (const vector<variableType>&);

	// Умножение матрицы на вектор слева
	template <class T> friend vector<T> operator* (const vector<T>&, Matrix<T>&);

	// Сложение матрицы с числом справа
	Matrix<variableType> operator+ (variableType);

	// Сложение матрицы с числом слева
	template <class T> friend Matrix<T> operator+ (T, Matrix<T>&);

	// перегрузка << для матриц
	template <class T> friend ostream& operator << (ostream&, Matrix<T>&);

	// перегрузка >> для матриц
	template <class T> friend istream& operator >> (istream&, Matrix<T>&);

	// Определитель через QR-разложение
	variableType LUDeterminant();

	// Обратная матрица через LU-разложение
	Matrix<variableType> LUInverse();

	// Обратная матрица через QR-разложение
	Matrix<variableType> QRInverse();

	// Норма матрицы
	variableType norm();

	Matrix<variableType> generateE(int);

	// Проверка матрицы на диагональное преобладание
	bool diagonal();

	// Обратная матрица методом Шульца (Newton-Schulz-Hotelling)
	Matrix<variableType> ShultzInverse(double = 0.001);

	// Перегрузка оператора [] для матриц
	vector<variableType>& operator[] (const int i);

	// Перегрузка оператора [] для константных матриц
	const vector<variableType>& operator[] (const int i) const;

	// Количество строк
	int rows;

	// Количество колонок
	int cols;

private:
	// Сама матрица представляется просто вектором векторов
	vector<vector<variableType>> arr;
};


// Далее следует функционал для нашего класса транспонированных векторов VectorT




template <class variableType>
class VectorT {
public:
	// Конструктор
	VectorT(int size = 1) : v(size) {}

	// Конструктор копирования
	VectorT(vector<variableType> v1) : v(v1) {}

	// Размер вектора VectorT
	int size();

	// Умножение транспонированного вектора на обычный - число
	variableType operator* (vector<variableType>&);

	// Сумма транспонированных векторов
	VectorT<variableType> operator+ (VectorT<variableType>&);

	// Разность транспонированных векторов
	VectorT<variableType> operator- (VectorT<variableType>&);

	// Оператор []
	variableType& operator[](int);
	const variableType& operator[](int i) const;

	// Сложение с числом справа
	VectorT<variableType> operator+ (variableType);

	// Сложение с числом слева
	template <class T> friend VectorT<T> operator+ (T, VectorT<T>&);

	// Вычитание числа
	VectorT<variableType> operator- (variableType);

	// Умножение на число справа
	VectorT<variableType> operator* (variableType);

	// Умножение на число слева
	template <class T> friend VectorT<T> operator* (T, VectorT<T>&);

	// Деление на число (возможно только справа)
	VectorT<variableType> operator/ (variableType);

	// Вывод транспонированного вектора
	template <class T> friend ostream& operator<< (ostream&, VectorT<T>&);

	// Считывание транспонированного вектора
	template <class T> friend istream& operator>> (istream&, VectorT<T>&);

private:
	// Поле, содержащее вектор
	vector<variableType> v;
};



// Далее следует ункционал для встроенных векторов std::vector



// Умножение вектора на число слева
template<typename variableType>
vector<variableType> operator* (variableType, const vector<variableType>&);

// Умножение вектора на число справа
template<typename variableType>
vector<variableType> operator* (const vector<variableType>&, variableType);

// Сложение вектора с числом слева
template<typename variableType>
vector<variableType> operator+ (variableType, const vector<variableType>&);

// Сложение вектора с числом справа
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, variableType);

// Вычитание числа из вектора (только справа возможно)
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, variableType);

// Сложение векторов
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, const vector<variableType>&);

// Вычитание векторов
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, const vector<variableType>&);

// Вычисление длины вектора
template<typename variableType>
variableType vectorLength(const vector<variableType>&);

// Скалярное произведение векторов
template<typename variableType>
variableType scalarOp(const vector<variableType>&, const vector<variableType>&);

// Транспонирование обычного вектора - получаем объект класса VectorT
template<typename variableType>
VectorT<variableType> transpose(vector<variableType>&);

// Тензорное умножение векторов
template<typename variableType>
Matrix<variableType> operator* (const vector<variableType>&, const VectorT<variableType>&);
#endif