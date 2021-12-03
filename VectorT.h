#ifndef VECTOR_T_H
#define VECTOR_T_H

#include <vector>
#include <iostream>

using namespace std;

class VectorT {
public:
	VectorT(int size = 1) : v(size) {}

	// Конструктор копирования
	VectorT(vector<double>& v1) : v(v1) {}

	// Конструктор присваивания копированием
	VectorT operator= (vector<double>& v1) {
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


VectorT transpose(vector<double>&);
#endif