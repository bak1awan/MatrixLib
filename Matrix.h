#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;

const double precision = 1e-5;
const int maxIteration = 40;

// ��������������� ���������� ������ VectorT
template <class variableType>
class VectorT;



// ����� ������� ���������� ������ ������ Matrix



template <class variableType>
class Matrix {
public:
	// �����������, �� ��������� ���������� ��������� �������
	Matrix(int nRow = 1, int nCol = 1) : rows(nRow), cols(nCol) {
		arr = vector<vector<variableType>>(nRow, vector<variableType>(nCol, 0));
	};

	// ���������������� �������
	Matrix<variableType> transpose() const;

	// ������� �����
	Matrix<variableType> operator- ();

	// ���������� ��������� + � ��������
	Matrix<variableType> operator+ (const Matrix<variableType>&);

	// ���������� ��������� - � ��������
	Matrix<variableType> operator- (const Matrix<variableType>&);

	// ���������� ��������� * � ��������
	Matrix<variableType> operator* (const Matrix<variableType>&);

	// ��������� ������� �� ����� ������
	Matrix<variableType> operator* (variableType);

	// ��������� ������� �� ����� �����
	template <class T> friend Matrix<T> operator* (T, Matrix<T>&);

	// ��������� ������� �� ������ ������
	vector<variableType> operator* (const vector<variableType>&);

	// ��������� ������� �� ������ �����
	template <class T> friend vector<T> operator* (const vector<T>&, Matrix<T>&);

	// �������� ������� � ������ ������
	Matrix<variableType> operator+ (variableType);

	// �������� ������� � ������ �����
	template <class T> friend Matrix<T> operator+ (T, Matrix<T>&);

	// ���������� << ��� ������
	template <class T> friend ostream& operator << (ostream&, Matrix<T>&);

	// ���������� >> ��� ������
	template <class T> friend istream& operator >> (istream&, Matrix<T>&);

	// ������������ ����� QR-����������
	variableType LUDeterminant();

	// �������� ������� ����� LU-����������
	Matrix<variableType> LUInverse();

	// �������� ������� ����� QR-����������
	Matrix<variableType> QRInverse();

	// ����� �������
	variableType norm();

	Matrix<variableType> generateE(int);

	// �������� ������� �� ������������ ������������
	bool diagonal();

	// �������� ������� ������� ������ (Newton-Schulz-Hotelling)
	Matrix<variableType> ShultzInverse(double = 0.001);

	// ���������� ��������� [] ��� ������
	vector<variableType>& operator[] (const int i);

	// ���������� ��������� [] ��� ����������� ������
	const vector<variableType>& operator[] (const int i) const;

	// ���������� �����
	int rows;

	// ���������� �������
	int cols;

private:
	// ���� ������� �������������� ������ �������� ��������
	vector<vector<variableType>> arr;
};


// ����� ������� ���������� ��� ������ ������ ����������������� �������� VectorT




template <class variableType>
class VectorT {
public:
	// �����������
	VectorT(int size = 1) : v(size) {}

	// ����������� �����������
	VectorT(vector<variableType> v1) : v(v1) {}

	// ������ ������� VectorT
	int size();

	// ��������� ������������������ ������� �� ������� - �����
	variableType operator* (vector<variableType>&);

	// ����� ����������������� ��������
	VectorT<variableType> operator+ (VectorT<variableType>&);

	// �������� ����������������� ��������
	VectorT<variableType> operator- (VectorT<variableType>&);

	// �������� []
	variableType& operator[](int);
	const variableType& operator[](int i) const;

	// �������� � ������ ������
	VectorT<variableType> operator+ (variableType);

	// �������� � ������ �����
	template <class T> friend VectorT<T> operator+ (T, VectorT<T>&);

	// ��������� �����
	VectorT<variableType> operator- (variableType);

	// ��������� �� ����� ������
	VectorT<variableType> operator* (variableType);

	// ��������� �� ����� �����
	template <class T> friend VectorT<T> operator* (T, VectorT<T>&);

	// ������� �� ����� (�������� ������ ������)
	VectorT<variableType> operator/ (variableType);

	// ����� ������������������ �������
	template <class T> friend ostream& operator<< (ostream&, VectorT<T>&);

	// ���������� ������������������ �������
	template <class T> friend istream& operator>> (istream&, VectorT<T>&);

private:
	// ����, ���������� ������
	vector<variableType> v;
};



// ����� ������� ��������� ��� ���������� �������� std::vector



// ��������� ������� �� ����� �����
template<typename variableType>
vector<variableType> operator* (variableType, const vector<variableType>&);

// ��������� ������� �� ����� ������
template<typename variableType>
vector<variableType> operator* (const vector<variableType>&, variableType);

// �������� ������� � ������ �����
template<typename variableType>
vector<variableType> operator+ (variableType, const vector<variableType>&);

// �������� ������� � ������ ������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, variableType);

// ��������� ����� �� ������� (������ ������ ��������)
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, variableType);

// �������� ��������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, const vector<variableType>&);

// ��������� ��������
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, const vector<variableType>&);

// ���������� ����� �������
template<typename variableType>
variableType vectorLength(const vector<variableType>&);

// ��������� ������������ ��������
template<typename variableType>
variableType scalarOp(const vector<variableType>&, const vector<variableType>&);

// ���������������� �������� ������� - �������� ������ ������ VectorT
template<typename variableType>
VectorT<variableType> transpose(vector<variableType>&);

// ��������� ��������� ��������
template<typename variableType>
Matrix<variableType> operator* (const vector<variableType>&, const VectorT<variableType>&);
#endif