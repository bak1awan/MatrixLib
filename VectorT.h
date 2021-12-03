#ifndef VECTOR_T_H
#define VECTOR_T_H

#include <vector>
#include <iostream>

using namespace std;

class VectorT {
public:
	VectorT(int size = 1) : v(size) {}

	// ����������� �����������
	VectorT(vector<double>& v1) : v(v1) {}

	// ����������� ������������ ������������
	VectorT operator= (vector<double>& v1) {
		this->v = v1;
	}

	int size();

	// ��������� ������������������ ������� �� ������� - �����
	double operator* (vector<double>&);

	// ����� ����������������� ��������
	VectorT operator+ (VectorT&);

	// �������� ����������������� ��������
	VectorT operator- (VectorT&);

	// �������� []
	double& operator[](int);
	const double& operator[](int i) const;

	// �������� � ������
	VectorT operator+ (double);
	friend VectorT operator+ (double, VectorT&);

	// ��������� �����
	VectorT operator- (double);

	// ��������� �� �����
	VectorT operator* (double);
	friend VectorT operator* (double, VectorT&);

	// ������� �� �����
	VectorT operator/ (double);

	// ����� ������������������ �������
	friend ostream& operator<< (ostream&, VectorT&);

	// ���������� ������������������ �������
	friend istream& operator>> (istream&, VectorT&);

private:
	// ����, ���������� ������
	vector<double> v;
};


VectorT transpose(vector<double>&);
#endif