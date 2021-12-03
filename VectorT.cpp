#include "VectorT.h"
const double precision = 1e-5;

double& VectorT::operator[](int i) {
	return v[i];
}

const double& VectorT::operator[](int i) const {
	return v[i];
}

int VectorT::size() {
	return this->v.size();
}

// ��������� ������������������ ������� �� ������� - �����
double VectorT::operator* (vector<double>& v1) {
	double result(0);
	for (int i = 0; i < v.size(); i++)
		result += v[i] * v1[i];
	return result;
}

// ����� ����������������� ��������
VectorT VectorT::operator+ (VectorT& v1) {
	VectorT result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] + v1[i];
	return result;
}

// �������� ����������������� ��������
VectorT VectorT::operator- (VectorT& v1) {
	VectorT result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] - v1[i];
	return result;
}

// �������� � ������
VectorT  VectorT::operator+ (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] + n;
	}
	return result;
}

VectorT operator+ (double n, VectorT& v) {
	return v + n;
}

// ��������� �����
VectorT VectorT::operator- (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] - n;
	}
	return result;
}

// ��������� �� �����
VectorT VectorT::operator* (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n * v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n * v[i];
	}
	return result;
}

VectorT operator* (double n, VectorT& v) {
	return v * n;
}

// ������� �� �����
VectorT VectorT::operator/ (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n / v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n / v[i];
	}
	return result;
}

// ����� ������������������ �������
ostream& operator<< (ostream& out, VectorT& vT) {
	for (double el : vT.v) {
		out << el << ' ';
	}
	return out;
}

// ���������� ������������������ �������
istream& operator>> (istream& in, VectorT& vT) {
	for (int i = 0; i < vT.size(); i++) {
		cout << " Input [" << i << "]" << " element : ";
		in >> vT[i];
	}
	return in;
}

VectorT transpose(vector<double>& v) {
	return VectorT(v);
}