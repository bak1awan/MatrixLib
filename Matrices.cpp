#include <iostream >
#include "Matrix.h"
#include <iomanip>
#include "Optimization.h"

using namespace std;

int main() {

	vector<double> v1 = { 1, 2, 3 };
	VectorT<double> v2 = transpose(v1);
	Matrix<double> E = v1 * v2;

	// Matrix<double> E = v1 * v2;

	/*
	vector<vector<Matrix<double>>> v;
	v.resize(3);
	vector<Matrix<double>>::iterator it;
	for (int i = 0; i < v.size(); i++) {
		it = v[i].begin();
		for (int j = 0; j < v.size(); j++) {
			it = v[i].emplace(it, 3, 3);
		}
	}
	return 0;
	*/

}