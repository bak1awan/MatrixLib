#include <iostream >
#include "Matrix.h"
#include <iomanip>
#include "Optimization.h"

using namespace std;

int main() {
	vector<vector<Matrix>> v;
	v.reserve(3);
	vector<Matrix> v1;
	v1.reserve(3);
	Matrix E = generateE(3);
	for (int i = 0; i < 3; i++) v1.push_back(E);
	for (int i = 0; i < 3; i++) v.push_back(v1);
	return 0;
}