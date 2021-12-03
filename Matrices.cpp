#include <iostream >
#include "Matrix.h"
#include <iomanip>
#include "Optimization.h"

using namespace std;

int main() {
	vector<vector<Matrix>> v;
	v.resize(3);
	vector<Matrix>::iterator it;
	for (int i = 0; i < v.size(); i++) {
		it = v[i].begin();
		for (int j = 0; j < v.size(); j++) {
			it = v[i].emplace(it, 3, 3);
		}
	}
	return 0;
}