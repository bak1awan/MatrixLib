#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include "Matrix.h"

// ������� ���� ����� QR-����������
template<typename variableType>
void QRSolution(const Matrix<variableType>&, vector<variableType>&, const vector<variableType>&);

// ������� ���� ����� LU-����������
template<typename variableType>
void LUSolution(const Matrix<variableType>&, vector<variableType>&, const vector<variableType>&);

template<typename variableType>
void GaussSeidelSolution(const Matrix<variableType>&, vector<variableType>&, const vector<variableType>&, variableType);

#endif
