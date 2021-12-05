#ifndef EIGEN_H
#define EIGEN_H

#include "Matrix.h"

// ����������� ����� ����� QR-����������
template <typename variableType>
void QREigen(const Matrix<variableType>&, vector<variableType>&);

// ����������� ����� ����� ����������� �����
template <typename variableType>
variableType RayleighEigen(const Matrix<variableType>&, variableType);

#endif