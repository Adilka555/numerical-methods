#include <stdio.h>
#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

using std::vector;

class Matrix {
private:
  unsigned m_rowSize;
  unsigned m_colSize;
  vector<vector<double>> m_matrix;
public:
  Matrix();
  Matrix(unsigned, unsigned, double);
  Matrix(const char *);
  //Matrix(const Matrix &);
  //~Matrix();
  Matrix operator+(Matrix &);
  Matrix operator-(Matrix &);
  Matrix operator*(Matrix &);
  Matrix operator+(double);
  Matrix operator-(double);
  Matrix operator*(double);
  Matrix operator/(double);
//  Matrix operator+(vector<double>);
  Matrix transpose();
  double& operator()(const unsigned &, const unsigned &);
  void print() const;
  void my_cin();
  unsigned getRows() const;
  unsigned getCols() const;
  vector<double>& operator[](const unsigned&);
};
