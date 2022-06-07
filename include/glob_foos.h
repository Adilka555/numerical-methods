#include <stdio.h>
#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <fstream>
#include "matrix.h"
using namespace std;
using std::vector;

vector<double> operator +(const vector<double>& , const vector<double>&);
void print_v(vector<double> );
vector<double> linspace(const double, const double, unsigned);
Matrix inver(Matrix );
double u_a(double, double);
double f(double, double);
double west(double);
double east(double);
double south(double);
double nord(double);
vector<double> mul(Matrix , vector<double> );
Matrix my_sweep(Matrix, Matrix);
double approx_m(unsigned, Matrix, Matrix, double);
