#include "glob_foos.h"
using namespace std;
using std::vector;
#include <tuple>


const double PI = 3.14159;

vector<double> operator + (const vector<double>& v1, const vector<double>& v2 )
{
  /*Do quick size check on vectors before proceeding*/
  vector<double> result(v1.size());
  for (unsigned int i = 0; i < result.size(); ++i)
  {
    result[i]=v1[i]+v2[i];
  }
  return result;
}

void print_v(vector<double> v) {
	for(unsigned i = 0; i < v.size(); i++)
		cout << v[i] << " ";
}

vector<double> linspace(const double start_in, const double end_in, unsigned num_in)
{
  vector<double> linspaced;

  double start = start_in;
  double end = end_in;
  double num = num_in;

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);
  return linspaced;
}


Matrix inver(Matrix a) {
	int n = a.getRows();
	Matrix e(n,n,0);
	for (int i = 0; i < n; i++) e(i, i) = 1;
	double eps = 1.0e-20;
	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < n; i++) {
			try {
				if (fabs(a(k, k) ) < eps) throw "Error!";
			}
			catch (const char* str) {cout << str << endl;}
			double tmp = a(i, k) / a(k, k);
			a(i, k) = 0;
			for (int j = k + 1; j < n; j++) a(i, j) -= tmp * a(k, j);
			for (int s = 0; s < n; s++) e(i, s) -= tmp * e(k, s);
		}
	}
	Matrix x(n,n,0);
	for (int k = 0; k < n; k++) {
		for (int i = n - 1; i >= 0; i--) {
			x(i, k) = e(i, k);
			double sum = 0 ;
			for (int j = n - 1; j > i; j--) sum += a(i, j) * x(j, k);
			x(i, k) = (x(i, k) - sum) / a(i, i);
		}
	}
	return x;
}




double u_a(double x, double y) {
	return sin(PI*x)*cos(PI*y);
//return x*x + y*y;
  //return exp(2*x)*sin(2*y);
}
double f(double x, double y) {
	return 2*PI*PI*sin(PI*x)*cos(PI*y);
//  return -4;
  //return 0;
}
double west(double y) {
	return u_a(0.0, y);
}
double east(double y) {
	return u_a(1.0, y);
}
double south(double x) {
	return u_a(x, 0.0);
}
double nord(double x) {
	return u_a(x, 1.0);
}

vector<double> mul(Matrix a, vector<double> b) {
	vector<double> res;
	unsigned n = a.getRows();
	Matrix res1(n,n,0);
	Matrix buf(n,n,0);
	for (unsigned i = 0; i < n; i++) {
			buf(i,0) = b[i];
	}
	res1 = a * buf;
	for (unsigned i = 0; i < n; i++) {
		res.push_back(res1(i,0));
	}
	return res;
}

Matrix my_sweep(Matrix c, Matrix f) {
	unsigned n = c.getRows();
	vector<Matrix> alpha;
	Matrix x(n+1,n+1,0);
	Matrix bet(n+1,n+1,0);
	alpha.push_back(inver(c));
	bet[0] = mul(alpha[0], f[0]);
	for (unsigned i = 0; i < n-1; i++) {
		alpha.push_back(inver(c - alpha[i]));
		bet[i+1] = mul(alpha[i+1], bet[i] + f[i+1]);
	}
	x[n-1] = bet[n-1];
	for (int i = n - 2; i >= 0; i--)
		x[i] = mul(alpha[i],x[i+1]+bet[i]);
	//x.print();
	return x;
}


double approx_m (unsigned n, Matrix u, Matrix f, double h2) {
	double a = 0.0;
	for (unsigned i = 1; i < n; i++)
		for (unsigned j = 1; j < n; j++) {
			double b = (u(i+1,j) - u(i,j)) / h2;
			b+= (u(i-1,j) - u(i,j)) / h2;
			b+= (u(i,j+1) - u(i,j)) / h2;
			b+= (u(i,j-1) - u(i,j)) / h2;
			b+= f(i,j);
			b = abs(b);
			a = max(a,b);
		}
	return a;
}
