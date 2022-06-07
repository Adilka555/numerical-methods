#include "glob_foos.h"
#include <string>
#include <iomanip>
#include <tuple>

const double PI = 3.14159;


auto seidel (Matrix u, Matrix F, double h2, Matrix au, Matrix &new_u) {
	unsigned n = u.getRows();
	double c, d, a = 0.0, b = 0.0, t;
	new_u = u;
		for (unsigned i = 1; i < n-1; i++) {
			for (unsigned j = 1; j < n-1; j++) {
				new_u(i,j) = (new_u(i-1,j) + u(i+1,j) + new_u(i,j-1) + u(i,j+1) + h2*F(i,j))/4;
				t = new_u(i,j);
				c = abs(t - u(i,j));
				d = abs(t - au(i,j));
				a = max(c,a);
				b = max(d,b);
			}
		}
	auto res1 = make_tuple(a,b);
return res1;
}

double m_norma(Matrix m) {
double max;
max = 0.0;
  for (unsigned i = 0; i < m.getRows(); i++) {
    double s = 0.0;
    for (unsigned j = 0; j < m.getCols(); j++) {
      	s += abs(m(i,j));
    }
    if (s > max)
      max = s;
  }
  return max;
}

double approx (unsigned n, Matrix u, Matrix f, double h2) {
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

int main() {
	unsigned n = 16;
	unsigned np1 = n + 1;
	double h = 1.0/n;
	double h2 = h*h;
	vector<double> x = linspace(0.0,1.0,np1);
	vector<double> y = linspace(0.0,1.0,np1);
	Matrix F(n,n,0);
	Matrix U(np1,np1,0);
	Matrix V(np1,np1,0);
	for (unsigned i = 1; i < n; i++)
		for (unsigned j = 0; j < n; j++)
			F(i,j) = f(x[i], y[j]);

	for (unsigned i = 0; i < n; i++) {
		V(i,0) = U(i,0) = south(x[i]);
		V(i,n) = U(i,n) = nord(x[i]);
		V(0,i) = U(0,i) = west(y[i]);
		V(n,i) = U(n,i) = east(y[i]);
	}
	Matrix au(np1,np1,0);
	double ea0 = 0.0;
	for (unsigned i = 0; i < np1; i++)
		for (unsigned j = 0; j < np1; j++) {
			double t = u_a(x[i], y[j]);
			au(i,j) = t;
			ea0 = max(ea0,abs(t - U(i,j)));
}
cout << "au - u0 : " << ea0 << endl;
double ea_apr = approx(n, au, F, h2);
double e0_apr = approx(n, U, F, h2);
cout << "E0 APR = " << e0_apr << endl;
double re_apr = ea_apr / e0_apr;
cout.precision(5);
cout << "ea_apr : " << ea_apr << " re_apr : " << re_apr << endl;
	int K = 3 * log(10);
	K/= (PI * PI * h * h);
	int K1 = ceil(K);
	double rs = cos(PI*h);
	double mr = rs / (1.0 - rs);
	Matrix u = U, v = V;
	double p_euk = 1.0;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << " k" << setw(13) << " F-AUk" << setw(15) << "F-AUk/F-AU0" << setw(15) << "Uk-aU" << setw(15) << "Uk-aU/U0-aU" << setw(10) << "Uk-Uk-1" << setw(15) << "pogr" << setw(15) << "rs_exp" << endl;
for(int k = 1; k <= K1; k++) {
auto res1 = seidel(u, F, h2, au, v);
	if(k % 100 == 0 || k == K1) {
			double c2 = approx(n, v, F, h2);
			double c3 = c2 / e0_apr;
			double c4 = get<1>(res1);
			double c5 = c4 / ea0;
			double c6 = get<0>(res1);
			double c7 = mr * c6;
			double rsk = get<0>(res1)/p_euk;
			double c8 = rsk;
			cout << k;
			cout << setw(15)  << c2;
			cout << setw(15)  << c3;
			cout << setw(15)  << c4;
			cout << setw(15)  << c5;
			cout << setw(15)  << c6;
			cout << setw(15)  << c7;
			cout << setw(15)  << c8 << endl;
	}
	u = v;
	p_euk = get<0>(res1);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	ofstream fout;
	fout.open("../output/seidel3.txt");
	for (unsigned i = 0; i < x.size() ; i++) {
		for (unsigned j = 0; j < x.size() ; j++) {
			fout << x[i] << " " << y[j] << " " << u(i,j) << '\n';
		}
	}
	return 0;
}
