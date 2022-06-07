#include "glob_foos.h"
#include <string>
#include <iomanip>
#include <tuple>
const double PI = 3.141592;

double approx_j (unsigned n, Matrix u, Matrix f, double h2) {
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

double m_norma(Matrix m) {
	double norma = 0.0;
	for(unsigned i = 0; i < m.getRows(); i++)
	  {
	     double temp=0.0;
	     for(unsigned j = 0; j < m.getRows(); j++)
	        temp+=abs(m(j,i));
	     if(temp > norma)
	         norma = temp;
	  }
	return norma;
}

auto Jacob (Matrix u, Matrix F, unsigned h2, Matrix &res, Matrix au) {
	unsigned n = u.getRows();
	double t;
	double c, d, a = 0.0, b = 0.0;
		for (unsigned i = 1; i < n-1; i++) {
			for (unsigned j = 1; j < n-1; j++) {
				res(i,j) = u(i-1,j);
				res(i,j) += u(i+1,j);
				res(i,j) += u(i,j-1);
				res(i,j) += u(i,j+1);
				res(i,j) += h2 * f(i,j);
				res(i,j) /= 4;
				t = res(i,j);
				c = abs(t - u(i,j));
				d = abs(t - au(i,j));
				a = max(c,a);
				b = max(d,b);
			}
		}
	auto res1 = make_tuple(a,b);
return res1;
}

int main() {
	unsigned n = 65;
	unsigned np1 = n + 1;
	double h = 1.0/n;
	double h2 = h*h;
	vector<double> x = linspace(0.0,1.0,np1);
	vector<double> y = linspace(0.0,1.0,np1);
	Matrix F(n,n,0);
	Matrix U(np1,np1,0);
	Matrix V(np1,np1,0);
	for (unsigned i = 1; i < n; i++)
		for (unsigned j = 1; j < n; j++)
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
double ea_apr = approx_j(n, au, F, h2);
double e0_apr = approx_j(n, U, F, h2);
cout << "e0_apr = " << e0_apr << endl;
double re_apr = ea_apr / e0_apr;
//cout.precision(5);
cout << "ea_apr : " << ea_apr << " re_apr : " << re_apr << endl;
int p = 2;
int K = ((1.25 * (2 * p * log(10) * n * n/(PI * PI + 0.5)) / 100) * 100);
cout << "K : " << K << endl;

double rs = cos(PI*h);
double mr = rs / (1.0 - rs);
Matrix u = U, v = V;
double p_euk = 1.0;

cout << " k " << setw(13) << " F-AUk" << setw(15) << "F-AUk/F-AU0" << setw(15) << "Uk-aU" << setw(15) << "Uk-aU/U0-aU" << setw(10) << "Uk-Uk-1 " << setw(15) << "pogr" << setw(15) << "rs_exp" << endl;
for(int k = 1; k <= K ; k++) {
auto res1 = Jacob(u, F, h2, v, au);
	if(k % 100 == 0 || k == K) {
		double c2 = approx_j(n, v, F, h2);
		double c3 = c2 / e0_apr;
		double c4 = get<1>(res1);
		double c5 = c4 / ea0;
		double c6 = get<0>(res1);
		double c7 = mr * c6;
		double rsk = get<0>(res1)/p_euk;
		double c8 = rsk;
		cout << k;
		cout << setw(15) << c2;
		cout << setw(15) << c3;
		cout << setw(15) << c4;
		cout << setw(15) << c5;
		cout << setw(15) << c6;
		cout << setw(15) << c7;
		cout << setw(15) << c8 << endl;
	}
	u = v;
	p_euk = get<0>(res1);
}
	ofstream fout;
	fout.open("../output/Jacob.txt");
	//cout << "sizes : " << x.size() << " " << y.size() << " " << V.getRows() << endl;
	for (unsigned i = 0; i < x.size() -1; i++) {
		for (unsigned j = 0; j < x.size()-1; j++) {
			fout << x[i] << " ," << y[j] << " ," << u(i,j) << '\n';
		}
	}
	return 0;
}
