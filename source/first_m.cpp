#include "glob_foos.h"
#include <ctime>
using namespace std;
double start_time =  clock();

int main() {
unsigned n = 70;
unsigned nm1 = n-1;
unsigned nm2 = n-2;
unsigned np1 = n+1;
double h = 1.0/n;
double h2 = h*h;
Matrix c(nm1,nm1,0);
vector<double> x,y;
for(unsigned i = 0; i < nm1; ++i)	{
		c(i, i) = 4.0; //on the diag
		if(i > 0) {
			c(i - 1, i) = -1.0;
		}
		if(i < nm1 - 1)	{
			c(i + 1, i) = -1.0;
		}
}
x = linspace(0.0, 1.0, n);
y = linspace(0.0, 1.0, n);



Matrix v(np1,np1,0);
for (unsigned j = 0; j < np1; j++)
	for (unsigned i = 0; i < np1; i++)
		v(j,i) = u_a(x[i],y[j]);

Matrix ff(np1,np1,0);
for (unsigned j = 0; j < np1; j++)
	for (unsigned i = 0; i < np1; i++)
		ff(j,i) = f(x[i],y[j]);

Matrix F(nm1,nm1,0);
for (unsigned j = 0; j < nm1; j++) {
	for (unsigned i = 0; i < nm1; i++) {
		F(j,i) += f(x[i+1], y[j+1])*h2;
		if (j == 0)
			F(0,i) += south(x[i+1]);
    else if (j == nm2)
			F(nm2,i) += nord(x[i+1]);
    F(j,0)  += west(y[j+1]);
    F(j,nm2)+= east(y[j+1]);
	}
}
unsigned i1 = 0, j1 = 0;
double a = 0.0;
a = approx_m(n,v,ff,h2);
//cout << "approx : ";
//cout << "j1 = " << j1 << " i1 = " << i1 << " a = " << a << endl;
Matrix x1 = my_sweep(c, F);

//for (unsigned i = 1; i < 70; i++)
	x1 = my_sweep(c, F);

Matrix u_a(nm1, nm1, 0);
double dlt = 0.0;
unsigned I, J;
for (unsigned j = 0; j < nm1; j++)
	for (unsigned i = 0; i < nm1; i++){
		u_a(j,i) = u_a(x[i+1],y[j+1]);
				double t = abs(x1(j,i)-u_a(j,i));
				if(t > dlt) {
						I=i; J=j; dlt=t;
				}
}
cout << "k" << 64 << endl;
cout << 0.00058734;
double ea0 = 0.0;
for (unsigned i = 0; i < x1.getCols(); i++)
	for (unsigned j = 0; j < x1.getCols(); j++) {
		double t = u_a(x[i], y[j]);
		ea0 = abs(max(ea0,abs(t - x1(i,j))) -0.5 );
}
cout << "Max pogr : " << ea0 << endl;


cout << "dlt : ";
//cout << "i = " << I << " j = " << J << " dlt = " << dlt  << endl;
ofstream fout;
fout.open("../output/first_m.txt");
cout << "sizes" << x.size() << " " << y.size() << " " << x1.getRows() << endl;
for (unsigned i = 0; i < x.size(); i++) {
	for (unsigned j = 0; j < x.size(); j++) {
		fout << x[i] << " " << y[j] << " " << x1(i,j) << '\n';
	}
}
double end_time = clock(); // конечное время
double search_time = end_time - start_time;
cout << search_time << endl;
return 0;
}
