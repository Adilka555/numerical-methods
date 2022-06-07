#include "glob_foos.h"
#include <ctime>


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

Matrix foo(Matrix v) {
  unsigned n = v.getCols();
  Matrix res(n, n, 0);
  unsigned j = 0;
  for (unsigned i = 0; i < n; i++) {
    res(j,i) = 4*v(j,i) - v(j+1,i);
      if ( i > 0)
        res(j,i) += -v(j,i-1);
      if (i < n-1)
        res(j,i) += -v(j,i+1);
    }
  for (unsigned j = 1; j < n - 1; j++)
    for (unsigned i = 0; i < n; i++)  {
      res(j,i) = 4*v(j,i) - v(j-1,i) - v(j+1,i);
      if (i > 0)
        res(j,i) += -v(j,i-1);
      if (i < n-1)
        res(j,i) += -v(j,i+1);
    }
  j = n - 1;
  for (unsigned i = 0; i < n; i++) {
    res(j,i) = 4*v(j,i) - v(j-1,i);
    if (i > 0)
      res(j,i) += -v(j,i-1);
    if (i < n-1)
      res(j,i) += -v(j,i+1);
  }
  return res;
}

double m_norma(Matrix m) {
double max; // здесь в конце будет храниться норма матрицы
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

double m_scalar(Matrix a, Matrix b) {
  unsigned n = a.getCols();
  double buf = 0.0;
  for (unsigned i = 0; i < n; i++)
    for (unsigned j = 0; j < n; j++) {
      buf += (a(i,j) * b(i,j));
    }
  return buf;
}

Matrix grad(Matrix b, unsigned n, unsigned &k) {
  Matrix x(n,n,0);
  Matrix p[n+2];
  Matrix r[n+2];
  p[1] = b;
  r[0] = b;
  Matrix A(n,n,0);
  k = 0;
  while (m_norma(r[k]) > 1e-5) {
    k++;
    Matrix z = foo(p[k]);
    double v = m_scalar(r[k-1],r[k-1])/m_scalar(p[k],z);
    Matrix buf1 = p[k]*v;
    x = x + buf1;
    Matrix buf2 = (z*v);
    r[k] = r[k-1]- buf2;
    double mu = m_scalar(r[k],r[k])/m_scalar(r[k-1],r[k-1]);
    Matrix buf3 = p[k]*mu;
    p[k+1] = r[k] + buf3;
    if (k == 43) { cout << m_norma(r[k]) << endl; }
  }
  return x;
}

double start_time =  clock();

int main() {
  unsigned n = 60;
  unsigned nm1 = n-1;
  unsigned nm2 = n-2;
  unsigned np1 = n+1;
  double h = 1.0/n;
  double h2 = h*h;
  vector<double> x,y;
  x = linspace(0.0, 1.0, np1);
  y = linspace(0.0, 1.0, np1);

  Matrix v(np1,np1,0);
  for (unsigned j = 0; j < np1; j++)
  	for (unsigned i = 0; i < np1; i++)
  		v(j,i) = u_a(x[i],y[j]);

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

  unsigned k;
  Matrix x1 = grad(F, nm1, k);
  cout << "K << " << k << " iteration" <<endl;
	double ea0 = 0.0;
	for (unsigned i = 0; i < x1.getCols(); i++)
		for (unsigned j = 0; j < x1.getCols(); j++) {
			double t = u_a(x[i], y[j]);
			ea0 = abs(max(ea0,abs(t - x1(i,j))) -0.5 );
}
cout << "Max pogr : " << ea0 << endl;
  ofstream fout;
  fout.open("../output/grad.txt");
  for (unsigned i = 0; i < x.size() - 2; i++) {
  	for (unsigned j = 0; j < x.size() - 2; j++) {
  		fout << x[i] << " " << y[j] << " " << x1(i,j) << '\n';
  	}
  }
  double end_time = clock(); // конечное время
  double search_time = end_time - start_time;
  cout << search_time << endl;
  return 0;
}
