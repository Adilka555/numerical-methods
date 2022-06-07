#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <cstring>
using namespace std;
int main() {
  int tap;
  cout << "select the method you want to use :" << endl;
  cout << "tap 1 = прогонка" << endl;
  cout << "tap 2 = Jakobi" << endl;
  cout << "tap 3 = Seidel" << endl;
  cout << "tap 4 = relax" << endl;
  cout << "tap 5 = Grafient" << endl;
  cout << "tap : ";
  int i, j;
  char stop[4];
  while (strcmp(stop, "stop") != 0){
    cin >> tap;
    if (tap == 1) {
      j = system("g++ -Wall -I ../include/ -o ../bin/first_m ../source/first_m.cpp ../source/glob_foos.cpp ../source/matrix.cpp");
      i = system("../bin/first_m");
      cout << "sys status : " << i << endl;
    }
    else if (tap == 2) {
      j = system("g++ -Wall -I ../include/ -o ../bin/jakob ../source/jakob.cpp ../source/glob_foos.cpp ../source/matrix.cpp");
      i = system("../bin/jakob");
      cout << "sys status : " << i << endl;
    }
    else if (tap == 3) {
      j = system("g++ -Wall -I ../include/ -o ../bin/seidel ../source/seidel.cpp ../source/glob_foos.cpp ../source/matrix.cpp");
      i = system("../bin/seidel");
      cout << "sys status : " << i << endl;
    }
    else if (tap == 4) {
      j = system("g++ -Wall -I ../include/ -o ../bin/relax ../source/relax.cpp ../source/glob_foos.cpp ../source/matrix.cpp");
      i = system("../bin/relax");
      cout << "sys status : " << i << endl;
    }
    else {
      j = system("g++ -Wall -I ../include/ -o ../bin/grad ../source/grad.cpp ../source/glob_foos.cpp ../source/matrix.cpp");
      i = system("../bin/grad");
      cout << "sys status : " << i << endl;
    }
    cout << "stop?" << endl;
    cin >> stop;
  }
  return 0;
}
