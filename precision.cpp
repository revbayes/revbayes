//#pragma GCC target("fpmath=sse,sse2") // Turns off excess precision
// #pragma GCC target("fpmath=387") // Turns on excess precision
 
#include <iostream>
#include <cstdlib>
#include <cfloat>
using namespace std;
 
int main() {
  cout << "This is compiled in mode "<< FLT_EVAL_METHOD << '\n';
  cout << "0 means no excess precision.\n";
  cout << "2 means there is excess precision.\n\n";
  
  cout << "The following test detects excess precision\n";
  cout << "0 if no excess precision, or 8e-17 if there is excess precision.\n";
  double a = atof("1.2345678");
  double b = a*a;
  cout << b - 1.52415765279683990130 << '\n';
  return 0;
}
