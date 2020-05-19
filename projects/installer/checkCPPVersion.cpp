#include<iostream>

int main() {
  if(__cplusplus >= 201103L) {
    std::cout << "1" << std::endl;
  } else {
    std::cout << "0" << std::endl;
  }
}
