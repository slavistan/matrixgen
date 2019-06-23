#include <iostream>

auto lambda = [](int a) { return a; };

int main()
{

  std::cout << "Hello World!" << std::endl;
  std::cout << "a = " << lambda(44) << std::endl;
}
