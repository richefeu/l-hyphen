#include <iostream>
#include "Lhyphen.hpp"

// This is the command line interface (CLI) for using l-hyphen
int main(int argc, char const *argv[]) {
  INIT_TIMERS();

  if (argc != 2) {
    std::cout << "usage: " << argv[0] << " <input_file>\n";
    return 0;
  }

  try {
    Lhyphen S;

    // All the inputs for a simulation is read from a file
    S.loadCONF(argv[1]);
    S.head();
    S.integrate();
  } catch (const std::exception& e) {
    std::cerr << "An error occurred: " << e.what() << std::endl;
    return 1;
  }

  PRINT_TIMERS("mpmbox");

  return 0;
}




