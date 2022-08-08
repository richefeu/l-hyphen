#include "Lhyphen.hpp"

// This is the command line interface (CLI) for using l-hyphen
int main(int argc, char const *argv[]) {

  if (argc != 2) {
    std::cout << "usage: " << argv[0] << "<input_file>\n";
    return 0;
  }

  Lhyphen S;

  // All the inputs for a simulation is read from a file
  S.loadCONF(argv[1]);
  S.head();
  S.run();

  return 0;
}
