#include "cellPrepro.hpp"

int main(int argc, char *argv[]) {

  if (argc == 2) {
    if (fileTool::fileExists(argv[1])) {
      CellPrepro P;
      P.prepro_commands(argv[1]);
    } else {
      std::cout << argv[0] << " -> pas trouvé\n";
    }
  } else {
    std::cout << "Usage : " << argv[0] << " command_filename\n";
  }

  return 0;
}
