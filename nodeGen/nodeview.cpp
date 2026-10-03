// nodeview : visualise un nodeFile (x y idCellule) sous forme de SVG.
//
// Les parois sont colorées selon qu'elles seront collées ou non par Lhyphen (un noeud d'une autre cellule à
// moins de barWidth + distGlue), ce qui fait apparaître les bords libres et les pré-fissures, et les barres
// courtes (raideur de flexion kr/l² élevée) sont mises en évidence.
//
// Usage : nodeview nodefile.txt [options]
//   -w <barWidth>          largeur des barres (défaut : distance minimale entre noeuds de cellules différentes)
//   -g <distGlue>          tolérance de collage (défaut 2e-7)
//   -r <ratio>             seuil des barres courtes, en fraction de la longueur moyenne (défaut 0.3)
//   -box <x0> <x1> <y0> <y1>  ne dessiner que cette fenêtre (zoom, ex. pointe de fissure)
//   -nodes                 dessiner les noeuds
//   -bare                  sans titre, statistiques ni légende (pour des figures)
//   -px <largeur>          largeur de l'image en pixels (défaut 1400)
//   -o <fichier.svg>       fichier de sortie (défaut : <nodefile>.svg)

#include "common.hpp"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage : nodeview nodefile.txt [-w barWidth] [-g distGlue] [-r ratio] [-box x0 x1 y0 y1] "
                 "[-nodes] [-bare] [-px largeur] [-o out.svg]\n";
    return 1;
  }
  std::string in = argv[1];
  std::string out = in.substr(0, in.find_last_of('.')) + ".svg";
  ViewOptions opt;
  for (int i = 2; i < argc; i++) {
    std::string a = argv[i];
    auto need = [&](int n) {
      if (i + n >= argc) {
        std::cerr << "Option " << a << " : valeur manquante\n";
        std::exit(1);
      }
    };
    if (a == "-w") {
      need(1);
      opt.barWidth = std::stod(argv[++i]);
    } else if (a == "-g") {
      need(1);
      opt.distGlue = std::stod(argv[++i]);
    } else if (a == "-r") {
      need(1);
      opt.shortRatio = std::stod(argv[++i]);
    } else if (a == "-box") {
      need(4);
      opt.useBox = true;
      opt.bx0 = std::stod(argv[++i]);
      opt.bx1 = std::stod(argv[++i]);
      opt.by0 = std::stod(argv[++i]);
      opt.by1 = std::stod(argv[++i]);
    } else if (a == "-bare") {
      opt.bare = true;
    } else if (a == "-nodes") {
      opt.showNodes = true;
    } else if (a == "-px") {
      need(1);
      opt.widthPx = std::stod(argv[++i]);
    } else if (a == "-o") {
      need(1);
      out = argv[++i];
    } else {
      std::cerr << "Option inconnue : " << a << "\n";
      return 1;
    }
  }

  std::vector<Polygon> cells = readNodeFile(in);
  if (cells.empty()) {
    return 1;
  }
  opt.title = "nodeview : " + in;
  MeshStats st = analyseAndWriteSVG(cells, opt, out);

  std::cout << std::setprecision(4);
  std::cout << in << "\n";
  std::cout << "  cellules : " << st.nCells << ", noeuds : " << st.nNodes << "\n";
  std::cout << "  barres   : l_moy = " << st.lMean << ", l_min = " << st.lMin << " (l_min/l_moy = " << st.lMin / st.lMean
            << "), " << st.nShort << " barres < " << opt.shortRatio << " l_moy\n";
  std::cout << "  côtés    : " << sidesHistogram(st) << "\n";
  std::cout << "  parois   : " << st.nGlued << " collées, " << st.nFree << " libres (barWidth = " << st.barWidth
            << ", distGlue = " << opt.distGlue << ")\n";
  std::cout << "  -> " << out << "\n";
  return 0;
}
