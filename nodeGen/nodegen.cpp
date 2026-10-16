// nodegen : génère un nodeFile (x y idCellule) pour Lhyphen::readNodeFile.
//
// Pavage de Voronoï d'un rectangle (germes tirés avec une distance minimale puis régularisés par Lloyd),
// traitement des arêtes trop courtes par transitions T1 ou étirement (jamais par fusion, qui créerait des
// sommets à 4 branches), éclatement des sommets à plus de 3 branches, puis décalage des parois de barWidth/2
// vers l'intérieur (même construction que cellPrepro) : deux cellules voisines sont exactement à barWidth
// l'une de l'autre. Seules les extrémités de la pré-fissure gardent un sommet à 4 branches.
//
// Pré-fissure nette : les germes proches de la ligne de fissure sont placés par paires symétriques par
// rapport à cette ligne, de sorte que la fissure suit exactement des arêtes du pavage. Les parois le long
// de la fissure sont décalées de opening/2 en plus : l'écart barWidth + opening empêche le collage.
//
// Usage : nodegen params.txt      (voir README.md et example/params.txt)

#include "common.hpp"

#include <filesystem>
#include <random>
#include <set>

// ======================================================================================================
// Paramètres
// ======================================================================================================

struct Params {
  double xmin{0.0}, ymin{0.0}, Lx{1e-2}, Ly{1e-2};
  double cellSize{3.5e-4}; // diamètre équivalent visé des cellules
  double barWidth{2e-6};
  int lloydIterations{10};
  double minEdgeFactor{0.2}; // longueur minimale des arêtes, en fraction de lRef = sqrt(Lx * Ly / nombre de cellules)
  double t1Margin{1.2};     // longueur de l'arête créée par une transition T1, en multiple du seuil
  unsigned seed{1};
  std::string lattice{"hex"};    // random : tirage aléatoire ; hex : réseau hexagonal perturbé
  double disorder{0.3};          // réseau hex : amplitude de la perturbation des germes (fraction du pas)
  int minSides{5};               // une transition T1 ne descend pas une cellule sous minSides côtés
  bool hasCrack{false};
  V c0, c1;               // pré-fissure : segment [c0, c1]
  double tipInterface{1.0}; // longueur d'interface collée prolongeant la fissure au-delà de la pointe (en cellules)
  double opening{2e-6};   // ouverture supplémentaire (écart entre parois = barWidth + opening)
  double distGlue{2e-7};  // tolérance de collage utilisée dans l'input (pour le contrôle et le SVG)
  std::string output{"nodefile.txt"};
  std::string svg{"nodefile.svg"};
  std::string input{"nodegen_input.txt"}; // lignes d'input l-hyphen correspondant au maillage
  // Valeurs facultatives reportées dans les lignes d'input (sinon lignes commentées avec des repères <...>)
  bool hasCellProps{false};
  double Kn{0}, Kr{0}, MzMax{0}, pInt{0};
  int gripRows{1};           // nombre de rangées de noeuds attrapées par chaque mors (défaut : la rangée du bord)
  bool hasGripRows{false};
  double gripHeight{-1.0};   // hauteur explicite des mors (prioritaire sur gripRows si gripRows n'est pas donné)
  double pullVelocity{3e-4}; // vitesse imposée aux mors (bas : -v, haut : +v)
  bool hasVelocity{false};   // sinon, en mode inputTemplate, modes et vitesses sont repris du modèle
  std::string inputTemplate;           // input l-hyphen existant servant de modèle
  std::string inputDeck{"input.txt"};  // input complet produit à partir du modèle
  bool hasGlueProps{false};
  double knCoh{0}, ktCoh{0}, Gc{0};
  double nodeMass{-1.0};     // pour estimer le pas de temps critique de flexion
};

static Params readParams(const std::string &name) {
  Params p;
  std::ifstream file(name);
  if (!file) {
    std::cerr << "Impossible d'ouvrir " << name << "\n";
    std::exit(1);
  }
  std::string line;
  while (std::getline(file, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) {
      line = line.substr(0, hash);
    }
    std::istringstream is(line);
    std::string k;
    if (!(is >> k)) {
      continue;
    }
    if (k == "xmin") is >> p.xmin;
    else if (k == "ymin") is >> p.ymin;
    else if (k == "Lx") is >> p.Lx;
    else if (k == "Ly") is >> p.Ly;
    else if (k == "cellSize") is >> p.cellSize;
    else if (k == "barWidth") is >> p.barWidth;
    else if (k == "lloydIterations") is >> p.lloydIterations;
    else if (k == "minEdgeFactor") is >> p.minEdgeFactor;
    else if (k == "minEdgeRatio") { // ancien mot-clé, rapporté à la longueur moyenne des arêtes (~0.62 lRef)
      double r = 0.0;
      is >> r;
      std::cerr << "minEdgeRatio n'est plus accepté : le seuil est maintenant rapporté à lRef = sqrt(Lx * Ly / nombre "
                   "de cellules), environ 1.6 fois la longueur moyenne des arêtes. Remplacer par minEdgeFactor "
                << std::setprecision(2) << 0.62 * r << " (équivalent de minEdgeRatio " << r << ").\n";
      std::exit(1);
    }
    else if (k == "t1Margin") is >> p.t1Margin;
    else if (k == "seed") is >> p.seed;
    else if (k == "lattice") is >> p.lattice;
    else if (k == "disorder") is >> p.disorder;
    else if (k == "minSides") is >> p.minSides;
    else if (k == "crack") {
      is >> p.c0.x >> p.c0.y >> p.c1.x >> p.c1.y;
      p.hasCrack = true;
    } else if (k == "opening") is >> p.opening;
    else if (k == "tipInterface") is >> p.tipInterface;
    else if (k == "distGlue") is >> p.distGlue;
    else if (k == "output") is >> p.output;
    else if (k == "svg") is >> p.svg;
    else if (k == "input") is >> p.input;
    else if (k == "cellProperties") {
      is >> p.Kn >> p.Kr >> p.MzMax >> p.pInt;
      p.hasCellProps = true;
    } else if (k == "grips") { // forme historique : hauteur et vitesse
      is >> p.gripHeight >> p.pullVelocity;
      p.hasVelocity = true;
    } else if (k == "gripRows") {
      is >> p.gripRows;
      p.hasGripRows = true;
    } else if (k == "gripHeight") is >> p.gripHeight;
    else if (k == "pullVelocity") {
      is >> p.pullVelocity;
      p.hasVelocity = true;
    } else if (k == "inputTemplate") {
      is >> p.inputTemplate;
      std::string out;
      if (is >> out) {
        p.inputDeck = out;
      }
      is.clear();
    }
    else if (k == "glueProperties") {
      is >> p.knCoh >> p.ktCoh >> p.Gc;
      p.hasGlueProps = true;
    } else if (k == "nodeMass") is >> p.nodeMass;
    else {
      std::cerr << "Mot-clé inconnu : " << k << "\n";
      std::exit(1);
    }
    if (is.fail()) {
      std::cerr << "Valeur invalide pour " << k << "\n";
      std::exit(1);
    }
  }
  return p;
}

// ======================================================================================================
// Germes
// ======================================================================================================

struct SeedGrid {
  double x0, y0, h;
  long nx, ny;
  std::vector<std::vector<int>> b;
  SeedGrid(double X0, double Y0, double Lx, double Ly, double H) : x0(X0), y0(Y0), h(H) {
    nx = std::max(1L, (long)std::ceil(Lx / h));
    ny = std::max(1L, (long)std::ceil(Ly / h));
    b.resize(nx * ny);
  }
  long ix(double x) const { return std::clamp((long)std::floor((x - x0) / h), 0L, nx - 1); }
  long iy(double y) const { return std::clamp((long)std::floor((y - y0) / h), 0L, ny - 1); }
  void add(int i, const V &p) { b[iy(p.y) * nx + ix(p.x)].push_back(i); }
};

// Repère de la fissure : origine c0, t le long de la fissure, n normale
struct CrackFrame {
  V o, t, n;
  double L;
  double along(const V &p) const { return dot(p - o, t); }
  double normal(const V &p) const { return dot(p - o, n); }
};

static bool inside(const Params &P, const V &q) {
  return q.x > P.xmin && q.x < P.xmin + P.Lx && q.y > P.ymin && q.y < P.ymin + P.Ly;
}

// Rend les germes symétriques par rapport à la ligne de fissure dans une bande autour de celle-ci :
// les germes du côté n < 0 sont supprimés, ceux du côté n >= 0 sont dupliqués par symétrie.
static std::vector<V> symmetrize(const std::vector<V> &seeds, const Params &P, const CrackFrame &F, double band,
                                 double ext, double rmin) {
  if (!P.hasCrack) {
    return seeds;
  }
  std::vector<V> out;
  for (V s : seeds) {
    double a = F.along(s), nn = F.normal(s);
    // la symétrie est prolongée au-delà de la pointe : la ligne de fissure s'y continue par une interface
    // collée entre deux cellules (la fissure ne bute pas sur une cellule)
    bool inBand = a >= -ext && a <= F.L + ext && std::fabs(nn) < band;
    if (!inBand) {
      out.push_back(s);
      continue;
    }
    if (nn < 0.0) {
      continue;
    }
    nn = std::max(nn, 0.5 * rmin); // une paire trop serrée donnerait des cellules très aplaties
    V up = F.o + F.t * a + F.n * nn;
    V down = F.o + F.t * a - F.n * nn;
    if (inside(P, up) && inside(P, down)) {
      out.push_back(up);
      out.push_back(down);
    }
  }
  return out;
}

// ======================================================================================================
// Voronoï par découpe de demi-plans (cellules convexes, découpées au rectangle)
// ======================================================================================================

// Garde la partie de P plus proche de p que de q
static Polygon clipHalfPlane(const Polygon &poly, const V &p, const V &q) {
  V m = (p + q) * 0.5, d = q - p;
  Polygon out;
  size_t n = poly.size();
  for (size_t i = 0; i < n; i++) {
    const V &A = poly[i];
    const V &B = poly[(i + 1) % n];
    double fa = dot(A - m, d), fb = dot(B - m, d);
    if (fa <= 0.0) {
      out.push_back(A);
    }
    if ((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0)) {
      out.push_back(A + (B - A) * (fa / (fa - fb)));
    }
  }
  return out;
}

static std::vector<Polygon> voronoi(const std::vector<V> &seeds, const Params &P, double h) {
  SeedGrid g(P.xmin, P.ymin, P.Lx, P.Ly, h);
  for (size_t i = 0; i < seeds.size(); i++) {
    g.add((int)i, seeds[i]);
  }
  Polygon rect = {V(P.xmin, P.ymin), V(P.xmin + P.Lx, P.ymin), V(P.xmin + P.Lx, P.ymin + P.Ly),
                  V(P.xmin, P.ymin + P.Ly)};
  std::vector<Polygon> cells(seeds.size());
  for (size_t i = 0; i < seeds.size(); i++) {
    const V &p = seeds[i];
    Polygon poly = rect;
    long cx = g.ix(p.x), cy = g.iy(p.y);
    for (long k = 0;; k++) {
      for (long jy = cy - k; jy <= cy + k; jy++) {
        for (long jx = cx - k; jx <= cx + k; jx++) {
          if (std::max(std::labs(jx - cx), std::labs(jy - cy)) != k || jx < 0 || jy < 0 || jx >= g.nx ||
              jy >= g.ny) {
            continue;
          }
          for (int j : g.b[jy * g.nx + jx]) {
            if (j != (int)i) {
              poly = clipHalfPlane(poly, p, seeds[j]);
            }
          }
        }
      }
      double R = 0.0;
      for (auto &q : poly) {
        R = std::max(R, norm(q - p));
      }
      if ((double)k * h >= 2.0 * R || (k > g.nx && k > g.ny)) {
        break;
      }
    }
    cells[i] = poly;
  }
  return cells;
}

// ======================================================================================================
// Topologie : sommets partagés, types, arêtes courtes et sommets à plus de 3 branches
// ======================================================================================================

enum : int { LEFT = 1, RIGHT = 2, BOTTOM = 4, TOP = 8, CRACK = 16 };

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage : nodegen params.txt\n";
    return 1;
  }
  Params P = readParams(argv[1]);

  const double cellArea = M_PI * P.cellSize * P.cellSize / 4.0;
  const double spacing = std::sqrt(cellArea);
  const size_t N = (size_t)std::llround(P.Lx * P.Ly / cellArea);
  const double rmin = 0.7 * spacing;
  const double band = 2.0 * spacing;
  // la symétrie s'étend de 'ext' au-delà des extrémités de la fissure : la ligne s'y prolonge par une interface
  // collée d'environ tipInterface cellules avant de rejoindre le réseau régulier
  const double ext = (std::max(0.0, P.tipInterface) + 0.5) * spacing;

  CrackFrame F;
  if (P.hasCrack) {
    F.o = P.c0;
    F.L = norm(P.c1 - P.c0);
    F.t = unit(P.c1 - P.c0);
    F.n = V(-F.t.y, F.t.x);
  }

  // 1. Germes
  std::mt19937_64 rng(P.seed);
  std::uniform_real_distribution<double> ux(P.xmin, P.xmin + P.Lx), uy(P.ymin, P.ymin + P.Ly);
  std::vector<V> seeds;
  if (P.lattice == "hex") {
    // Réseau hexagonal (cellules de même aire que cellSize), orienté selon la fissure : les rangées sont
    // parallèles à la fissure, symétriques par rapport à elle, et un sommet tombe sur la pointe.
    const double a = std::sqrt(2.0 * cellArea / std::sqrt(3.0));
    const double dy = 0.5 * std::sqrt(3.0) * a;
    V o = P.hasCrack ? F.o : V(P.xmin, P.ymin);
    V t = P.hasCrack ? F.t : V(1.0, 0.0);
    V n = P.hasCrack ? F.n : V(0.0, 1.0);
    double L = P.hasCrack ? F.L : 0.0;
    double amin = 1e300, amax = -1e300, nmin = 1e300, nmax = -1e300;
    for (V c : {V(P.xmin, P.ymin), V(P.xmin + P.Lx, P.ymin), V(P.xmin + P.Lx, P.ymin + P.Ly), V(P.xmin, P.ymin + P.Ly)}) {
      amin = std::min(amin, dot(c - o, t));
      amax = std::max(amax, dot(c - o, t));
      nmin = std::min(nmin, dot(c - o, n));
      nmax = std::max(nmax, dot(c - o, n));
    }
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (long k = (long)std::floor(nmin / dy) - 2; k <= (long)std::ceil(nmax / dy) + 2; k++) {
      double nk = ((double)k + 0.5) * dy;
      for (long i = (long)std::floor((amin - L) / a) - 2; i <= (long)std::ceil((amax - L) / a) + 2; i++) {
        // Le long de la fissure (et jusqu'à 'band' au-delà de ses extrémités), les rangées sont symétriques
        // par rapport à elle ; ailleurs le réseau hexagonal est régulier, pour ne pas prolonger la fissure par
        // une interface rectiligne qui serait un chemin de rupture privilégié.
        double a0 = L + ((double)i + 0.5) * a;
        bool mirrored = P.hasCrack && a0 >= -ext && a0 <= L + ext;
        long parity = mirrored ? ((k >= 0) ? k : -k - 1) % 2 : ((k % 2) + 2) % 2;
        double ak = a0 + (double)parity * 0.5 * a; // sommets du réseau sur la ligne en L + i*a
        double r = P.disorder * a * std::sqrt(u01(rng)), th = 2.0 * M_PI * u01(rng);
        V q = o + t * (ak + r * std::cos(th)) + n * (nk + r * std::sin(th));
        if (inside(P, q)) {
          seeds.push_back(q);
        }
      }
    }
  } else {
    // Tirage avec distance minimale rmin (puis complément aléatoire si besoin)
    SeedGrid g(P.xmin, P.ymin, P.Lx, P.Ly, rmin);
    size_t attempts = 0;
    while (seeds.size() < N && attempts < 200 * N) {
      attempts++;
      V q(ux(rng), uy(rng));
      bool ok = true;
      long cx = g.ix(q.x), cy = g.iy(q.y);
      for (long jy = cy - 1; jy <= cy + 1 && ok; jy++) {
        for (long jx = cx - 1; jx <= cx + 1 && ok; jx++) {
          if (jx < 0 || jy < 0 || jx >= g.nx || jy >= g.ny) {
            continue;
          }
          for (int j : g.b[jy * g.nx + jx]) {
            if (norm(seeds[j] - q) < rmin) {
              ok = false;
              break;
            }
          }
        }
      }
      if (ok) {
        g.add((int)seeds.size(), q);
        seeds.push_back(q);
      }
    }
    while (seeds.size() < N) {
      seeds.emplace_back(ux(rng), uy(rng));
    }
  }
  seeds = symmetrize(seeds, P, F, band, ext, rmin);

  // 2. Régularisation de Lloyd (la symétrie autour de la fissure est rétablie à chaque itération)
  std::vector<Polygon> vor;
  for (int it = 0; it < P.lloydIterations; it++) {
    vor = voronoi(seeds, P, spacing);
    for (size_t i = 0; i < seeds.size(); i++) {
      if (vor[i].size() >= 3) {
        seeds[i] = centroid(vor[i]);
      }
    }
    seeds = symmetrize(seeds, P, F, band, ext, rmin);
  }
  vor = voronoi(seeds, P, spacing);

  // 3. Sommets partagés (fusion à tolérance près des sommets calculés indépendamment par chaque cellule)
  const double tolV = 1e-8 * spacing;
  std::vector<V> verts;
  std::vector<std::vector<int>> cellV;
  {
    std::unordered_map<long long, std::vector<int>> hash;
    const double hb = 1e-6 * spacing;
    auto key = [](long long i, long long j) { return i * 1000003LL + j; };
    for (auto &poly : vor) {
      std::vector<int> ids;
      for (auto &q : poly) {
        long long i = (long long)std::floor(q.x / hb), j = (long long)std::floor(q.y / hb);
        int found = -1;
        for (long long di = -1; di <= 1 && found < 0; di++) {
          for (long long dj = -1; dj <= 1 && found < 0; dj++) {
            auto it = hash.find(key(i + di, j + dj));
            if (it == hash.end()) {
              continue;
            }
            for (int v : it->second) {
              if (norm(verts[v] - q) < tolV) {
                found = v;
                break;
              }
            }
          }
        }
        if (found < 0) {
          found = (int)verts.size();
          verts.push_back(q);
          hash[key(i, j)].push_back(found);
        }
        if (ids.empty() || ids.back() != found) {
          ids.push_back(found);
        }
      }
      while (ids.size() > 1 && ids.front() == ids.back()) {
        ids.pop_back();
      }
      cellV.push_back(ids);
    }
  }

  // Types des sommets (bords du domaine, ligne de fissure)
  const double tolLine = 1e-6 * spacing;
  std::vector<int> flags(verts.size(), 0);
  for (size_t v = 0; v < verts.size(); v++) {
    const V &q = verts[v];
    if (std::fabs(q.x - P.xmin) < tolLine) flags[v] |= LEFT;
    if (std::fabs(q.x - (P.xmin + P.Lx)) < tolLine) flags[v] |= RIGHT;
    if (std::fabs(q.y - P.ymin) < tolLine) flags[v] |= BOTTOM;
    if (std::fabs(q.y - (P.ymin + P.Ly)) < tolLine) flags[v] |= TOP;
    if (P.hasCrack && std::fabs(F.normal(q)) < tolLine && F.along(q) > -band - tolLine &&
        F.along(q) < F.L + band + tolLine) {
      flags[v] |= CRACK;
    }
  }

  // 4. Arêtes courtes et sommets à plus de 3 branches. Une jonction de plus de 3 parois est mécaniquement
  //    indésirable : les arêtes courtes ne sont donc jamais fusionnées (la fusion créerait un sommet à 4
  //    branches). Elles sont traitées par une transition T1 (l'arête est remplacée par une arête de longueur
  //    t1Margin * seuil, orientée transversalement : ses deux cellules cessent d'être voisines, les deux autres
  //    le deviennent) si ses deux cellules gardent au moins minSides côtés, sinon par étirement (étape 4bis).
  //    Les sommets à 4 branches ou plus (sommets quasi confondus, ligne de symétrie prolongeant la fissure
  //    au-delà de ses extrémités) sont éclatés en sommets à 3 branches ; sur la ligne prolongée, l'interface
  //    collée devient ainsi un zigzag. Seules les extrémités de la fissure gardent 4 branches.
  std::vector<V> pos = verts;
  std::vector<int> rflags = flags;
  std::vector<std::vector<int>> poly = cellV;
  // Longueur caractéristique : côté du carré de même aire qu'une cellule, lRef = sqrt(Lx * Ly / nombre de cellules).
  // Elle ne dépend que de la densité de cellules (la longueur moyenne des arêtes, environ 0.62 lRef pour des
  // hexagones, dépend aussi de la forme des cellules et donc du traitement des arêtes courtes).
  size_t nCellsTiling = 0;
  for (auto &ids : poly) {
    if (ids.size() >= 3) {
      nCellsTiling++;
    }
  }
  const double lRef = std::sqrt(P.Lx * P.Ly / (double)std::max<size_t>(1, nCellsTiling));
  const double lCrit = P.minEdgeFactor * lRef;
  // longueur des arêtes créées par T1 ou éclatement (sans seuil, minEdgeFactor 0 : 0.06 lRef, environ 0.1 fois la
  // longueur moyenne des arêtes)
  const double lT1 = std::max(1.0, P.t1Margin) * (lCrit > 0.0 ? lCrit : 0.06 * lRef);

  // Extrémités de la fissure : sommets de la ligne les plus proches de c0 et de c1. La fissure s'ouvre entre
  // eux ; leurs sommets (lèvres libres) ne sont pas des jonctions à 4 branches au sens mécanique.
  double startAlong = 0.0, tipAlong = 0.0;
  auto findCrackEnds = [&]() {
    double dS = 1e300, dT = 1e300;
    for (auto &ids : poly) {
      for (int r : ids) {
        if (!(rflags[r] & CRACK)) {
          continue;
        }
        double a = F.along(pos[r]);
        if (std::fabs(a) < dS) {
          dS = std::fabs(a);
          startAlong = a;
        }
        if (std::fabs(a - F.L) < dT) {
          dT = std::fabs(a - F.L);
          tipAlong = a;
        }
      }
    }
  };
  if (P.hasCrack) {
    findCrackEnds();
  }
  auto onOpenCrack = [&](int r) {
    const double tol = 1e-6 * spacing;
    if (!P.hasCrack || !(rflags[r] & CRACK)) {
      return false;
    }
    double a = F.along(pos[r]);
    return a >= startAlong - tol && a <= tipAlong + tol;
  };

  std::vector<std::set<int>> vCells(pos.size());
  for (size_t c = 0; c < poly.size(); c++) {
    for (int r : poly[c]) {
      vCells[r].insert((int)c);
    }
  }
  auto center = [&](int c) {
    V m;
    for (int r : poly[c]) {
      m += pos[r];
    }
    return m * (1.0 / (double)poly[c].size());
  };
  auto sortCell = [&](int c) {
    V m = center(c);
    std::sort(poly[c].begin(), poly[c].end(), [&](int a, int b) {
      return std::atan2(pos[a].y - m.y, pos[a].x - m.x) < std::atan2(pos[b].y - m.y, pos[b].x - m.x);
    });
  };
  auto convex = [&](int c) {
    auto &ids = poly[c];
    size_t n = ids.size();
    double area2 = 0.0;
    for (size_t k = 0; k < n; k++) {
      area2 += cross(pos[ids[k]], pos[ids[(k + 1) % n]]);
    }
    double sgn = area2 >= 0.0 ? 1.0 : -1.0;
    for (size_t k = 0; k < n; k++) {
      V e1 = pos[ids[(k + 1) % n]] - pos[ids[k]];
      V e2 = pos[ids[(k + 2) % n]] - pos[ids[(k + 1) % n]];
      if (sgn * cross(e1, e2) < -1e-12 * lRef * lRef || norm(e1) < 1e-3 * lCrit) {
        return false;
      }
    }
    return true;
  };
  // arêtes issues du sommet r toutes au moins égales au seuil
  auto longEdgesAt = [&](int r) {
    for (int c : vCells[r]) {
      auto &ids = poly[c];
      size_t n = ids.size();
      for (size_t k = 0; k < n; k++) {
        if (ids[k] != r) {
          continue;
        }
        if (norm(pos[ids[(k + 1) % n]] - pos[r]) < lCrit * (1.0 - 1e-9) ||
            norm(pos[ids[(k + n - 1) % n]] - pos[r]) < lCrit * (1.0 - 1e-9)) {
          return false;
        }
      }
    }
    return true;
  };
  // cellules autour du sommet r, dans l'ordre angulaire
  auto cellsAround = [&](int r) {
    std::vector<int> cyc(vCells[r].begin(), vCells[r].end());
    std::vector<double> ang(poly.size());
    for (int c : cyc) {
      V d = center(c) - pos[r];
      ang[c] = std::atan2(d.y, d.x);
    }
    std::sort(cyc.begin(), cyc.end(), [&](int a, int b) { return ang[a] < ang[b]; });
    return cyc;
  };
  // Sauvegarde locale (cellules et sommets touchés) pour annuler une opération refusée
  struct Snapshot {
    std::map<int, std::vector<int>> polys;
    std::map<int, std::set<int>> vc;
    std::map<int, V> p;
    std::map<int, int> fl;
  };
  auto save = [&](const std::vector<int> &cs, const std::vector<int> &rs) {
    Snapshot S;
    for (int c : cs) {
      S.polys[c] = poly[c];
    }
    for (int r : rs) {
      S.vc[r] = vCells[r];
      S.p[r] = pos[r];
      S.fl[r] = rflags[r];
    }
    return S;
  };
  auto restore = [&](const Snapshot &S) {
    for (auto &kv : S.polys) {
      poly[kv.first] = kv.second;
    }
    for (auto &kv : S.vc) {
      vCells[kv.first] = kv.second;
    }
    for (auto &kv : S.p) {
      pos[kv.first] = kv.second;
    }
    for (auto &kv : S.fl) {
      rflags[kv.first] = kv.second;
    }
  };
  // Éclatement du sommet v (cellules cyc dans l'ordre angulaire) : v garde cyc[i], cyc[i+1], cyc[i+2] et se
  // déplace vers cyc[i+1] ; w reçoit les autres cellules plus cyc[i] et cyc[i+2], qui deviennent voisines par
  // la nouvelle arête [v, w] de longueur len.
  auto split = [&](int v, int w, const std::vector<int> &cyc, size_t i, double len) {
    size_t k = cyc.size();
    int ca = cyc[i], cm = cyc[(i + 1) % k], cb = cyc[(i + 2) % k];
    V p = pos[v];
    V u = unit(center(cm) - p);
    pos[v] = p + u * (0.5 * len);
    pos[w] = p - u * (0.5 * len);
    rflags[v] = rflags[w] = 0;
    vCells[w].clear();
    for (int c : cyc) {
      if (c == ca || c == cm || c == cb) {
        continue;
      }
      std::replace(poly[c].begin(), poly[c].end(), v, w);
      vCells[v].erase(c);
      vCells[w].insert(c);
    }
    for (int c : {ca, cb}) {
      poly[c].push_back(w);
      vCells[w].insert(c);
    }
    for (int c : cyc) {
      sortCell(c);
    }
  };

  // Éclate le sommet v, puis les sommets créés, jusqu'à ce qu'ils aient tous 3 branches. Chaque éclatement
  // préfère rendre voisines les cellules qui ont le moins de côtés et refuse une cellule non convexe ; il
  // respecte si possible le seuil sur les arêtes créées (sinon l'étirement complétera). Renvoie le nombre
  // d'éclatements, ou -1 en cas d'échec (l'état est alors partiellement modifié). noise : part aléatoire du
  // choix (petite : départage des égalités ; grande : essais successifs d'autres éclatements).
  std::uniform_real_distribution<double> uTie(0.0, 1.0);
  auto explode = [&](int v0, double noise) {
    int nDone = 0;
    std::vector<int> todo = {v0};
    while (!todo.empty()) {
      int v = todo.back();
      todo.pop_back();
      if (vCells[v].size() < 4) {
        continue;
      }
      std::vector<int> cyc = cellsAround(v);
      size_t k = cyc.size();
      std::vector<std::pair<double, size_t>> order;
      for (size_t i = 0; i < k; i++) {
        order.push_back({(double)(poly[cyc[i]].size() + poly[cyc[(i + 2) % k]].size()) + noise * uTie(rng), i});
      }
      std::sort(order.begin(), order.end());
      int w = (int)pos.size();
      pos.emplace_back();
      rflags.push_back(0);
      vCells.emplace_back();
      // passe 0 : arêtes au seuil respectées ; passes 1 à 3 : seulement convexe, avec une arête créée de plus en
      // plus courte (petites arêtes voisines), que l'étirement allongera ensuite
      int done = -1;
      for (int pass = 0; pass < 4 && done < 0; pass++) {
        double len = pass <= 1 ? lT1 : (pass == 2 ? 0.5 * lT1 : 0.25 * lT1);
        for (auto &o : order) {
          Snapshot S = save(cyc, {v, w});
          split(v, w, cyc, o.second, len);
          bool ok = std::all_of(cyc.begin(), cyc.end(), [&](int c) { return convex(c); });
          if (ok && (pass >= 1 || (longEdgesAt(v) && longEdgesAt(w)))) {
            done = pass;
            break;
          }
          restore(S);
        }
      }
      if (done < 0) {
        pos.pop_back();
        rflags.pop_back();
        vCells.pop_back();
        return -1;
      }
      nDone++;
      todo.push_back(v);
      todo.push_back(w);
    }
    return nDone;
  };

  // 4a. Sommets à 4 branches ou plus (hors extrémités et lèvres de la fissure). Un échec est retenté après
  //     chaque tour de traitement des arêtes courtes.
  size_t nSplit = 0, nSplitFailed = 0;
  auto explodeAll = [&]() {
    nSplitFailed = 0;
    for (int v = 0, nv = (int)pos.size(); v < nv; v++) {
      if (vCells[v].size() >= 4 && !onOpenCrack(v) && (rflags[v] & ~CRACK) == 0) {
        int n = explode(v, 0.5);
        if (n < 0) {
          nSplitFailed++;
        } else {
          nSplit += (size_t)n;
        }
      }
    }
  };
  explodeAll();

  // arêtes sous le seuil, de la plus courte à la plus longue
  auto collectShort = [&]() {
    std::set<std::pair<int, int>> E;
    for (auto &ids : poly) {
      for (size_t k = 0; k < ids.size(); k++) {
        int a = ids[k], b = ids[(k + 1) % ids.size()];
        E.insert({std::min(a, b), std::max(a, b)});
      }
    }
    std::vector<std::pair<double, std::pair<int, int>>> shortEdges;
    for (auto &e : E) {
      double l = norm(pos[e.second] - pos[e.first]);
      if (l < lCrit * (1.0 - 1e-9)) {
        shortEdges.push_back({l, e});
      }
    }
    std::sort(shortEdges.begin(), shortEdges.end());
    return shortEdges;
  };
  // nombre d'arêtes sous le seuil et longueur minimale des arêtes des cellules cs (hors arête [ea, eb])
  auto localEdges = [&](const std::vector<int> &cs, int ea, int eb) {
    size_t nShort = 0;
    double lMin = std::numeric_limits<double>::max();
    for (int c : cs) {
      auto &ids = poly[c];
      for (size_t k = 0; k < ids.size(); k++) {
        int a = ids[k], b = ids[(k + 1) % ids.size()];
        if ((a == ea && b == eb) || (a == eb && b == ea)) {
          nShort++;
          continue;
        }
        double l = norm(pos[b] - pos[a]);
        lMin = std::min(lMin, l);
        if (l < lCrit * (1.0 - 1e-9)) {
          nShort++;
        }
      }
    }
    return std::make_pair(nShort, lMin);
  };

  // 4b. Arêtes courtes intérieures : transition T1 si les deux cellules qui perdent un côté en gardent au moins
  //     minSides et restent convexes, et si la transition réduit le nombre d'arêtes courtes de ces cellules sans
  //     raccourcir leur plus petite arête (des arêtes courtes voisines sont traitées l'une après l'autre).
  //     Les transitions T1, l'étirement (4bis), la contraction des amas (4c) et l'éclatement des sommets restants
  //     (4a) alternent jusqu'à ce que le nombre de défauts ne diminue plus.
  size_t nT1 = 0, nShortInitial = 0, nStretched = 0, nCluster = 0;
  size_t nDefectsPrev = std::numeric_limits<size_t>::max();
  for (int round = 0; round < 8; round++) {
    for (int pass = 0; pass < 10; pass++) {
      std::vector<std::pair<double, std::pair<int, int>>> shortEdges = collectShort();
      if (round == 0 && pass == 0) {
        nShortInitial = shortEdges.size();
      }
      size_t done = 0;
      for (auto &se : shortEdges) {
        int a = se.second.first, b = se.second.second;
        if (norm(pos[b] - pos[a]) >= lCrit * (1.0 - 1e-9) || rflags[a] || rflags[b] || vCells[a].size() != 3 ||
            vCells[b].size() != 3) {
          continue;
        }
        std::vector<int> LR, U, D;
        for (int c : vCells[a]) {
          (vCells[b].count(c) ? LR : U).push_back(c);
        }
        for (int c : vCells[b]) {
          if (!vCells[a].count(c)) {
            D.push_back(c);
          }
        }
        if (LR.size() != 2 || U.size() != 1 || D.size() != 1 || U[0] == D[0]) {
          continue;
        }
        if ((int)poly[LR[0]].size() - 1 < std::max(3, P.minSides) ||
            (int)poly[LR[1]].size() - 1 < std::max(3, P.minSides)) {
          continue;
        }
        std::vector<int> four = {LR[0], LR[1], U[0], D[0]};
        auto before = localEdges(four, a, b);
        Snapshot S = save(four, {a, b});
        // contraction de [a, b] en a (sommet provisoire à 4 branches), puis éclatement dans l'autre sens
        pos[a] = (pos[a] + pos[b]) * 0.5;
        for (int c : LR) {
          poly[c].erase(std::remove(poly[c].begin(), poly[c].end(), b), poly[c].end());
        }
        std::replace(poly[D[0]].begin(), poly[D[0]].end(), b, a);
        vCells[a] = {LR[0], LR[1], U[0], D[0]};
        vCells[b].clear();
        std::vector<int> cyc = cellsAround(a);
        size_t iU = std::find(cyc.begin(), cyc.end(), U[0]) - cyc.begin();
        bool ok = cyc.size() == 4 && cyc[(iU + 2) % 4] == D[0];
        if (ok) {
          split(a, b, cyc, iU, lT1);
          auto after = localEdges(four, -1, -1);
          ok = std::all_of(four.begin(), four.end(), [&](int c) { return convex(c); }) && after.first < before.first &&
               after.second >= std::min(lCrit, before.second) * (1.0 - 1e-9);
        }
        if (!ok) {
          restore(S);
          continue;
        }
        nT1++;
        done++;
      }
      if (done == 0) {
        break;
      }
    }

    // 4bis. Arêtes encore trop courtes (T1 refusée ou impossible) : leurs deux sommets sont écartés
    //       le long de l'arête jusqu'à 1.02 fois le seuil, ce qui conserve le nombre de côtés. Un sommet de bord ne se
    //       déplace que le long du bord, un sommet de la fissure le long de la fissure, un coin reste fixe ; un
    //       déplacement qui rendrait une cellule voisine non convexe est annulé.
    {
      std::vector<std::vector<int>> rootCells(pos.size());
      for (size_t c = 0; c < poly.size(); c++) {
        for (int r : poly[c]) {
          rootCells[r].push_back((int)c);
        }
      }
      auto allowedMove = [&](int r, V d) {
        int f = rflags[r];
        int nb = 0;
        for (int b : {LEFT, RIGHT, BOTTOM, TOP, CRACK}) {
          if (f & b) {
            nb++;
          }
        }
        if (nb >= 2) {
          return V(0.0, 0.0);
        }
        if (f & (LEFT | RIGHT)) {
          return V(0.0, d.y);
        }
        if (f & (BOTTOM | TOP)) {
          return V(d.x, 0.0);
        }
        if (f & CRACK) {
          return F.t * dot(d, F.t);
        }
        return d;
      };
      for (int pass = 0; pass < 20; pass++) {
        std::set<std::pair<int, int>> E;
        for (auto &ids : poly) {
          for (size_t k = 0; k < ids.size(); k++) {
            int a = ids[k], b = ids[(k + 1) % ids.size()];
            E.insert({std::min(a, b), std::max(a, b)});
          }
        }
        size_t moved = 0;
        for (auto &e : E) {
          int a = e.first, b = e.second;
          double l = norm(pos[b] - pos[a]);
          if (l >= lCrit * (1.0 - 1e-9)) {
            continue;
          }
          V d = unit(pos[b] - pos[a]);
          double deficit = 1.02 * lCrit - l; // marge : deux arêtes courtes voisines ne se raccourcissent plus en alternance
          V ma = allowedMove(a, d * (-0.5 * deficit)), mb = allowedMove(b, d * (0.5 * deficit));
          // si un sommet est bloqué, l'autre fait tout le chemin (projeté sur sa ligne)
          if (norm(ma) < 1e-15) {
            mb = allowedMove(b, d * deficit);
          } else if (norm(mb) < 1e-15) {
            ma = allowedMove(a, d * (-deficit));
          }
          V pa = pos[a], pb = pos[b];
          pos[a] += ma;
          pos[b] += mb;
          bool ok = true;
          for (int r : {a, b}) {
            for (int c : rootCells[r]) {
              if (!convex(c)) {
                ok = false;
              }
            }
          }
          if (!ok) {
            pos[a] = pa;
            pos[b] = pb;
            } else if (norm(ma) + norm(mb) > 0.0) {
            moved++;
            if (pass == 0) {
              nStretched++;
            }
          }
        }
        if (moved == 0) {
          break;
        }
      }
    }
    // 4c. Amas d'arêtes courtes voisines (aucune T1 isolée acceptable, étirement bloqué) : l'amas est contracté
    //     en un sommet, aussitôt éclaté en sommets à 3 branches. Accepté s'il réduit le nombre d'arêtes courtes
    //     sans faire passer une cellule sous minSides côtés (ni sous son nombre de côtés initial s'il est plus
    //     petit).
    for (auto &se : collectShort()) {
      int a = se.second.first, b = se.second.second;
      if (norm(pos[b] - pos[a]) >= lCrit * (1.0 - 1e-9)) {
        continue;
      }
      std::vector<int> Q = {a};
      std::set<int> inQ = {a};
      bool okQ = true;
      for (size_t qi = 0; qi < Q.size() && okQ; qi++) {
        int r = Q[qi];
        if (rflags[r] || vCells[r].size() != 3) {
          okQ = false;
          break;
        }
        for (int c : vCells[r]) {
          auto &ids = poly[c];
          size_t n = ids.size();
          for (size_t k = 0; k < n; k++) {
            if (ids[k] != r) {
              continue;
            }
            for (int q : {ids[(k + 1) % n], ids[(k + n - 1) % n]}) {
              if (!inQ.count(q) && norm(pos[q] - pos[r]) < 1.1 * lCrit) { // voisines tout juste au seuil comprises
                inQ.insert(q);
                Q.push_back(q);
              }
            }
          }
        }
      }
      if (!okQ || Q.size() > 6) {
        continue;
      }
      std::set<int> cs;
      for (int r : Q) {
        cs.insert(vCells[r].begin(), vCells[r].end());
      }
      std::vector<int> cl(cs.begin(), cs.end());
      auto before = localEdges(cl, -1, -1);
      std::map<int, size_t> sidesBefore;
      for (int c : cl) {
        sidesBefore[c] = poly[c].size();
      }
      const size_t nPos = pos.size();
      Snapshot S = save(cl, Q);
      bool ok = false;
      for (int attempt = 0; attempt < 30 && !ok; attempt++) {
        if (attempt > 0) {
          restore(S);
          pos.resize(nPos);
          rflags.resize(nPos);
          vCells.resize(nPos);
        }
        int r0 = Q[0];
        V m;
        for (int r : Q) {
          m += pos[r];
        }
        pos[r0] = m * (1.0 / (double)Q.size());
        ok = true;
        for (int c : cl) {
          std::vector<int> ids;
          for (int r : poly[c]) {
            int rr = inQ.count(r) ? r0 : r;
            if (ids.empty() || ids.back() != rr) {
              ids.push_back(rr);
            }
          }
          while (ids.size() > 1 && ids.front() == ids.back()) {
            ids.pop_back();
          }
          if (ids.size() < 3 || std::count(ids.begin(), ids.end(), r0) > 1) {
            ok = false;
          }
          poly[c] = ids;
        }
        if (!ok) {
          break; // contraction impossible : inutile d'essayer d'autres éclatements
        }
        for (int r : Q) {
          vCells[r].clear();
        }
        vCells[r0] = cs;
        ok = explode(r0, attempt == 0 ? 0.5 : 4.0) >= 0;
        if (ok) {
          auto after = localEdges(cl, -1, -1);
          ok = std::all_of(cl.begin(), cl.end(), [&](int c) { return convex(c); }) && after.first < before.first;
          for (int c : cl) {
            if (poly[c].size() < std::min(sidesBefore[c], (size_t)std::max(3, P.minSides))) {
              ok = false;
            }
          }
        }
      }
      if (!ok) {
        restore(S);
        pos.resize(nPos);
        rflags.resize(nPos);
        vCells.resize(nPos);
        continue;
      }
      nCluster++;
    }

    if (nSplitFailed > 0) {
      explodeAll();
    }
    size_t nDefects = collectShort().size() + nSplitFailed; // arêtes courtes + sommets non éclatés
    if (nDefects == 0 || nDefects >= nDefectsPrev) {
      break;
    }
    nDefectsPrev = nDefects;
  }

  // 5. Décalage des parois vers l'intérieur : barWidth/2, plus opening/2 le long de la pré-fissure
  // La fissure s'ouvre entre les sommets de la ligne les plus proches de c0 et de c1 : au-delà de la pointe,
  // la ligne se poursuit par une interface collée entre deux cellules.
  if (P.hasCrack) {
    findCrackEnds(); // l'étirement a pu déplacer les extrémités le long de la fissure
  }
  auto isCrackEdge = [&](int a, int b) {
    if (!P.hasCrack || !(rflags[a] & CRACK) || !(rflags[b] & CRACK)) {
      return false;
    }
    const double tol = 1e-6 * spacing;
    double aa = F.along(pos[a]), ab = F.along(pos[b]);
    return std::min(aa, ab) >= startAlong - tol && std::max(aa, ab) <= tipAlong + tol;
  };
  std::vector<Polygon> cells;
  std::vector<std::vector<int>> cellRoots; // sommet du pavage de chaque noeud (même ordre que cells)
  size_t nCrackEdges = 0, nBad = 0;
  for (size_t c = 0; c < poly.size(); c++) {
    auto &ids = poly[c];
    size_t n = ids.size();
    if (n < 3) {
      continue;
    }
    Polygon Pv;
    for (int r : ids) {
      Pv.push_back(pos[r]);
    }
    if (signedArea(Pv) < 0.0) {
      std::reverse(Pv.begin(), Pv.end());
      std::reverse(ids.begin(), ids.end());
    }
    std::vector<double> hEdge(n);
    for (size_t k = 0; k < n; k++) {
      bool ce = isCrackEdge(ids[k], ids[(k + 1) % n]);
      hEdge[k] = 0.5 * P.barWidth + (ce ? 0.5 * P.opening : 0.0);
      if (ce) {
        nCrackEdges++;
      }
    }
    Polygon out(n);
    for (size_t k = 0; k < n; k++) {
      size_t kp = (k + n - 1) % n;
      const V &Pk = Pv[k];
      V t1 = unit(Pk - Pv[kp]), t2 = unit(Pv[(k + 1) % n] - Pk);
      V n1(-t1.y, t1.x), n2(-t2.y, t2.x); // normales intérieures (polygone direct)
      double h1 = hEdge[kp], h2 = hEdge[k];
      double cr = cross(n1, n2);
      if (std::fabs(cr) < 1e-9) {
        out[k] = Pk + n1 * std::min(h1, h2); // côtés alignés (pointe de fissure)
      } else {
        double c1 = dot(n1, Pk) + h1, c2 = dot(n2, Pk) + h2;
        out[k] = V((c1 * n2.y - n1.y * c2) / cr, (n1.x * c2 - n2.x * c1) / cr);
      }
    }
    for (size_t k = 0; k < n; k++) {
      if (dot(out[(k + 1) % n] - out[k], Pv[(k + 1) % n] - Pv[k]) <= 0.0) {
        nBad++;
        break;
      }
    }
    cells.push_back(out);
    cellRoots.push_back(ids);
  }

  writeNodeFile(P.output, cells);

  // Contrôle final : nombre de branches de chaque sommet du pavage (arêtes issues du sommet, bords compris) et
  // arêtes sous le seuil. Seuls les sommets de la fissure ouverte (lèvres libres et extrémités) peuvent avoir
  // 4 branches.
  std::map<int, size_t> degreeHist;
  std::vector<V> badJunctions;
  size_t nShortFinal = 0, nCrackJunctions = 0;
  {
    std::set<std::pair<int, int>> E;
    for (auto &ids : poly) {
      for (size_t k = 0; k < ids.size() && ids.size() >= 3; k++) {
        int a = ids[k], b = ids[(k + 1) % ids.size()];
        E.insert({std::min(a, b), std::max(a, b)});
      }
    }
    std::map<int, int> degree;
    for (auto &e : E) {
      degree[e.first]++;
      degree[e.second]++;
      if (norm(pos[e.second] - pos[e.first]) < lCrit * (1.0 - 1e-6)) {
        nShortFinal++;
      }
    }
    for (auto &kv : degree) {
      degreeHist[kv.second]++;
      if (kv.second >= 4) {
        if (onOpenCrack(kv.first)) {
          nCrackJunctions++;
        } else {
          badJunctions.push_back(pos[kv.first]);
        }
      }
    }
  }

  // Mors : hauteur des boîtes de contrôle. Les rangées de noeuds sont définies topologiquement : la rangée 1 est
  // formée des sommets du bord, la rangée k des sommets à k-1 barres du bord. Une boîte ne sélectionne que par
  // la hauteur : elle s'arrête à mi-distance entre le noeud le plus haut des gripRows premières rangées et le noeud
  // suivant. La rangée 1 est ainsi attrapée exactement ; au-delà, les rangées d'un pavage de Voronoï ne sont pas
  // horizontales et la boîte peut attraper en plus quelques noeuds des rangées suivantes (comptés dans nExtra).
  struct Grip {
    double h{0.0};
    size_t nNodes{0}, nExtra{0};
    bool byRows{false};
  };
  const bool gripsByRows = P.hasGripRows || !(P.gripHeight > 0.0);
  const int nRows = std::max(1, P.gripRows);
  auto gripBox = [&](int sideFlag, bool bottom) {
    const double Yb = bottom ? P.ymin : P.ymin + P.Ly;
    auto depth = [&](const V &q) { return bottom ? q.y - Yb : Yb - q.y; };
    Grip g;
    std::vector<std::vector<char>> inRows(cells.size());
    if (gripsByRows) {
      g.byRows = true;
      std::map<int, std::set<int>> adj;
      for (auto &ids : cellRoots) {
        for (size_t k = 0; k < ids.size(); k++) {
          int a = ids[k], b = ids[(k + 1) % ids.size()];
          adj[a].insert(b);
          adj[b].insert(a);
        }
      }
      std::map<int, int> dist;
      std::vector<int> front;
      for (auto &kv : adj) {
        if (rflags[kv.first] & sideFlag) {
          dist[kv.first] = 0;
          front.push_back(kv.first);
        }
      }
      for (int d = 1; d < nRows && !front.empty(); d++) {
        std::vector<int> next;
        for (int u : front) {
          for (int v : adj[u]) {
            if (!dist.count(v)) {
              dist[v] = d;
              next.push_back(v);
            }
          }
        }
        front = next;
      }
      double hIn = 0.0, hNext = std::numeric_limits<double>::max();
      for (size_t c = 0; c < cells.size(); c++) {
        inRows[c].assign(cells[c].size(), 0);
        for (size_t k = 0; k < cells[c].size(); k++) {
          if (dist.count(cellRoots[c][k])) {
            inRows[c][k] = 1;
            hIn = std::max(hIn, depth(cells[c][k]));
          }
        }
      }
      for (size_t c = 0; c < cells.size(); c++) {
        for (size_t k = 0; k < cells[c].size(); k++) {
          double dq = depth(cells[c][k]);
          if (!inRows[c][k] && dq > hIn) {
            hNext = std::min(hNext, dq);
          }
        }
      }
      g.h = (hNext < std::numeric_limits<double>::max()) ? 0.5 * (hIn + hNext) : hIn + P.barWidth;
    } else {
      g.h = P.gripHeight;
    }
    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t k = 0; k < cells[c].size(); k++) {
        if (depth(cells[c][k]) <= g.h) {
          g.nNodes++;
          if (g.byRows && !inRows[c][k]) {
            g.nExtra++;
          }
        }
      }
    }
    return g;
  };
  const Grip gBot = gripBox(BOTTOM, true), gTop = gripBox(TOP, false);
  auto gripText = [&](const Grip &g) {
    std::ostringstream o;
    o << g.nNodes << " noeuds";
    if (g.byRows) {
      o << " (" << nRows << (nRows > 1 ? " rangées" : " rangée");
      if (g.nExtra > 0) {
        o << " + " << g.nExtra << " noeuds des rangées suivantes";
      }
      o << ")";
    } else {
      o << " (hauteur imposée " << std::setprecision(4) << g.h << ")";
    }
    return o.str();
  };

  ViewOptions vo;
  vo.barWidth = P.barWidth;
  vo.distGlue = P.distGlue;
  // barres courtes : celles dont l'arête du pavage (avant décalage des parois, qui les raccourcit d'environ barWidth)
  // est sous le seuil, comme dans le contrôle final ; une arête étirée tout juste au seuil n'est pas signalée
  vo.shortLength = lCrit;
  vo.shortBars.resize(cells.size());
  for (size_t c = 0; c < cells.size(); c++) {
    auto &ids = cellRoots[c];
    size_t n = ids.size();
    vo.shortBars[c].assign(n, 0);
    for (size_t k = 0; k < n; k++) {
      vo.shortBars[c][k] = norm(pos[ids[(k + 1) % n]] - pos[ids[k]]) < lCrit * (1.0 - 1e-6);
    }
  }
  vo.title = "nodegen : " + P.output;
  if (P.hasCrack) {
    vo.guides.push_back({P.c0, P.c1});
  }
  MeshStats st = analyseAndWriteSVG(cells, vo, P.svg);

  std::cout << std::setprecision(4);
  std::cout << "nodegen\n";
  std::cout << "  domaine           : [" << P.xmin << ", " << P.xmin + P.Lx << "] x [" << P.ymin << ", "
            << P.ymin + P.Ly << "]\n";
  std::cout << "  cellules          : " << st.nCells << " (visé " << N << "), noeuds : " << st.nNodes << "\n";
  std::cout << "  arêtes courtes    : " << nShortInitial << " sous le seuil " << P.minEdgeFactor << " x lRef = " << lCrit
            << " (lRef = sqrt(Lx Ly / " << nCellsTiling << ") = " << lRef << ")\n";
  std::cout << "                      " << nT1 << " T1 (longueur " << lT1 << "), " << nStretched << " étirées, " << nCluster
            << " amas contractés puis éclatés\n";
  std::cout << "  sommets éclatés   : " << nSplit << " (4 branches ou plus)"
            << (nSplitFailed ? ", " + std::to_string(nSplitFailed) + " échecs" : std::string()) << "\n";
  std::cout << "  branches/sommet   : ";
  for (auto &kv : degreeHist) {
    std::cout << kv.first << ": " << kv.second << "  ";
  }
  std::cout << "(" << nCrackJunctions << " sur la fissure ouverte, extrémités comprises)\n";
  std::cout << "  barres            : l_moy = " << st.lMean << ", l_min = " << st.lMin
            << " (l_min/l_moy = " << st.lMin / st.lMean << ", l_min/lRef = " << st.lMin / lRef << ")\n";
  {
    // cellules intérieures : aucun sommet sur le bord du domaine (les cellules de bord sont coupées par le rectangle)
    MeshStats inner;
    for (auto &ids : poly) {
      if (ids.size() < 3) {
        continue;
      }
      bool onBorder = std::any_of(ids.begin(), ids.end(), [&](int r) { return rflags[r] & (LEFT | RIGHT | BOTTOM | TOP); });
      if (!onBorder) {
        inner.sides[ids.size()]++;
        inner.nCells++;
      }
    }
    std::cout << "  cellules par nombre de côtés (intérieures : sans sommet sur le bord du domaine) :\n"
              << sidesTable({{"toutes", st}, {"intérieures", inner}}, "      ");
  }
  std::cout << "  parois            : " << st.nGlued << " collées, " << st.nFree << " libres (bords + fissure)\n";
  if (P.hasCrack) {
    V tip = F.o + F.t * tipAlong;
    std::cout << "  pré-fissure       : " << nCrackEdges << " parois ouvertes (barWidth + " << P.opening
              << "), pointe effective en (" << tip.x << ", " << tip.y << ") au lieu de (" << P.c1.x << ", "
              << P.c1.y << ")\n";
  }
  if (nBad > 0) {
    std::cout << "  ATTENTION         : " << nBad << " cellule(s) retournée(s) par le décalage des parois\n";
  }
  if (!badJunctions.empty()) {
    std::cout << "  ERREUR            : " << badJunctions.size() << " sommet(s) à 4 branches ou plus, ex. ("
              << badJunctions[0].x << ", " << badJunctions[0].y << ")\n";
  }
  if (nShortFinal > 0) {
    std::cout << "  ERREUR            : " << nShortFinal << " arête(s) encore sous le seuil " << lCrit << "\n";
  }
  std::cout << "  fichiers          : " << P.output << ", " << P.svg << ", " << P.input << "\n";
  // 6. Lignes d'input l-hyphen
  {
    const double X0 = P.xmin, X1 = P.xmin + P.Lx, Y0 = P.ymin, Y1 = P.ymin + P.Ly;
    const size_t nBottom = gBot.nNodes, nTop = gTop.nNodes;

    std::ofstream in(P.input);
    in << std::setprecision(15);
    auto num = [](double v) { // valeurs transmises à l-hyphen : précision complète
      std::ostringstream o;
      o << std::setprecision(15) << v;
      return o.str();
    };
    auto sh = [](double v) { // valeurs des commentaires
      std::ostringstream o;
      o << std::setprecision(4) << v;
      return o.str();
    };
    std::string cm = P.hasCellProps ? "" : "# ";
    in << "# ===== Lignes d'input générées par nodegen (" << argv[1] << ") =====\n";
    in << "# Domaine [" << X0 << ", " << X1 << "] x [" << Y0 << ", " << Y1 << "], " << st.nCells << " cellules, "
       << st.nNodes << " noeuds\n";
    in << "# Barres : l_moy = " << sh(st.lMean) << ", l_min = " << sh(st.lMin) << " (l_min/l_moy = "
       << sh(st.lMin / st.lMean) << ") : cleanShortBars est inutile\n";
    if (P.hasCrack) {
      V tip = F.o + F.t * tipAlong;
      in << "# Pré-fissure : de (" << P.c0.x << ", " << P.c0.y << ") à la pointe effective (" << sh(tip.x) << ", "
         << sh(tip.y) << "), lèvres écartées de barWidth + opening = " << P.barWidth + P.opening << "\n";
    }
    in << "# Les lignes commentées contiennent des valeurs <...> à compléter (ou à fournir dans le fichier de\n"
          "# paramètres de nodegen : cellProperties, glueProperties).\n\n";

    in << "# --- Échantillon (section SAMPLE)\n";
    in << "#             fichier        barWidth   Kn   Kr   Mz_max   p_int\n";
    in << cm << "readNodeFile  " << P.output << "  " << P.barWidth << "  "
       << (P.hasCellProps ? num(P.Kn) + "  " + num(P.Kr) + "  " + num(P.MzMax) + "  " + num(P.pInt)
                          : std::string("<Kn>  <Kr>  <Mz_max>  <p_int>"))
       << "\n";
    if (P.hasCellProps && P.nodeMass > 0.0 && P.Kr > 0.0) {
      double dtFlex = st.lMin * std::sqrt(P.nodeMass / P.Kr);
      in << "# setNodeMasses " << P.nodeMass << "  ->  dt_crit(flexion) = l_min * sqrt(m/Kr) = " << sh(dtFlex)
         << " s (dt <= " << sh(dtFlex / 10.0) << " conseillé)\n";
    }
    in << "# (masses, amortissements, cellContent... ici)\n\n";

    in << "# --- Chargement : mors bas et haut, libres en x (mode 1 = force nulle), vitesse imposée en y (mode 0)\n";
    in << "#                    xmin  xmax  ymin  ymax  xmode xvalue ymode yvalue\n";
    in << "setNodeControlInBox  " << X0 << "  " << X1 << "  " << Y0 << "  " << num(Y0 + gBot.h) << "  1 0.0  0 "
       << -P.pullVelocity << "    # bas  : " << gripText(gBot) << "\n";
    in << "setNodeControlInBox  " << X0 << "  " << X1 << "  " << num(Y1 - gTop.h) << "  " << Y1 << "  1 0.0  0 "
       << P.pullVelocity << "    # haut : " << gripText(gTop) << "\n";
    in << "captureNodes  bottom.txt  " << X0 << "  " << X1 << "  " << Y0 << "  " << num(Y0 + gBot.h) << "\n";
    in << "captureNodes  top.txt     " << X0 << "  " << X1 << "  " << num(Y1 - gTop.h) << "  " << Y1 << "\n\n";

    in << "# --- Collage (après les contrôles) : distGcGlue doit rester inférieur à opening = " << P.opening << "\n";
    in << "distGcGlue " << P.distGlue << "\n";
    in << (P.hasGlueProps ? "" : "# ") << "setGcGlueSameProperties  "
       << (P.hasGlueProps ? num(P.knCoh) + "  " + num(P.ktCoh) + "  " + num(P.Gc)
                          : std::string("<kn_coh>  <kt_coh>  <Gc>"))
       << "\n";

    if (P.inputTemplate.empty()) {
    std::cout << "  lignes d'input    : " << P.input << "\n";
    std::cout << "  mors bas          : " << gripText(gBot) << "\n";
    std::cout << "  mors haut         : " << gripText(gTop) << "\n";
    if (nBottom == 0 || nTop == 0) {
      std::cout << "  ATTENTION         : un mors ne contient aucun noeud (augmenter la hauteur des mors : grips)\n";
    }
    if (P.hasCrack && (std::min(P.c0.y, P.c1.y) < Y0 + gBot.h || std::max(P.c0.y, P.c1.y) > Y1 - gTop.h)) {
      std::cout << "  ATTENTION         : la pré-fissure entre dans un mors\n";
    }
    if (!P.hasCellProps) {
      std::cout << "  (readNodeFile commenté : fournir cellProperties Kn Kr Mz_max p_int pour l'activer)\n";
    }
    }
  }
  // 7. Input complet à partir d'un modèle : les lignes liées au maillage sont remplacées, le reste est recopié
  if (!P.inputTemplate.empty()) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (fs::exists(P.inputDeck) && fs::equivalent(P.inputTemplate, P.inputDeck, ec)) {
      std::cerr << "inputTemplate : le fichier produit (" << P.inputDeck << ") ne peut pas être le modèle lui-même\n";
      return 1;
    }
    std::ifstream tf(P.inputTemplate);
    if (!tf) {
      std::cerr << "inputTemplate : impossible d'ouvrir " << P.inputTemplate << "\n";
      return 1;
    }
    std::vector<std::string> lines;
    for (std::string l; std::getline(tf, l);) {
      lines.push_back(l);
    }
    auto tokens = [](const std::string &l) {
      std::vector<std::string> t;
      std::istringstream is(l);
      for (std::string w; is >> w;) {
        if (w[0] == '#') {
          break;
        }
        t.push_back(w);
      }
      return t;
    };
    auto toNum = [](const std::string &w, double &v) {
      try {
        size_t pos = 0;
        v = std::stod(w, &pos);
        return pos == w.size();
      } catch (...) {
        return false;
      }
    };
    auto num = [](double v) {
      std::ostringstream o;
      o << std::setprecision(15) << v;
      return o.str();
    };
    const std::string tag = "    # [nodegen]";

    // repérage des lignes du modèle
    std::vector<size_t> ctrl, capt;
    for (size_t i = 0; i < lines.size(); i++) {
      auto t = tokens(lines[i]);
      if (t.empty() || t[0][0] == '/' || t[0][0] == '!') {
        continue;
      }
      if (t[0] == "setNodeControlInBox" && t.size() >= 9) ctrl.push_back(i);
      if (t[0] == "captureNodes" && t.size() >= 6) capt.push_back(i);
    }
    // bas = boîte de ymin le plus petit, haut = la plus grande
    auto ordered = [&](const std::vector<size_t> &idx, int yminTok) {
      std::vector<size_t> o = idx;
      if (o.size() == 2) {
        double a = 0, b = 0;
        toNum(tokens(lines[o[0]])[yminTok], a);
        toNum(tokens(lines[o[1]])[yminTok], b);
        if (a > b) {
          std::swap(o[0], o[1]);
        }
      }
      return o;
    };
    ctrl = ordered(ctrl, 3);
    capt = ordered(capt, 4);
    const bool doCtrl = ctrl.size() == 2, doCapt = capt.size() == 2;

    const double X0 = P.xmin, X1 = P.xmin + P.Lx, Y0 = P.ymin, Y1 = P.ymin + P.Ly;
    const double hBot = gBot.h, hTop = gTop.h; // hauteurs des mors (gripRows / gripHeight), pas celles du modèle
    // chemin du nodeFile relatif au répertoire de l'input produit (run est lancé dans ce répertoire)
    fs::path deckDir = fs::absolute(P.inputDeck).parent_path();
    std::string nodePath = fs::relative(fs::absolute(P.output), deckDir, ec).generic_string();
    if (ec || nodePath.empty()) {
      nodePath = P.output;
    }

    std::ofstream out(P.inputDeck);
    out << "# Input généré par nodegen à partir du modèle " << P.inputTemplate << " et de " << argv[1] << "\n";
    out << "# Les lignes modifiées par nodegen sont marquées [nodegen] ; les autres sont celles du modèle.\n";
    std::vector<std::string> notes;
    double tplMass = -1.0, tplKr = -1.0, tplDt = -1.0;
    for (size_t i = 0; i < lines.size(); i++) {
      auto t = tokens(lines[i]);
      std::string k = t.empty() ? "" : t[0];
      if (k == "readNodeFile" && t.size() >= 7) {
        std::string props = P.hasCellProps ? num(P.Kn) + "  " + num(P.Kr) + "  " + num(P.MzMax) + "  " + num(P.pInt)
                                           : t[3] + "  " + t[4] + "  " + t[5] + "  " + t[6];
        if (!P.hasCellProps) {
          toNum(t[4], tplKr);
        }
        out << "readNodeFile  " << nodePath << "  " << P.barWidth << "  " << props << tag << "\n";
      } else if (k == "cleanShortBars") {
        out << "# " << lines[i] << tag << " inutile : l_min/l_moy = " << std::setprecision(3) << st.lMin / st.lMean
            << "\n";
      } else if (k == "setNodeControlInBox" && doCtrl && (i == ctrl[0] || i == ctrl[1])) {
        bool bottom = (i == ctrl[0]);
        double ya = bottom ? Y0 : Y1 - hTop, yb = bottom ? Y0 + hBot : Y1;
        std::string modes = P.hasVelocity ? std::string("1 0.0  0 ") + num(bottom ? -P.pullVelocity : P.pullVelocity)
                                       : t[5] + " " + t[6] + "  " + t[7] + " " + t[8];
        out << "setNodeControlInBox  " << num(X0) << "  " << num(X1) << "  " << num(ya) << "  " << num(yb) << "  "
            << modes << tag << (bottom ? " mors bas" : " mors haut") << "\n";
      } else if (k == "captureNodes" && doCapt && (i == capt[0] || i == capt[1])) {
        bool bottom = (i == capt[0]);
        double ya = bottom ? Y0 : Y1 - hTop, yb = bottom ? Y0 + hBot : Y1;
        out << "captureNodes  " << t[1] << "  " << num(X0) << "  " << num(X1) << "  " << num(ya) << "  " << num(yb)
            << tag << "\n";
      } else if (k == "setGcGlueSameProperties" && P.hasGlueProps) {
        out << "setGcGlueSameProperties  " << num(P.knCoh) << "  " << num(P.ktCoh) << "  " << num(P.Gc) << tag
            << "\n";
      } else {
        if ((k == "distGcGlue" || k == "GcGlue" || k == "glue") && t.size() >= 2) {
          double d;
          if (toNum(t[1], d) && d >= P.opening) {
            notes.push_back(k + " " + t[1] + " >= opening = " + num(P.opening) +
                            " : les lèvres de la pré-fissure seraient collées");
          }
        }
        if (k == "setNodeMasses" && t.size() >= 2) toNum(t[1], tplMass);
        if (k == "dt" && t.size() >= 2) toNum(t[1], tplDt);
        out << lines[i] << "\n";
      }
    }
    if (!doCtrl) {
      notes.push_back(std::to_string(ctrl.size()) +
                      " setNodeControlInBox dans le modèle (2 attendus : bas et haut) : lignes recopiées sans modification");
    }
    if (!doCapt) {
      notes.push_back(std::to_string(capt.size()) +
                      " captureNodes dans le modèle (2 attendus : bas et haut) : lignes recopiées sans modification");
    }
    std::cout << "  input complet     : " << P.inputDeck << " (modèle " << P.inputTemplate << ", fragments dans "
              << P.input << ")\n";
    if (doCtrl) {
      std::cout << "  mors bas          : " << gripText(gBot) << "\n";
      std::cout << "  mors haut         : " << gripText(gTop) << "\n";
    }
    double m = P.nodeMass > 0.0 ? P.nodeMass : tplMass;
    double kr = P.hasCellProps ? P.Kr : tplKr;
    if (m > 0.0 && kr > 0.0) {
      double dtFlex = st.lMin * std::sqrt(m / kr);
      std::cout << "                      dt_crit(flexion) = " << std::setprecision(4) << dtFlex;
      if (tplDt > 0.0) {
        std::cout << ", dt du modèle = " << tplDt << " (rapport " << dtFlex / tplDt << ")";
      }
      std::cout << "\n";
    }
    for (auto &n : notes) {
      std::cout << "  ATTENTION         : " << n << "\n";
    }
  }

  if (!badJunctions.empty() || nShortFinal > 0) {
    std::cerr << "nodegen : maillage non conforme (sommets à 4 branches ou arêtes courtes), voir ci-dessus\n";
    return 2;
  }
  return 0;
}
