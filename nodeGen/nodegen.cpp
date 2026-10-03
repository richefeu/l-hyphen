// nodegen : génère un nodeFile (x y idCellule) pour Lhyphen::readNodeFile.
//
// Pavage de Voronoï d'un rectangle (germes tirés avec une distance minimale puis régularisés par Lloyd),
// fusion des arêtes trop courtes, puis décalage des parois de barWidth/2 vers l'intérieur (même
// construction que cellPrepro) : deux cellules voisines sont exactement à barWidth l'une de l'autre.
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
  double minEdgeRatio{0.3}; // fusion des arêtes de Voronoï plus courtes que minEdgeRatio * longueur moyenne
  unsigned seed{1};
  std::string lattice{"hex"};    // random : tirage aléatoire ; hex : réseau hexagonal perturbé
  double disorder{0.3};          // réseau hex : amplitude de la perturbation des germes (fraction du pas)
  int minSides{5};               // la fusion des arêtes courtes ne descend pas une cellule sous minSides côtés
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
  double gripHeight{-1.0};   // hauteur des mors haut et bas (< 0 : une taille de cellule)
  double pullVelocity{3e-4}; // vitesse imposée aux mors (bas : -v, haut : +v)
  bool hasGrips{false};
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
    else if (k == "minEdgeRatio") is >> p.minEdgeRatio;
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
    } else if (k == "grips") {
      is >> p.gripHeight >> p.pullVelocity;
      p.hasGrips = true;
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
// Topologie : sommets partagés, types et fusion des arêtes courtes
// ======================================================================================================

enum : int { LEFT = 1, RIGHT = 2, BOTTOM = 4, TOP = 8, CRACK = 16 };

static int priority(int flags) {
  int nb = 0;
  for (int b : {LEFT, RIGHT, BOTTOM, TOP, CRACK}) {
    if (flags & b) {
      nb++;
    }
  }
  return nb >= 2 ? 3 : (nb == 1 ? 2 : 1);
}

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

  // 4. Fusion des arêtes courtes (union-find sur les sommets, une cellule garde au moins 3 sommets)
  std::vector<std::vector<int>> vertCells(verts.size());
  std::map<std::pair<int, int>, double> edges;
  for (size_t c = 0; c < cellV.size(); c++) {
    auto &ids = cellV[c];
    for (size_t k = 0; k < ids.size(); k++) {
      vertCells[ids[k]].push_back((int)c);
      int a = ids[k], b = ids[(k + 1) % ids.size()];
      edges[{std::min(a, b), std::max(a, b)}] = norm(verts[a] - verts[b]);
    }
  }
  double eMean = 0.0;
  for (auto &e : edges) {
    eMean += e.second;
  }
  eMean /= (double)edges.size();

  std::vector<int> parent(verts.size());
  std::vector<std::vector<int>> members(verts.size());
  for (size_t v = 0; v < verts.size(); v++) {
    parent[v] = (int)v;
    members[v] = {(int)v};
  }
  auto find = [&](int v) {
    while (parent[v] != v) {
      parent[v] = parent[parent[v]];
      v = parent[v];
    }
    return v;
  };
  auto distinctInCell = [&](int c) {
    std::set<int> r;
    for (int v : cellV[c]) {
      r.insert(find(v));
    }
    return r.size();
  };
  std::vector<std::pair<double, std::pair<int, int>>> shortEdges;
  for (auto &e : edges) {
    if (e.second < P.minEdgeRatio * eMean) {
      shortEdges.push_back({e.second, e.first});
    }
  }
  std::sort(shortEdges.begin(), shortEdges.end());
  size_t nMerged = 0, nRefused = 0;
  for (auto &se : shortEdges) {
    int ra = find(se.second.first), rb = find(se.second.second);
    if (ra == rb) {
      continue;
    }
    std::set<int> ca, both;
    for (int v : members[ra]) {
      ca.insert(vertCells[v].begin(), vertCells[v].end());
    }
    for (int v : members[rb]) {
      for (int c : vertCells[v]) {
        if (ca.count(c)) {
          both.insert(c);
        }
      }
    }
    bool ok = true;
    for (int c : both) {
      if ((int)distinctInCell(c) <= std::max(3, P.minSides)) {
        ok = false;
        break;
      }
    }
    if (!ok) {
      nRefused++;
      continue;
    }
    if (members[ra].size() < members[rb].size()) {
      std::swap(ra, rb);
    }
    parent[rb] = ra;
    members[ra].insert(members[ra].end(), members[rb].begin(), members[rb].end());
    members[rb].clear();
    nMerged++;
  }

  // Position et type des sommets fusionnés : les sommets contraints (coins, bords, fissure) imposent leur position
  std::vector<V> pos = verts;
  std::vector<int> rflags = flags;
  for (size_t r = 0; r < verts.size(); r++) {
    if (members[r].size() < 2) {
      continue;
    }
    int best = 0;
    for (int v : members[r]) {
      best = std::max(best, priority(flags[v]));
    }
    std::vector<int> top;
    for (int v : members[r]) {
      if (priority(flags[v]) == best) {
        top.push_back(v);
      }
    }
    bool sameType = std::all_of(top.begin(), top.end(), [&](int v) { return flags[v] == flags[top[0]]; });
    V m;
    if (sameType) {
      for (int v : top) {
        m += verts[v];
      }
      m = m * (1.0 / (double)top.size());
    } else {
      m = verts[top[0]];
    }
    pos[r] = m;
    rflags[r] = flags[top[0]];
  }

  // Polygones du pavage après fusion
  std::vector<std::vector<int>> poly(cellV.size());
  for (size_t c = 0; c < cellV.size(); c++) {
    for (int v : cellV[c]) {
      int r = find(v);
      if (poly[c].empty() || poly[c].back() != r) {
        poly[c].push_back(r);
      }
    }
    while (poly[c].size() > 1 && poly[c].front() == poly[c].back()) {
      poly[c].pop_back();
    }
  }

  // 4bis. Arêtes encore trop courtes (fusion refusée pour respecter minSides) : leurs deux sommets sont écartés
  //       le long de l'arête jusqu'à la longueur seuil, ce qui conserve le nombre de côtés. Un sommet de bord ne se
  //       déplace que le long du bord, un sommet de la fissure le long de la fissure, un coin reste fixe ; un
  //       déplacement qui rendrait une cellule voisine non convexe est annulé.
  size_t nStretched = 0, nStillShort = 0;
  {
    const double lCrit = P.minEdgeRatio * eMean;
    std::vector<std::vector<int>> rootCells(verts.size());
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
        if (sgn * cross(e1, e2) < -1e-12 * eMean * eMean || norm(e1) < 1e-3 * lCrit) {
          return false;
        }
      }
      return true;
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
      nStillShort = 0;
      for (auto &e : E) {
        int a = e.first, b = e.second;
        double l = norm(pos[b] - pos[a]);
        if (l >= lCrit * (1.0 - 1e-9)) {
          continue;
        }
        V d = unit(pos[b] - pos[a]);
        double deficit = lCrit - l;
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
          nStillShort++;
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

  // 5. Décalage des parois vers l'intérieur : barWidth/2, plus opening/2 le long de la pré-fissure
  // La fissure s'ouvre entre les sommets de la ligne les plus proches de c0 et de c1 : au-delà de la pointe,
  // la ligne se poursuit par une interface collée entre deux cellules.
  double startAlong = 0.0, tipAlong = 0.0;
  if (P.hasCrack) {
    double dS = 1e300, dT = 1e300;
    for (size_t c = 0; c < poly.size(); c++) {
      for (int r : poly[c]) {
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
  }

  writeNodeFile(P.output, cells);

  ViewOptions vo;
  vo.barWidth = P.barWidth;
  vo.distGlue = P.distGlue;
  vo.shortRatio = 0.3;
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
  std::cout << "  arêtes fusionnées : " << nMerged << " (" << nRefused << " refusées : minSides), seuil "
            << P.minEdgeRatio << " x " << eMean << "\n";
  std::cout << "  arêtes étirées    : " << nStretched << " (fusion refusée pour garder >= " << P.minSides
            << " côtés)" << (nStillShort ? ", " + std::to_string(nStillShort) + " encore trop courtes" : std::string())
            << "\n";
  std::cout << "  barres            : l_moy = " << st.lMean << ", l_min = " << st.lMin
            << " (l_min/l_moy = " << st.lMin / st.lMean << ")\n";
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
    std::cout << "  côtés             : " << sidesHistogram(st) << "\n";
    std::cout << "  côtés (intérieur) : " << sidesHistogram(inner) << "(" << inner.nCells << " cellules)\n";
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
  std::cout << "  fichiers          : " << P.output << ", " << P.svg << ", " << P.input << "\n";
  // 6. Lignes d'input l-hyphen
  {
    const double X0 = P.xmin, X1 = P.xmin + P.Lx, Y0 = P.ymin, Y1 = P.ymin + P.Ly;
    const double grip = P.gripHeight > 0.0 ? P.gripHeight : P.cellSize;
    auto countIn = [&](double ya, double yb) {
      size_t nb = 0;
      for (auto &c : cells) {
        for (auto &q : c) {
          if (q.x >= X0 && q.x <= X1 && q.y >= ya && q.y <= yb) {
            nb++;
          }
        }
      }
      return nb;
    };
    size_t nBottom = countIn(Y0, Y0 + grip), nTop = countIn(Y1 - grip, Y1);

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
    in << "setNodeControlInBox  " << X0 << "  " << X1 << "  " << Y0 << "  " << num(Y0 + grip) << "  1 0.0  0 "
       << -P.pullVelocity << "    # bas  : " << nBottom << " noeuds\n";
    in << "setNodeControlInBox  " << X0 << "  " << X1 << "  " << num(Y1 - grip) << "  " << Y1 << "  1 0.0  0 "
       << P.pullVelocity << "    # haut : " << nTop << " noeuds\n";
    in << "captureNodes  bottom.txt  " << X0 << "  " << X1 << "  " << Y0 << "  " << num(Y0 + grip) << "\n";
    in << "captureNodes  top.txt     " << X0 << "  " << X1 << "  " << num(Y1 - grip) << "  " << Y1 << "\n\n";

    in << "# --- Collage (après les contrôles) : distGcGlue doit rester inférieur à opening = " << P.opening << "\n";
    in << "distGcGlue " << P.distGlue << "\n";
    in << (P.hasGlueProps ? "" : "# ") << "setGcGlueSameProperties  "
       << (P.hasGlueProps ? num(P.knCoh) + "  " + num(P.ktCoh) + "  " + num(P.Gc)
                          : std::string("<kn_coh>  <kt_coh>  <Gc>"))
       << "\n";

    if (P.inputTemplate.empty()) {
    std::cout << "  lignes d'input    : " << P.input << " (mors de " << grip << " : " << nBottom << " noeuds en bas, "
              << nTop << " en haut)\n";
    if (nBottom == 0 || nTop == 0) {
      std::cout << "  ATTENTION         : un mors ne contient aucun noeud (augmenter la hauteur des mors : grips)\n";
    }
    if (P.hasCrack && (std::min(P.c0.y, P.c1.y) < Y0 + grip || std::max(P.c0.y, P.c1.y) > Y1 - grip)) {
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
    double hBot = P.gripHeight > 0.0 ? P.gripHeight : P.cellSize, hTop = hBot;
    if (!P.hasGrips && doCtrl) { // hauteur des mors reprise du modèle
      auto tb = tokens(lines[ctrl[0]]), tt = tokens(lines[ctrl[1]]);
      double a, b;
      if (toNum(tb[3], a) && toNum(tb[4], b) && b > a) hBot = b - a;
      if (toNum(tt[3], a) && toNum(tt[4], b) && b > a) hTop = b - a;
    }
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
        std::string modes = P.hasGrips ? std::string("1 0.0  0 ") + num(bottom ? -P.pullVelocity : P.pullVelocity)
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
    auto countIn = [&](double ya, double yb) {
      size_t nb = 0;
      for (auto &c : cells) {
        for (auto &q : c) {
          if (q.x >= X0 && q.x <= X1 && q.y >= ya && q.y <= yb) {
            nb++;
          }
        }
      }
      return nb;
    };
    std::cout << "  input complet     : " << P.inputDeck << " (modèle " << P.inputTemplate << ", fragments dans "
              << P.input << ")\n";
    if (doCtrl) {
      std::cout << "                      mors : " << countIn(Y0, Y0 + hBot) << " noeuds en bas (hauteur " << hBot
                << "), " << countIn(Y1 - hTop, Y1) << " en haut (hauteur " << hTop << ")\n";
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

  return 0;
}
