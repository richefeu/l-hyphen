// Outils communs à nodegen et nodeview : géométrie 2D, lecture/écriture des nodeFiles
// (format 'x y idCellule' lu par Lhyphen::readNodeFile) et rendu SVG d'analyse.
#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ======================================================================================================
// Géométrie
// ======================================================================================================

struct V {
  double x{0.0}, y{0.0};
  V() = default;
  V(double X, double Y) : x(X), y(Y) {}
  V operator+(const V &o) const { return V(x + o.x, y + o.y); }
  V operator-(const V &o) const { return V(x - o.x, y - o.y); }
  V operator*(double k) const { return V(x * k, y * k); }
  V &operator+=(const V &o) {
    x += o.x;
    y += o.y;
    return *this;
  }
};

inline double dot(const V &a, const V &b) { return a.x * b.x + a.y * b.y; }
inline double cross(const V &a, const V &b) { return a.x * b.y - a.y * b.x; }
inline double norm(const V &a) { return std::sqrt(dot(a, a)); }
inline V unit(const V &a) { return a * (1.0 / norm(a)); }

using Polygon = std::vector<V>;

inline double signedArea(const Polygon &p) {
  double a = 0.0;
  for (size_t i = 0; i < p.size(); i++) {
    a += cross(p[i], p[(i + 1) % p.size()]);
  }
  return 0.5 * a;
}

inline V centroid(const Polygon &p) {
  double a = 0.0;
  V c;
  for (size_t i = 0; i < p.size(); i++) {
    const V &u = p[i];
    const V &w = p[(i + 1) % p.size()];
    double k = cross(u, w);
    a += k;
    c += (u + w) * k;
  }
  if (std::fabs(a) < 1e-300) {
    V m;
    for (auto &q : p) {
      m += q;
    }
    return m * (1.0 / (double)p.size());
  }
  return c * (1.0 / (3.0 * a));
}

// Distance d'un point au segment [a, b]
inline double distPointSegment(const V &p, const V &a, const V &b) {
  V ab = b - a;
  double l2 = dot(ab, ab);
  double s = (l2 > 0.0) ? std::clamp(dot(p - a, ab) / l2, 0.0, 1.0) : 0.0;
  return norm(a + ab * s - p);
}

// ======================================================================================================
// nodeFile
// ======================================================================================================

// Lit un nodeFile comme Lhyphen::readNodeFile : regroupement par identifiant puis tri angulaire
// des noeuds autour de leur barycentre (les cellules sont supposées convexes).
inline std::vector<Polygon> readNodeFile(const std::string &name) {
  std::ifstream file(name);
  if (!file) {
    std::cerr << "Impossible d'ouvrir " << name << "\n";
    return {};
  }
  std::map<long, Polygon> byId;
  double x, y;
  long id;
  while (file >> x >> y >> id) {
    byId[id].emplace_back(x, y);
  }
  std::vector<Polygon> cells;
  for (auto &kv : byId) {
    Polygon p = kv.second;
    V c;
    for (auto &q : p) {
      c += q;
    }
    c = c * (1.0 / (double)p.size());
    std::stable_sort(p.begin(), p.end(), [&](const V &a, const V &b) {
      return std::atan2(a.y - c.y, a.x - c.x) < std::atan2(b.y - c.y, b.x - c.x);
    });
    cells.push_back(p);
  }
  return cells;
}

inline void writeNodeFile(const std::string &name, const std::vector<Polygon> &cells) {
  std::ofstream file(name);
  file << std::setprecision(17);
  for (size_t c = 0; c < cells.size(); c++) {
    for (auto &q : cells[c]) {
      file << q.x << ' ' << q.y << ' ' << c << '\n';
    }
  }
}

// ======================================================================================================
// Analyse et rendu SVG
// ======================================================================================================

struct ViewOptions {
  double barWidth{-1.0};   // < 0 : estimée (distance minimale entre noeuds de cellules différentes)
  double distGlue{2e-7};   // tolérance de collage (distGcGlue / glue dans l'input)
  double shortRatio{0.3};  // barres plus courtes que shortRatio * longueur moyenne mises en évidence
  bool showNodes{false};   // dessiner les noeuds
  bool useBox{false};      // ne dessiner qu'une fenêtre [bx0, bx1] x [by0, by1]
  double bx0{0}, bx1{0}, by0{0}, by1{0};
  double widthPx{1400.0};  // largeur de l'image
  bool bare{false};        // sans titre, statistiques ni légende (figures)
  std::string title;
  std::vector<std::pair<V, V>> guides; // segments de repère (ex. ligne de pré-fissure), en pointillés
};

struct MeshStats {
  size_t nCells{0}, nNodes{0}, nBars{0};
  size_t nGlued{0}, nFree{0}, nShort{0};
  double lMean{0.0}, lMin{0.0}, barWidth{0.0};
  std::map<size_t, size_t> sides; // nombre de cellules par nombre de côtés
};

// Répartition des cellules par nombre de côtés, ex. "4: 2.1%  5: 30.2%  6: 60.0%  7: 7.7%"
inline std::string sidesHistogram(const MeshStats &st) {
  std::ostringstream os;
  os << std::fixed << std::setprecision(1);
  for (auto &kv : st.sides) {
    os << kv.first << ": " << 100.0 * (double)kv.second / (double)st.nCells << "%  ";
  }
  return os.str();
}

// Distance minimale entre deux noeuds de cellules différentes (= barWidth pour un maillage de cellPrepro/nodegen)
inline double minInterCellNodeDistance(const std::vector<Polygon> &cells, double h) {
  std::unordered_map<long long, std::vector<std::pair<int, V>>> grid;
  auto key = [](long long i, long long j) { return i * 1000003LL + j; };
  for (size_t c = 0; c < cells.size(); c++) {
    for (auto &q : cells[c]) {
      grid[key((long long)std::floor(q.x / h), (long long)std::floor(q.y / h))].push_back({(int)c, q});
    }
  }
  double dmin = std::numeric_limits<double>::max();
  for (size_t c = 0; c < cells.size(); c++) {
    for (auto &q : cells[c]) {
      long long i = (long long)std::floor(q.x / h), j = (long long)std::floor(q.y / h);
      for (long long di = -1; di <= 1; di++) {
        for (long long dj = -1; dj <= 1; dj++) {
          auto it = grid.find(key(i + di, j + dj));
          if (it == grid.end()) {
            continue;
          }
          for (auto &e : it->second) {
            if (e.first != (int)c) {
              dmin = std::min(dmin, norm(e.second - q));
            }
          }
        }
      }
    }
  }
  return dmin;
}

// Analyse le maillage (barres collées ou libres au sens de Lhyphen::glue, barres courtes) et écrit un SVG.
// Une barre est « collée » si un noeud d'une autre cellule est à moins de barWidth + distGlue du segment
// (même critère que la création des liens cohésifs) ; sinon elle est libre (bord, pré-fissure, trou).
inline MeshStats analyseAndWriteSVG(const std::vector<Polygon> &cells, ViewOptions opt, const std::string &svgFile) {
  MeshStats st;
  st.nCells = cells.size();
  double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
  double lSum = 0.0;
  st.lMin = std::numeric_limits<double>::max();
  for (auto &p : cells) {
    st.nNodes += p.size();
    st.sides[p.size()]++;
    for (size_t i = 0; i < p.size(); i++) {
      double l = norm(p[(i + 1) % p.size()] - p[i]);
      lSum += l;
      st.lMin = std::min(st.lMin, l);
      st.nBars++;
      xmin = std::min(xmin, p[i].x);
      xmax = std::max(xmax, p[i].x);
      ymin = std::min(ymin, p[i].y);
      ymax = std::max(ymax, p[i].y);
    }
  }
  if (st.nBars == 0) {
    return st;
  }
  st.lMean = lSum / (double)st.nBars;
  if (opt.barWidth < 0.0) {
    opt.barWidth = minInterCellNodeDistance(cells, st.lMean);
  }
  st.barWidth = opt.barWidth;

  // Test de collage d'une paroi : son milieu est-il à moins de barWidth + distGlue d'une barre d'une autre
  // cellule ? (les extrémités d'une barre sont toujours proches des noeuds voisins au sommet commun, elles ne
  // permettent pas de distinguer une paroi partagée d'une paroi libre)
  const double reach = opt.barWidth + opt.distGlue;
  struct Seg {
    V a, b;
    int cell;
    int cls; // 0 collée, 1 libre
    bool isShort;
  };
  std::vector<Seg> segs;
  double lMax = 0.0;
  for (size_t c = 0; c < cells.size(); c++) {
    const Polygon &p = cells[c];
    for (size_t i = 0; i < p.size(); i++) {
      V a = p[i], b = p[(i + 1) % p.size()];
      lMax = std::max(lMax, norm(b - a));
      segs.push_back({a, b, (int)c, 1, norm(b - a) < opt.shortRatio * st.lMean});
    }
  }
  const double h = st.lMean;
  std::unordered_map<long long, std::vector<int>> grid; // barres indexées par leur milieu
  auto key = [](long long i, long long j) { return i * 1000003LL + j; };
  for (size_t s = 0; s < segs.size(); s++) {
    V m = (segs[s].a + segs[s].b) * 0.5;
    grid[key((long long)std::floor(m.x / h), (long long)std::floor(m.y / h))].push_back((int)s);
  }
  const long long rng = (long long)std::ceil((0.5 * lMax + reach) / h);
  for (auto &sg : segs) {
    V m = (sg.a + sg.b) * 0.5;
    long long i = (long long)std::floor(m.x / h), j = (long long)std::floor(m.y / h);
    for (long long di = -rng; di <= rng && sg.cls == 1; di++) {
      for (long long dj = -rng; dj <= rng && sg.cls == 1; dj++) {
        auto it = grid.find(key(i + di, j + dj));
        if (it == grid.end()) {
          continue;
        }
        for (int o : it->second) {
          if (segs[o].cell != sg.cell && distPointSegment(m, segs[o].a, segs[o].b) < reach) {
            sg.cls = 0;
            break;
          }
        }
      }
    }
    (sg.cls == 0 ? st.nGlued : st.nFree)++;
    if (sg.isShort) {
      st.nShort++;
    }
  }

  // Fenêtre et transformation monde -> pixels (y vers le haut)
  double X0 = opt.useBox ? opt.bx0 : xmin, X1 = opt.useBox ? opt.bx1 : xmax;
  double Y0 = opt.useBox ? opt.by0 : ymin, Y1 = opt.useBox ? opt.by1 : ymax;
  const double margin = opt.bare ? 4.0 : 20.0, header = opt.bare ? 0.0 : 90.0;
  double sc = (opt.widthPx - 2.0 * margin) / (X1 - X0);
  double heightPx = (Y1 - Y0) * sc + 2.0 * margin + header;
  auto px = [&](const V &q) { return V(margin + (q.x - X0) * sc, header + margin + (Y1 - q.y) * sc); };
  auto inBox = [&](const V &a, const V &b) {
    if (!opt.useBox) {
      return true;
    }
    return !(std::max(a.x, b.x) < X0 || std::min(a.x, b.x) > X1 || std::max(a.y, b.y) < Y0 ||
             std::min(a.y, b.y) > Y1);
  };
  // épaisseurs : à l'échelle si la fenêtre est petite (zoom), sinon valeurs minimales lisibles
  double wGlued = std::max(0.5, 0.5 * opt.barWidth * sc);
  double wFree = std::max(1.6, 1.0 * opt.barWidth * sc);
  double wShort = std::max(2.6, 1.5 * opt.barWidth * sc);

  std::ofstream f(svgFile);
  f << std::fixed << std::setprecision(2);
  f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << opt.widthPx << "\" height=\"" << heightPx
    << "\" viewBox=\"0 0 " << opt.widthPx << ' ' << heightPx << "\" font-family=\"Helvetica, Arial, sans-serif\">\n";
  f << "<rect width=\"100%\" height=\"100%\" fill=\"#fcfcfb\"/>\n";
  f << "<defs><clipPath id=\"win\"><rect x=\"" << margin << "\" y=\"" << header + margin << "\" width=\""
    << (X1 - X0) * sc << "\" height=\"" << (Y1 - Y0) * sc << "\"/></clipPath></defs>\n";
  if (!opt.bare) {
    f << "<text x=\"" << margin << "\" y=\"26\" font-size=\"18\" fill=\"#1f1e1b\">" << opt.title << "</text>\n";
    std::ostringstream info;
    info << std::setprecision(3) << st.nCells << " cellules, " << st.nNodes << " noeuds  |  l_moy = " << st.lMean
         << ", l_min = " << st.lMin << " (l_min/l_moy = " << st.lMin / st.lMean << ")  |  barWidth = " << opt.barWidth
         << ", distGlue = " << opt.distGlue;
    f << "<text x=\"" << margin << "\" y=\"48\" font-size=\"12.5\" fill=\"#3b3a35\">" << info.str() << "</text>\n";
    f << "<text x=\"" << margin << "\" y=\"66\" font-size=\"12.5\" fill=\"#3b3a35\">cellules par nombre de côtés  "
      << sidesHistogram(st) << "</text>\n";
    // légende
    double lx = margin, ly = 84;
    auto legend = [&](const char *col, double w, const std::string &txt) {
      f << "<line x1=\"" << lx << "\" y1=\"" << ly << "\" x2=\"" << lx + 22 << "\" y2=\"" << ly << "\" stroke=\"" << col
        << "\" stroke-width=\"" << std::min(w, 4.0) << "\" stroke-linecap=\"round\"/>";
      f << "<text x=\"" << lx + 28 << "\" y=\"" << ly + 4 << "\" font-size=\"12\" fill=\"#3b3a35\">" << txt
        << "</text>\n";
      lx += 30 + 7.0 * (double)txt.size();
    };
    legend("#9a9990", 1.0, "paroi collée (" + std::to_string(st.nGlued) + ")");
    legend("#eb6834", 2.0, "paroi libre : bord / pré-fissure (" + std::to_string(st.nFree) + ")");
    {
      std::ostringstream s;
      s << "barre &lt; " << opt.shortRatio << " l_moy (" << st.nShort << ")";
      legend("#2a78d6", 3.0, s.str());
    }
  }

  f << "<g clip-path=\"url(#win)\" fill=\"none\" stroke-linecap=\"round\">\n";
  auto path = [&](const char *col, double w, auto pred) {
    f << "<path stroke=\"" << col << "\" stroke-width=\"" << w << "\" d=\"";
    for (auto &s : segs) {
      if (pred(s) && inBox(s.a, s.b)) {
        V A = px(s.a), B = px(s.b);
        f << 'M' << A.x << ' ' << A.y << 'L' << B.x << ' ' << B.y;
      }
    }
    f << "\"/>\n";
  };
  path("#9a9990", wGlued, [](const Seg &s) { return s.cls == 0; });
  path("#eb6834", wFree, [](const Seg &s) { return s.cls == 1; });
  path("#2a78d6", wShort, [](const Seg &s) { return s.isShort; });
  for (auto &g : opt.guides) {
    V A = px(g.first), B = px(g.second);
    f << "<line x1=\"" << A.x << "\" y1=\"" << A.y << "\" x2=\"" << B.x << "\" y2=\"" << B.y
      << "\" stroke=\"#1f1e1b\" stroke-width=\"0.8\" stroke-dasharray=\"4 4\"/>\n";
  }
  if (opt.showNodes) {
    double r = std::max(0.8, 0.5 * opt.barWidth * sc);
    f << "<g fill=\"#3b3a35\" stroke=\"none\">";
    for (auto &p : cells) {
      for (auto &q : p) {
        if (!opt.useBox || (q.x >= X0 && q.x <= X1 && q.y >= Y0 && q.y <= Y1)) {
          V Q = px(q);
          f << "<circle cx=\"" << Q.x << "\" cy=\"" << Q.y << "\" r=\"" << r << "\"/>";
        }
      }
    }
    f << "</g>\n";
  }
  f << "</g>\n</svg>\n";
  return st;
}
