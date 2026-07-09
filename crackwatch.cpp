#include "Lhyphen.hpp"
#include <iostream>

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <queue>

#include "fileTool.hpp"

Lhyphen Conf;

struct BreakEvent {
  double time;
  size_t a_ci, a_cj, a_in, a_jn;
  size_t b_ci, b_cj, b_in, b_jn;
  double nrj;
};
std::vector<BreakEvent> breakEvents;

void readBreakHistory(const char *fname) {
  breakEvents.clear();
  if (!fileTool::fileExists(fname)) {
    return;
  }
  std::ifstream file(fname);
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    BreakEvent ev;
    if (iss >> ev.time >> ev.a_ci >> ev.a_cj >> ev.a_in >> ev.a_jn >> ev.b_ci >> ev.b_cj >> ev.b_in >> ev.b_jn >>
        ev.nrj) {
      breakEvents.push_back(ev);
    }
  }
  std::cout << "Read " << fname << " : " << breakEvents.size() << " break events" << std::endl;
}

// Renvoie true et remplit pos si (ci, cj, in, jn) est indexable dans le conf courant
static bool crackContactPos(size_t ci, size_t cj, size_t in, size_t jn, vec2r &pos) {
  if (ci >= Conf.cells.size() || cj >= Conf.cells.size()) {
    return false;
  }
  if (in >= Conf.cells[ci].nodes.size() || jn >= Conf.cells[cj].nodes.size()) {
    return false;
  }
  Conf.getPosition(ci, cj, in, jn, pos);
  return true;
}

struct Connection {
  size_t i, j;
};

struct CrackNode {
  vec2r pos;
  int order; // nb de micro fissures connectées au noeud
};

struct Crack {
  std::vector<CrackNode> nodes;
  std::vector<Connection> Conns;
};

void saveCrackSVG(const Crack &crack, const char *fname, const std::vector<size_t> *mainPathNodes = nullptr,
                  int CanvasWidth = 800) {
  if (crack.nodes.empty()) {
    std::cout << "saveCrackSVG : graphe vide, rien à écrire" << std::endl;
    return;
  }

  // Boite qui englobe les noeuds
  double xmin = crack.nodes[0].pos.x, xmax = crack.nodes[0].pos.x;
  double ymin = crack.nodes[0].pos.y, ymax = crack.nodes[0].pos.y;
  for (size_t k = 1; k < crack.nodes.size(); k++) {
    xmin = std::min(xmin, crack.nodes[k].pos.x);
    xmax = std::max(xmax, crack.nodes[k].pos.x);
    ymin = std::min(ymin, crack.nodes[k].pos.y);
    ymax = std::max(ymax, crack.nodes[k].pos.y);
  }
  // Garde-fou
  double spanx = xmax - xmin;
  double spany = ymax - ymin;
  double span = std::max(std::max(spanx, spany), 1.0e-12);

  std::ofstream ofs(fname);
  SVGfile svg(ofs);

  int CanvasHeight = (int)(CanvasWidth * (spany > 0.0 ? spany : span) / (spanx > 0.0 ? spanx : span));
  svg.begin(CanvasWidth, CanvasHeight);
  viewZone vz(0, 0, CanvasWidth, CanvasHeight);
  double delta = 0.05 * span; // bordure
  vz.adjustRange(xmin - delta, xmax + delta, ymin - delta, ymax + delta);

  // Connexions (segments de fissure).
  for (size_t c = 0; c < crack.Conns.size(); c++) {
    const vec2r &A = crack.nodes[crack.Conns[c].i].pos;
    const vec2r &B = crack.nodes[crack.Conns[c].j].pos;
    svg.line(vz, A.x, A.y, B.x, B.y, "stroke:red;stroke-width:2;stroke-linecap:round");
  }

  // Chemin principal (par-dessus les segments), tracé comme une polyligne épaisse orange.
  if (mainPathNodes != nullptr && mainPathNodes->size() >= 2) {
    for (size_t k = 1; k < mainPathNodes->size(); k++) {
      const vec2r &A = crack.nodes[(*mainPathNodes)[k - 1]].pos;
      const vec2r &B = crack.nodes[(*mainPathNodes)[k]].pos;
      svg.line(vz, A.x, A.y, B.x, B.y, "stroke:orange;stroke-width:4;stroke-linecap:round");
    }
  }

  // Noeuds du graphe.
  for (size_t k = 0; k < crack.nodes.size(); k++) {
    double cx = crack.nodes[k].pos.x * vz.scalex + vz.x0;
    double cy = crack.nodes[k].pos.y * vz.scaley + vz.y0;

    // std::cout << crack.nodes[k].order << std::endl;

    if (crack.nodes[k].order == 1)
      svg.circle(cx, cy, 3.0, "stroke:none;fill:black");
    else if (crack.nodes[k].order == 2)
      svg.circle(cx, cy, 3.0, "stroke:none;fill:green");
    else if (crack.nodes[k].order == 3)
      svg.circle(cx, cy, 3.0, "stroke:none;fill:blue");
    else if (crack.nodes[k].order > 3)
      svg.circle(cx, cy, 3.0, "stroke:none;fill:red");
  }

  svg.end();
  std::cout << "save file: " << fname << std::endl;
}

// Ajoute un point au graphe en fusionnant si distance < sqrt(tol2)
size_t addNode(Crack &crack, const vec2r &p, double tol2) {
  for (size_t k = 0; k < crack.nodes.size(); k++) {

    if (norm2(crack.nodes[k].pos - p) < tol2) {
      crack.nodes[k].order += 1;
      return k;
    }
  }
  CrackNode N;
  N.pos = p;
  N.order = 1;
  crack.nodes.push_back(N);
  return crack.nodes.size() - 1;
}

// Scinde un graphe de fissures en ses composantes connexes (chaque Crack renvoyé est
// entièrement connecté). Union-find sur les Conns, puis on recopie noeuds et connexions en
// réindexant localement chaque composante.
std::vector<Crack> splitCracks(const Crack &crack) {
  size_t n = crack.nodes.size();

  std::vector<size_t> parent(n);
  for (size_t i = 0; i < n; i++) {
    parent[i] = i;
  }
  std::function<size_t(size_t)> find = [&](size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]]; // compression de chemin
      x = parent[x];
    }
    return x;
  };
  for (size_t c = 0; c < crack.Conns.size(); c++) {
    parent[find(crack.Conns[c].i)] = find(crack.Conns[c].j);
  }

  // root -> index du Crack dans le résultat ; nodeMap : ancien indice global -> nouvel indice local.
  std::map<size_t, size_t> rootToCrack;
  std::vector<size_t> nodeMap(n);
  std::vector<Crack> result;
  for (size_t i = 0; i < n; i++) {
    size_t r = find(i);
    auto it = rootToCrack.find(r);
    if (it == rootToCrack.end()) {
      it = rootToCrack.emplace(r, result.size()).first;
      result.push_back(Crack());
    }
    Crack &dst = result[it->second];
    nodeMap[i] = dst.nodes.size();
    dst.nodes.push_back(crack.nodes[i]);
  }
  for (size_t c = 0; c < crack.Conns.size(); c++) {
    size_t ci = rootToCrack[find(crack.Conns[c].i)];
    result[ci].Conns.push_back({nodeMap[crack.Conns[c].i], nodeMap[crack.Conns[c].j]});
  }
  return result;
}

// Chemin principal d'une fissure connectée : parmi toutes les paires d'extrémités
// (noeuds d'order 1), on calcule le plus court chemin (Dijkstra pondéré par la longueur
// géométrique des segments) et on retient le plus long de ces plus courts chemins. Renvoie la
// liste ordonnée des noeuds du chemin et met sa longueur dans `length`.
std::vector<size_t> mainPath(const Crack &crack, double &length) {
  length = 0.0;
  size_t n = crack.nodes.size();

  // Liste d'adjacence pondérée.
  std::vector<std::vector<std::pair<size_t, double>>> adj(n);
  for (size_t c = 0; c < crack.Conns.size(); c++) {
    size_t i = crack.Conns[c].i, j = crack.Conns[c].j;
    double w = norm(crack.nodes[i].pos - crack.nodes[j].pos);
    adj[i].push_back({j, w});
    adj[j].push_back({i, w});
  }

  // Extrémités = noeuds d'order 1.
  std::vector<size_t> ends;
  for (size_t i = 0; i < n; i++) {
    if (crack.nodes[i].order == 1) {
      ends.push_back(i);
    }
  }

  const double INF = std::numeric_limits<double>::infinity();
  const size_t NONE = std::numeric_limits<size_t>::max();
  std::vector<size_t> best;

  for (size_t s : ends) {
    std::vector<double> dist(n, INF);
    std::vector<size_t> prev(n, NONE);
    dist[s] = 0.0;
    std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>,
                        std::greater<std::pair<double, size_t>>>
        pq;
    pq.push({0.0, s});
    while (!pq.empty()) {
      double d = pq.top().first;
      size_t u = pq.top().second;
      pq.pop();
      if (d > dist[u]) {
        continue;
      }
      for (size_t e = 0; e < adj[u].size(); e++) {
        size_t v = adj[u][e].first;
        double nd = d + adj[u][e].second;
        if (nd < dist[v]) {
          dist[v] = nd;
          prev[v] = u;
          pq.push({nd, v});
        }
      }
    }

    // On ne teste chaque paire qu'une fois (t > s) et on garde le plus long plus-court-chemin.
    for (size_t t : ends) {
      if (t <= s || dist[t] == INF) {
        continue;
      }
      if (dist[t] > length) {
        length = dist[t];
        best.clear();
        for (size_t v = t; v != NONE; v = prev[v]) {
          best.push_back(v);
        }
        std::reverse(best.begin(), best.end()); // de s vers t
      }
    }
  }
  return best;
}

int main(/*int argc, const char *argv[]*/) {

  Conf.loadCONF("conf0");
  readBreakHistory("breakHistory.txt");

  // Algo naif appliqué à toutes les ruptures
  double Ltot = 0.0;
  double Wcum = 0.0;
  std::ofstream file("out.txt");
  for (size_t i = 0 ; i < breakEvents.size() ; i++) {
    vec2r pos1, pos2;
    crackContactPos(breakEvents[i].a_ci, breakEvents[i].a_cj, breakEvents[i].a_in, breakEvents[i].a_jn, pos1);
    crackContactPos(breakEvents[i].b_ci, breakEvents[i].b_cj, breakEvents[i].b_in, breakEvents[i].b_jn, pos2);
    vec2r Lvec = pos2 - pos1;
    double L = Lvec.length();
    Ltot += L;
    Wcum += breakEvents[i].nrj;
    file << breakEvents[i].time << " " << Ltot << " " << Wcum << std::endl;
  }

  // construction de la structure de crack entier
  struct Segment {
    vec2r a, b;
    size_t ev; // indice de l'évènement de rupture d'origine (pour journalisation)
  };

  // Premier passage : on construit un segment par rupture indexable (une interface rompue relie ses
  // deux points de contact A et B). On garde aussi la longueur de chaque segment.
  std::vector<Segment> allSegs;
  std::vector<double> lens;
  for (size_t i = 0; i < breakEvents.size(); i++) {
    Segment seg;
    seg.ev = i;
    if (!crackContactPos(breakEvents[i].a_ci, breakEvents[i].a_cj, breakEvents[i].a_in, breakEvents[i].a_jn, seg.a)) {
      continue;
    }
    if (!crackContactPos(breakEvents[i].b_ci, breakEvents[i].b_cj, breakEvents[i].b_in, breakEvents[i].b_jn, seg.b)) {
      continue;
    }
    allSegs.push_back(seg);
    lens.push_back(norm(seg.b - seg.a));
  }

  // Longueur de référence robuste = médiane des longueurs de segments. Contrairement au minimum,
  // elle n'est pas empoisonnée par quelques ruptures dégénérées (interfaces dont les deux points de
  // contact sont quasi confondus)
  double medLen = 0.0;
  if (!lens.empty()) {
    std::vector<double> tmp = lens;
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
    medLen = tmp[tmp.size() / 2];
  }

  // Une rupture est jugée « problématique » (dégénérée) si son segment est beaucoup plus court que la
  // médiane : l'interface rompue n'a alors pas d'étendue physique et parasiterait la reconstruction du
  // graphe (empoisonnement de la tolérance de fusion + composantes fantômes). On l'écarte et on la
  // journalise.
  const double degRatio = 0.1; // seuil : 10 % de la médiane
  double degThreshold = std::max(medLen * degRatio, 1.0e-9);

  // Second passage : on écarte les ruptures problématiques et on calcule minLen sur les segments
  // conservés (sans le biais du tout premier segment qu'avait l'ancien calcul).
  std::vector<Segment> segments;
  double minLen = -1.0;
  int nDropped = 0;
  for (size_t s = 0; s < allSegs.size(); s++) {
    double len = lens[s];
    if (len < degThreshold) {
      const BreakEvent &e = breakEvents[allSegs[s].ev];
      std::cout << "  rupture problematique ecartee : t=" << e.time << "  cells " << e.a_ci << "-" << e.a_cj
                << "  len=" << len << " (< " << degThreshold << ")" << std::endl;
      nDropped++;
      continue;
    }
    if (minLen < 0.0 || len < minLen) {
      minLen = len;
    }
    segments.push_back(allSegs[s]);
  }
  std::cout << nDropped << " rupture(s) problematique(s) ecartee(s) / " << allSegs.size()
            << "  (mediane des longueurs = " << medLen << ")" << std::endl;

  double tol = (minLen > 0.0) ? 0.9 * minLen : 1.0e-6;
  double tol2 = tol * tol;

  // Second passage
  Crack crack;
  for (size_t i = 0; i < segments.size(); i++) {
    size_t na = addNode(crack, segments[i].a, tol2);
    size_t nb = addNode(crack, segments[i].b, tol2);
    if (na != nb) {
      crack.Conns.push_back({na, nb});
    }
  }
  for (size_t i = 0; i < crack.nodes.size(); i++) {
    if (crack.nodes[i].order > 3)
      crack.nodes[i].order = 3;
  }

  std::cout << "Crack : " << crack.nodes.size() << " nodes, " << crack.Conns.size() << " connections" << std::endl;
  saveCrackSVG(crack, "crack.svg");

  // Scinder le graphe global en fissures individuelles (composantes connexes).
  std::vector<Crack> cracks = splitCracks(crack);
  std::cout << "Nb de fissures (composantes connexes) : " << cracks.size() << std::endl;

  // Pour chaque fissure : chemin principal (plus long des plus courts chemins entre
  // extrémités) et sauvegarde SVG avec ce chemin surligné.
  for (size_t c = 0; c < cracks.size(); c++) {
    double L = 0.0;
    std::vector<size_t> path = mainPath(cracks[c], L);

    std::cout << "  fissure " << c << " : " << cracks[c].nodes.size() << " nodes, " << cracks[c].Conns.size()
              << " connections, chemin principal = " << path.size() << " noeuds, longueur = " << L << std::endl;

    char fname[256];
    snprintf(fname, sizeof(fname), "crack_%zu.svg", c);
    saveCrackSVG(cracks[c], fname, path.empty() ? nullptr : &path);
  }

  return 0;
}
