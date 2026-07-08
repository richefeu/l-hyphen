//  Copyright or © or Copr. l-hyphen
//
//  This software is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "Lhyphen.hpp"
#include <iostream>

#include "fileTool.hpp"

Lhyphen Conf;

/// Un évènement de rupture lu depuis breakHistory.txt.
/// On ne stocke que les références topologiques des deux côtés de l'interface (côté A = le contact,
/// côté B = son frère). L'interface rompue est reconstruite à l'affichage (cf. drawCrackPath) via
/// Conf.getPosition pour chaque côté, puis tracée comme un trait épais entre les deux points de
/// contact — comme le mode 'g' le fait pour les interfaces encore collées. Les positions étant lues
/// dans le conf affiché, le trait suit les cellules même si elles se sont déplacées depuis la rupture.
struct BreakEvent {
  double time;
  size_t a_ci, a_cj, a_in, a_jn; // côté A (ce contact)
  size_t b_ci, b_cj, b_in, b_jn; // côté B (le frère ; = côté A s'il n'y a pas de frère)
  double nrj;
};
std::vector<BreakEvent> breakEvents; ///< chargé une fois à l'ouverture du premier conf


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


// Renvoie true et remplit pos si (ci, cj, in, jn) est indexable dans le conf courant (garde-fou au cas
// où breakHistory.txt ne correspondrait pas au jeu de conf chargé), puis délègue à Conf.getPosition.
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



void saveCrackSVG(const Crack &crack, const char *fname, int CanvasWidth = 800) {
  if (crack.nodes.empty()) {
    std::cout << "saveCrackSVG : graphe vide, rien à écrire" << std::endl;
    return;
  }

  // Boîte englobante des noeuds.
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

  // Noeuds du graphe.
  for (size_t k = 0; k < crack.nodes.size(); k++) {
    double cx = crack.nodes[k].pos.x * vz.scalex + vz.x0;
    double cy = crack.nodes[k].pos.y * vz.scaley + vz.y0;
    
    std::cout << crack.nodes[k].order << std::endl;
    
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


// Ajoute un point au graphe en fusionnant 
size_t addNode(Crack &crack, const vec2r &p, double tol2){
  //std::cout << "tytytytytytytytyty" << std::endl;
  for (size_t k = 0; k < crack.nodes.size(); k++) {
    
    if (norm2(crack.nodes[k].pos - p) < tol2) {
      crack.nodes[k].order += 1;
      //std::cout << norm2(crack.nodes[k].pos - p) << "  " << tol2 << std::endl;  
      return k;
    } 
  }
  CrackNode N;
  N.pos = p; N.order = 1;
  crack.nodes.push_back(N);
  return crack.nodes.size() - 1;
}


int main(/*int argc, const char *argv[]*/) {

  Conf.loadCONF("conf0");
  readBreakHistory("breakHistory.txt");
  
  
  // Algo naif
  /*
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
  */
  
  // construction de la structure de crack entier
  

  struct Segment {
    vec2r a, b;
  };
  std::vector<Segment> segments;
  double minLen = -1.0;
  for (size_t i = 0; i < breakEvents.size(); i++) {
    Segment seg;
    if (!crackContactPos(breakEvents[i].a_ci, breakEvents[i].a_cj, breakEvents[i].a_in, breakEvents[i].a_jn, seg.a)) {
      continue;
    }
    if (!crackContactPos(breakEvents[i].b_ci, breakEvents[i].b_cj, breakEvents[i].b_in, breakEvents[i].b_jn, seg.b)) {
      continue;
    }
    double len = norm(seg.b - seg.a);
    if (minLen < 0.0 || (len < minLen && len > 1.0e-6)) {
      minLen = len;
    }
    segments.push_back(seg);
  }

  double tol = (minLen > 0.0) ? 0.9 * minLen : 1.0e-6;
  //double tol = 2*8.41e-07;
  double tol2 = tol * tol;

  // Second passage 
  Crack crack;
  for (size_t i = 0; i < segments.size(); i++) {
    size_t na = addNode(crack, segments[i].a, tol2);
    size_t nb = addNode(crack, segments[i].b, tol2);
    if (na != nb) {
      crack.Conns.push_back({na, nb});
      //crack.nodes[na].order += 1;
      //crack.nodes[nb].order += 1;
    }
  }
  for (size_t i = 0; i < crack.nodes.size(); i++){
    if (crack.nodes[i].order > 3) crack.nodes[i].order = 3; 
  }

  // il faudrait scinder le Crack en plusieurs Cracks (un crack est forcément entièrement connecté)




  std::cout << "Crack : " << crack.nodes.size() << " nodes, " << crack.Conns.size() << " connections" << std::endl;

  saveCrackSVG(crack, "crack.svg");

  return 0;
}
