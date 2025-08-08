#ifndef CELLPREPRO_HPP
#define CELLPREPRO_HPP

#include "ColorTable.hpp"
#include "convexHull.hpp"
#include "fastSort3.hpp"
#include "fileTool.hpp"
#include "kwParser.hpp"
#include "vec2.hpp"

#include "./delaunator.hpp"

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

// Compilation
// JYD : g++ -std=c++11 colorCells.cpp -I /Users/jyd/Documents/Codes/richefeu-toofus-6c99da3
// VR : g++ -O3 -std=c++11 colorCells.cpp -I /Users/vrichefeu/Documents/devel/TOOFUS
//
// Convertit jpg en pgm
// convert ./figure.jpg -compress none pgm:- | tr -s '\012\015' ' '  > figure.pgm

struct Cell {
  double x;
  double y;
  double Req;
  std::vector<int> stack_pixels;
  std::vector<int> neighbors;
  std::vector<vec2r> convex_hull;
  std::vector<vec2r> shifted;
};

class CellPrepro {
private:
  kwParser parser;
  int nb_fill_colors;
  std::string neighbor_strategy;
  double max_wall_width;

public:
  std::vector<std::vector<int>> pixels;
  int lx, ly;
  int maxPGMGreyLevel;
  int MinPixelsPerCell;
  int repair_cell_flag;

  std::vector<Cell> cells;

  CellPrepro() {
    init_command_parser();
    neighbor_strategy = "disks";
    max_wall_width = 50.0;
    nb_fill_colors = 0;
    repair_cell_flag = 0;
  }

private:
  void init_command_parser() {

    parser.kwMap["read_pgm"] = __DO__(file) {
      std::string pgm_file_name;
      file >> pgm_file_name;
      read_pgm(pgm_file_name.c_str());
    };

    parser.kwMap["write_pgm"] = __DO__(file) {
      std::string pgm_file_name;
      file >> pgm_file_name;
      write_pgm(pgm_file_name.c_str());
    };

    parser.kwMap["make_histo"] = __DO__() { make_histo(); };

    parser.kwMap["segmentation"] = __DO__(file) {
      int threshold;
      file >> threshold;
      segmentation(threshold);
    };

    parser.kwMap["fill_cells"] = __DO__() { fill_cells(); };

    parser.kwMap["create_cell_dataset"] = __DO__() { create_cell_dataset(); };
    parser.kwMap["save_cell_dataset"] = __DO__(file) {
      std::string filename;
      file >> filename;
      save_cell_dataset(filename.c_str());
    };

    parser.kwMap["add_centers_image"] = __DO__(file) {
      int size;
      file >> size;
      add_centers_image(size);
    };

    parser.kwMap["clean_cell_data"] = __DO__(file) {
      int minPixNumber;
      file >> minPixNumber;
      clean_cell_data(minPixNumber);
    };

    parser.kwMap["neighbor_strategy"] = __DO__(file) { file >> neighbor_strategy; };
    parser.kwMap["max_wall_width"] = __DO__(file) { file >> max_wall_width; };

    parser.kwMap["build_cells"] = __DO__() { build_cells(); };

    parser.kwMap["create_lhyphen_input"] = __DO__(file) {
      double pixSize, barWidth;
      file >> pixSize >> barWidth;
      create_lhyphen_input(pixSize, barWidth);
    };

    parser.kwMap["repair_cell"] = __DO__(file) { file >> repair_cell_flag; };
  }

  void prepro_commands(const char *name) {
    std::ifstream file(name);
    if (!file.is_open()) {
      std::cout << "Cannot read file " << name << '\n';
      return;
    }

    parser.parse(file);
  }

  // Read a PGM P2 file (grey levels)
  void read_pgm(const char *filename) {
    std::cout << "read_pgm " << filename << '\n';
    std::ifstream ifile(filename);
    if (!ifile) {
      std::cout << "Error! Could not find: " << std::endl << "> " << filename << std::endl;
      exit(0);
    }

    if (ifile) {
      std::string tmp;
      ifile >> tmp;
      std::cout << "Reading a " << tmp << " file" << std::endl;
      ifile >> lx >> ly;
      std::cout << "lx = " << lx << '\n';
      std::cout << "ly = " << ly << '\n';
      pixels.resize(lx);
      for (int i = 0; i < lx; i++) {
        pixels[i].resize(ly);
      }

      ifile >> maxPGMGreyLevel;
      std::cout << "Max PGM grey levels: " << maxPGMGreyLevel << std::endl;
      int greyLevel;
      for (int y = 0; y < ly; y++) {
        for (int x = 0; x < lx; x++) {
          ifile >> greyLevel;
          pixels[x][y] = greyLevel;
        }
      }
    }
    ifile.close();
  }

public:
  // Save a PGM P3 file (rgb colors)
  void write_pgm(const std::string &output_file) {
    std::cout << "write_pgm " << output_file << '\n';
    std::ofstream myfile(output_file, std::ios::out | std::ios::trunc);

    ColorTable CT;
    CT.setTableID(20); // reach integer index corresponds to a random color
    CT.setSize(nb_fill_colors + 1);
    CT.setMinMax(0.0, CT.getSize());
    CT.Rebuild();
    colorRGBA col;

    if (myfile) {
      myfile << "P3\n";
      myfile << lx << ' ' << ly << '\n';
      myfile << maxPGMGreyLevel << '\n';

      int r, v, b;

      for (int y = 0; y < ly; y++) {
        for (int x = 0; x < lx; x++) {
          if (pixels[x][y] == 0) { // paroi
            r = v = b = 0;
          } else if (pixels[x][y] > 0) { // couleur interne cellule
            CT.getRGB(pixels[x][y], &col);
            r = col.r;
            v = col.g;
            b = col.b;
          } else { // cellule sans floodfill
            r = v = b = maxPGMGreyLevel;
          }
          myfile << r << ' ' << v << ' ' << b << ' ';
          if (x % 5 == 0) {
            myfile << std::endl; // Les lignes ne doivent pas dépasser 70 caracteres
          }
        }
      }
    }
  }

  void make_histo() {
    std::cout << "make_histo\n";
    std::vector<int> count(maxPGMGreyLevel, 0);
    for (int y = 0; y < ly; y++) {
      for (int x = 0; x < lx; x++) {
        count[pixels[x][y]] += 1;
      }
    }
    std::ofstream histo("histo.txt");
    for (size_t i = 0; i < count.size(); i++) {
      histo << i << ' ' << count[i] << '\n';
    }
  }

  void segmentation(int threshold) {
    std::cout << "segmentation " << threshold << '\n';
    int greyLevel;
    for (int y = 0; y < ly; y++) {
      for (int x = 0; x < lx; x++) {
        greyLevel = pixels[x][y];
        if (greyLevel < threshold) {
          pixels[x][y] = 0; // C'est une paroi
        } else {
          pixels[x][y] = -1; // c'est du vide
        }
      }
    }
  }

  void fill_cells() {
    std::cout << "fill_cells\n";
    int newcol = 1;
    for (int x = 1; x < lx - 1; x++) {
      for (int y = 1; y < ly - 1; y++) {
        if (pixels[x][y] == -1) {
          if (floodFill8Stack(x, y, newcol, -1))
            newcol++;
        }
      }
    }
    nb_fill_colors = newcol - 1;
  }

  void add_centers_image(int size) {
    std::cout << "add_centers_image " << size << '\n';
    for (size_t i = 0; i < cells.size(); i++) {
      int x = cells[i].x;
      int y = cells[i].y;
      for (int xx = x - size; xx <= x + size; xx++) {
        for (int yy = y - size; yy <= y + size; yy++) {
          if (xx >= 0 && xx < lx && yy >= 0 && yy < ly) {
            pixels[xx][yy] = -2;
          }
        }
      }
    }
  }

  void clean_cell_data(int minPixNumber) {
    std::cout << "clean_cell_data " << minPixNumber << '\n';
    std::vector<Cell> clean_cells;
    for (size_t i = 0; i < cells.size(); i++) {
      if (cells[i].stack_pixels.size() > (size_t)minPixNumber) {
        // cells.erase(cells.begin() + i);
        clean_cells.push_back(cells[i]);
      }
    }
    cells.resize(clean_cells.size());
    for (size_t i = 0; i < clean_cells.size(); i++) {
      cells[i] = clean_cells[i];
    }
  }

  void create_cell_dataset() {
    std::cout << "create_cell_dataset\n";
    cells.clear();
    cells.resize(nb_fill_colors);
    for (int y = 0; y < ly; y++) {
      for (int x = 0; x < lx; x++) {
        int icell = pixels[x][y] - 1;
        if (icell < 0 || icell >= (int)cells.size())
          continue;
        cells[icell].x += x;
        cells[icell].y += y;
        cells[icell].stack_pixels.push_back(ly * x + y);
      }
    }

    for (size_t i = 0; i < cells.size(); i++) {
      double nbPix = (double)cells[i].stack_pixels.size();
      cells[i].x /= nbPix;
      cells[i].y /= nbPix;
      cells[i].Req = sqrt(nbPix / 3.14159);
    }
  }

  void save_cell_dataset(const char *filename) {
    std::cout << "save_cell_dataset " << filename << '\n';
    std::ofstream file(filename);
    file << "#id x y nbPixels\n";
    if (file) {
      for (size_t i = 0; i < cells.size(); i++) {
        file << i << ' ' << cells[i].x << ' ' << cells[i].y << ' ' << cells[i].stack_pixels.size() << '\n';
      }
    }
  }

  void build_cells() {
    std::cout << "build_cells\n";

    if (neighbor_strategy == "delaunay") {
      find_neighbors_delaunay();
    } else if (neighbor_strategy == "disks") {
      find_neighbors_disks();
    }

    for (size_t i = 0; i < cells.size(); i++) {
      std::vector<vec2r> tmp;
      tmp.push_back(vec2r(0.0, 0.0));
      tmp.push_back(vec2r(lx, 0.0));
      tmp.push_back(vec2r(lx, ly));
      tmp.push_back(vec2r(0.0, ly));
      cells[i].convex_hull = convexHull(tmp);
    }

    std::ofstream pp("cuts.txt");
    for (size_t i = 0; i < cells.size(); i++) {
      for (size_t k = 0; k < cells[i].neighbors.size(); k++) {
        size_t j = cells[i].neighbors[k];

        // vec2r pos_plan(0.5 * (cells[i].x + cells[j].x), 0.5 * (cells[i].y + cells[j].y));
        vec2r pos_plan = find_wall_pos(i, j);

        vec2r normal(cells[i].x - cells[j].x, cells[i].y - cells[j].y);

        pp << pos_plan << '\n';
        normal.normalize();
        // std::cout << "normal = " << normal << '\n';
        cut_polyg_with_plan(cells[i].convex_hull, pos_plan, normal);
      }
    }

    std::ofstream file("voro.txt");
    for (size_t i = 0; i < cells.size(); i++) {
      if (cells[i].convex_hull.empty())
        continue;
      for (size_t k = 0; k < cells[i].convex_hull.size(); k++) {
        file << cells[i].convex_hull[k] << '\n';
      }
      file << cells[i].convex_hull[0] << "\n\n";
    }
  }

  void place(double Px, double Py, double w, double n1x, double n1y, double n2x, double n2y, double &newPx,
             double &newPy) {
    double A1 = (Px + 0.5 * w * n1x) * n1x + (Py + 0.5 * w * n1y) * n1y;
    double A2 = (Px + 0.5 * w * n2x) * n2x + (Py + 0.5 * w * n2y) * n2y;
    double cross = n1x * n2y - n1y * n2x;
    newPx = (n2y * A1 - n1y * A2) / cross;
    newPy = (n1x * A2 - n2x * A1) / cross;
  }

  bool repair_cell(size_t i) {
    if (cells[i].shifted.size() < 4)
      return true;

    int nbVertices = cells[i].shifted.size();
    for (int k0 = 0; k0 < nbVertices; k0++) {
      int k1 = k0 + 1;
      if (k1 >= nbVertices) {
        k1 -= nbVertices;
      }
      int k2 = k1 + 1;
      if (k2 >= nbVertices) {
        k2 -= nbVertices;
      }
      int k3 = k2 + 1;
      if (k3 >= nbVertices) {
        k3 -= nbVertices;
      }

      vec2r A = cells[i].shifted[k0];
      vec2r B = cells[i].shifted[k1];
      vec2r C = cells[i].shifted[k2];
      vec2r D = cells[i].shifted[k3];
      double alphaAB = 0.0;
      if (intersectSeg(A, B, C, D, alphaAB)) {
        cells[i].shifted[k2] = A + alphaAB * (B - A);
        cells[i].shifted.erase(cells[i].shifted.begin() + k1);
        return false;
      }
    }
    return true;
  }

  void create_lhyphen_input(double pixSize, double barWidth) {
    std::ofstream file("sample.txt");

    for (size_t i = 0; i < cells.size(); i++) {
      if (cells[i].convex_hull.empty()) {
        continue;
      }

      int nbVertices = (int)cells[i].convex_hull.size();
      for (int k = 0; k < nbVertices; k++) {
        int knext = k + 1;
        if (knext >= nbVertices) {
          knext = 0;
        }
        int kprev = k - 1;
        if (kprev < 0) {
          kprev = nbVertices - 1;
        }

        vec2r t = cells[i].convex_hull[knext] - cells[i].convex_hull[k];
        t.y *= -1.0;
        t.normalize();
        vec2r nnext(t.y, -t.x);
        t = cells[i].convex_hull[k] - cells[i].convex_hull[kprev];
        t.y *= -1.0;
        t.normalize();
        vec2r nprev(t.y, -t.x);

        double Sx, Sy;
        double Px = pixSize * cells[i].convex_hull[k].x;
        double Py = pixSize * (ly - cells[i].convex_hull[k].y);
        place(Px, Py, barWidth, nprev.x, nprev.y, nnext.x, nnext.y, Sx, Sy);
        cells[i].shifted.emplace_back(Sx, Sy);
      }
    }

    //
    if (repair_cell_flag > 0) {
      for (size_t i = 0; i < cells.size(); i++) {
        while (!repair_cell(i)) {
          // ?????
        }
      }
    }

    //
    for (size_t i = 0; i < cells.size(); i++) {
      if (cells[i].convex_hull.empty()) {
        continue;
      }
      for (int k = 0; k < (int)cells[i].shifted.size(); k++) {
        file << cells[i].shifted[k].x << ' ' << cells[i].shifted[k].y << '\n';
      }
      file << cells[i].shifted[0].x << ' ' << cells[i].shifted[0].y << "\n\n";
    }

    std::ofstream nodefile("nodeFile.txt");
    int ID = 0;
    for (size_t i = 0; i < cells.size(); i++) {
      if (cells[i].convex_hull.empty()) {
        continue;
      }

      for (int k = 0; k < (int)cells[i].shifted.size(); k++) {
        nodefile << cells[i].shifted[k].x << ' ' << cells[i].shifted[k].y << ' ' << ID << '\n';
      }
      ID++;
    }
  }

private:
  void find_neighbors_disks() {
    for (size_t i = 0; i < cells.size(); i++) {
      for (size_t j = i + 1; j < cells.size(); j++) {
        double dx = cells[j].x - cells[i].x;
        double dy = cells[j].y - cells[i].y;
        double dist = sqrt(dx * dx + dy * dy) - (cells[i].Req + cells[j].Req);
        if (dist <= max_wall_width) {
          cells[i].neighbors.push_back(j);
          cells[j].neighbors.push_back(i);
        }
      }
    }
  }

  void find_neighbors_delaunay() {
    std::vector<double> coords;
    for (size_t i = 0; i < cells.size(); i++) {
      coords.push_back(cells[i].x);
      coords.push_back(cells[i].y);
    }

    // triangulation happens here
    delaunator::Delaunator d(coords);

    std::set<std::pair<int, int>> set_close_cells;

    std::ofstream tt("triang.txt");

    for (size_t i = 0; i < d.triangles.size(); i += 3) {
      int i0 = d.triangles[i];
      int i1 = d.triangles[i + 1];
      int i2 = d.triangles[i + 2];
      fastSort3(i0, i1, i2);
      // std::cout << i0 << ' ' << i1 << ' ' << i2 << '\n';

      tt << cells[i0].x << ' ' << cells[i0].y << '\n';
      tt << cells[i1].x << ' ' << cells[i1].y << '\n';
      tt << cells[i2].x << ' ' << cells[i2].y << '\n';
      tt << cells[i0].x << ' ' << cells[i0].y << '\n';
      tt << '\n';

      set_close_cells.insert(std::make_pair(i0, i1));
      set_close_cells.insert(std::make_pair(i1, i2));
      set_close_cells.insert(std::make_pair(i0, i2));
    }

    std::ofstream ll("lines.txt");
    for (auto p : set_close_cells) {
      int i0 = p.first;
      int i1 = p.second;
      ll << cells[i0].x << ' ' << cells[i0].y << '\n';
      ll << cells[i1].x << ' ' << cells[i1].y << '\n';
      ll << '\n';

      cells[i0].neighbors.push_back(i1);
      cells[i1].neighbors.push_back(i0);
    }
  }

  vec2r find_wall_pos(size_t i, size_t j) {

    vec2r branch(cells[j].x - cells[i].x, cells[j].y - cells[i].y);
    vec2r u = branch;
    double len = u.normalize();

    double iprojmax = 0.0;
    for (size_t k = 0; k < cells[i].stack_pixels.size(); k++) {
      int p = cells[i].stack_pixels[k];
      double x = (double)(p / ly);
      double y = (double)(p % ly);
      vec2r a(x - cells[i].x, y - cells[i].y);
      double proj = a * u;
      if (proj > iprojmax) {
        iprojmax = proj;
      }
    }
    if (iprojmax > len) {
      iprojmax = 0.5 * len;
    }

    double jprojmin = len;
    for (size_t k = 0; k < cells[j].stack_pixels.size(); k++) {
      int p = cells[j].stack_pixels[k];
      double x = (double)(p / ly);
      double y = (double)(p % ly);
      vec2r a(x - cells[i].x, y - cells[i].y);
      double proj = a * u;
      if (proj < jprojmin) {
        jprojmin = proj;
      }
    }
    if (jprojmin < 0.0) {
      jprojmin = 0.5 * len;
    }

    vec2r mid(cells[i].x, cells[i].y);
    mid += 0.5 * (iprojmax + jprojmin) * u;
    return mid;
  }

  // alpha = 0.0 means intersect in A
  // alpha = 1.0 means intersect in B
  // alpha in [0.0; 1.0] -> intersect segment AB
  bool intersect(const vec2r &A, const vec2r &B, const vec2r &pl_pos, const vec2r &pl_normal, double &alpha) {
    vec2r PA = A - pl_pos;
    vec2r u = B - A; // do not normalize!
    double un = u * pl_normal;
    if (fabs(un) < 1e-20) {
      return false;
    }
    double PAn = PA * pl_normal;
    alpha = -PAn / un;
    return true;
  }

  bool intersectSeg(const vec2r &A, const vec2r &B, const vec2r &C, const vec2r &D, double &alphaAB) {
    vec2r t = C - D;
    vec2r n(t.y, -t.x);
    n.normalize();
    intersect(A, B, C, n, alphaAB);
    if (alphaAB < 0.0 || alphaAB > 1.0) {
      return false;
    }
    t = B - A;
    n.set(t.y, -t.x);
    n.normalize();
    double alphaCD = 0.0;
    intersect(C, D, A, n, alphaCD);
    if (alphaCD < 0.0 || alphaCD > 1.0) {
      return false;
    }

    return true;
  }

  void cut_polyg_with_plan(std::vector<vec2r> &polyg, vec2r &pos_plan, vec2r &normal) {
    if (polyg.size() < 3) {
      polyg.clear();
      return;
    }
    std::vector<vec2r> hull0 = convexHull(polyg);
    std::vector<vec2r> lstPts;
    for (size_t i0 = 0; i0 < hull0.size(); i0++) {
      size_t i1 = i0 + 1;
      if (i1 == hull0.size()) {
        i1 = 0;
      }
      double alpha;
      if (!intersect(hull0[i0], hull0[i1], pos_plan, normal, alpha)) {
        return;
      }
      if (alpha > 0.0 && alpha < 1.0) {
        vec2r newPt = hull0[i0] + alpha * (hull0[i1] - hull0[i0]);
        lstPts.push_back(newPt);
      }
    }

    for (size_t i = 0; i < hull0.size(); i++) {
      vec2r a = hull0[i] - pos_plan;
      double sgn = a * normal;
      if (sgn >= 0.0) {
        lstPts.push_back(hull0[i]);
      }
    }
    if (lstPts.size() < 3) {
      return;
    }
    polyg.clear();
    polyg = convexHull(lstPts);
  }

  //
  // Floodfill algorithm (adapted from http://lodev.org/cgtutor/floodfill.html)
  //

  std::vector<int> stack_pixels;

  bool pop(int &x, int &y) {
    if (!stack_pixels.empty()) {
      int p = stack_pixels.back();
      x = p / ly;
      y = p % ly;
      stack_pixels.pop_back();
      return true;
    } else {
      return false;
    }
  }

  void push(int x, int y) { stack_pixels.push_back(ly * x + y); }

  // 8-way floodfill using our own stack routines
  const int dx[8] = {0, 1, 1, 1, 0, -1, -1, -1}; // relative neighbor x coordinates
  const int dy[8] = {-1, -1, 0, 1, 1, 1, 0, -1}; // relative neighbor y coordinates

  bool floodFill8Stack(int x, int y, int newColor, int oldColor) {
    if (newColor == oldColor) {
      return false;
    } // if you don't do this: infinite loop!
    
    // empty stack
    stack_pixels.clear();

    push(x, y);

    while (pop(x, y)) {
      pixels[x][y] = newColor;
      for (int i = 0; i < 8; i++) {
        int nx = x + dx[i];
        int ny = y + dy[i];
        if (nx >= 0 && nx < lx && ny >= 0 && ny < ly && pixels[nx][ny] == oldColor) {
          push(nx, ny);
        }
      }
    }
    return true; // floodfill has painted a zone
  }
};

#endif /* end of include guard: CELLPREPRO_HPP */
