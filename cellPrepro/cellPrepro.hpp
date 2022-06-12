#ifndef CELLPREPRO_HPP
#define CELLPREPRO_HPP

#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "kwParser.hpp"

#include "./delaunator.hpp"

#include <fstream>
#include <iostream>
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
  // double a; // area
  std::vector<int> stack_pixels;
};

class CellPrepro {
private:
  kwParser parser;
  int nb_fill_colors;

public:
  std::vector<std::vector<int>> pixels;
  int lx, ly;
  int maxPGMGreyLevel;
  int MinPixelsPerCell;

  std::vector<Cell> cells;

  CellPrepro() {
    init_command_parser();

    nb_fill_colors = 0;
  }

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

    parser.kwMap["build_cells"] = __DO__() { build_cells(); };
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
      cells[i].x /= (double)cells[i].stack_pixels.size();
      cells[i].y /= (double)cells[i].stack_pixels.size();
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
    
    std::vector<double> coords;
    for (size_t i = 0; i < cells.size(); i++) {
      coords.push_back(cells[i].x);
      coords.push_back(cells[i].y);
    }

    // triangulation happens here
    delaunator::Delaunator d(coords);

    for(size_t i = 0; i < d.triangles.size(); i+=3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            d.coords[2 * d.triangles[i]],        //tx0
            d.coords[2 * d.triangles[i] + 1],    //ty0
            d.coords[2 * d.triangles[i + 1]],    //tx1
            d.coords[2 * d.triangles[i + 1] + 1],//ty1
            d.coords[2 * d.triangles[i + 2]],    //tx2
            d.coords[2 * d.triangles[i + 2] + 1] //ty2
        );
    }
    
    
  }

private:
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
    if (newColor == oldColor)
      return false; // if you don't do this: infinite loop!
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
