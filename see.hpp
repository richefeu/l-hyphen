#ifndef SEE_HPP
#define SEE_HPP

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "Lhyphen.hpp"
#include "null_size_t.hpp"

//#include "ColorTable.hpp"

Lhyphen Conf;
int confNum = 1;

struct AABB {
  double xmin,xmax, ymin, ymax;
};

AABB worldBox;

int main_window;

// flags
int show_cells = 1;

// arrow sizes
double arrowSize = 0.0005;
double arrowAngle = 0.7;
double vScale = 10.0;


// window sizes
int width = 800;
int height = 600;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int display_mode = 0;  // sample or slice rotation
int mouse_start[2];

// drawing functions
void drawBar(size_t ci, size_t in, size_t jn, double radius);
void drawCells();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void printHelp();
void fit_view();
bool fileExists(const char* fileName);
void try_to_readConf(int num, Lhyphen& conf, int& OKNum);

#endif /* end of include guard: SEE_HPP */
