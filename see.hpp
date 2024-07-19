#ifndef SEE_HPP
#define SEE_HPP

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "Lhyphen.hpp"
#include "null_size_t.hpp"

#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "triangulatePolygon.hpp"

Lhyphen Conf;
int confNum = 1;

struct AABB_see {
  double xmin, xmax, ymin, ymax;
};

AABB_see worldBox;

ColorTable BarRedTable;
ColorTable BarBlueTable;
ColorTable NodeRedTable;
ColorTable NodeBlueTable;

int main_window;

// flags
int show_cells = 1;
int show_glue_points = 0;
int show_bar_colors = 0;
int show_inter_cells_forces = 0;
int show_pressure = 0;

// arrow sizes
double arrowSize = 0.15;
double arrowAngle = 0.35;
double vScale = 0.05;

// window sizes
int width = 800;
int height = 600;
float wh_ratio = (float)width / (float)height;

// Miscellaneous global variables
enum class MouseMode { NOTHING, ROTATION, ZOOM, PAN };
MouseMode mouse_mode = MouseMode::NOTHING;

int display_mode = 0; // sample or slice rotation
int mouse_start[2];

// drawing functions
void drawCircle(double xc, double yc, double radius, int nbDiv = 18);
void drawBar(size_t ci, size_t in, size_t jn, double radius, color4f &BarColor, color4f &NodeColor);
void drawCells();
void drawGluePoints();
void arrow(double x0, double y0, double x1, double y1);
void drawForceActionReaction(const Neighbor &InterIt, vec2r &pos);
void drawForces();
void drawPressure();
//void drawColorTable(ColorTable &CT, const char * title);

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
bool fileExists(const char *fileName);
void try_to_readConf(int num, Lhyphen &conf, int &OKNum);

#endif /* end of include guard: SEE_HPP */
