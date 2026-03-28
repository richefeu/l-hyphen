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

#pragma once

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "Lhyphen.hpp"
#include "null_size_t.hpp"

#include "AABB.hpp"
#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "glTools.hpp"
#include "message.hpp"
#include "triangulatePolygon.hpp"

Lhyphen Conf;
int confNum = 1;

AABB worldBox;

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
int show_nodes = 1;
int show_control_boxes = 0;
int show_background = 1;

// arrow sizes
double arrowSize = 0.15;
double arrowAngle = 0.35;
double vScale = 0.05;

// window sizes
int width = 800;
int height = 600;
float wh_ratio = (float)width / (float)height;
glTextZone textZone(1, &width, &height);

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
void drawControlBoxes();

void updateTextLine();

// Callback functions
void keyboard(GLFWwindow *window, int key, int scancode, int action, int mods);
void mouse_button(GLFWwindow *window, int button, int action, int mods);
void cursor_pos(GLFWwindow *window, double xpos, double ypos);
void display(GLFWwindow *window);
void reshape(GLFWwindow *window, int width, int height);

// Helper functions
void printHelp();
void fit_view(GLFWwindow *window);
void try_to_readConf(int num, Lhyphen &conf, int &OKNum);
