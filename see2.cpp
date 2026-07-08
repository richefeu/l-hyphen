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

#include "see2.hpp"

void printHelp() {
  std::cout << std::endl;
  std::cout << "Commandes:" << std::endl;
  std::cout << "a           show/hide control area boxes" << std::endl;
  std::cout << "b           colorize the cell bars" << std::endl;
  std::cout << "c           show/hide the cells" << std::endl;
  std::cout << "f           show/hide the forces" << std::endl;
  std::cout << "g           show/hide the glue points" << std::endl;
  std::cout << "r           show/hide the crack path (broken links up to current time)" << std::endl;
  std::cout << "h           print this help" << std::endl;
  std::cout << "n           show/hide cell contours" << std::endl;
  std::cout << "v           show/hide nodes (points)" << std::endl;
  std::cout << "p           show/hide pressure" << std::endl;
  std::cout << "q           quit" << std::endl;
  std::cout << "s/S         force-chain lines thinner/thicker" << std::endl;
  std::cout << "t/T         force filter lower/higher (show more/fewer chains)" << std::endl;
  std::cout << "z/Z         zoom in/out" << std::endl;
  std::cout << "->          load next configuration file" << std::endl;
  std::cout << "<-          load previous configuration file" << std::endl;
  std::cout << "Shift+<-    jump to conf 0" << std::endl;
  std::cout << "=           fit the view" << std::endl;
  std::cout << "x           save screenshot.png" << std::endl;
  std::cout << "Shift+x     batch screenshots of all confs" << std::endl;
  std::cout << "space       save options to see2-options.toml" << std::endl;
  std::cout << "Shift+space reload options from see2-options.toml" << std::endl;
  std::cout << std::endl;
}

void keyboard(GLFWwindow *window, int key, int /*scancode*/, int action, int mods) {
  if (action != GLFW_PRESS)
    return;

  switch (key) {

  case GLFW_KEY_Q: { // a
    show_control_boxes = 1 - show_control_boxes;
    textZone.addLine("show_control_boxes = %d", show_control_boxes);
  } break;

  case GLFW_KEY_B: {
    show_bar_colors = 1 - show_bar_colors;
    textZone.addLine("show_bar_colors = %d", show_bar_colors);
  } break;

  case GLFW_KEY_C: {
    show_cells = 1 - show_cells;
    textZone.addLine("show_cells = %d", show_cells);
  } break;

  case GLFW_KEY_F: {
    show_inter_cells_forces = 1 - show_inter_cells_forces;
    textZone.addLine("show_forces = %d", show_inter_cells_forces);
  } break;

  case GLFW_KEY_G: {
    show_glue_points = 1 - show_glue_points;
    textZone.addLine("show_glue = %d", show_glue_points);
  } break;

  case GLFW_KEY_R: {
    show_crack_path = 1 - show_crack_path;
    textZone.addLine("show_crack_path = %d (%zu events)", show_crack_path, breakEvents.size());
  } break;

  case GLFW_KEY_H: {
    printHelp();
  } break;

  case GLFW_KEY_SPACE: {
    if (mods == GLFW_MOD_SHIFT) {
      readTomlOptions();
      glfwSetWindowSize(window, width, height);
      reshape(window, width, height);
      std::cout << "Loaded 'see2-options.toml'" << std::endl;
      textZone.addLine("Loaded 'see2-options.toml'");
    } else {
      saveTomlOptions();
      std::cout << "Saved 'see2-options.toml'" << std::endl;
      textZone.addLine("Saved 'see2-options.toml'");
    }
  } break;

  case GLFW_KEY_N: {
    show_contours = 1 - show_contours;
    textZone.addLine("show_contours = %d", show_contours);
  } break;

  case GLFW_KEY_V: {
    show_nodes = 1 - show_nodes;
    textZone.addLine("show_nodes = %d", show_nodes);
  } break;

  case GLFW_KEY_P: {
    show_pressure = 1 - show_pressure;
    textZone.addLine("show_pressure = %d", show_pressure);
  } break;

  case GLFW_KEY_A: { // q
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  } break;

  case GLFW_KEY_S: {
    if (mods == GLFW_MOD_SHIFT) {
      fnWidthFactor *= 1.1;
    } else {
      fnWidthFactor *= 0.9;
    }
    textZone.addLine("fnWidthFactor = %g", fnWidthFactor);
  } break;

  case GLFW_KEY_T: {
    if (mods == GLFW_MOD_SHIFT) {
      forceFilter *= 1.2; // seuil plus haut -> ne garde que les chaînes les plus fortes
    } else {
      forceFilter *= 0.8; // seuil plus bas -> montre plus de contacts
    }
    if (forceFilter < 0.0) forceFilter = 0.0;
    textZone.addLine("forceFilter = %g x mean|fn|", forceFilter);
  } break;

  case GLFW_KEY_W: { // z
    if (mods == GLFW_MOD_SHIFT) {
      double dy = worldBox.max.y - worldBox.min.y;
      double ddy = -0.2 * dy;
      double ddx = -0.2 * dy;
      worldBox.min.x -= ddx;
      worldBox.max.x += ddx;
      worldBox.min.y -= ddy;
      worldBox.max.y += ddy;
    } else {
      double dy = worldBox.max.y - worldBox.min.y;
      double ddy = 0.2 * dy;
      double ddx = 0.2 * dy;
      worldBox.min.x -= ddx;
      worldBox.max.x += ddx;
      worldBox.min.y -= ddy;
      worldBox.max.y += ddy;
    }
    reshape(window, width, height);
  } break;

  case GLFW_KEY_X: {
    if (mods == GLFW_MOD_SHIFT) {
      do {
        char filename[256];
        snprintf(filename, 256, "screenshot%d.png", confNum);
        display(window);
        captureScreenshot(filename);
        std::cout << filename << " saved" << std::endl;
        updateTextLine();
      } while (try_to_readConf(confNum + 1, Conf, confNum));
    } else {
      display(window);
      char filename[256];
      snprintf(filename, 256, "screenshot%d.png", confNum);
      captureScreenshot(filename);
      std::cout << filename << " saved" << std::endl;
      textZone.addLine("%s saved", filename);
    }
  } break;

  case GLFW_KEY_LEFT: {
    if (mods == GLFW_MOD_SHIFT) {
      try_to_readConf(0, Conf, confNum);
    } else if (confNum > 0) {
      try_to_readConf(confNum - 1, Conf, confNum);
    }
    updateTextLine();
  } break;

  case GLFW_KEY_RIGHT: {
    try_to_readConf(confNum + 1, Conf, confNum);
    updateTextLine();
  } break;

  case GLFW_KEY_SLASH: { // '='
    fit_view(window);
  } break;

  case GLFW_KEY_UP: {
    textZone.increase_nbLine();
  } break;

  case GLFW_KEY_DOWN: {
    textZone.decrease_nbLine();
  } break;
  };

  needsRedraw = true;
}

void updateTextLine() {
  textZone.addLine("Conf%d, time = %g", confNum, Conf.t);
  if (g_window) {
    char title[256];
    snprintf(title, 256, "see2 — conf%d  t = %g", confNum, Conf.t);
    glfwSetWindowTitle(g_window, title);
  }
}

void mouse_button(GLFWwindow *window, int button, int action, int mods) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);

  if (action == GLFW_RELEASE) {
    mouse_mode = MouseMode::NOTHING;
  } else if (action == GLFW_PRESS) {
    mouse_start[0] = static_cast<int>(x);
    mouse_start[1] = static_cast<int>(y);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
      if (mods == GLFW_MOD_SHIFT) {
        mouse_mode = MouseMode::PAN;
      } else {
        mouse_mode = MouseMode::ROTATION;
      }
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
      mouse_mode = MouseMode::ZOOM;
    }
  }

  needsRedraw = true;
}

void cursor_pos(GLFWwindow *window, double xpos, double ypos) {
  if (mouse_mode == MouseMode::NOTHING) {
    return;
  }

  double dx = (xpos - mouse_start[0]) / static_cast<double>(width);
  double dy = (ypos - mouse_start[1]) / static_cast<double>(height);

  switch (mouse_mode) {
  case MouseMode::ZOOM: {
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    double ddx = (worldBox.max.x - worldBox.min.x) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x += ddx;
    worldBox.min.y -= ddy;
    worldBox.max.y += ddy;
  } break;

  case MouseMode::PAN: {
    double Lx = worldBox.max.x - worldBox.min.x;
    double Ly = worldBox.max.y - worldBox.min.y;
    double L = (Lx + Ly);
    double ddx = L * dx;
    double ddy = L * dy;

    worldBox.min.x -= ddx;
    worldBox.max.x -= ddx;
    worldBox.min.y += ddy;
    worldBox.max.y += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = static_cast<int>(xpos);
  mouse_start[1] = static_cast<int>(ypos);

  reshape(window, width, height);
  needsRedraw = true;
}

void display(GLFWwindow *window) {
  glTools::clearBackground(show_background, bottom_r, bottom_g, bottom_b, top_r, top_g, top_b);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (show_pressure) {
    drawPressure();
  }
  if (show_cells) {
    drawCells();
  }
  if (show_glue_points) {
    drawGluePoints();
  }
  if (show_crack_path) {
    drawCrackPath();
  }
  if (show_inter_cells_forces) {
    drawForces();
  }

  if (show_control_boxes) {
    drawControlBoxes();
  }

  textZone.draw();

  glFlush();
  glfwSwapBuffers(window);
}

void fit_view(GLFWwindow *window) {
  worldBox.min.x = Conf.xmin;
  worldBox.max.x = Conf.xmax;
  worldBox.min.y = Conf.ymin;
  worldBox.max.y = Conf.ymax;
  reshape(window, width, height);
}

void reshape(GLFWwindow *window, int w, int h) {
  if (window) {
    glfwGetFramebufferSize(window, &w, &h);
  }

  width = w;
  height = h;

  double left = worldBox.min.x;
  double right = worldBox.max.x;
  double bottom = worldBox.min.y;
  double top = worldBox.max.y;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;

  if (worldW > worldH) {
    worldH = worldW * ((GLfloat)height / (GLfloat)width);
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * ((GLfloat)width / (GLfloat)height);
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);
}

void drawCircle(double xc, double yc, double radius, int nbDiv) {
  glBegin(GL_POLYGON);
  double da = 2.0 * M_PI / (double)nbDiv;
  for (double a = 0.0; a < 2.0 * M_PI; a += da) {
    glVertex2d(xc + radius * cos(a), yc + radius * sin(a));
  }
  glEnd();
}

/**
 * Draws a bar between two nodes in a given cell.
 *
 * @param ci        The index of the cell
 * @param i         The index of the first node
 * @param j         The index of the second node
 * @param radius    The radius of the bar
 * @param BarColor  The color of the bar
 * @param NodeColor The color of the nodes
 */
void drawBar(size_t ci, size_t i, size_t j, double radius, color4f &BarColor, color4f &NodeColor) {
  if (i == null_size_t || j == null_size_t) {
    return;
  }

  double xi = Conf.cells[ci].nodes[i].pos.x;
  double yi = Conf.cells[ci].nodes[i].pos.y;
  double xj = Conf.cells[ci].nodes[j].pos.x;
  double yj = Conf.cells[ci].nodes[j].pos.y;

  double nxij = xj - xi;
  double nyij = yj - yi;
  double nij = sqrt(nxij * nxij + nyij * nyij);
  nxij /= nij;
  nyij /= nij;
  double txij = -nyij;
  double tyij = nxij;

  glColor4f(BarColor.r, BarColor.g, BarColor.b, 1.0f);
  glLineWidth(2.0f);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // draw the bar inclined rectangle
  glBegin(GL_POLYGON);
  glVertex2d(xi - radius * txij, yi - radius * tyij);
  glVertex2d(xj - radius * txij, yj - radius * tyij);
  glVertex2d(xj + radius * txij, yj + radius * tyij);
  glVertex2d(xi + radius * txij, yi + radius * tyij);
  glEnd();

  glColor4f(NodeColor.r, NodeColor.g, NodeColor.b, 1.0f);

  glBegin(GL_POLYGON);
  for (double a = 0.0; a < 2.0 * M_PI; a += M_PI / 18.0) {
    glVertex2d(xi + radius * cos(a), yi + radius * sin(a));
  }
  glEnd();

  if (Conf.cells[ci].nodes[j].nextNode == null_size_t || Conf.cells[ci].nodes[j].nextNode == 1) {
    glBegin(GL_POLYGON);
    for (double a = 0.0; a < 2.0 * M_PI; a += M_PI / 18.0) {
      glVertex2d(xj + radius * cos(a), yj + radius * sin(a));
    }
    glEnd();
  }
}

/**
 * Draws the pressure of the cells in the Conf object.
 */
void drawPressure() {
  double pmax = -1.0e20;
  double pmin = 1.0e20;
  for (size_t i = 0; i < Conf.cells.size(); ++i) {
    if (Conf.cells[i].close == false) {
      continue;
    }
    if (Conf.cells[i].p_int < pmin) {
      pmin = Conf.cells[i].p_int;
    }
    if (Conf.cells[i].p_int > pmax) {
      pmax = Conf.cells[i].p_int;
    }
  }

  color4f col;
  ColorTable pTable;
  pTable.setTableID(16);
  pTable.setSwap(true);
  pTable.setMinMax((float)pmin, (float)pmax);
  std::cout << "pmin = " << pmin << '\n';
  std::cout << "pmax = " << pmax << '\n';

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  for (size_t i = 0; i < Conf.cells.size(); ++i) {
    if (Conf.cells[i].close == false)
      continue;

    pTable.getColor4f((float)Conf.cells[i].p_int, &col);
    glColor3f(col.r, col.g, col.b);

    std::vector<vec2r> contour;
    for (size_t n = 0; n < Conf.cells[i].nodes.size(); ++n) {
      contour.push_back(Conf.cells[i].nodes[n].pos);
    }
    std::vector<int> result;
    TriangulatePolygon::Process(contour, result);
    glBegin(GL_TRIANGLES);
    for (size_t s = 0; s < result.size(); s += 3) {
      vec2r p0 = contour[result[s]];
      vec2r p1 = contour[result[s + 1]];
      vec2r p2 = contour[result[s + 2]];

      glVertex2d(p0.x, p0.y);
      glVertex2d(p1.x, p1.y);
      glVertex2d(p2.x, p2.y);
    }
    glEnd();
  }
}

/**
 * Draws the cells on the screen.
 */
void drawCells() {
  glLineWidth(2.0f);

  color4f BarColor, NodeColor;
  BarColor.r = NodeColor.r = 0.4f;
  BarColor.g = NodeColor.g = 0.8f;
  BarColor.b = NodeColor.b = 1.0f;
  BarColor.a = NodeColor.a = 1.0f;

  if (show_bar_colors) {
    double fnredmax = 0.0;
    double fnbluemax = 0.0;
    for (size_t i = 0; i < Conf.cells.size(); ++i) {
      for (size_t b = 0; b < Conf.cells[i].bars.size(); ++b) {
        if (Conf.cells[i].bars[b].fn >= 0.0) {
          fnredmax = (fnredmax > Conf.cells[i].bars[b].fn) ? fnredmax : Conf.cells[i].bars[b].fn;
        } else {
          fnbluemax = std::max(fnbluemax, -Conf.cells[i].bars[b].fn);
        }
      }
    }
    BarRedTable.setMinMax(0.0f, (float)fnredmax);
    BarBlueTable.setMinMax(0.0f, (float)fnbluemax);
  }

  for (size_t i = 0; i < Conf.cells.size(); ++i) {
    for (size_t b = 0; b < Conf.cells[i].bars.size(); ++b) {

      if (show_bar_colors) {
        if (Conf.cells[i].bars[b].fn > 0.0) {
          BarRedTable.getColor4f((float)(Conf.cells[i].bars[b].fn), &BarColor);
        } else if (Conf.cells[i].bars[b].fn < 0.0) {
          BarBlueTable.getColor4f((float)(-Conf.cells[i].bars[b].fn), &BarColor);
        } else {
          BarColor.r = 0.0f;
          BarColor.g = 1.0f;
          BarColor.b = 0.0f;
          BarColor.a = 1.0f;
        }
      }

      size_t in = Conf.cells[i].bars[b].i;
      size_t jn = Conf.cells[i].bars[b].j;
      drawBar(i, in, jn, Conf.cells[i].radius, BarColor, NodeColor);
    }

    if (show_contours) {
      glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

      // Draw a CCW arc from a_start to a_end (sweeping CCW, expanding range if needed)
      auto drawArcCCW = [&](double xn, double yn, double r, double a_start, double a_end) {
        while (a_end < a_start) a_end += 2.0 * M_PI;
        int nSteps = std::max(2, (int)std::ceil((a_end - a_start) / (M_PI / 18.0)));
        double da = (a_end - a_start) / nSteps;
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k <= nSteps; k++) {
          double a = a_start + k * da;
          glVertex2d(xn + r * std::cos(a), yn + r * std::sin(a));
        }
        glEnd();
      };

      bool isClosed = Conf.cells[i].close;

      for (size_t n = 0; n < Conf.cells[i].nodes.size(); ++n) {
        size_t pv = Conf.cells[i].nodes[n].prevNode;
        size_t nx = Conf.cells[i].nodes[n].nextNode;

        double xn = Conf.cells[i].nodes[n].pos.x;
        double yn = Conf.cells[i].nodes[n].pos.y;
        double r  = Conf.cells[i].radius;

        if (isClosed) {
          // Exterior contour only (both prev and next must exist)
          if (pv == null_size_t || nx == null_size_t) continue;

          double d1x = xn - Conf.cells[i].nodes[pv].pos.x;
          double d1y = yn - Conf.cells[i].nodes[pv].pos.y;
          double theta1 = std::atan2(d1y, d1x);

          double d2x = Conf.cells[i].nodes[nx].pos.x - xn;
          double d2y = Conf.cells[i].nodes[nx].pos.y - yn;
          double theta2 = std::atan2(d2y, d2x);

          double delta = theta2 - theta1;
          while (delta >  M_PI) delta -= 2.0 * M_PI;
          while (delta < -M_PI) delta += 2.0 * M_PI;

          if (delta >= 0.0) {
            // Convex node: exterior arc (right side, CCW)
            drawArcCCW(xn, yn, r, theta1 - M_PI * 0.5, theta2 - M_PI * 0.5);
          } else {
            // Concave node: straight segment connecting exterior tangent endpoints
            glBegin(GL_LINES);
            glVertex2d(xn + r * std::cos(theta1 - M_PI * 0.5), yn + r * std::sin(theta1 - M_PI * 0.5));
            glVertex2d(xn + r * std::cos(theta2 - M_PI * 0.5), yn + r * std::sin(theta2 - M_PI * 0.5));
            glEnd();
          }
        } else {
          // Open cell: draw both sides of the bar chain
          if (pv != null_size_t && nx != null_size_t) {
            // Interior node
            double d1x = xn - Conf.cells[i].nodes[pv].pos.x;
            double d1y = yn - Conf.cells[i].nodes[pv].pos.y;
            double theta1 = std::atan2(d1y, d1x);

            double d2x = Conf.cells[i].nodes[nx].pos.x - xn;
            double d2y = Conf.cells[i].nodes[nx].pos.y - yn;
            double theta2 = std::atan2(d2y, d2x);

            double delta = theta2 - theta1;
            while (delta >  M_PI) delta -= 2.0 * M_PI;
            while (delta < -M_PI) delta += 2.0 * M_PI;

            if (delta > 0.0) {
              // CCW (left) turn: exterior arc on right side, straight segment on concave left side
              drawArcCCW(xn, yn, r, theta1 - M_PI * 0.5, theta2 - M_PI * 0.5);
              glBegin(GL_LINES);
              glVertex2d(xn + r * std::cos(theta1 + M_PI * 0.5), yn + r * std::sin(theta1 + M_PI * 0.5));
              glVertex2d(xn + r * std::cos(theta2 + M_PI * 0.5), yn + r * std::sin(theta2 + M_PI * 0.5));
              glEnd();
            } else if (delta < 0.0) {
              // CW (right) turn: exterior arc on left side, straight segment on concave right side
              drawArcCCW(xn, yn, r, theta2 + M_PI * 0.5, theta1 + M_PI * 0.5);
              glBegin(GL_LINES);
              glVertex2d(xn + r * std::cos(theta1 - M_PI * 0.5), yn + r * std::sin(theta1 - M_PI * 0.5));
              glVertex2d(xn + r * std::cos(theta2 - M_PI * 0.5), yn + r * std::sin(theta2 - M_PI * 0.5));
              glEnd();
            }
            // delta == 0: straight, tangent lines meet, nothing to close
          } else if (pv == null_size_t && nx != null_size_t) {
            // Start node: semicircle cap (back end)
            double d2x = Conf.cells[i].nodes[nx].pos.x - xn;
            double d2y = Conf.cells[i].nodes[nx].pos.y - yn;
            double theta2 = std::atan2(d2y, d2x);
            // Cap sweeps from left side (+pi/2) CCW by pi to opposite left (-pi/2 = theta2-pi/2+pi)
            drawArcCCW(xn, yn, r, theta2 + M_PI * 0.5, theta2 + M_PI * 1.5);
          } else if (nx == null_size_t && pv != null_size_t) {
            // End node: semicircle cap (front end)
            double d1x = xn - Conf.cells[i].nodes[pv].pos.x;
            double d1y = yn - Conf.cells[i].nodes[pv].pos.y;
            double theta1 = std::atan2(d1y, d1x);
            // Cap sweeps from right side (-pi/2) CCW by pi to left side (+pi/2)
            drawArcCCW(xn, yn, r, theta1 - M_PI * 0.5, theta1 + M_PI * 0.5);
          }
        }
      }

      // Tangent lines along bars
      glBegin(GL_LINES);
      for (size_t b = 0; b < Conf.cells[i].bars.size(); ++b) {
        double xi = Conf.cells[i].nodes[Conf.cells[i].bars[b].i].pos.x;
        double yi = Conf.cells[i].nodes[Conf.cells[i].bars[b].i].pos.y;
        double xj = Conf.cells[i].nodes[Conf.cells[i].bars[b].j].pos.x;
        double yj = Conf.cells[i].nodes[Conf.cells[i].bars[b].j].pos.y;

        double nxij = xj - xi;
        double nyij = yj - yi;
        double nij = sqrt(nxij * nxij + nyij * nyij);
        nxij /= nij;
        nyij /= nij;
        double txij = -nyij;
        double tyij = nxij;
        double rad = Conf.cells[i].radius;

        // Right side (exterior for closed cells)
        glVertex2d(xj - rad * txij, yj - rad * tyij);
        glVertex2d(xi - rad * txij, yi - rad * tyij);

        if (!isClosed) {
          // Left side for open cells
          glVertex2d(xj + rad * txij, yj + rad * tyij);
          glVertex2d(xi + rad * txij, yi + rad * tyij);
        }
      }
      glEnd();
    }

    if (show_nodes) {
      glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
      glPointSize(4.0f);
      glBegin(GL_POINTS);
      for (size_t n = 0; n < Conf.cells[i].nodes.size(); ++n) {
        glVertex2d(Conf.cells[i].nodes[n].pos.x, Conf.cells[i].nodes[n].pos.y);
      }
      glEnd();
    }
  }
}

/**
 * Draws the glue points in the OpenGL context.
 */
void drawGluePoints() {

  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(1.0);
  double rp = Conf.cells[0].radius * 0.4;
  int rnb = 4;

  vec2r pos;
  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = Conf.cells[ci].neighbors.begin();
         InterIt != Conf.cells[ci].neighbors.end(); ++InterIt) {
      if (InterIt->glueState > 0) {
        size_t cj = InterIt->jc;
        size_t in = InterIt->in;
        size_t jn = InterIt->jn;
        int type = Conf.getPosition(ci, cj, in, jn, pos);
        if (type == 0) { // ne peut pas arriver normalement
          glColor3f(0.5f, 0.5f, 0.5f);
          rnb = 22;
        } else if (type == 1) {
          glColor3f(1.f, 0.f, 0.f);
          rnb = 12;
        } else if (type == 2) {
          glColor3f(1.f, 0.f, 0.f);
          rnb = 12;
        } else if (type == 3) {
          glColor3f(1.f, 0.6f, 0.f);
          rnb = 12;
        }

        drawCircle(pos.x, pos.y, rp, rnb);
      }
    }
  }

  // trace un trait entre les points frères
  glBegin(GL_LINES);
  glLineWidth(6.0f);
  glColor3f(0.f, 0.f, 1.0f);
  vec2r pos1, pos2;
  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = Conf.cells[ci].neighbors.begin();
         InterIt != Conf.cells[ci].neighbors.end(); ++InterIt) {
      if (InterIt->glueState > 0 && InterIt->brother != nullptr) {

        Conf.getPosition(InterIt->ic, InterIt->jc, InterIt->in, InterIt->jn, pos1);
        Conf.getPosition(InterIt->brother->ic, InterIt->brother->jc, InterIt->brother->in, InterIt->brother->jn, pos2);

        glVertex2d(pos1.x, pos1.y);
        glVertex2d(pos2.x, pos2.y);
      }
    }
  }
  glEnd();
}

// Charge les évènements de rupture depuis breakHistory.txt (si présent), une seule fois.
// Les lignes de commentaire (#) sont ignorées. Format attendu (cf. Lhyphen::recordBreakEvent) :
//   time a_ci a_cj a_in a_jn b_ci b_cj b_in b_jn released_NRJ
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

// Trace les interfaces rompues entre l'instant initial et le temps du conf affiché (Conf.t), comme le
// mode 'g' (drawGluePoints) le fait pour les interfaces encore collées : un trait épais entre les deux
// points de contact (côté A et son frère côté B). Les points sont reconstruits à partir des positions
// COURANTES, donc les traits suivent les cellules même si elles se sont déplacées depuis la rupture.
void drawCrackPath() {
  if (breakEvents.empty()) {
    return;
  }

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glLineWidth(6.0f);
  glColor3f(0.85f, 0.0f, 0.0f);

  glBegin(GL_LINES);
  vec2r pos1, pos2;
  for (const BreakEvent &ev : breakEvents) {
    if (ev.time > Conf.t) {
      continue; // rupture postérieure au conf affiché
    }
    if (!crackContactPos(ev.a_ci, ev.a_cj, ev.a_in, ev.a_jn, pos1)) {
      continue;
      //pos1.reset();
    }
    if (!crackContactPos(ev.b_ci, ev.b_cj, ev.b_in, ev.b_jn, pos2)) {
      continue;
      //pos2.reset();
    }
    glVertex2d(pos1.x, pos1.y);
    glVertex2d(pos2.x, pos2.y);
  }
  glEnd();
}

/**
 * Draws an arrow from point (x0, y0) to point (x1, y1) using OpenGL.
 *
 * @param x0  x-coordinate of the starting point
 * @param y0  y-coordinate of the starting point
 * @param x1  x-coordinate of the ending point
 * @param y1  y-coordinate of the ending point
 */
void arrow(double x0, double y0, double x1, double y1) {

  double nx = x1 - x0;
  double ny = y1 - y0;
  double len = sqrt(nx * nx + ny * ny);
  if (len == 0.0)
    return;
  nx /= len;
  ny /= len;

  glVertex2d(x0, y0);
  glVertex2d(x1, y1);

  double c = cos(arrowAngle);
  double s = sin(arrowAngle);

  double ex = c * nx - s * ny;
  double ey = s * nx + c * ny;
  glVertex2d(x1, y1);
  glVertex2d(x1 - arrowSize * ex, y1 - arrowSize * ey);

  ex = c * nx + s * ny;
  ey = -s * nx + c * ny;
  glVertex2d(x1, y1);
  glVertex2d(x1 - arrowSize * ex, y1 - arrowSize * ey);
}

/**
 * Draws the inter-cell force network ("force chains").
 *
 * For each contact between cells ci and cj, a poly-line center_i -> contact point -> center_j is
 * drawn, with thickness proportional to |fn| (normalized by the mean contact force, so the scaling
 * adapts to the absolute force units) and colored by sign: red = compression, blue = traction
 * (cohesion included). Magnitude is encoded ONLY by thickness; the length is geometric. This reveals
 * the load-bearing skeleton of the assembly, the standard granular force-chain picture.
 *
 * 's'/'S' tune the live thickness factor fnWidthFactor.
 */
void drawForces() {
  if (Conf.cells.empty()) return;

  // centres de cellules à jour (positions courantes du conf affiché)
  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    Conf.cells[ci].CellCenter();
  }

  // Statistiques sur les contacts actifs : moyenne (référence du filtre) et max (échelle d'épaisseur).
  double sum_f = 0.0, max_f = 0.0;
  size_t n_f = 0;
  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    for (const Neighbor &Inter : Conf.cells[ci].neighbors) {
      double f = std::fabs(Inter.fn + Inter.fn_coh);
      if (f > 0.0) {
        sum_f += f;
        if (f > max_f) max_f = f;
        n_f++;
      }
    }
  }
  if (n_f == 0 || max_f == 0.0) return;
  double mean_f = sum_f / (double)n_f;
  double threshold = forceFilter * mean_f; // ne garde que les chaînes porteuses (|fn| >= seuil)

  const float max_lw = 10.0f;
  const float min_lw = 0.4f;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  vec2r pc;
  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    for (const Neighbor &Inter : Conf.cells[ci].neighbors) {
      double fn_total = Inter.fn + Inter.fn_coh;
      if (std::fabs(fn_total) < threshold) continue; // filtre des forces faibles
      size_t cj = Inter.jc;
      if (cj >= Conf.cells.size()) continue;

      // point de contact sur l'interface (reconstruit depuis les positions courantes)
      Conf.getPosition(ci, cj, Inter.in, Inter.jn, pc);

      // rouge = compression (fn > 0), bleu = traction (fn < 0, cohésion comprise)
      if (fn_total > 0.0) glColor4f(0.85f, 0.15f, 0.15f, 0.9f);
      else                glColor4f(0.15f, 0.35f, 0.95f, 0.9f);

      float lw = (float)(max_lw * std::fabs(fn_total) / max_f * fnWidthFactor);
      lw = std::max(min_lw, std::min(max_lw, lw));
      glLineWidth(lw);

      // chaîne de force : centre_i -> contact -> centre_j
      glBegin(GL_LINE_STRIP);
      glVertex2d(Conf.cells[ci].center.x, Conf.cells[ci].center.y);
      glVertex2d(pc.x, pc.y);
      glVertex2d(Conf.cells[cj].center.x, Conf.cells[cj].center.y);
      glEnd();
    }
  }

  glLineWidth(1.0f);
}

void drawControlBoxes() {
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(2.0f);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  for (size_t i = 0; i < Conf.controlBoxAreas.size(); i++) {
    GLfloat dy = (GLfloat)Conf.controlBoxAreas[i].ymax - (GLfloat)Conf.controlBoxAreas[i].ymin;
    glText::print((GLfloat)Conf.controlBoxAreas[i].xmin, (GLfloat)Conf.controlBoxAreas[i].ymin + 1.1 * dy, 0.0f,
                  "%sx = %g, %sy = %g", (Conf.controlBoxAreas[i].xmode == VELOCITY_CONTROL) ? "V" : "F",
                  Conf.controlBoxAreas[i].xvalue, (Conf.controlBoxAreas[i].ymode == VELOCITY_CONTROL) ? "V" : "F",
                  Conf.controlBoxAreas[i].yvalue);
    glBegin(GL_LINE_LOOP);
    glVertex2d(Conf.controlBoxAreas[i].xmin, Conf.controlBoxAreas[i].ymin);
    glVertex2d(Conf.controlBoxAreas[i].xmax, Conf.controlBoxAreas[i].ymin);
    glVertex2d(Conf.controlBoxAreas[i].xmax, Conf.controlBoxAreas[i].ymax);
    glVertex2d(Conf.controlBoxAreas[i].xmin, Conf.controlBoxAreas[i].ymax);
    glEnd();
  }
}

bool try_to_readConf(int num, Lhyphen &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    OKNum = num;
    CF.loadCONF(file_name);
    return true;
  }
  std::cout << file_name << " does not exist" << std::endl;
  return false;
}

void captureScreenshot(const char *filename) {
  unsigned char *pixels  = new unsigned char[width * height * 3];
  unsigned char *flipped = new unsigned char[width * height * 3];
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  for (int y = 0; y < height; y++) {
    memcpy(flipped + (height - 1 - y) * width * 3, pixels + y * width * 3, width * 3);
  }
  stbi_write_png(filename, width, height, 3, flipped, width * 3);
  delete[] pixels;
  delete[] flipped;
}

// =====================================================================
// Main function
// =====================================================================

void readTomlOptions() {
  if (!fileTool::fileExists("see2-options.toml")) {
    saveTomlOptions();
    return;
  }

  toml::table tbl = toml::parse_file("see2-options.toml");

  if (tbl.contains("display")) {
    show_cells              = tbl["display"]["show_cells"].value_or(show_cells);
    show_glue_points        = tbl["display"]["show_glue_points"].value_or(show_glue_points);
    show_bar_colors         = tbl["display"]["show_bar_colors"].value_or(show_bar_colors);
    show_inter_cells_forces = tbl["display"]["show_inter_cells_forces"].value_or(show_inter_cells_forces);
    show_pressure           = tbl["display"]["show_pressure"].value_or(show_pressure);
    show_contours           = tbl["display"]["show_contours"].value_or(show_contours);
    show_nodes              = tbl["display"]["show_nodes"].value_or(show_nodes);
    show_control_boxes      = tbl["display"]["show_control_boxes"].value_or(show_control_boxes);
    show_background         = tbl["display"]["show_background"].value_or(show_background);
    show_crack_path         = tbl["display"]["show_crack_path"].value_or(show_crack_path);

    if (tbl["display"].as_table()->contains("bottom_color")) {
      auto arr = tbl["display"]["bottom_color"].as_array();
      if (arr && arr->size() >= 3) {
        bottom_r = (*arr)[0].value_or(bottom_r);
        bottom_g = (*arr)[1].value_or(bottom_g);
        bottom_b = (*arr)[2].value_or(bottom_b);
      }
    }
    if (tbl["display"].as_table()->contains("top_color")) {
      auto arr = tbl["display"]["top_color"].as_array();
      if (arr && arr->size() >= 3) {
        top_r = (*arr)[0].value_or(top_r);
        top_g = (*arr)[1].value_or(top_g);
        top_b = (*arr)[2].value_or(top_b);
      }
    }
  }

  if (tbl.contains("window")) {
    width  = tbl["window"]["width"].value_or(width);
    height = tbl["window"]["height"].value_or(height);
  }

  if (tbl.contains("view")) {
    fit_at_loading = tbl["view"]["fit_at_loading"].value_or(fit_at_loading);
    worldBox.min.x = tbl["view"]["xmin"].value_or(worldBox.min.x);
    worldBox.max.x = tbl["view"]["xmax"].value_or(worldBox.max.x);
    worldBox.min.y = tbl["view"]["ymin"].value_or(worldBox.min.y);
    worldBox.max.y = tbl["view"]["ymax"].value_or(worldBox.max.y);
  }

  if (tbl.contains("arrows")) {
    vScale         = tbl["arrows"]["vScale"].value_or(vScale);
    fnWidthFactor  = tbl["arrows"]["fnWidthFactor"].value_or(fnWidthFactor);
    forceFilter    = tbl["arrows"]["forceFilter"].value_or(forceFilter);
  }
}

void saveTomlOptions() {
  // clang-format off
  auto tbl = toml::table{
    {"display", toml::table{
      {"show_cells",              show_cells},
      {"show_glue_points",        show_glue_points},
      {"show_bar_colors",         show_bar_colors},
      {"show_inter_cells_forces", show_inter_cells_forces},
      {"show_pressure",           show_pressure},
      {"show_contours",           show_contours},
      {"show_nodes",              show_nodes},
      {"show_control_boxes",      show_control_boxes},
      {"show_background",         show_background},
      {"show_crack_path",         show_crack_path},
      {"bottom_color", toml::array{bottom_r, bottom_g, bottom_b}},
      {"top_color",    toml::array{top_r,    top_g,    top_b}},
    }},
    {"window", toml::table{
      {"width",  width},
      {"height", height},
    }},
    {"view", toml::table{
      {"fit_at_loading", fit_at_loading},
      {"xmin", worldBox.min.x},
      {"xmax", worldBox.max.x},
      {"ymin", worldBox.min.y},
      {"ymax", worldBox.max.y},
    }},
    {"arrows", toml::table{
      {"vScale",        vScale},
      {"fnWidthFactor", fnWidthFactor},
      {"forceFilter",   forceFilter},
    }},
  };
  // clang-format on

  std::ofstream file("see2-options.toml");
  file << tbl << "\n";
}

int main(int argc, char *argv[]) {

  if (argc == 1) {
    confNum = 0;
    try_to_readConf(confNum, Conf, confNum);
  } else if (argc == 2) {
    if (fileTool::fileExists(argv[1])) {
      std::cout << "Read " << argv[1] << std::endl;
      Conf.loadCONF(argv[1]);
    } else {
      confNum = atoi(argv[1]);
      try_to_readConf(confNum, Conf, confNum);
    }
  }

  Conf.findDisplayArea(1.15);

  // breakHistory.txt est lu une seule fois ici (à l'ouverture du premier conf) ; les évènements
  // sont ensuite filtrés par le temps du conf affiché dans drawCrackPath().
  readBreakHistory();

  readTomlOptions();

  // init color tables
  BarRedTable.setSize(128);
  BarRedTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {255, 0, 0, 255}});
  // BarRedTable.savePpm("BarRedTable.ppm");
  BarBlueTable.setSize(128);
  BarBlueTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {0, 0, 255, 255}});
  // BarBlueTable.savePpm("BarBlueTable.ppm");

  NodeRedTable.setSize(128);
  NodeRedTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {255, 0, 0, 255}});
  NodeBlueTable.setSize(128);
  NodeBlueTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {0, 0, 255, 255}});

  // ==== Init GLFW and create window
  if (!glfwInit()) {
    fprintf(stderr, "Failed to initialize GLFW\n");
    return -1;
  }

  // Désactive le support Retina (force une fenêtre en résolution 1x)
  glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);

  GLFWwindow *window = glfwCreateWindow(width, height, "see2", NULL, NULL);
  if (!window) {
    fprintf(stderr, "Failed to create GLFW window\n");
    glfwTerminate();
    return -1;
  }
  g_window = window;

  // ==== Register callbacks
  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, keyboard);
  glfwSetMouseButtonCallback(window, mouse_button);
  glfwSetCursorPosCallback(window, cursor_pos);
  glfwSetFramebufferSizeCallback(window, reshape);

  mouse_mode = MouseMode::NOTHING;

  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Other initialisations
  glText::init();
  updateTextLine();

  // ==== mainloop
  if (fit_at_loading) fit_view(window);
  updateTextLine();

  while (!glfwWindowShouldClose(window)) {
    glfwWaitEvents();
    if (needsRedraw) {
      reshape(window, width, height);
      display(window);
      needsRedraw = false;
    }
  }

  glfwTerminate();

  return 0;
}
