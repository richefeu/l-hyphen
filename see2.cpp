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
  std::cout << "e           show/hide the nodal velocity arrows" << std::endl;
  std::cout << "g           show/hide the glue points" << std::endl;
  std::cout << "r           show/hide the crack path (broken links up to current time)" << std::endl;
  std::cout << "d           show/hide the background gradient" << std::endl;
  std::cout << "i           show/hide the state panel (HUD)" << std::endl;
  std::cout << "k           cell strain colors: off / eps_v / eps_q / eps_xx / eps_yy / eps_xy" << std::endl;
  std::cout << "j           show/hide principal strain directions (thick major, thin minor; red tension, blue compression)" << std::endl;
  std::cout << "l           cell stress colors: off / sig_m / sig_q / p / p+sig_m / sig_xx / sig_yy / sig_xy" << std::endl;
  std::cout << "m           show/hide principal stress directions (thick major, thin minor; red tension, blue compression)" << std::endl;
  std::cout << "u/U         principal strain/stress lines shorter/longer (eScale)" << std::endl;
  std::cout << "o           use the displayed conf as strain reference (Shift+o: conf0)" << std::endl;
  std::cout << "h           show/hide the on-screen help (and print it here)" << std::endl;
  std::cout << "n           show/hide cell contours" << std::endl;
  std::cout << "v           show/hide nodes (points)" << std::endl;
  std::cout << "p           show/hide pressure" << std::endl;
  std::cout << "q           quit" << std::endl;
  std::cout << "s/S         force-chain lines thinner/thicker" << std::endl;
  std::cout << "t/T         force filter lower/higher (show more/fewer chains)" << std::endl;
  std::cout << "y/Y         velocity arrows shorter/longer (vScale)" << std::endl;
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
    show_help = 1 - show_help;
    printHelp();
  } break;

  case GLFW_KEY_K: {
    show_strain = (show_strain + 1) % nbStrainModes;
    if (show_strain) {
      show_stress = 0; // un seul remplissage à la fois
    }
    textZone.addLine("show_strain = %s (ref. conf%d)", strainModeNames[show_strain], refConfNum);
  } break;

  case GLFW_KEY_J: {
    show_strain_dirs = 1 - show_strain_dirs;
    if (show_strain_dirs) {
      show_stress_dirs = 0; // un seul jeu de traits à la fois
    }
    textZone.addLine("show_strain_dirs = %d (ref. conf%d)", show_strain_dirs, refConfNum);
  } break;

  case GLFW_KEY_L: {
    show_stress = (show_stress + 1) % nbStressModes;
    if (show_stress) {
      show_strain = 0;
    }
    textZone.addLine("show_stress = %s", stressModeNames[show_stress]);
  } break;

  case GLFW_KEY_SEMICOLON: { // m
    show_stress_dirs = 1 - show_stress_dirs;
    if (show_stress_dirs) {
      show_strain_dirs = 0;
    }
    textZone.addLine("show_stress_dirs = %d", show_stress_dirs);
  } break;

  case GLFW_KEY_U: {
    if (mods == GLFW_MOD_SHIFT) {
      eScale *= 1.25;
    } else {
      eScale *= 0.8;
    }
    textZone.addLine("eScale = %g", eScale);
  } break;

  case GLFW_KEY_O: {
    refConfNum = (mods == GLFW_MOD_SHIFT) ? 0 : confNum;
    textZone.addLine("strain reference = conf%d", refConfNum);
  } break;

  case GLFW_KEY_I: {
    show_hud = 1 - show_hud;
    textZone.addLine("show_hud = %d", show_hud);
  } break;

  case GLFW_KEY_D: {
    show_background = 1 - show_background;
    textZone.addLine("show_background = %d", show_background);
  } break;

  case GLFW_KEY_E: {
    show_velocities = 1 - show_velocities;
    textZone.addLine("show_velocities = %d", show_velocities);
  } break;

  case GLFW_KEY_Y: {
    if (mods == GLFW_MOD_SHIFT) {
      vScale *= 1.25;
    } else {
      vScale *= 0.8;
    }
    textZone.addLine("vScale = %g", vScale);
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

  if (show_strain || show_strain_dirs) {
    computeStrains();
  }
  if ((show_stress && show_stress != 3) || show_stress_dirs) {
    computeStresses();
  }
  // champ affiché en couleur (déformation ou contrainte, un seul à la fois) et borne de son échelle
  CellScalarField colorField;
  double colorBound = 0.0;
  if (show_strain) {
    colorField = strainField(show_strain);
    colorBound = fieldColorBound(colorField, strainColorMax);
  } else if (show_stress) {
    colorField = stressField(show_stress);
    colorBound = fieldColorBound(colorField, stressColorMax);
  }

  if (show_pressure) {
    drawPressure();
  }
  if (show_strain || show_stress) {
    drawCellScalars(colorField, colorBound);
  }
  if (show_cells) {
    drawCells();
  }
  lastDirsName = nullptr;
  if (show_strain_dirs) {
    lastDirsMax = drawTensorDirections(cellStrains, strainDirsMax);
    lastDirsName = "strain";
  }
  if (show_stress_dirs) {
    lastDirsMax = drawTensorDirections(cellStresses, stressDirsMax);
    lastDirsName = "stress";
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
  if (show_velocities) {
    drawVelocities();
  }

  if (show_control_boxes) {
    drawControlBoxes();
  }

  if (!snapshotMode) {
    textZone.draw();
  }

  lastFieldName = show_strain ? strainModeNames[show_strain] : (show_stress ? stressModeNames[show_stress] : nullptr);
  lastFieldBound = colorBound;
  lastFieldDivergent = colorField.divergent;
  if (!show_colorbar) {
    // barre de couleur masquée (images assemblées avec une barre commune)
  } else if (show_strain) {
    char extra[64];
    snprintf(extra, 64, "  ref conf%d", refConfNum);
    drawTensorColorBar(strainModeNames[show_strain], colorField.divergent ? -colorBound : 0.0, colorBound,
                       strainColorMax <= 0.0, extra);
  } else if (show_stress) {
    drawTensorColorBar(stressModeNames[show_stress], colorField.divergent ? -colorBound : 0.0, colorBound,
                       stressColorMax <= 0.0, "");
  }

  if (show_hud) {
    drawHUD();
  }
  if (show_help) {
    drawHelpOverlay();
  }

  glFlush();
  if (!snapshotMode) { // en mode instantané, l'image est lue dans le tampon arrière avant tout échange
    glfwSwapBuffers(window);
  }
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

// Callback de redimensionnement du framebuffer.
// On redessine ici même plutôt que de se contenter de lever needsRedraw : sous macOS,
// le glissement de la bordure enferme l'application dans une boucle d'événements Cocoa
// modale d'où glfwWaitEvents() ne revient pas, donc la boucle principale ne tournerait
// qu'au relâchement de la souris.
void framebuffer_size(GLFWwindow *window, int w, int h) {
  if (w <= 0 || h <= 0) {
    return; // fenêtre réduite : reshape() diviserait par zéro
  }
  reshape(window, w, h);
  display(window);
  needsRedraw = false;
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
// Trace une flèche (hampe + 2 barbules) sous forme de segments.
// À appeler entre glBegin(GL_LINES) et glEnd().
// Les barbules mesurent arrowSize * (longueur de la flèche) : elles restent donc
// proportionnées quelle que soit la longueur, contrairement à une taille absolue.
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
  double barb = arrowSize * len;

  double ex = c * nx - s * ny;
  double ey = s * nx + c * ny;
  glVertex2d(x1, y1);
  glVertex2d(x1 - barb * ex, y1 - barb * ey);

  ex = c * nx + s * ny;
  ey = -s * nx + c * ny;
  glVertex2d(x1, y1);
  glVertex2d(x1 - barb * ex, y1 - barb * ey);
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

// Flèches de vitesse aux noeuds.
// Comme drawForces(), on normalise par la vitesse max : la flèche la plus longue mesure
// vScale * (rayon moyen des cellules). L'affichage est donc lisible quel que soit l'ordre
// de grandeur des vitesses, et indépendant du zoom. vScale se règle avec 'y'/'Y'.
void drawVelocities() {
  double max_v = 0.0;
  double sum_radius = 0.0;
  size_t n_cell = Conf.cells.size();
  if (n_cell == 0) return;

  for (size_t ci = 0; ci < n_cell; ++ci) {
    sum_radius += Conf.cells[ci].radius;
    for (size_t in = 0; in < Conf.cells[ci].nodes.size(); ++in) {
      double v = Conf.cells[ci].nodes[in].vel.length();
      if (v > max_v) max_v = v;
    }
  }
  if (max_v == 0.0) return;

  double refLen = sum_radius / (double)n_cell;
  double scale  = vScale * refLen / max_v;

  glColor4f(0.10f, 0.20f, 0.85f, 0.9f);
  glLineWidth(1.5f);

  glBegin(GL_LINES);
  for (size_t ci = 0; ci < n_cell; ++ci) {
    for (size_t in = 0; in < Conf.cells[ci].nodes.size(); ++in) {
      const Node &N = Conf.cells[ci].nodes[in];
      arrow(N.pos.x, N.pos.y, N.pos.x + scale * N.vel.x, N.pos.y + scale * N.vel.y);
    }
  }
  glEnd();

  glLineWidth(1.0f);
}

// =====================================================================
// Tenseurs de déformation et de contrainte des cellules fermées
// =====================================================================

// Valeurs et directions principales du tenseur symétrique [[a, b], [b, d]] (l1 >= l2)
static void symEigen(double a, double b, double d, double &l1, double &l2, vec2r &u1, vec2r &u2) {
  double m = 0.5 * (a + d);
  double r = sqrt(0.25 * (a - d) * (a - d) + b * b);
  l1 = m + r;
  l2 = m - r;
  double theta = 0.5 * atan2(2.0 * b, a - d);
  u1.set(cos(theta), sin(theta));
  u2.set(-sin(theta), cos(theta));
}

// Centre (moyenne des noeuds) et rayon moyen (distance moyenne des noeuds au centre) d'une cellule
static void cellCenterAndRadius(const Cell &C, vec2r &center, double &R) {
  center.reset();
  for (size_t n = 0; n < C.nodes.size(); n++) {
    center += C.nodes[n].pos;
  }
  center /= (double)C.nodes.size();
  R = 0.0;
  for (size_t n = 0; n < C.nodes.size(); n++) {
    R += (C.nodes[n].pos - center).length();
  }
  R /= (double)C.nodes.size();
}

// Charge la conf de référence refConfNum dans RefConf (si ce n'est pas déjà fait)
bool loadRefConf() {
  if (loadedRefConfNum == refConfNum) {
    return true;
  }
  char file_name[256];
  snprintf(file_name, 256, "conf%d", refConfNum);
  if (!fileTool::fileExists(file_name)) {
    std::cout << "Reference " << file_name << " does not exist" << std::endl;
    textZone.addLine("reference %s does not exist", file_name);
    loadedRefConfNum = -1;
    return false;
  }
  std::cout << "Read reference " << file_name << std::endl;
  RefConf.loadCONF(file_name);
  loadedRefConfNum = refConfNum;
  return true;
}

// Pour chaque cellule fermée, gradient de transformation moyen F (moindres carrés) entre la conf de
// référence et la conf affichée : x_k = F X_k, positions des noeuds relatives au centre de la cellule,
// F = (sum x_k X_k^T) (sum X_k X_k^T)^-1. Les déformations principales sont celles de Hencky,
// e_i = 0.5 ln(b_i), où les b_i sont les valeurs propres de B = F F^T ; elles ne dépendent pas de la
// rotation de la cellule. Les vecteurs propres de B donnent les directions dans la conf affichée.
void computeStrains() {
  cellStrains.assign(Conf.cells.size(), CellTensor());
  if (!loadRefConf()) {
    return;
  }

  for (size_t c = 0; c < Conf.cells.size(); c++) {
    CellTensor &S = cellStrains[c];
    if (c >= RefConf.cells.size() || !Conf.cells[c].close) {
      continue;
    }
    const std::vector<Node> &cur = Conf.cells[c].nodes;
    const std::vector<Node> &ref = RefConf.cells[c].nodes;
    size_t nn = cur.size();
    if (nn < 3 || ref.size() != nn) {
      continue;
    }

    vec2r xc, Xc;
    double R;
    cellCenterAndRadius(Conf.cells[c], xc, R);
    for (size_t n = 0; n < nn; n++) {
      Xc += ref[n].pos;
    }
    Xc /= (double)nn;

    // A = sum x X^T, G = sum X X^T
    double A11 = 0.0, A12 = 0.0, A21 = 0.0, A22 = 0.0;
    double G11 = 0.0, G12 = 0.0, G22 = 0.0;
    for (size_t n = 0; n < nn; n++) {
      vec2r x = cur[n].pos - xc;
      vec2r X = ref[n].pos - Xc;
      A11 += x.x * X.x;
      A12 += x.x * X.y;
      A21 += x.y * X.x;
      A22 += x.y * X.y;
      G11 += X.x * X.x;
      G12 += X.x * X.y;
      G22 += X.y * X.y;
    }
    double detG = G11 * G22 - G12 * G12;
    if (detG <= 0.0) {
      continue;
    }
    // G^-1
    double I11 = G22 / detG, I12 = -G12 / detG, I22 = G11 / detG;
    // F = A G^-1
    double F11 = A11 * I11 + A12 * I12;
    double F12 = A11 * I12 + A12 * I22;
    double F21 = A21 * I11 + A22 * I12;
    double F22 = A21 * I12 + A22 * I22;
    if (F11 * F22 - F12 * F21 <= 0.0) {
      continue; // cellule retournée
    }
    // B = F F^T (symétrique)
    double b1, b2;
    symEigen(F11 * F11 + F12 * F12, F11 * F21 + F12 * F22, F21 * F21 + F22 * F22, b1, b2, S.u1, S.u2);
    if (b2 <= 0.0) {
      continue;
    }

    S.ok = true;
    S.center = xc;
    S.R = R;
    S.v1 = 0.5 * log(b1);
    S.v2 = 0.5 * log(b2);
    // composantes du tenseur de Hencky sum v_i u_i (x) u_i
    S.xx = S.v1 * S.u1.x * S.u1.x + S.v2 * S.u2.x * S.u2.x;
    S.yy = S.v1 * S.u1.y * S.u1.y + S.v2 * S.u2.y * S.u2.y;
    S.xy = S.v1 * S.u1.x * S.u1.y + S.v2 * S.u2.x * S.u2.y;
  }
}

// Contrainte moyenne de chaque cellule fermée (Love-Weber), à partir des seules forces d'interaction
// avec les autres cellules (contact et cohésion, sans les efforts internes de la cellule ni sa
// pression) : sigma = (1/A) sym( sum (x_c - centre) (x) f_c ), A étant la surface actuelle de la cellule.
// Une interaction (ci, in) / (cj, jn) exerce f = (fn + fn_coh) n + (ft + ft_coh) T sur ci et -f sur cj,
// toutes deux au point de contact x_c (Lhyphen::getPosition). Les forces visqueuses, non sauvegardées
// dans les conf, ne sont pas prises en compte. Les cellules ayant des noeuds pilotés (mors) ne sont pas
// affichées : la réaction du contrôle leur manque, leur tenseur serait incomplet.
void computeStresses() {
  size_t nc = Conf.cells.size();
  cellStresses.assign(nc, CellTensor());

  std::vector<vec2r> center(nc);
  std::vector<double> radius(nc, 0.0);
  for (size_t c = 0; c < nc; c++) {
    if (!Conf.cells[c].nodes.empty()) {
      cellCenterAndRadius(Conf.cells[c], center[c], radius[c]);
    }
  }

  std::vector<double> M11(nc, 0.0), M12(nc, 0.0), M21(nc, 0.0), M22(nc, 0.0); // sum x (x) f
  vec2r pc;
  for (size_t ci = 0; ci < nc; ci++) {
    for (const Neighbor &Inter : Conf.cells[ci].neighbors) {
      size_t cj = Inter.jc;
      if (cj >= nc || cj == ci) {
        continue;
      }
      double fn = Inter.fn + Inter.fn_coh;
      double ft = Inter.ft + Inter.ft_coh;
      if (fn == 0.0 && ft == 0.0) {
        continue;
      }
      vec2r T(-Inter.n.y, Inter.n.x);
      vec2r f = fn * Inter.n + ft * T;
      Conf.getPosition(ci, cj, Inter.in, Inter.jn, pc);

      vec2r xi = pc - center[ci];
      M11[ci] += xi.x * f.x;
      M12[ci] += xi.x * f.y;
      M21[ci] += xi.y * f.x;
      M22[ci] += xi.y * f.y;

      vec2r xj = pc - center[cj];
      M11[cj] -= xj.x * f.x;
      M12[cj] -= xj.x * f.y;
      M21[cj] -= xj.y * f.x;
      M22[cj] -= xj.y * f.y;
    }
  }

  for (size_t c = 0; c < nc; c++) {
    const Cell &C = Conf.cells[c];
    if (!C.close || C.nodes.size() < 3) {
      continue;
    }
    bool controlled = false;
    for (size_t n = 0; n < C.nodes.size(); n++) {
      if (C.nodes[n].ictrl != null_size_t) {
        controlled = true;
        break;
      }
    }
    if (controlled) {
      continue;
    }
    // surface actuelle (formule du lacet)
    double A = 0.0;
    for (size_t n = 0; n < C.nodes.size(); n++) {
      const vec2r &p = C.nodes[n].pos;
      const vec2r &q = C.nodes[(n + 1) % C.nodes.size()].pos;
      A += p.x * q.y - p.y * q.x;
    }
    A = 0.5 * fabs(A);
    if (A <= 0.0) {
      continue;
    }
    CellTensor &S = cellStresses[c];
    S.xx = M11[c] / A;
    S.yy = M22[c] / A;
    S.xy = 0.5 * (M12[c] + M21[c]) / A;
    symEigen(S.xx, S.xy, S.yy, S.v1, S.v2, S.u1, S.u2);
    S.ok = true;
    S.center = center[c];
    S.R = radius[c];
  }
}

// Champ de couleur de la déformation (voir strainModeNames) : 1 = eps_v = v1 + v2, 2 = eps_q = v1 - v2,
// 3/4/5 = eps_xx, eps_yy, eps_xy (composantes du tenseur de Hencky)
CellScalarField strainField(int mode) {
  CellScalarField F;
  F.divergent = (mode != 2);
  F.value.assign(cellStrains.size(), 0.0);
  F.ok.assign(cellStrains.size(), 0);
  for (size_t c = 0; c < cellStrains.size(); c++) {
    const CellTensor &T = cellStrains[c];
    if (!T.ok) {
      continue;
    }
    double v = 0.0;
    switch (mode) {
    case 1: v = T.v1 + T.v2; break;
    case 2: v = T.v1 - T.v2; break;
    case 3: v = T.xx; break;
    case 4: v = T.yy; break;
    case 5: v = T.xy; break;
    default: break;
    }
    F.ok[c] = 1;
    F.value[c] = v;
  }
  return F;
}

// Champ de couleur de la contrainte (voir stressModeNames) : 1 = sig_m = (v1 + v2)/2, 2 = sig_q = v1 - v2,
// 3 = p (pression interne, toutes les cellules fermées), 4 = p + sig_m, 5/6/7 = sig_xx, sig_yy, sig_xy
CellScalarField stressField(int mode) {
  CellScalarField F;
  F.divergent = (mode != 2);
  F.value.assign(Conf.cells.size(), 0.0);
  F.ok.assign(Conf.cells.size(), 0);
  for (size_t c = 0; c < Conf.cells.size(); c++) {
    if (mode == 3) {
      if (Conf.cells[c].close) {
        F.ok[c] = 1;
        F.value[c] = Conf.cells[c].p_int;
      }
      continue;
    }
    if (c >= cellStresses.size() || !cellStresses[c].ok) {
      continue;
    }
    const CellTensor &T = cellStresses[c];
    double sm = 0.5 * (T.v1 + T.v2);
    double v = 0.0;
    switch (mode) {
    case 1: v = sm; break;
    case 2: v = T.v1 - T.v2; break;
    case 4: v = Conf.cells[c].p_int + sm; break;
    case 5: v = T.xx; break;
    case 6: v = T.yy; break;
    case 7: v = T.xy; break;
    default: break;
    }
    F.ok[c] = 1;
    F.value[c] = v;
  }
  return F;
}

// Borne de l'échelle de couleur : fixedMax si > 0, sinon max de |valeur| sur les cellules coloriées
double fieldColorBound(const CellScalarField &field, double fixedMax) {
  if (fixedMax > 0.0) {
    return fixedMax;
  }
  double vmax = 0.0;
  for (size_t c = 0; c < field.value.size(); c++) {
    if (field.ok[c]) {
      vmax = std::max(vmax, fabs(field.value[c]));
    }
  }
  return (vmax > 0.0) ? vmax : 1.0e-12;
}

// Remplissage des cellules : échelle bleu-blanc-rouge symétrique [-vmax, vmax] pour un champ signé,
// blanc-jaune-rouge [0, vmax] sinon
void drawCellScalars(const CellScalarField &field, double vmax) {
  ColorTable &table = field.divergent ? TensorSphTable : TensorDevTable;
  if (field.divergent) {
    table.setMinMax((float)(-vmax), (float)vmax);
  } else {
    table.setMinMax(0.0f, (float)vmax);
  }

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  color4f col;
  for (size_t i = 0; i < Conf.cells.size(); ++i) {
    if (i >= field.value.size() || !field.ok[i]) {
      continue;
    }
    table.getColor4f((float)field.value[i], &col);
    glColor3f(col.r, col.g, col.b);

    std::vector<vec2r> contour;
    for (size_t n = 0; n < Conf.cells[i].nodes.size(); ++n) {
      contour.push_back(Conf.cells[i].nodes[n].pos);
    }
    std::vector<int> result;
    TriangulatePolygon::Process(contour, result);
    glBegin(GL_TRIANGLES);
    for (size_t s = 0; s < result.size(); s += 3) {
      glVertex2d(contour[result[s]].x, contour[result[s]].y);
      glVertex2d(contour[result[s + 1]].x, contour[result[s + 1]].y);
      glVertex2d(contour[result[s + 2]].x, contour[result[s + 2]].y);
    }
    glEnd();
  }
}

// Directions principales : un trait centré sur la cellule pour chaque direction (épais = majeure v1,
// fin = mineure v2), de demi-longueur eScale * R * |v_i| / vmax. Par défaut vmax = max|v| de la conf : le plus
// grand trait mesure alors eScale rayons de sa cellule ; les longueurs sont comparables entre cellules de même
// taille. Avec fixedMax > 0 (strain_dirsMax, stress_dirsMax), vmax = fixedMax : les longueurs sont aussi
// comparables d'une conf ou d'un calcul à l'autre. Retourne vmax.
// Couleur selon le signe (convention tension positive) : rouge = tension (v_i > 0),
// bleu = compression (v_i < 0).
double drawTensorDirections(const std::vector<CellTensor> &tensors, double fixedMax) {
  double vmax = 0.0;
  for (size_t c = 0; c < tensors.size(); c++) {
    if (tensors[c].ok) {
      vmax = std::max(vmax, std::max(fabs(tensors[c].v1), fabs(tensors[c].v2)));
    }
  }
  if (fixedMax > 0.0) {
    vmax = fixedMax;
  }
  if (vmax == 0.0) {
    return 0.0;
  }

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // un passage par direction (l'épaisseur ne peut pas changer entre glBegin et glEnd)
  for (int dir = 1; dir <= 2; dir++) {
    glLineWidth((dir == 1) ? 4.0f : 1.5f);
    glBegin(GL_LINES);
    for (size_t c = 0; c < tensors.size(); c++) {
      const CellTensor &S = tensors[c];
      if (!S.ok) {
        continue;
      }
      double v = (dir == 1) ? S.v1 : S.v2;
      const vec2r &u = (dir == 1) ? S.u1 : S.u2;
      double l = eScale * S.R * fabs(v) / vmax;
      if (v > 0.0) {
        glColor4f(0.85f, 0.05f, 0.05f, 1.0f); // tension
      } else {
        glColor4f(0.05f, 0.15f, 0.90f, 1.0f); // compression
      }
      glVertex2d(S.center.x - l * u.x, S.center.y - l * u.y);
      glVertex2d(S.center.x + l * u.x, S.center.y + l * u.y);
    }
    glEnd();
  }
  glLineWidth(1.0f);
  return vmax;
}

// Barre de couleur (coin bas-droit) ; vmin < 0 pour la partie sphérique (échelle divergente)
void drawTensorColorBar(const char *name, double vmin, double vmax, bool autoBound, const char *extra) {
  ColorTable &table = (vmin < 0.0) ? TensorSphTable : TensorDevTable;
  table.setMinMax((float)vmin, (float)vmax);

  const int glyphW = 10; // largeur approximative d'un caractère de glText
  const int barW = 260;
  const int barH = 14;
  const int x0 = width - barW - 30;
  const int y0 = 30;
  const int nSeg = 64;

  char title[128], smin[32], smax[32];
  snprintf(title, 128, "%s%s%s", name, extra, autoBound ? "  (auto)" : "");
  snprintf(smin, 32, "%.3g", vmin);
  snprintf(smax, 32, "%.3g", vmax);

  switch2D::go(width, height);

  glColor4f(0.0f, 0.0f, 0.0f, 0.3f);
  glBegin(GL_QUADS);
  glVertex2i(x0 - 10, y0 - 22);
  glVertex2i(x0 + barW + 10, y0 - 22);
  glVertex2i(x0 + barW + 10, y0 + barH + 26);
  glVertex2i(x0 - 10, y0 + barH + 26);
  glEnd();

  color4f col;
  glBegin(GL_QUADS);
  for (int k = 0; k < nSeg; k++) {
    double v = vmin + (vmax - vmin) * (k + 0.5) / (double)nSeg;
    table.getColor4f((float)v, &col);
    glColor3f(col.r, col.g, col.b);
    int xa = x0 + (barW * k) / nSeg;
    int xb = x0 + (barW * (k + 1)) / nSeg;
    glVertex2i(xa, y0);
    glVertex2i(xb, y0);
    glVertex2i(xb, y0 + barH);
    glVertex2i(xa, y0 + barH);
  }
  glEnd();

  glColor3f(0.95f, 0.95f, 0.95f);
  glText::print(x0, y0 + barH + 8, "%s", title);
  glText::print(x0, y0 - 16, "%s", smin);
  glText::print(x0 + barW - glyphW * (int)strlen(smax), y0 - 16, "%s", smax);

  switch2D::back();
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

// Panneau d'état permanent (coin haut-gauche) : liste des toggles et leur état.
void drawHUD() {
  struct Item {
    const char *key;
    const char *name;
    int state;
  };
  // NB: labels de touches en AZERTY, comme dans printHelp()
  Item items[] = {
      {"c", "cells", show_cells},
      {"n", "contours", show_contours},
      {"v", "nodes", show_nodes},
      {"p", "pressure", show_pressure},
      {"f", "forces", show_inter_cells_forces},
      {"e", "velocities", show_velocities},
      {"g", "glue", show_glue_points},
      {"r", "crack path", show_crack_path},
      {"k", "strain", show_strain},
      {"j", "strain dirs", show_strain_dirs},
      {"l", "stress", show_stress},
      {"m", "stress dirs", show_stress_dirs},
      {"b", "bar colors", show_bar_colors},
      {"a", "ctrl boxes", show_control_boxes},
      {"d", "background", show_background},
  };
  const int n = (int)(sizeof(items) / sizeof(items[0]));

  const int nInfo  = 6; // conf/temps + les réglages continus + la référence des déformations
  const int glyphH = 13;
  const int lineH  = 16;
  const int padX   = 8;
  const int padY   = 8;
  const int panelW = 250;
  const int panelH = padY * 2 + (n + 1 + 1 + nInfo) * lineH; // titre + toggles + séparateur + infos

  const int x0 = 6;
  const int y1 = height - 6;   // bord haut du panneau
  const int y0 = y1 - panelH;  // bord bas

  switch2D::go(width, height);

  // Fond translucide
  glColor4f(0.0f, 0.0f, 0.0f, 0.3f);
  glBegin(GL_QUADS);
  glVertex2i(x0, y0);
  glVertex2i(x0 + panelW, y0);
  glVertex2i(x0 + panelW, y1);
  glVertex2i(x0, y1);
  glEnd();

  int ty = y1 - padY - glyphH;
  glColor3f(0.85f, 0.85f, 0.90f);
  glText::print(x0 + padX, ty, "display    i:hud  h:help");
  ty -= lineH;

  for (int i = 0; i < n; ++i) {
    if (items[i].state) {
      glColor3f(0.35f, 0.95f, 0.45f);
    } else {
      //glColor3f(0.55f, 0.55f, 0.55f);
      glColor3f(0.95f, 0.35f, 0.35f);
    }
    glText::print(x0 + padX, ty, "%s  %-11s %s", items[i].key, items[i].name, items[i].state ? "ON" : "off");
    ty -= lineH;
  }

  // Séparateur
  ty -= lineH / 2;
  glColor4f(0.7f, 0.7f, 0.75f, 0.5f);
  glBegin(GL_LINES);
  glVertex2i(x0 + padX, ty + 4);
  glVertex2i(x0 + panelW - padX, ty + 4);
  glEnd();
  ty -= lineH / 2;

  // Valeurs courantes des réglages continus
  glColor3f(0.85f, 0.85f, 0.90f);
  glText::print(x0 + padX, ty, "conf %-4d t = %.4g", confNum, Conf.t);
  ty -= lineH;
  glText::print(x0 + padX, ty, "s/S fn width  %.3g", fnWidthFactor);
  ty -= lineH;
  glText::print(x0 + padX, ty, "t/T filter    %.3g", forceFilter);
  ty -= lineH;
  glText::print(x0 + padX, ty, "y/Y vel scale %.3g", vScale);
  ty -= lineH;
  glText::print(x0 + padX, ty, "u/U tens scale %.3g", eScale);
  ty -= lineH;
  glText::print(x0 + padX, ty, "o   eps ref   conf%d", refConfNum);

  switch2D::back();
}

// Overlay d'aide (touche 'h') : rappel de tous les raccourcis clavier.
void drawHelpOverlay() {
  const char *lines[] = {
      "Keyboard shortcuts",
      "",
      "a           show/hide control area boxes",
      "b           colorize the cell bars",
      "c           show/hide the cells",
      "f           show/hide the forces",
      "e           show/hide the nodal velocity arrows",
      "g           show/hide the glue points",
      "r           show/hide the crack path",
      "n           show/hide cell contours",
      "v           show/hide nodes (points)",
      "p           show/hide pressure",
      "d           show/hide the background gradient",
      "i           show/hide this state panel (HUD)",
      "k           strain colors (eps_v, eps_q, eps_xx, ...)",
      "j           principal strain directions",
      "l           stress colors (sig_m, sig_q, p, ...)",
      "m           principal stress directions",
      "u/U         tensor lines shorter/longer (eScale)",
      "o           strain reference = this conf (Shift: 0)",
      "h           show/hide this help",
      "q           quit",
      "s/S         force-chain lines thinner/thicker",
      "t/T         force filter lower/higher",
      "y/Y         velocity arrows shorter/longer (vScale)",
      "z/Z         zoom in/out",
      "->  / <-    load next / previous configuration",
      "Shift+<-    jump to conf 0",
      "=           fit the view",
      "x           save screenshot.png",
      "Shift+x     batch screenshots of all confs",
      "space       save options to see2-options.toml",
      "Shift+space reload options from see2-options.toml",
  };
  const int n = (int)(sizeof(lines) / sizeof(lines[0]));

  const int glyphH  = 13;
  const int glyphW  = 10;
  const int lineH   = 16;
  const int padX    = 14;
  const int padY    = 12;
  const int panelW  = 50 * glyphW + 2 * padX;
  const int panelH  = 2 * padY + n * lineH;

  const int x0 = (width - panelW) / 2;
  const int y0 = (height - panelH) / 2;
  const int y1 = y0 + panelH;

  switch2D::go(width, height);

  // Fond translucide + liseré
  glColor4f(0.05f, 0.05f, 0.08f, 0.82f);
  glBegin(GL_QUADS);
  glVertex2i(x0, y0);
  glVertex2i(x0 + panelW, y0);
  glVertex2i(x0 + panelW, y1);
  glVertex2i(x0, y1);
  glEnd();
  glColor4f(0.6f, 0.6f, 0.7f, 0.9f);
  glBegin(GL_LINE_LOOP);
  glVertex2i(x0, y0);
  glVertex2i(x0 + panelW, y0);
  glVertex2i(x0 + panelW, y1);
  glVertex2i(x0, y1);
  glEnd();

  int ty = y1 - padY - glyphH;
  for (int i = 0; i < n; ++i) {
    if (i == 0) {
      glColor3f(0.95f, 0.85f, 0.35f);
    } else {
      glColor3f(0.92f, 0.92f, 0.92f);
    }
    glText::print(x0 + padX, ty, "%s", lines[i]);
    ty -= lineH;
  }

  switch2D::back();
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
  glPixelStorei(GL_PACK_ALIGNMENT, 1); // lignes non alignées sur 4 octets si la largeur n'est pas multiple de 4
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  for (int y = 0; y < height; y++) {
    memcpy(flipped + (height - 1 - y) * width * 3, pixels + y * width * 3, width * 3);
  }
  stbi_write_png(filename, width, height, 3, flipped, width * 3);
  delete[] pixels;
  delete[] flipped;
}

// Rendu hors écran (--snapshot) : la conf est dessinée dans un framebuffer de la taille exacte demandée
// (multiéchantillonné puis résolu), indépendamment de la fenêtre, qui reste cachée, et de la résolution de
// l'écran (sous macOS, une fenêtre cachée annonce un framebuffer Retina qu'elle n'a pas).
bool renderOffscreen(GLFWwindow *window, int w, int h, const char *filename) {
  GLint maxSamples = 0;
  glGetIntegerv(GL_MAX_SAMPLES_EXT, &maxSamples);
  GLint samples = std::min(4, (int)maxSamples);

  GLuint fboMS = 0, rbMS = 0, fbo = 0, rb = 0;
  glGenFramebuffersEXT(1, &fboMS);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fboMS);
  glGenRenderbuffersEXT(1, &rbMS);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbMS);
  glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, samples, GL_RGBA8, w, h);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, rbMS);
  if (glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT) != GL_FRAMEBUFFER_COMPLETE_EXT) {
    std::cerr << "see2 : framebuffer hors écran incomplet" << std::endl;
    return false;
  }

  reshape(nullptr, w, h); // fixe width, height, la vue et la projection
  display(window);        // en mode instantané : pas d'échange de tampons

  // résolution du multiéchantillonnage dans un framebuffer simple, puis lecture
  glGenFramebuffersEXT(1, &fbo);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
  glGenRenderbuffersEXT(1, &rb);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rb);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA8, w, h);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, rb);
  glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, fboMS);
  glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, fbo);
  glBlitFramebufferEXT(0, 0, w, h, 0, 0, w, h, GL_COLOR_BUFFER_BIT, GL_NEAREST);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  captureScreenshot(filename);

  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
  glDeleteRenderbuffersEXT(1, &rb);
  glDeleteFramebuffersEXT(1, &fbo);
  glDeleteRenderbuffersEXT(1, &rbMS);
  glDeleteFramebuffersEXT(1, &fboMS);
  return true;
}

// =====================================================================
// Main function
// =====================================================================

void readTomlOptions() {
  if (!fileTool::fileExists(optionsFile.c_str())) {
    if (!snapshotMode) { // en mode instantané, on n'écrit rien dans le répertoire du calcul
      saveTomlOptions();
    }
    return;
  }

  toml::table tbl = toml::parse_file(optionsFile);

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
    show_velocities         = tbl["display"]["show_velocities"].value_or(show_velocities);
    show_hud                = tbl["display"]["show_hud"].value_or(show_hud);
    show_strain             = tbl["display"]["show_strain"].value_or(show_strain);
    show_strain_dirs        = tbl["display"]["show_strain_dirs"].value_or(show_strain_dirs);
    show_stress             = tbl["display"]["show_stress"].value_or(show_stress);
    show_stress_dirs        = tbl["display"]["show_stress_dirs"].value_or(show_stress_dirs);
    if (show_strain < 0 || show_strain >= nbStrainModes) show_strain = 0;
    if (show_stress < 0 || show_stress >= nbStressModes) show_stress = 0;

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
    arrowSize      = tbl["arrows"]["arrowSize"].value_or(arrowSize);
    arrowAngle     = tbl["arrows"]["arrowAngle"].value_or(arrowAngle);
    fnWidthFactor  = tbl["arrows"]["fnWidthFactor"].value_or(fnWidthFactor);
    forceFilter    = tbl["arrows"]["forceFilter"].value_or(forceFilter);
  }

  if (tbl.contains("strain")) {
    refConfNum     = tbl["strain"]["refConf"].value_or(refConfNum);
    eScale         = tbl["strain"]["eScale"].value_or(eScale);
    strainColorMax = tbl["strain"]["colorMax"].value_or(strainColorMax);
  }

  if (tbl.contains("stress")) {
    stressColorMax = tbl["stress"]["colorMax"].value_or(stressColorMax);
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
      {"show_velocities",         show_velocities},
      {"show_hud",                show_hud},
      {"show_strain",             show_strain},
      {"show_strain_dirs",        show_strain_dirs},
      {"show_stress",             show_stress},
      {"show_stress_dirs",        show_stress_dirs},
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
      {"arrowSize",     arrowSize},
      {"arrowAngle",    arrowAngle},
      {"fnWidthFactor", fnWidthFactor},
      {"forceFilter",   forceFilter},
    }},
    {"strain", toml::table{
      {"refConf",  refConfNum},
      {"eScale",   eScale},
      {"colorMax", strainColorMax},
    }},
    {"stress", toml::table{
      {"colorMax", stressColorMax},
    }},
  };
  // clang-format on

  std::ofstream file(optionsFile);
  file << tbl << "\n";
}

// ---------------------------------------------------------------------
// Ligne de commande
// ---------------------------------------------------------------------

void printUsage() {
  std::cout << "Usage : see2 [conf] [options]\n"
               "  conf                 numéro N (fichier confN), nom de fichier, ou « last » (dernier conf) ; défaut : 0\n"
               "  --options FICHIER    fichier d'options (défaut : see2-options.toml)\n"
               "  --set NOM=VALEUR     surcharge une option après lecture du fichier (répétable), ex. :\n"
               "                       show_crack_path=1  show_hud=0  show_strain=eps_yy  strain_colorMax=0.01\n"
               "  --size LxH           taille de la fenêtre / de l'image en pixels\n"
               "  --snapshot IMAGE.png rend la conf dans IMAGE.png (fenêtre cachée) puis quitte\n"
               "  --list-options       liste les noms utilisables avec --set\n"
               "  -h, --help           cette aide\n";
}

// Options modifiables par --set (mêmes noms que dans see2-options.toml)
struct SetOption {
  const char *name;
  int *i;
  double *d;
  const char *help;
};

std::vector<SetOption> setOptions() {
  return {
      {"show_cells", &show_cells, nullptr, "cellules (0/1)"},
      {"show_contours", &show_contours, nullptr, "contours des cellules (0/1)"},
      {"show_nodes", &show_nodes, nullptr, "noeuds (0/1)"},
      {"show_pressure", &show_pressure, nullptr, "pression interne (0/1)"},
      {"show_inter_cells_forces", &show_inter_cells_forces, nullptr, "forces entre cellules (0/1)"},
      {"show_velocities", &show_velocities, nullptr, "vitesses (0/1)"},
      {"show_glue_points", &show_glue_points, nullptr, "points de colle (0/1)"},
      {"show_crack_path", &show_crack_path, nullptr, "chemin de fissure, liens rompus (0/1)"},
      {"show_strain", &show_strain, nullptr, "déformation : off eps_v eps_q eps_xx eps_yy eps_xy (nom ou 0-5)"},
      {"show_strain_dirs", &show_strain_dirs, nullptr, "directions principales de déformation (0/1)"},
      {"show_stress", &show_stress, nullptr, "contrainte : off sig_m sig_q p p+sig_m sig_xx sig_yy sig_xy (nom ou 0-7)"},
      {"show_stress_dirs", &show_stress_dirs, nullptr, "directions principales de contrainte (0/1)"},
      {"show_bar_colors", &show_bar_colors, nullptr, "couleur des barres selon l'effort (0/1)"},
      {"show_control_boxes", &show_control_boxes, nullptr, "boîtes de contrôle (0/1)"},
      {"show_background", &show_background, nullptr, "fond en dégradé (0/1)"},
      {"show_hud", &show_hud, nullptr, "panneau d'état des options (0/1)"},
      {"show_colorbar", &show_colorbar, nullptr, "barre de couleur du champ affiché (0/1)"},
      {"fit_at_loading", &fit_at_loading, nullptr, "cadrage automatique sur l'échantillon (0/1)"},
      {"refConf", &refConfNum, nullptr, "conf de référence des déformations"},
      {"eScale", nullptr, &eScale, "longueur des traits de directions principales"},
      {"strain_colorMax", nullptr, &strainColorMax, "borne de l'échelle de déformation (0 = automatique)"},
      {"strain_dirsMax", nullptr, &strainDirsMax, "déformation principale du plus grand trait (0 = max de la conf)"},
      {"stress_dirsMax", nullptr, &stressDirsMax, "contrainte principale du plus grand trait (0 = max de la conf)"},
      {"stress_colorMax", nullptr, &stressColorMax, "borne de l'échelle de contrainte (0 = automatique)"},
      {"xmin", nullptr, &worldBox.min.x, "fenêtre affichée (avec fit_at_loading=0)"},
      {"xmax", nullptr, &worldBox.max.x, ""},
      {"ymin", nullptr, &worldBox.min.y, ""},
      {"ymax", nullptr, &worldBox.max.y, ""},
  };
}

bool applySetOption(const std::string &arg) {
  auto eq = arg.find('=');
  if (eq == std::string::npos) {
    std::cerr << "see2 : --set attend NOM=VALEUR (reçu « " << arg << " »)" << std::endl;
    return false;
  }
  std::string name = arg.substr(0, eq), value = arg.substr(eq + 1);
  // les champs de déformation et de contrainte peuvent être donnés par leur nom
  if (name == "show_strain" || name == "show_stress") {
    const char **names = (name == "show_strain") ? strainModeNames : stressModeNames;
    int nb = (name == "show_strain") ? nbStrainModes : nbStressModes;
    for (int k = 0; k < nb; k++) {
      if (value == names[k]) {
        value = std::to_string(k);
      }
    }
  }
  for (auto &o : setOptions()) {
    if (name != o.name) {
      continue;
    }
    try {
      size_t pos = 0;
      if (o.i) {
        *o.i = std::stoi(value, &pos);
      } else {
        *o.d = std::stod(value, &pos);
      }
      if (pos != value.size()) {
        throw std::invalid_argument(value);
      }
    } catch (...) {
      std::cerr << "see2 : valeur invalide pour " << name << " : « " << value << " »" << std::endl;
      return false;
    }
    if (show_strain < 0 || show_strain >= nbStrainModes) show_strain = 0;
    if (show_stress < 0 || show_stress >= nbStressModes) show_stress = 0;
    return true;
  }
  std::cerr << "see2 : option inconnue pour --set : « " << name << " » (voir --list-options)" << std::endl;
  return false;
}

int main(int argc, char *argv[]) {

  // --- ligne de commande
  std::string confArg, snapshotFile;
  std::vector<std::string> sets;
  int cliWidth = 0, cliHeight = 0;
  for (int a = 1; a < argc; a++) {
    std::string arg = argv[a];
    auto need = [&]() {
      if (a + 1 >= argc) {
        std::cerr << "see2 : valeur manquante après " << arg << std::endl;
        exit(1);
      }
      return std::string(argv[++a]);
    };
    if (arg == "-h" || arg == "--help") {
      printUsage();
      return 0;
    } else if (arg == "--list-options") {
      for (auto &o : setOptions()) {
        printf("  %-24s %s\n", o.name, o.help);
      }
      return 0;
    } else if (arg == "--options") {
      optionsFile = need();
    } else if (arg == "--set") {
      sets.push_back(need());
    } else if (arg == "--snapshot") {
      snapshotFile = need();
      snapshotMode = true;
    } else if (arg == "--size") {
      std::string v = need();
      if (sscanf(v.c_str(), "%dx%d", &cliWidth, &cliHeight) != 2 || cliWidth <= 0 || cliHeight <= 0) {
        std::cerr << "see2 : --size attend LxH, par exemple 800x600" << std::endl;
        return 1;
      }
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "see2 : option inconnue « " << arg << " »" << std::endl;
      printUsage();
      return 1;
    } else {
      confArg = arg;
    }
  }

  // --- conf à afficher
  if (confArg.empty()) {
    confNum = 0;
    try_to_readConf(confNum, Conf, confNum);
  } else if (confArg == "last") {
    int n = 0;
    char name[256];
    do {
      snprintf(name, 256, "conf%d", n + 1);
    } while (fileTool::fileExists(name) && ++n);
    if (!try_to_readConf(n, Conf, confNum)) {
      return 1;
    }
  } else if (fileTool::fileExists(confArg.c_str())) {
    std::cout << "Read " << confArg << std::endl;
    Conf.loadCONF(confArg.c_str());
  } else {
    char *end = nullptr;
    long n = strtol(confArg.c_str(), &end, 10);
    if (*end != '\0' || !try_to_readConf((int)n, Conf, confNum)) {
      std::cerr << "see2 : conf introuvable : " << confArg << std::endl;
      return 1;
    }
  }

  Conf.findDisplayArea(1.15);

  // breakHistory.txt est lu une seule fois ici (à l'ouverture du premier conf) ; les évènements
  // sont ensuite filtrés par le temps du conf affiché dans drawCrackPath().
  readBreakHistory();

  readTomlOptions();
  for (auto &a : sets) {
    if (!applySetOption(a)) {
      return 1;
    }
  }
  if (cliWidth > 0) {
    width = cliWidth;
    height = cliHeight;
  }

  // init color tables
  BarRedTable.setSize(128);
  BarRedTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {255, 0, 0, 255}});
  // BarRedTable.savePpm("BarRedTable.ppm");
  BarBlueTable.setSize(128);
  BarBlueTable.rebuild_interp_rgba({0, 127}, {{0, 255, 0, 255}, {0, 0, 255, 255}});
  // BarBlueTable.savePpm("BarBlueTable.ppm");

  TensorSphTable.setSize(128);
  TensorSphTable.rebuild_interp_rgba({0, 63, 127}, {{40, 60, 200, 255}, {255, 255, 255, 255}, {200, 30, 30, 255}});
  TensorDevTable.setSize(128);
  TensorDevTable.rebuild_interp_rgba({0, 63, 127}, {{255, 255, 255, 255}, {255, 210, 60, 255}, {200, 30, 30, 255}});

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

  // Anti-aliasing (MSAA 4x) pour des contours et chaînes de force plus nets
  glfwWindowHint(GLFW_SAMPLES, 4);

  const int imageWidth = width, imageHeight = height; // taille demandée (reshape la remplace par celle du framebuffer)
  if (snapshotMode) {
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE); // fenêtre cachée : rendu hors écran puis lecture de l'image
  }

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
  glfwSetFramebufferSizeCallback(window, framebuffer_size);
  glfwSetWindowRefreshCallback(window, display); // ré-exposition (dé-minimisation, passage au 1er plan)

  mouse_mode = MouseMode::NOTHING;

  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Anti-aliasing : MSAA + lissage des lignes/points
  glEnable(GL_MULTISAMPLE);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  // ==== Other initialisations
  glText::init();
  updateTextLine();

  // ==== mainloop
  if (fit_at_loading) fit_view(window);
  updateTextLine();

  if (snapshotMode) {
    if (!renderOffscreen(window, imageWidth, imageHeight, snapshotFile.c_str())) {
      glfwTerminate();
      return 1;
    }
    std::cout << "snapshot " << snapshotFile << " : conf" << confNum << ", t = " << Conf.t << ", " << width << "x"
              << height << std::endl;
    if (lastFieldName) {
      std::cout << "field " << lastFieldName << " bound " << lastFieldBound << " divergent " << lastFieldDivergent
                << std::endl;
    }
    if (lastDirsName) {
      std::cout << "dirs " << lastDirsName << " max " << lastDirsMax << std::endl;
    }
    glfwTerminate();
    return 0;
  }

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
