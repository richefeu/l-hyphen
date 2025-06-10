#include "see.hpp"

void printHelp() {
  std::cout << std::endl;
  std::cout << "Commandes:" << std::endl;
  std::cout << "b         colorize the cell bars" << std::endl;
  std::cout << "c         show the cells" << std::endl;
  std::cout << "f         show the forces" << std::endl;
  std::cout << "g         show the glue points" << std::endl;
  std::cout << "h         show this help" << std::endl;
  std::cout << "p         show pressure" << std::endl;
  std::cout << "q         quit" << std::endl;
  std::cout << "s/S       scale the force length" << std::endl;
  std::cout << "z/Z       zoom in or out" << std::endl;
  std::cout << "+         load next configuration file" << std::endl;
  std::cout << "-         load previous configuration file" << std::endl;
  std::cout << "=         fit the view" << std::endl;
  std::cout << std::endl;
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case 'b': {
    show_bar_colors = 1 - show_bar_colors;
  } break;

  case 'c': {
    show_cells = 1 - show_cells;
  } break;

  case 'f': {
    show_inter_cells_forces = 1 - show_inter_cells_forces;
  } break;

  case 'g': {
    show_glue_points = 1 - show_glue_points;
  } break;

  case 'h': {
    printHelp();
  } break;

  case 'p': {
    show_pressure = 1 - show_pressure;
  } break;

  case 'q': {
    exit(0);
  } break;

  case 's': {
    vScale *= 0.9;
  } break;
  case 'S': {
    vScale *= 1.1;
  } break;

  case 'z': {
    double dy = worldBox.ymax - worldBox.ymin;
    double ddy = 0.2 * dy;
    double ddx = 0.2 * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax += ddx;
    worldBox.ymin -= ddy;
    worldBox.ymax += ddy;
    reshape(width, height);
  } break;
  case 'Z': {
    double dy = worldBox.ymax - worldBox.ymin;
    double ddy = -0.2 * dy;
    double ddx = -0.2 * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax += ddx;
    worldBox.ymin -= ddy;
    worldBox.ymax += ddy;
    reshape(width, height);
  } break;

  case '-': {
    if (confNum > 0)
      try_to_readConf(confNum - 1, Conf, confNum);
  } break;

  case '+': {
    try_to_readConf(confNum + 1, Conf, confNum);
  } break;

  case '=': {
    fit_view();
    reshape(width, height);
  } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = MouseMode::NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON:
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        mouse_mode = MouseMode::PAN;
      else
        mouse_mode = MouseMode::ROTATION;
      break;
    case GLUT_MIDDLE_BUTTON:
      mouse_mode = MouseMode::ZOOM;
      break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == MouseMode::NOTHING)
    return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case MouseMode::ZOOM: {
    double ddy = (worldBox.ymax - worldBox.ymin) * dy;
    double ddx = (worldBox.xmax - worldBox.xmin) * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax += ddx;
    worldBox.ymin -= ddy;
    worldBox.ymax += ddy;
  } break;

  case MouseMode::PAN: {
    double ddx = (worldBox.xmax - worldBox.xmin) * dx;
    double ddy = (worldBox.ymax - worldBox.ymin) * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax -= ddx;
    worldBox.ymin += ddy;
    worldBox.ymax += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  glutPostRedisplay();
  //display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_pressure)
    drawPressure();
  if (show_cells)
    drawCells();
  if (show_glue_points)
    drawGluePoints();
  if (show_inter_cells_forces)
    drawForces();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  worldBox.xmin = Conf.xmin;
  worldBox.xmax = Conf.xmax;
  worldBox.ymin = Conf.ymin;
  worldBox.ymax = Conf.ymax;
}

void reshape(int w, int h) {
  width = w;
  height = h;
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  double left = worldBox.xmin;
  double right = worldBox.xmax;
  double bottom = worldBox.ymin;
  double top = worldBox.ymax;
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
  if (worldW >= worldH) {
    worldH = worldW / aspect;
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * aspect;
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
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

  if (Conf.cells[ci].nodes[j].nextNode == null_size_t) {
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
    if (Conf.cells[i].close == false)
      continue;
    if (Conf.cells[i].p_int < pmin)
      pmin = Conf.cells[i].p_int;
    if (Conf.cells[i].p_int > pmax)
      pmax = Conf.cells[i].p_int;
  }

  color4f col;
  ColorTable pTable;
  pTable.setTableID(16);
  pTable.setSwap(true);
  pTable.setMinMax((float)pmin, (float)pmax);
  // std::cout << "pmax = " << pmax << '\n';
  std::cout << "pmin = " << pmin << '\n';
  std::cout << "pmax = " << pmax << '\n';
  // pTable.Rebuild();

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  for (size_t i = 1; i < Conf.cells.size(); ++i) {
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
  // glShadeModel(GL_SMOOTH);

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
          rnb = 12;
        } else if (type == 1) {
          glColor3f(1.f, 0.f, 0.f);
          rnb = 12;
        } else if (type == 2) {
          glColor3f(1.f, 0.f, 0.f);
          rnb = 12;
        } else if (type == 3) {
          glColor3f(1.f, 0.6f, 0.f);
          rnb = 4;
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
 * Draws the force, action, and reaction on a given position.
 *
 * @param InterIt  Neighbor object representing the interaction
 * @param pos      position where the force, action, and reaction will be drawn
 */
void drawForceActionReaction(const Neighbor &InterIt, vec2r &pos) {
  vec2r T(-InterIt.n.y, InterIt.n.x);
  vec2r vf = InterIt.n * InterIt.fn + T * InterIt.ft;
  vf *= vScale;
  if (InterIt.fn > 0.0) {
    arrow(pos.x - vf.x, pos.y - vf.y, pos.x, pos.y);
    arrow(pos.x + vf.x, pos.y + vf.y, pos.x, pos.y);
  } else {
    arrow(pos.x, pos.y, pos.x - vf.x, pos.y - vf.y);
    arrow(pos.x, pos.y, pos.x + vf.x, pos.y + vf.y);
  }
}

/**
 * Draws the forces between cells.
 */
void drawForces() {
  glLineWidth(2.0f);
  glShadeModel(GL_SMOOTH);

  glBegin(GL_LINES);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  for (size_t ci = 0; ci < Conf.cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = Conf.cells[ci].neighbors.begin();
         InterIt != Conf.cells[ci].neighbors.end(); ++InterIt) {

      size_t cj = InterIt->jc;
      size_t in = InterIt->in;
      size_t jn = InterIt->jn;

      if (Conf.cells[cj].nodes[jn].nextNode == null_size_t) {
        vec2r b = Conf.cells[cj].nodes[jn].pos - Conf.cells[ci].nodes[in].pos;
        b *= (Conf.cells[ci].radius / Conf.cells[ci].radius + Conf.cells[cj].radius);
        vec2r pos = Conf.cells[cj].nodes[jn].pos + b;
        drawForceActionReaction(*InterIt, pos);
      } else {
        size_t jnext = Conf.cells[cj].nodes[jn].nextNode;
        // jnext est le numéro du noeud à la fin de la barre dans la cellule cj (jn c'est le début)
        vec2r b = Conf.cells[ci].nodes[in].pos - Conf.cells[cj].nodes[jn].pos;
        vec2r u = Conf.cells[cj].nodes[jnext].pos - Conf.cells[cj].nodes[jn].pos;
        double u_length = u.normalize();
        double proj = b * u;

        double fact = Conf.cells[cj].radius / (Conf.cells[ci].radius + Conf.cells[cj].radius);
        if (proj <= 0.0) {
          vec2r ut = Conf.cells[ci].nodes[in].pos - Conf.cells[cj].nodes[jn].pos;
          vec2r pos = Conf.cells[cj].nodes[jn].pos + fact * ut;
          drawForceActionReaction(*InterIt, pos);
        } else if (proj >= u_length) {
          vec2r ut = Conf.cells[ci].nodes[in].pos - Conf.cells[cj].nodes[jnext].pos;
          vec2r pos = Conf.cells[cj].nodes[jn].pos + u_length * u + fact * ut;
          drawForceActionReaction(*InterIt, pos);
        } else {
          vec2r pos = Conf.cells[cj].nodes[jn].pos + proj * u;
          vec2r ut = Conf.cells[ci].nodes[in].pos - (Conf.cells[cj].nodes[jn].pos + proj * u);
          pos += fact * ut;
          drawForceActionReaction(*InterIt, pos);
        }
      }
    }
  }
  glEnd();
}

void try_to_readConf(int num, Lhyphen &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    OKNum = num;
    CF.loadCONF(file_name);
  } else
    std::cout << file_name << " does not exist" << std::endl;
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char *argv[]) {

  if (argc == 1) {
    confNum = 0;
  } else if (argc == 2) {
    confNum = atoi(argv[1]);
  }

  try_to_readConf(confNum, Conf, confNum);
  Conf.findDisplayArea(1.15);

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

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  mouse_mode = MouseMode::NOTHING;

  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
