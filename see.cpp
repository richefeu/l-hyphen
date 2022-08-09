#include "see.hpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next configuration file" << endl;
  cout << "-         load previous configuration file" << endl;
  cout << "=         fit the view" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void printInfo() {
  using namespace std;

  // cout << "Reference Conf = " << refConfNum << "\n";
  // cout << "Current Conf = " << confNum << "\n";
  //
  // cout << "Box:\n";
  // cout << "top left = (" << Conf.Box.x[0] << ", " << Conf.Box.y[0] << ")\n";
  // cout << "top right = (" << Conf.Box.x[1] << ", " << Conf.Box.y[1] << ")\n";
  // cout << "bottom right = (" << Conf.Box.x[2] << ", " << Conf.Box.y[2] << ")\n";
  // cout << "bottom left = (" << Conf.Box.x[3] << ", " << Conf.Box.y[3] << ")\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case 'i': {
    printInfo();
  } break;
   
  case 'q': {
    exit(0);
  } break;

  case '-': {
    std::cout << "Current Configuration: ";
    if (confNum > 0)
      try_to_readConf(confNum - 1, Conf, confNum);
  } break;

  case '+': {
    std::cout << "Current Configuration: ";
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
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON:
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        mouse_mode = PAN;
      else
        mouse_mode = ROTATION;
      break;
    case GLUT_MIDDLE_BUTTON:
      mouse_mode = ZOOM;
      break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING)
    return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case ZOOM: {
    double ddy = (worldBox.ymax - worldBox.ymin) * dy;
    double ddx = (worldBox.xmax - worldBox.xmin) * dy;
    worldBox.xmin -= ddx;
    worldBox.xmax += ddx;
    worldBox.ymin -= ddy;
    worldBox.ymax += ddy;
  } break;

  case PAN: {
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
  display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  drawCells();

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



void drawBar(size_t ci, size_t i, size_t j, double radius) {
  if (i == null_size_t || j == null_size_t)
    return;

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

  glBegin(GL_POLYGON);
  glVertex2f(xi - radius * txij, yi - radius * tyij);
  glVertex2f(xj - radius * txij, yj - radius * tyij);
  glVertex2f(xj + radius * txij, yj + radius * tyij);
  glVertex2f(xi + radius * txij, yi + radius * tyij);
  glEnd();

  glBegin(GL_POLYGON);
  for (double a = 0.0; a < 2.0 * M_PI; a += M_PI / 18.0) {
    glVertex2f(xi + radius * cos(a), yi + radius * sin(a));
  }
  glEnd();

  if (Conf.cells[ci].nodes[j].nextNode == null_size_t) {
    glBegin(GL_POLYGON);
    for (double a = 0.0; a < 2.0 * M_PI; a += M_PI / 18.0) {
      glVertex2f(xj + radius * cos(a), yj + radius * sin(a));
    }
    glEnd();
  }
}

void drawCells() {
  glLineWidth(4.0f);
  glColor4f(0.8f, 0.8f, 0.9f, 1.0f);

  for (size_t i = 0; i < Conf.cells.size(); ++i) {
    for (size_t n = 0; n < Conf.cells[i].nodes.size(); ++n) {
      size_t nextNode = Conf.cells[i].nodes[n].nextNode;
      drawBar(i, n, nextNode, Conf.cells[i].radius);
    }
  }
}

/// Robust and portable function to test if a file exists
bool fileExists(const char *fileName) {
  std::fstream fin;
  fin.open(fileName, std::ios::in);
  if (fin.is_open()) {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}

void try_to_readConf(int num, Lhyphen &CF, int &OKNum) {
  char file_name[256];
  sprintf(file_name, "conf%d", num);
  if (fileExists(file_name)) {
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

  std::cout << "Current Configuration: ";
  try_to_readConf(confNum, Conf, confNum);
  Conf.findDisplayArea(1.1);

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


  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
