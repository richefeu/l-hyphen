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
#include <sstream>

#include "Lhyphen.hpp"
#include "null_size_t.hpp"

#include "AABB.hpp"
#include "ColorTable.hpp"
#include "fileTool.hpp"
#include "glTools.hpp"
#include "message.hpp"
#include "triangulatePolygon.hpp"

#include "toofus/toofus-gate/toml++/toml.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "toofus/toofus-gate/stb/stb_image_write.h"

Lhyphen Conf;
int confNum = 1;

AABB worldBox;

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

ColorTable BarRedTable;
ColorTable BarBlueTable;
ColorTable NodeRedTable;
ColorTable NodeBlueTable;

int main_window;

// global window pointer (for updateTextLine title update)
GLFWwindow *g_window{nullptr};

// redraw flag (avoids busy-loop CPU waste)
bool needsRedraw{true};

// flags
int show_cells = 1;
int show_glue_points = 0;
int show_bar_colors = 0;
int show_inter_cells_forces = 0;
int show_pressure = 0;
int show_contours = 1;
int show_nodes = 0;
int show_control_boxes = 0;
int show_background = 0;
int show_crack_path = 0; // trace les liens rompus jusqu'au temps du conf affiché
int show_velocities = 0; // flèches de vitesse aux noeuds ('e'), échelle vScale ('y'/'Y')
int show_strain = 0;      // couleur des cellules fermées selon la déformation ('k') : voir strainModeNames
int show_strain_dirs = 0; // directions principales de déformation ('j')
int show_stress = 0;      // couleur des cellules fermées selon la contrainte ('l') : voir stressModeNames
int show_stress_dirs = 0; // directions principales de contrainte ('m')
int show_hud = 0;        // panneau d'état permanent des toggles ('i')
int show_help = 0;       // overlay d'aide des raccourcis clavier ('h')

// arrow/force sizes
double arrowSize = 0.25;  // longueur des barbules, en fraction de la longueur de la flèche
double arrowAngle = 0.35; // demi-angle des barbules [rad]
double vScale = 1.0;      // longueur de la plus grande flèche de vitesse, en rayons de cellule moyens ('y'/'Y')
double fnWidthFactor = 1.0; // multiplicateur d'épaisseur des chaînes de force ('s'/'S')
double forceFilter = 1.0;   // seuil du filtre = forceFilter * |fn|_moyen : ne montre que les chaînes porteuses ('t'/'T')
double eScale = 0.8;         // demi-longueur du plus grand trait de direction principale, en rayons de cellule ('u'/'U')
double strainColorMax = 0.0; // borne de l'échelle de couleur des déformations (0 = automatique, pour chaque conf)
double stressColorMax = 0.0; // borne de l'échelle de couleur des contraintes (0 = automatique, pour chaque conf)

/// Tenseur symétrique 2D d'une cellule fermée, par ses valeurs et directions principales
/// (v1 >= v2, u1 et u2 unitaires, dans la conf affichée). Convention tension positive.
/// - déformation : v_i = ln(lambda_i) (Hencky), à partir du gradient de transformation moyen F de la
///   cellule entre la conf de référence (RefConf) et la conf affichée ;
///   eps_v = v1 + v2 = ln(J) et eps_q = v1 - v2.
/// - contrainte : sigma = (1/A) sum (x_c - centre) (x) f_c, sur les forces de contact et de cohésion
///   exercées par les autres cellules (Love-Weber) ; sig_m = (v1 + v2)/2 et sig_q = v1 - v2.
struct CellTensor {
  bool ok{false}; // false si cellule ouverte, ou tenseur non calculable
  vec2r center;   // centre (moyenne des noeuds) dans la conf affichée
  double R{0.0};  // rayon moyen (distance moyenne des noeuds au centre)
  double v1{0.0}, v2{0.0};
  vec2r u1, u2;
  double xx{0.0}, yy{0.0}, xy{0.0}; // composantes dans le repère (x, y)
};
std::vector<CellTensor> cellStrains;
std::vector<CellTensor> cellStresses;
Lhyphen RefConf;           // configuration de référence pour les déformations
int refConfNum = 0;        // son numéro ('o' : conf affichée, Shift+o : conf0)
int loadedRefConfNum = -1; // numéro de la conf de référence effectivement chargée (-1 : aucune)
ColorTable TensorSphTable; // divergente bleu-blanc-rouge pour les grandeurs signées (eps_v, sig_m, p, sig_xx...)
ColorTable TensorDevTable; // blanc-jaune-rouge pour les grandeurs positives (eps_q, sig_q)

// Modes de couleur ('k' pour la déformation, 'l' pour la contrainte), 0 = rien.
// p est la pression interne p_int de la cellule (positive quand elle tend à la faire gonfler).
const char *strainModeNames[] = {"off", "eps_v", "eps_q", "eps_xx", "eps_yy", "eps_xy"};
const char *stressModeNames[] = {"off", "sig_m", "sig_q", "p", "p+sig_m", "sig_xx", "sig_yy", "sig_xy"};
const int nbStrainModes = 6;
const int nbStressModes = 8;

/// Un scalaire par cellule à afficher en couleur
struct CellScalarField {
  std::vector<double> value;
  std::vector<char> ok;  // la cellule est-elle coloriée ?
  bool divergent{true};  // échelle symétrique bleu-blanc-rouge (sinon 0..max, blanc-jaune-rouge)
};

// window sizes
int width = 800;
int height = 600;
float wh_ratio = (float)width / (float)height;
glTextZone textZone(1, &width, &height);
int fit_at_loading{1};

// background gradient colors (RGB 0-255)
int bottom_r{135}, bottom_g{206}, bottom_b{250};
int top_r{255}, top_g{255}, top_b{255};

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
void drawCrackPath();
void readBreakHistory(const char *fname = "breakHistory.txt");
void arrow(double x0, double y0, double x1, double y1);

void drawForces();
void drawVelocities();
void drawPressure();
void drawControlBoxes();
bool loadRefConf();
void computeStrains();
void computeStresses();
CellScalarField strainField(int mode);
CellScalarField stressField(int mode);
double fieldColorBound(const CellScalarField &field, double fixedMax);
void drawCellScalars(const CellScalarField &field, double vmax);
void drawTensorDirections(const std::vector<CellTensor> &tensors);
void drawTensorColorBar(const char *name, double vmin, double vmax, bool autoBound, const char *extra);
void drawHUD();
void drawHelpOverlay();

void updateTextLine();

// Callback functions
void keyboard(GLFWwindow *window, int key, int scancode, int action, int mods);
void mouse_button(GLFWwindow *window, int button, int action, int mods);
void cursor_pos(GLFWwindow *window, double xpos, double ypos);
void display(GLFWwindow *window);
void reshape(GLFWwindow *window, int width, int height);
void framebuffer_size(GLFWwindow *window, int width, int height);

// Screenshot
void captureScreenshot(const char *filename);

// Option file
void readTomlOptions();
void saveTomlOptions();

// Helper functions
void printHelp();
void fit_view(GLFWwindow *window);
bool try_to_readConf(int num, Lhyphen &conf, int &OKNum);
