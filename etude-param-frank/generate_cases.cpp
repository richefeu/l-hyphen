// generate_cases : crée les répertoires de calcul d'une étude paramétrique l-hyphen.
//
// Lit un fichier de plan (plan.txt) :
//
//   model   input_model.txt            # modèle d'input (lignes « define NOM valeur »)
//   copy    echantillon/nodefile.txt   # fichier(s) copié(s) dans chaque cas (plusieurs lignes possibles)
//   output  cases                      # répertoire des cas
//   XI      0.1 0.5 1 5 10             # paramètre : NOM suivi de ses valeurs
//   ZETA    0.1 0.5 1 5 10
//
// Toutes les combinaisons des valeurs (produit cartésien) sont générées. Pour chaque cas, un répertoire
// <output>/XI0.1_ZETA0.5/ reçoit les fichiers à copier et un input.txt : copie du modèle dans laquelle seule la
// valeur des lignes « define NOM ... » des paramètres est remplacée (le commentaire de fin de ligne est gardé).
// Le fichier <output>/cases.txt liste les cas (numéro, répertoire, valeurs) pour le lancement et le
// dépouillement ; <output>/plan.txt garde une copie du plan utilisé.
//
// Un cas déjà existant n'est pas modifié (ses résultats sont préservés), sauf avec --force, qui réécrit
// input.txt et les fichiers copiés sans toucher aux résultats.
//
// Usage : generate_cases [plan.txt] [--force] [--dry-run]

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Plan {
  std::string model;
  std::vector<std::string> copies;
  std::string output{"cases"};
  std::vector<std::pair<std::string, std::vector<std::string>>> params; // dans l'ordre du fichier
};

static std::string stripComment(const std::string &line) {
  auto p = line.find('#');
  return p == std::string::npos ? line : line.substr(0, p);
}

static Plan readPlan(const std::string &name) {
  std::ifstream f(name);
  if (!f) {
    throw std::runtime_error("impossible d'ouvrir le plan " + name);
  }
  Plan P;
  std::string line;
  int ln = 0;
  while (std::getline(f, line)) {
    ln++;
    std::istringstream is(stripComment(line));
    std::string key;
    if (!(is >> key)) {
      continue;
    }
    std::vector<std::string> vals;
    for (std::string v; is >> v;) {
      vals.push_back(v);
    }
    auto need1 = [&]() {
      if (vals.size() != 1) {
        throw std::runtime_error(name + ":" + std::to_string(ln) + " : « " + key + " » attend une seule valeur");
      }
      return vals[0];
    };
    if (key == "model") {
      P.model = need1();
    } else if (key == "copy") {
      if (vals.empty()) {
        throw std::runtime_error(name + ":" + std::to_string(ln) + " : « copy » sans fichier");
      }
      P.copies.insert(P.copies.end(), vals.begin(), vals.end());
    } else if (key == "output") {
      P.output = need1();
    } else {
      if (vals.empty()) {
        throw std::runtime_error(name + ":" + std::to_string(ln) + " : le paramètre « " + key + " » n'a pas de valeur");
      }
      for (auto &p : P.params) {
        if (p.first == key) {
          throw std::runtime_error(name + ":" + std::to_string(ln) + " : paramètre « " + key + " » défini deux fois");
        }
      }
      P.params.push_back({key, vals});
    }
  }
  if (P.model.empty()) {
    throw std::runtime_error(name + " : mot-clé « model » manquant");
  }
  if (P.params.empty()) {
    throw std::runtime_error(name + " : aucun paramètre");
  }
  return P;
}

// Remplace la valeur de « define NOM valeur [# commentaire] » en gardant l'indentation et le commentaire.
// Retourne le nombre de lignes modifiées.
static int setDefine(std::vector<std::string> &lines, const std::string &name, const std::string &value) {
  int n = 0;
  for (auto &l : lines) {
    std::istringstream is(stripComment(l));
    std::string kw, nm, old;
    if (!(is >> kw >> nm >> old) || kw != "define" || nm != name) {
      continue;
    }
    auto c = l.find('#');
    std::string comment = (c == std::string::npos) ? "" : "      " + l.substr(c);
    l = "define " + name + "  " + value + comment;
    n++;
  }
  return n;
}

// Nom de répertoire lisible : XI0.1_ZETA0.5 (caractères gênants remplacés)
static std::string caseName(const std::vector<std::pair<std::string, std::string>> &combo) {
  std::string s;
  for (auto &kv : combo) {
    if (!s.empty()) {
      s += "_";
    }
    s += kv.first + kv.second;
  }
  for (auto &ch : s) {
    if (ch == '/' || ch == ' ' || ch == '$' || ch == '*' || ch == '(' || ch == ')' || ch == '\\') {
      ch = '-';
    }
  }
  return s;
}

int main(int argc, char *argv[]) {
  std::string planFile = "plan.txt";
  bool force = false, dryRun = false;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--force") {
      force = true;
    } else if (a == "--dry-run") {
      dryRun = true;
    } else if (a == "-h" || a == "--help") {
      std::cout << "Usage : generate_cases [plan.txt] [--force] [--dry-run]\n";
      return 0;
    } else {
      planFile = a;
    }
  }

  try {
    Plan P = readPlan(planFile);
    // les chemins du plan sont relatifs au répertoire du plan
    fs::path base = fs::absolute(planFile).parent_path();
    auto resolve = [&](const std::string &p) { return fs::path(p).is_absolute() ? fs::path(p) : base / p; };

    std::ifstream mf(resolve(P.model));
    if (!mf) {
      throw std::runtime_error("impossible d'ouvrir le modèle " + P.model);
    }
    std::vector<std::string> model;
    for (std::string l; std::getline(mf, l);) {
      model.push_back(l);
    }
    for (auto &p : P.params) { // chaque paramètre doit exister dans le modèle
      std::vector<std::string> tmp = model;
      if (setDefine(tmp, p.first, "0") != 1) {
        throw std::runtime_error("le modèle " + P.model + " doit contenir exactement une ligne « define " + p.first +
                                 " ... »");
      }
    }
    for (auto &c : P.copies) {
      if (!fs::exists(resolve(c))) {
        throw std::runtime_error("fichier à copier introuvable : " + c);
      }
    }

    // produit cartésien des valeurs
    size_t nCases = 1;
    for (auto &p : P.params) {
      nCases *= p.second.size();
    }
    fs::path out = resolve(P.output);
    if (!dryRun) {
      fs::create_directories(out);
    }

    std::ostringstream index;
    index << "# cas générés par generate_cases à partir de " << planFile << "\n# id  répertoire";
    for (auto &p : P.params) {
      index << "  " << p.first;
    }
    index << "\n";

    size_t nNew = 0, nKept = 0, nForced = 0;
    std::vector<size_t> idx(P.params.size(), 0);
    for (size_t id = 1; id <= nCases; id++) {
      std::vector<std::pair<std::string, std::string>> combo;
      for (size_t k = 0; k < P.params.size(); k++) {
        combo.push_back({P.params[k].first, P.params[k].second[idx[k]]});
      }
      std::string name = caseName(combo);
      fs::path dir = out / name;
      index << std::setw(4) << id << "  " << name;
      for (auto &kv : combo) {
        index << "  " << kv.second;
      }
      index << "\n";

      bool exists = fs::exists(dir / "input.txt");
      if (exists && !force) {
        nKept++;
      } else if (!dryRun) {
        fs::create_directories(dir);
        std::vector<std::string> lines = model;
        for (auto &kv : combo) {
          setDefine(lines, kv.first, kv.second);
        }
        std::ofstream in(dir / "input.txt");
        in << "# Cas " << name << " : généré par generate_cases (plan " << planFile << ", modèle " << P.model
           << ")\n";
        for (auto &l : lines) {
          in << l << "\n";
        }
        for (auto &c : P.copies) {
          fs::copy_file(resolve(c), dir / fs::path(c).filename(), fs::copy_options::overwrite_existing);
        }
        (exists ? nForced : nNew)++;
      } else {
        (exists ? nForced : nNew)++;
      }

      // combinaison suivante (le dernier paramètre varie le plus vite)
      for (size_t k = P.params.size(); k-- > 0;) {
        if (++idx[k] < P.params[k].second.size()) {
          break;
        }
        idx[k] = 0;
      }
    }

    if (!dryRun) {
      std::ofstream(out / "cases.txt") << index.str();
      fs::copy_file(resolve(planFile), out / "plan.txt", fs::copy_options::overwrite_existing);
    }
    std::cout << (dryRun ? "[simulation] " : "") << nCases << " cas dans " << P.output << "/ : " << nNew << " créés, "
              << nForced << " réécrits (--force), " << nKept << " existants conservés\n";
    if (dryRun) {
      std::cout << index.str();
    } else {
      std::cout << "liste des cas : " << P.output << "/cases.txt\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "generate_cases : " << e.what() << "\n";
    return 1;
  }
  return 0;
}
