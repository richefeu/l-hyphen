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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

class Lhyphen;
template <typename T> class exprParser;
using ExpressionParser = exprParser<double>; // même alias que dans Lhyphen.hpp

/// Un événement est vérifié à chaque pas de temps de la boucle de Lhyphen::integrate.
/// La méthode check met la variable 'active' à true quand l'événement doit se déclencher ;
/// la méthode action fait alors quelque chose (sauvegarder un fichier, modifier un paramètre,
/// arrêter la simulation, etc.). Un événement ne se déclenche qu'une seule fois (done = true ensuite).
///
/// Pour ajouter un nouveau type d'événement : dériver de Event, puis l'ajouter dans Event::create.
///
class Event {
public:
  bool active{false}; ///< mis à true par check quand l'événement doit se déclencher
  bool done{false};   ///< l'action a déjà été faite

  virtual ~Event() = default;

  virtual std::string name() const = 0;                         ///< mot-clé utilisé dans les fichiers d'entrée
  virtual void read(std::istream &is, ExpressionParser *ep) = 0; ///< lecture des paramètres (après le mot-clé)
  virtual void write(std::ostream &os) const = 0;               ///< écriture des paramètres (pour saveCONF)
  virtual void check(Lhyphen *lh) = 0;                          ///< met active à true si nécessaire
  virtual void action(Lhyphen *lh) = 0;                         ///< ce que fait l'événement

  static std::unique_ptr<Event> create(const std::string &name); ///< nullptr si le nom est inconnu
};

/// Sauvegarde une configuration (conf<iconf>) quand le temps atteint tEvent
///
/// event saveConfAtTime <tEvent>
///
class SaveConfAtTime : public Event {
public:
  double tEvent{0.0};

  std::string name() const override { return "saveConfAtTime"; }
  void read(std::istream &is, ExpressionParser *ep) override;
  void write(std::ostream &os) const override;
  void check(Lhyphen *lh) override;
  void action(Lhyphen *lh) override;
};

/// Sauvegarde une configuration (conf<iconf>) quand la longueur cumulée d'interfaces rompues
/// dépasse brokenLength (avec 0, on capture la toute première rupture)
///
/// event saveConfAtBrokenLength <brokenLength>
///
class SaveConfAtBrokenLength : public Event {
public:
  double brokenLength{0.0};

  std::string name() const override { return "saveConfAtBrokenLength"; }
  void read(std::istream &is, ExpressionParser *ep) override;
  void write(std::ostream &os) const override;
  void check(Lhyphen *lh) override;
  void action(Lhyphen *lh) override;
};

/// Sauvegarde une configuration puis arrête la simulation quand la longueur cumulée
/// d'interfaces rompues dépasse brokenLength
///
/// event stopAtBrokenLength <brokenLength>
///
class StopAtBrokenLength : public Event {
public:
  double brokenLength{0.0};

  std::string name() const override { return "stopAtBrokenLength"; }
  void read(std::istream &is, ExpressionParser *ep) override;
  void write(std::ostream &os) const override;
  void check(Lhyphen *lh) override;
  void action(Lhyphen *lh) override;
};

/// Arrête la simulation (après avoir sauvegardé une configuration) un certain temps 'delay' après
/// avoir détecté une chute de contrainte de 'dropPercent' %. La contrainte est mesurée par la force
/// de réaction (composante x ou y) sommée sur les noeuds pilotés par le control numéro 'ictrl'
/// (numérotés dans l'ordre des setNodeControl, setCellControl et setNodeControlInBox, à partir de 0).
/// La chute est comptée par rapport au pic : |F| <= (1 - dropPercent/100) * max|F|.
/// La force est lissée par une moyenne glissante exponentielle de temps caractéristique tau
/// (tau = 0 : pas de lissage), car la réaction présente des pics brefs dus aux chocs.
/// Le pic n'est suivi qu'à partir du temps tStart, pour ne pas être trompé par le régime
/// transitoire du début de chargement.
///
/// event stopAfterStressDrop <ictrl> <x|y> <dropPercent> <delay> <tStart> <tau>
///
class StopAfterStressDrop : public Event {
public:
  size_t ictrl{0};
  char component{'y'};
  double dropPercent{0.0};
  double delay{0.0};
  double tStart{0.0};
  double tau{0.0};

  double Fsmooth{0.0}; ///< force lissée
  double peak{0.0};    ///< max de |Fsmooth| atteint jusqu'ici
  bool dropDetected{false};
  double tDrop{0.0};  ///< temps auquel la chute a été détectée

  std::string name() const override { return "stopAfterStressDrop"; }
  void read(std::istream &is, ExpressionParser *ep) override;
  void write(std::ostream &os) const override;
  void check(Lhyphen *lh) override;
  void action(Lhyphen *lh) override;

private:
  bool nodesFound{false};
  std::vector<std::pair<size_t, size_t>> nodes; ///< (cellule, noeud) pilotés par ictrl
  double force(Lhyphen *lh);
};
