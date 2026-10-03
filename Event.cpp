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

#include "Event.hpp"
#include "Lhyphen.hpp"

std::unique_ptr<Event> Event::create(const std::string &name) {
  if (name == "saveConfAtTime") {
    return std::make_unique<SaveConfAtTime>();
  } else if (name == "saveConfAtBrokenLength") {
    return std::make_unique<SaveConfAtBrokenLength>();
  } else if (name == "stopAtBrokenLength") {
    return std::make_unique<StopAtBrokenLength>();
  } else if (name == "stopAfterStressDrop") {
    return std::make_unique<StopAfterStressDrop>();
  }
  return nullptr;
}

// ------------------------------------------------------------------------------------------------------

void SaveConfAtTime::read(std::istream &is, ExpressionParser *ep) { ep->getValue(is, tEvent); }

void SaveConfAtTime::write(std::ostream &os) const { os << tEvent; }

void SaveConfAtTime::check(Lhyphen *lh) {
  if (lh->t >= tEvent) {
    active = true;
  }
}

void SaveConfAtTime::action(Lhyphen *lh) {
  std::cout << "* Event saveConfAtTime (t = " << lh->t << ")" << std::endl;
  lh->saveCONF(lh->iconf);
  lh->iconf++;
}

// ------------------------------------------------------------------------------------------------------

void SaveConfAtBrokenLength::read(std::istream &is, ExpressionParser *ep) { ep->getValue(is, brokenLength); }

void SaveConfAtBrokenLength::write(std::ostream &os) const { os << brokenLength; }

void SaveConfAtBrokenLength::check(Lhyphen *lh) {
  if (lh->cumulatedL > brokenLength) {
    active = true;
  }
}

void SaveConfAtBrokenLength::action(Lhyphen *lh) {
  std::cout << "* Event saveConfAtBrokenLength (t = " << lh->t << ", broken length = " << lh->cumulatedL << ")"
            << std::endl;
  lh->saveCONF(lh->iconf);
  lh->iconf++;
}

// ------------------------------------------------------------------------------------------------------

void StopAtBrokenLength::read(std::istream &is, ExpressionParser *ep) { ep->getValue(is, brokenLength); }

void StopAtBrokenLength::write(std::ostream &os) const { os << brokenLength; }

void StopAtBrokenLength::check(Lhyphen *lh) {
  if (lh->cumulatedL > brokenLength) {
    active = true;
  }
}

void StopAtBrokenLength::action(Lhyphen *lh) {
  std::cout << "* Event stopAtBrokenLength (t = " << lh->t << ", broken length = " << lh->cumulatedL
            << "): end of simulation" << std::endl;
  lh->saveCONF(lh->iconf);
  lh->iconf++;
  lh->stopRequested = true;
}

// ------------------------------------------------------------------------------------------------------

void StopAfterStressDrop::read(std::istream &is, ExpressionParser *ep) {
  ep->getValue(is, ictrl);
  is >> component;
  if (component != 'x' && component != 'y') {
    std::cout << "@StopAfterStressDrop::read, component must be x or y (y is used)" << std::endl;
    component = 'y';
  }
  ep->getValue(is, dropPercent);
  ep->getValue(is, delay);
  ep->getValue(is, tStart);
  ep->getValue(is, tau);
}

void StopAfterStressDrop::write(std::ostream &os) const {
  os << ictrl << ' ' << component << ' ' << dropPercent << ' ' << delay << ' ' << tStart << ' ' << tau;
}

// Somme de la composante choisie des forces sur les noeuds pilotés par ictrl
// (la liste des noeuds est construite au premier appel)
double StopAfterStressDrop::force(Lhyphen *lh) {
  if (!nodesFound) {
    for (size_t c = 0; c < lh->cells.size(); c++) {
      for (size_t n = 0; n < lh->cells[c].nodes.size(); n++) {
        if (lh->cells[c].nodes[n].ictrl == ictrl) {
          nodes.push_back(std::make_pair(c, n));
        }
      }
    }
    nodesFound = true;
    if (nodes.empty()) {
      std::cout << "@StopAfterStressDrop, no node is driven by the control " << ictrl << std::endl;
    }
  }

  double F = 0.0;
  for (size_t i = 0; i < nodes.size(); i++) {
    const vec2r &f = lh->cells[nodes[i].first].nodes[nodes[i].second].force;
    F += (component == 'x') ? f.x : f.y;
  }
  return F;
}

void StopAfterStressDrop::check(Lhyphen *lh) {
  if (!dropDetected) {
    // le lissage commence dès le début, pour être établi à tStart
    double F = force(lh);
    if (tau > 0.0) {
      Fsmooth += (F - Fsmooth) * std::min(1.0, lh->dt / tau);
    } else {
      Fsmooth = F;
    }
  }

  if (!dropDetected && lh->t >= tStart) {
    double F = fabs(Fsmooth);
    if (F > peak) {
      peak = F;
    } else if (peak > 0.0 && F <= (1.0 - 0.01 * dropPercent) * peak) {
      dropDetected = true;
      tDrop = lh->t;
      std::cout << "* Event stopAfterStressDrop: stress drop detected at t = " << tDrop << " (|F| = " << F
                << ", peak = " << peak << ", smoothed values), stop at t = " << tDrop + delay << std::endl;
    }
  }

  if (dropDetected && lh->t >= tDrop + delay) {
    active = true;
  }
}

void StopAfterStressDrop::action(Lhyphen *lh) {
  std::cout << "* Event stopAfterStressDrop (t = " << lh->t << "): end of simulation" << std::endl;
  lh->saveCONF(lh->iconf);
  lh->iconf++;
  lh->stopRequested = true;
}
