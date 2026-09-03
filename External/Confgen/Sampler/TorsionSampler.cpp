//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Sampler/TorsionSampler.h"
#include "Sampler/TorsionSymmetries.h"

#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/RDLog.h>

#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>

namespace RDKit {
namespace {

double normalizeAngle(double angle) {
  angle = std::fmod(angle, 360.0);
  if (angle < 0.0) {
    angle += 360.0;
  }
  return angle;
}

double angleDistance(double a, double b) {
  double diff = std::abs(normalizeAngle(a) - normalizeAngle(b));
  return std::min(diff, 360.0 - diff);
}

bool addUniqueAngle(std::vector<double> &angles, double angle,
                    double minSeparation) {
  angle = normalizeAngle(angle);
  for (const auto existing : angles) {
    if (angleDistance(existing, angle) < minSeparation) {
      return false;
    }
  }
  angles.push_back(angle);
  return true;
}

double calcETKDGTorsionEnergy(const std::vector<double> &forceConstants,
                              const std::vector<int> &signs, double cosPhi) {
  if (forceConstants.size() < 6 || signs.size() < 6) {
    return std::numeric_limits<double>::max();
  }
  const double cosPhi2 = cosPhi * cosPhi;
  const double cosPhi3 = cosPhi * cosPhi2;
  const double cosPhi4 = cosPhi * cosPhi3;
  const double cosPhi5 = cosPhi * cosPhi4;
  const double cosPhi6 = cosPhi * cosPhi5;

  const double cos2Phi = 2.0 * cosPhi2 - 1.0;
  const double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
  const double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
  const double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
  const double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;

  return forceConstants[0] * (1.0 + signs[0] * cosPhi) +
         forceConstants[1] * (1.0 + signs[1] * cos2Phi) +
         forceConstants[2] * (1.0 + signs[2] * cos3Phi) +
         forceConstants[3] * (1.0 + signs[3] * cos4Phi) +
         forceConstants[4] * (1.0 + signs[4] * cos5Phi) +
         forceConstants[5] * (1.0 + signs[5] * cos6Phi);
}

std::vector<double> uniformAngles(double stepDeg) {
  std::vector<double> angles;
  if (stepDeg <= 0.0) {
    angles.push_back(0.0);
    return angles;
  }
  for (double a = 0.0; a < 360.0 - 1e-9; a += stepDeg) {
    angles.push_back(a);
  }
  return angles;
}

std::string xmlAttr(const std::string &tag, const std::string &name) {
  const std::string needle = name + "=\"";
  const auto start = tag.find(needle);
  if (start == std::string::npos) {
    return "";
  }
  const auto valueStart = start + needle.size();
  const auto valueEnd = tag.find('"', valueStart);
  if (valueEnd == std::string::npos) {
    return "";
  }
  return tag.substr(valueStart, valueEnd - valueStart);
}

bool isMappedCentralBond(const MatchVectType &match,
                         const std::array<unsigned int, 4> &queryAtoms,
                         unsigned int j, unsigned int k,
                         std::array<unsigned int, 4> &libAtoms) {
  std::array<int, 4> atomMap = {{-1, -1, -1, -1}};
  for (const auto &pair : match) {
    for (unsigned int idx = 0; idx < queryAtoms.size(); ++idx) {
      if (static_cast<unsigned int>(pair.first) == queryAtoms[idx]) {
        atomMap[idx] = pair.second;
      }
    }
  }
  for (const auto mapped : atomMap) {
    if (mapped < 0) {
      return false;
    }
  }

  if (atomMap[1] == static_cast<int>(j) && atomMap[2] == static_cast<int>(k)) {
    for (unsigned int idx = 0; idx < atomMap.size(); ++idx) {
      libAtoms[idx] = static_cast<unsigned int>(atomMap[idx]);
    }
    return true;
  }
  if (atomMap[1] == static_cast<int>(k) && atomMap[2] == static_cast<int>(j)) {
    libAtoms[0] = static_cast<unsigned int>(atomMap[3]);
    libAtoms[1] = static_cast<unsigned int>(atomMap[2]);
    libAtoms[2] = static_cast<unsigned int>(atomMap[1]);
    libAtoms[3] = static_cast<unsigned int>(atomMap[0]);
    return true;
  }
  return false;
}

//! Unordered element-pair key for the central-bond discrimination index (bond
//! is matched in both orientations at verify time, so the pair is
//! order-independent).
inline int centralElementPairKey(int a, int b) {
  return a < b ? a * 256 + b : b * 256 + a;
}

//! ANCHORED central-bond match.  For speed, force bond j-k to always be
//!  included in the search via SubstructMatchParameters::extraAtomCheck
template <typename F>
bool matchCentralBondAnchored(const ROMol &mol, const ROMol &query,
                              const std::array<unsigned int, 4> &queryAtoms,
                              unsigned int j, unsigned int k, bool collectAll,
                              F &&onMatch) {
  const unsigned int cJ = queryAtoms[1], cK = queryAtoms[2];
  SubstructMatchParameters ps;
  ps.uniquify = true;
  ps.maxMatches = collectAll ? 1024u : 8u;
  ps.extraAtomCheck = [cJ, cK, j, k, &query](const Atom &qa, const Atom &ma) {
    // If we are in a recusive query, disable the check
    if (&qa.getOwningMol() != &query) return true;
    const unsigned int qi = qa.getIdx();
    if (qi == cJ || qi == cK) {
      const unsigned int mi = ma.getIdx();
      return mi == j || mi == k;  // central atoms pinned to the known bond
    }
    return true;  // not a central atom: no extra constraint
  };
  bool any = false;
  std::array<unsigned int, 4> ruleAtoms{};
  for (const auto &m : SubstructMatch(mol, query, ps)) {
    if (!isMappedCentralBond(m, queryAtoms, j, k, ruleAtoms)) continue;
    any = true;
    onMatch(ruleAtoms);
    if (!collectAll) break;
  }
  return any;
}

//! Sample torsions around a sample atom switching if
//!  the reference geometry is planar: aromatic / conjugated), sp2/aromatic ->
//!  180.
double getAngleOffsetDeg(const ROMol &mol, unsigned int atomIdx) {
  const Atom *a = mol.getAtomWithIdx(atomIdx);
  const auto hyb = a->getHybridization();
  if (hyb == Atom::SP3) {
    // A planar sp3 nitrogen (amide / aromatic-adjacent) behaves like sp2 ->
    // 180.
    if (a->getIsAromatic()) return 180.0;
    if (a->getAtomicNum() == 7) {
      for (const auto nbr : mol.atomNeighbors(a)) {
        if (nbr->getIsAromatic() || nbr->getHybridization() == Atom::SP2) {
          return 180.0;
        }
      }
    }
    return 120.0;
  }
  if (hyb == Atom::SP2) return 180.0;
  return 0.0;
}

}  // namespace

std::vector<double> UniformTorsionSampler::getAngles(const ROMol &,
                                                     unsigned int, unsigned int,
                                                     unsigned int, unsigned int,
                                                     bool) const {
  return uniformAngles(d_stepDeg);
}

std::vector<double> TorsionSampler::foldSymmetricAngles(
    const ROMol &mol, unsigned int j, unsigned int k,
    std::vector<double> angles, double minSeparation, bool useSmarts) {
  const unsigned int rotSym = useSmarts
                                  ? torsionRotationalSymmetryBySmarts(mol, j, k)
                                  : torsionRotationalSymmetryByRanks(mol, j, k);
  if (rotSym > 1) {
    angles = foldAnglesByTorsionSymmetry(rotSym, angles, minSeparation);
  }
  return angles;
}

std::vector<double> ETKDGTorsionSampler::getAngles(
    const ROMol &mol, unsigned int i, unsigned int j, unsigned int k,
    unsigned int l, bool /*addBasinAngle*/) const {
  const double scanStep = d_scanStepDeg > 0.0 ? d_scanStepDeg : 10.0;
  const double minSeparation = std::max(10.0, scanStep * 1.5);
  //! remove redundant and similar torsions
  auto folded = [&](std::vector<double> a) {
    return foldBySymmetry
               ? foldSymmetricAngles(mol, j, k, std::move(a), minSeparation)
               : a;
  };
  std::vector<double> fallback = uniformAngles(d_fallbackStepDeg);
  if (!mol.getNumConformers()) {
    return folded(fallback);
  }

  const Bond *bond = mol.getBondBetweenAtoms(j, k);
  if (!bond) {
    return folded(fallback);
  }

  ForceFields::CrystalFF::CrystalFFDetails details;
  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const ForceFields::CrystalFF::ExpTorsionAngle *>>
      torsionBonds;
  ForceFields::CrystalFF::getExperimentalTorsions(
      mol, details, torsionBonds, true, d_useSmallRingTorsions,
      d_useMacrocycleTorsions, false, d_version, false);

  const auto matchIt = std::find_if(
      torsionBonds.begin(), torsionBonds.end(),
      [&](const auto &match) { return std::get<0>(match) == bond->getIdx(); });
  if (matchIt == torsionBonds.end()) {
    return folded(fallback);
  }

  auto ruleAtoms = std::get<1>(*matchIt);
  const auto *params = std::get<2>(*matchIt);
  if (!params || ruleAtoms.size() != 4) {
    return folded(fallback);
  }
  if (ruleAtoms[1] == k && ruleAtoms[2] == j) {
    std::reverse(ruleAtoms.begin(), ruleAtoms.end());
  }
  if (ruleAtoms[1] != j || ruleAtoms[2] != k) {
    return folded(fallback);
  }

  const Conformer &conf = mol.getConformer();
  const double currentTarget = MolTransforms::getDihedralDeg(conf, i, j, k, l);
  const double currentAngle = MolTransforms::getDihedralDeg(
      conf, ruleAtoms[0], ruleAtoms[1], ruleAtoms[2], ruleAtoms[3]);
  const double frameOffset = currentTarget - currentAngle;

  struct Candidate {
    double energy;
    double angle;
  };
  std::vector<Candidate> candidates;
  for (double angle = 0.0; angle < 360.0 - 1e-9; angle += scanStep) {
    const double rad = angle * M_PI / 180.0;
    candidates.push_back(Candidate{
        calcETKDGTorsionEnergy(params->V, params->signs, std::cos(rad)),
        normalizeAngle(angle + frameOffset)});
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &a, const Candidate &b) {
              if (a.energy != b.energy) {
                return a.energy < b.energy;
              }
              return a.angle < b.angle;
            });

  // Priority order: the library's own angles best-energy-first, then the
  // uniform grid as a
  //  sampling backstop.
  std::vector<double> ordered;
  ordered.reserve(candidates.size() + fallback.size());
  for (const auto &candidate : candidates) {
    ordered.push_back(candidate.angle);
  }
  ordered.insert(ordered.end(), fallback.begin(), fallback.end());

  // Symmetry fold angles
  ordered = folded(std::move(ordered));

  std::vector<double> angles;
  // 0 == no limit
  const unsigned int maxAngles =
      d_maxAngles ? d_maxAngles : std::numeric_limits<unsigned int>::max();
  for (const double angle : ordered) {
    if (angles.size() >= maxAngles) {
      break;
    }
    addUniqueAngle(angles, angle, minSeparation);
  }
  return angles;
}

struct TorsionLibrarySampler::Rule {
  struct Angle {
    double value{0.0};
    double tolerance1{0.0};
    double tolerance2{0.0};
    double score{0.0};
  };

  std::string pattern;
  std::shared_ptr<ROMol> query;
  std::array<unsigned int, 4> queryAtoms{{0, 0, 0, 0}};
  std::vector<Angle> angles;
};

//! Holds the rule table
//!  see acquireRuleTable.
struct TorsionLibrarySampler::RuleTable {
  std::vector<Rule> rules;
  //! For unambiguous queries
  std::unordered_map<int, std::vector<std::size_t>> unambiguousRules;
  //! For ambiguous queries, i.e. OR, [C,c], wildcards etc.
  std::vector<std::size_t> ambiguousRules;
  TorsionSamplerStatus status{TorsionSamplerStatus::Ok};
};

const char *torsionSamplerStatusMessage(TorsionSamplerStatus status) {
  switch (status) {
    case TorsionSamplerStatus::Ok:
      return "ok";
    case TorsionSamplerStatus::NoPathGiven:
      return "no library path was given";
    case TorsionSamplerStatus::FileUnreadable:
      return "could not open library file";
    case TorsionSamplerStatus::NoRulesParsed:
      return "failed to parse rules from the library";
  }
  return "unknown";
}

TorsionSamplerStatus TorsionLibrarySampler::status() const {
  return d_table ? d_table->status : TorsionSamplerStatus::NoRulesParsed;
}

std::shared_ptr<const TorsionLibrarySampler::RuleTable>
TorsionLibrarySampler::acquireRuleTable(const std::string &torsionLibraryPath) {
  // Return a rule table, creating if necessary.
  //  This is a thread safe function.
  static std::mutex cacheMutex;
  static std::map<std::string, std::shared_ptr<const RuleTable>> cache;
  {
    const std::lock_guard<std::mutex> lock(cacheMutex);
    const auto it = cache.find(torsionLibraryPath);
    if (it != cache.end()) {
      return it->second;
    }
  }

  auto table = std::make_shared<RuleTable>();
  auto *rules = &table->rules;
  bool readable = false;
  auto parseFile = [&](const std::string &path) {
    if (path.empty()) {
      table->status = TorsionSamplerStatus::NoPathGiven;
      return;
    }
    std::ifstream in(path);
    if (!in.good()) {
      table->status = TorsionSamplerStatus::FileUnreadable;
      return;
    }
    readable = true;
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string xml = ss.str();

    size_t pos = 0;
    RDLog::LogStateSetter blocker;
    while ((pos = xml.find("<rule", pos)) != std::string::npos) {
      const auto tagEnd = xml.find('>', pos);
      if (tagEnd == std::string::npos) {
        break;
      }
      const auto end = xml.find("</rule>", tagEnd);
      if (end == std::string::npos) {
        break;
      }
      const std::string tag = xml.substr(pos, tagEnd - pos + 1);
      const std::string pattern = xmlAttr(tag, "pattern");
      const std::string block = xml.substr(tagEnd + 1, end - tagEnd - 1);
      pos = end + 7;
      if (pattern.empty()) {
        continue;
      }

      std::unique_ptr<ROMol> query(SmartsToMol(pattern));
      if (!query) {
        continue;
      }

      std::array<unsigned int, 4> queryAtoms{{0, 0, 0, 0}};
      std::array<bool, 4> found{{false, false, false, false}};
      for (const auto atom : query->atoms()) {
        const int mapNum = atom->getAtomMapNum();
        if (mapNum >= 1 && mapNum <= 4) {
          queryAtoms[mapNum - 1] = atom->getIdx();
          found[mapNum - 1] = true;
        }
      }
      if (!found[0] || !found[1] || !found[2] || !found[3]) {
        continue;
      }

      Rule rule;
      rule.pattern = pattern;
      rule.query.reset(query.release());
      rule.queryAtoms = queryAtoms;

      size_t anglePos = 0;
      while ((anglePos = block.find("<angle", anglePos)) != std::string::npos) {
        const auto angleEnd = block.find("/>", anglePos);
        if (angleEnd == std::string::npos) {
          break;
        }
        const std::string angleTag =
            block.substr(anglePos, angleEnd - anglePos + 2);
        anglePos = angleEnd + 2;
        const std::string value = xmlAttr(angleTag, "value");
        if (value.empty()) {
          continue;
        }
        Rule::Angle angle;
        angle.value = std::atof(value.c_str());
        angle.score = std::atof(xmlAttr(angleTag, "score").c_str());
        angle.tolerance1 = std::atof(xmlAttr(angleTag, "tolerance1").c_str());
        angle.tolerance2 = std::atof(xmlAttr(angleTag, "tolerance2").c_str());
        rule.angles.push_back(angle);
      }
      if (rule.angles.empty()) {
        continue;
      }
      std::sort(rule.angles.begin(), rule.angles.end(),
                [](const Rule::Angle &a, const Rule::Angle &b) {
                  if (a.score != b.score) {
                    return a.score > b.score;
                  }
                  return a.value < b.value;
                });
      rules->push_back(std::move(rule));
    }
  };

  parseFile(torsionLibraryPath);

  // Seperate into ambiguous and non-ambiguous rules
  //  un-ambiguous rules can be run faster via an extraAtomCheck in the query.
  for (std::size_t ri = 0; ri < rules->size(); ++ri) {
    const Rule &r = (*rules)[ri];
    const int eA = r.query->getAtomWithIdx(r.queryAtoms[1])->getAtomicNum();
    const int eB = r.query->getAtomWithIdx(r.queryAtoms[2])->getAtomicNum();
    if (eA == 0 || eB == 0) {
      table->ambiguousRules.push_back(ri);
    } else {
      table->unambiguousRules[centralElementPairKey(eA, eB)].push_back(ri);
    }
  }

  if (readable && rules->empty()) {
    table->status = TorsionSamplerStatus::NoRulesParsed;
  }

  std::shared_ptr<const RuleTable> shared = std::move(table);
  const std::lock_guard<std::mutex> lock(cacheMutex);
  // Another thread may have won the race; keep whichever landed first so all
  // samplers on this path share one table.
  return cache.emplace(torsionLibraryPath, std::move(shared)).first->second;
}

std::vector<double> TorsionLibrarySampler::getAngles(
    const ROMol &mol, unsigned int i, unsigned int j, unsigned int k,
    unsigned int l, bool addBasinAngle) const {
  std::vector<double> fallback = uniformAngles(d_fallbackStepDeg);
  if (!mol.getNumConformers()) {
    return fallback;
  }
  if (!d_table || d_table->rules.empty()) {
    return fallback;
  }

  const Conformer &conf = mol.getConformer();
  const double currentTarget = MolTransforms::getDihedralDeg(conf, i, j, k, l);
  const unsigned int maxAngles =
      d_maxAngles ? d_maxAngles : std::numeric_limits<unsigned int>::max();
  const double minSeparation = 10.0;

  // if mostSpecificMatch: Match among ALL rules whose SMARTS matches.
  // i.e. pick the one with the most query atoms
  const bool firstMatch = !options.mostSpecificMatch;
  const Rule *bestRule = nullptr;
  std::array<unsigned int, 4> bestRuleAtoms{{0, 0, 0, 0}};
  unsigned int bestSpecificity = 0;

  // Prune all rules at once via j-k
  const int eJ = mol.getAtomWithIdx(j)->getAtomicNum();
  const int eK = mol.getAtomWithIdx(k)->getAtomicNum();
  std::vector<std::size_t> candidates;
  const auto &ambiguousRules = d_table->ambiguousRules;
  const auto bit =
      d_table->unambiguousRules.find(centralElementPairKey(eJ, eK));
  if (bit == d_table->unambiguousRules.end()) {
    candidates = ambiguousRules;
  } else if (ambiguousRules.empty()) {
    candidates = bit->second;
  } else {
    candidates.reserve(bit->second.size() + ambiguousRules.size());
    std::merge(bit->second.begin(), bit->second.end(), ambiguousRules.begin(),
               ambiguousRules.end(), std::back_inserter(candidates));
  }

  for (const std::size_t ri : candidates) {
    const Rule &rule = d_table->rules[ri];
    std::array<unsigned int, 4> ruleAtoms{};
    if (!matchCentralBondAnchored(
            mol, *rule.query, rule.queryAtoms, j, k,
            /*collectAll=*/false,
            [&](const std::array<unsigned int, 4> &la) { ruleAtoms = la; })) {
      continue;
    }
    const unsigned int spec = rule.query->getNumAtoms();
    if (spec > bestSpecificity) {
      bestSpecificity = spec;
      bestRule = &rule;
      bestRuleAtoms = ruleAtoms;
    }
    if (bestRule && firstMatch) {
      break;  // legacy: take the first matching rule in file order
    }
  }

  // By default return the matched rule's preferred angles plus a tolerance
  // depending on the basin/well size
  const bool doSym = !options.lean && options.symmetryFolding;
  const bool doOffset = !options.lean && options.terminalOffsets;
  // addBasinAngle=false forces the bare base angles (no tolerance/basin)
  //  this is a coarse sample option.
  const bool doTol = !options.lean && options.toleranceRanges && addBasinAngle;
  static const double kMinTol = 30.0;
  if (bestRule) {
    const double currentAngle =
        MolTransforms::getDihedralDeg(conf, bestRuleAtoms[0], bestRuleAtoms[1],
                                      bestRuleAtoms[2], bestRuleAtoms[3]);
    const double frameOffset = currentTarget - currentAngle;

    // Build the angle set in the torsion library  frame then map to the
    // driver's i-j-k-l frame with frameOffset at the end.
    // (1) base preferred angles (score-sorted when the table was built).
    std::vector<double> angles;
    for (const auto &entry : bestRule->angles) {
      addUniqueAngle(angles, entry.value, minSeparation);
    }
    // (2) +/- tolerance samples for wide wells.  The Hamburg library gives each
    // angle a tight (tolerance1) and a wide (tolerance2) half-width; sample the
    // well edges at both so a bioactive dihedral sitting off the peak but
    // inside the well is reachable.
    if (doTol) {
      for (const auto &entry : bestRule->angles) {
        if (entry.tolerance1 >= kMinTol) {
          addUniqueAngle(angles, entry.value + entry.tolerance1, minSeparation);
          addUniqueAngle(angles, entry.value - entry.tolerance1, minSeparation);
        }
        if (options.wideToleranceRanges && entry.tolerance2 >= kMinTol) {
          addUniqueAngle(angles, entry.value + entry.tolerance2, minSeparation);
          addUniqueAngle(angles, entry.value - entry.tolerance2, minSeparation);
        }
      }
    }
    // (3) ambiguous-terminal-reference offsets: if this rule maps onto the
    // central bond with more than one distinct terminal atom on a side, sample
    // each angle +/- offset about that side's central atom (j for map-1 side, k
    // for map-4 side).
    //   n.b. we can extend symmetry perhaps by canonical ranking ( note: this
    //   is done in below)
    if (doOffset) {
      std::set<unsigned int> ref1, ref4;
      matchCentralBondAnchored(
          mol, *bestRule->query, bestRule->queryAtoms, j, k,
          /*collectAll=*/true, [&](const std::array<unsigned int, 4> &la) {
            ref1.insert(la[0]);
            ref4.insert(la[3]);
          });
      const double off1 = ref1.size() > 1 ? getAngleOffsetDeg(mol, j) : 0.0;
      const double off4 = ref4.size() > 1 ? getAngleOffsetDeg(mol, k) : 0.0;
      if (off1 > 0.0) {
        const std::vector<double> base = angles;
        for (const double a : base) {
          addUniqueAngle(angles, a + off1, minSeparation);
          addUniqueAngle(angles, a - off1, minSeparation);
        }
      }
      if (off4 > 0.0) {
        const std::vector<double> base = angles;
        for (const double a : base) {
          addUniqueAngle(angles, a + off4, minSeparation);
          addUniqueAngle(angles, a - off4, minSeparation);
        }
      }
    }
    // (4) fold rotationally-equivalent angles together for symmetric
    // substituents.
    const unsigned int rotSym =
        options.symmetryUseSmarts ? torsionRotationalSymmetryBySmarts(mol, j, k)
                                  : torsionRotationalSymmetryByRanks(mol, j, k);
    if (doSym && rotSym > 1) {
      angles = foldSymmetricAngles(mol, j, k, std::move(angles), minSeparation,
                                   options.symmetryUseSmarts);
    }

    // DIAGNOSTIC (options.dump): print the matched rule + symmetry +
    // library-frame angles.
    if (options.dump) {
      std::string line = "[torlib] bond " + std::to_string(j) + "-" +
                         std::to_string(k) + " rule=" + bestRule->pattern +
                         " symC" + std::to_string(rotSym) + " angles=";
      std::vector<double> sorted = angles;
      std::sort(sorted.begin(), sorted.end());
      for (double a : sorted) line += std::to_string(a) + " ";
      BOOST_LOG(rdWarningLog) << line << "\n";
      std::fflush(stderr);
    }

    // Map to the driver frame and limit by maxAngles.
    std::vector<double> finalAngles;
    for (const double a : angles) {
      if (angles.size() >= maxAngles) break;
      addUniqueAngle(finalAngles, a + frameOffset, minSeparation);
    }
    if (!finalAngles.empty()) {
      return finalAngles;
    }
  }
  // If we got here, no library rule matched so use the fallback grid
  unsigned int fbSym = 1;
  if (doSym || options.dump) {
    fbSym = options.symmetryUseSmarts
                ? torsionRotationalSymmetryBySmarts(mol, j, k)
                : torsionRotationalSymmetryByRanks(mol, j, k);
  }
  if (doSym && fbSym > 1) {
    fallback = foldSymmetricAngles(mol, j, k, std::move(fallback),
                                   minSeparation, options.symmetryUseSmarts);
  }
  if (options.dump) {
    BOOST_LOG(rdWarningLog)
        << "[torlib] bond " << j << "-" << k << " rule=NONE symC" << fbSym
        << " angles=FALLBACK" << d_fallbackStepDeg << " n=" << fallback.size()
        << "\n";
    std::fflush(stderr);
  }
  return fallback;
}

std::vector<double> CompositeTorsionSampler::getAngles(
    const ROMol &mol, unsigned int i, unsigned int j, unsigned int k,
    unsigned int l, bool addBasinAngle) const {
  std::vector<double> merged;
  for (const auto &sampler : d_samplers) {
    if (!sampler) {
      continue;
    }
    const auto angles = sampler->getAngles(mol, i, j, k, l, addBasinAngle);
    merged.insert(merged.end(), angles.begin(), angles.end());
  }
  // Deduplicate angles within 1 degree of each other.
  std::sort(merged.begin(), merged.end());
  std::vector<double> unique;
  for (const double a : merged) {
    if (unique.empty() || std::abs(a - unique.back()) > 1.0) {
      unique.push_back(a);
    }
  }
  return unique;
}

TorsionSamplerStatus CompositeTorsionSampler::status() const {
  for (const auto &sampler : d_samplers) {
    if (sampler && !sampler->isValid()) {
      return sampler->status();
    }
  }
  return TorsionSamplerStatus::Ok;
}

// ---- SMIRNOFF (OpenFF .offxml) proper-torsion sampler
// -----------------------------------

namespace {
//! One cosine term of a SMIRNOFF proper torsion (idivf omitted: it only scales
//! k, not the minima locations we care about).
struct SmirnoffTerm {
  int periodicity{1};
  double phase{0.0};  // degrees
  double k{0.0};      // kcal/mol
};

//! Energy minima (preferred dihedral angles, degrees) of a SMIRNOFF cosine
//! torsion E(theta) = sum_n k_n (1 + cos(n theta - phase_n)).  Numerically
//! scans the circle and returns the local minima; an (almost) flat profile
//! yields none
std::vector<double> smirnoffTorsionMinima(
    const std::vector<SmirnoffTerm> &terms, double stepDeg,
    double minSeparation) {
  std::vector<double> minima;
  const int n = std::max(2, static_cast<int>(std::lround(360.0 / stepDeg)));
  std::vector<double> e(n);
  double emin = std::numeric_limits<double>::max();
  double emax = -std::numeric_limits<double>::max();
  for (int t = 0; t < n; ++t) {
    const double theta = t * 360.0 / n;
    double en = 0.0;
    for (const auto &tm : terms) {
      en +=
          tm.k *
          (1.0 + std::cos((tm.periodicity * theta - tm.phase) * M_PI / 180.0));
    }
    e[t] = en;
    emin = std::min(emin, en);
    emax = std::max(emax, en);
  }
  if (emax - emin < 1e-6) {
    return minima;  // flat torsion: no preference
  }
  for (int t = 0; t < n; ++t) {
    const double prev = e[(t - 1 + n) % n];
    const double next = e[(t + 1) % n];
    // strict on one side so a plateau contributes a single representative
    if (e[t] <= prev && e[t] < next) {
      addUniqueAngle(minima, t * 360.0 / n, minSeparation);
    }
  }
  return minima;
}
}  // namespace

struct SmirnoffTorsionSampler::Proper {
  std::string smirks;
  std::shared_ptr<ROMol> query;
  std::array<unsigned int, 4> queryAtoms{{0, 0, 0, 0}};
  std::vector<SmirnoffTerm> terms;
};

//! See TorsionLibrarySampler::RuleTable -- same shape for SMIRNOFF propers.
struct SmirnoffTorsionSampler::ProperTable {
  std::vector<Proper> propers;
  std::unordered_map<int, std::vector<std::size_t>> unambiguousRules;
  std::vector<std::size_t> ambiguousRules;
  TorsionSamplerStatus status{TorsionSamplerStatus::Ok};
};

TorsionSamplerStatus SmirnoffTorsionSampler::status() const {
  return d_table ? d_table->status : TorsionSamplerStatus::NoRulesParsed;
}

std::shared_ptr<const SmirnoffTorsionSampler::ProperTable>
SmirnoffTorsionSampler::acquireProperTable(const std::string &offxmlPath) {
  static std::mutex cacheMutex;
  static std::map<std::string, std::shared_ptr<const ProperTable>> cache;
  {
    const std::lock_guard<std::mutex> lock(cacheMutex);
    const auto it = cache.find(offxmlPath);
    if (it != cache.end()) {
      return it->second;
    }
  }

  auto table = std::make_shared<ProperTable>();
  auto *propers = &table->propers;
  std::ifstream in(offxmlPath);
  if (offxmlPath.empty()) {
    table->status = TorsionSamplerStatus::NoPathGiven;
  } else if (!in.good()) {
    table->status = TorsionSamplerStatus::FileUnreadable;
  }
  if (in.good()) {
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string xml = ss.str();
    RDLog::LogStateSetter blocker;

    // Restrict to the <ProperTorsions> ... </ProperTorsions> section.
    size_t sect = xml.find("<ProperTorsions");
    const size_t sectEnd = sect == std::string::npos
                               ? std::string::npos
                               : xml.find("</ProperTorsions>", sect);
    size_t pos =
        sect == std::string::npos ? std::string::npos : xml.find('>', sect);
    while (pos != std::string::npos && pos < sectEnd) {
      pos = xml.find("<Proper", pos);
      if (pos == std::string::npos || pos >= sectEnd) {
        break;
      }
      const char after = pos + 7 < xml.size() ? xml[pos + 7] : '\0';
      if (after != ' ' && after != '\t' && after != '\n' && after != '\r') {
        pos += 7;  // e.g. a stray "<ProperTorsions" -- skip
        continue;
      }
      const size_t tagEnd = xml.find('>', pos);
      if (tagEnd == std::string::npos) {
        break;
      }
      const std::string tag = xml.substr(pos, tagEnd - pos + 1);
      pos = tagEnd + 1;

      const std::string smirks = xmlAttr(tag, "smirks");
      if (smirks.empty()) {
        continue;
      }
      std::unique_ptr<ROMol> query(SmartsToMol(smirks));
      if (!query) {
        continue;
      }
      std::array<unsigned int, 4> queryAtoms{{0, 0, 0, 0}};
      std::array<bool, 4> found{{false, false, false, false}};
      for (const auto atom : query->atoms()) {
        const int mapNum = atom->getAtomMapNum();
        if (mapNum >= 1 && mapNum <= 4) {
          queryAtoms[mapNum - 1] = atom->getIdx();
          found[mapNum - 1] = true;
        }
      }
      if (!found[0] || !found[1] || !found[2] || !found[3]) {
        continue;
      }
      Proper p;
      p.smirks = smirks;
      p.query.reset(query.release());
      p.queryAtoms = queryAtoms;
      for (int term = 1; term <= 6; ++term) {
        const std::string ks = xmlAttr(tag, "k" + std::to_string(term));
        if (ks.empty()) {
          break;
        }
        SmirnoffTerm t;
        t.k = std::atof(ks.c_str());  // "0.11 * kilocalorie..." -> 0.11
        t.periodicity = std::atoi(
            xmlAttr(tag, "periodicity" + std::to_string(term)).c_str());
        t.phase =
            std::atof(xmlAttr(tag, "phase" + std::to_string(term)).c_str());
        if (t.periodicity <= 0) {
          continue;
        }
        p.terms.push_back(t);
      }
      if (p.terms.empty()) {
        continue;
      }
      propers->push_back(std::move(p));
    }
  }
  // Central-element discrimination index; see
  // TorsionLibrarySampler::acquireRuleTable for the ordering contract.
  for (std::size_t ri = 0; ri < propers->size(); ++ri) {
    const Proper &p = (*propers)[ri];
    const int eA = p.query->getAtomWithIdx(p.queryAtoms[1])->getAtomicNum();
    const int eB = p.query->getAtomWithIdx(p.queryAtoms[2])->getAtomicNum();
    if (eA == 0 || eB == 0) {
      table->ambiguousRules.push_back(ri);
    } else {
      table->unambiguousRules[centralElementPairKey(eA, eB)].push_back(ri);
    }
  }

  if (table->status == TorsionSamplerStatus::Ok && propers->empty()) {
    table->status = TorsionSamplerStatus::NoRulesParsed;
  }

  std::shared_ptr<const ProperTable> shared = std::move(table);
  const std::lock_guard<std::mutex> lock(cacheMutex);
  return cache.emplace(offxmlPath, std::move(shared)).first->second;
}

std::vector<double> SmirnoffTorsionSampler::getAngles(
    const ROMol &mol, unsigned int i, unsigned int j, unsigned int k,
    unsigned int l, bool /*addBasinAngle*/) const {
  std::vector<double> fallback = uniformAngles(d_fallbackStepDeg);
  if (!mol.getNumConformers()) {
    return fallback;
  }
  if (!d_table || d_table->propers.empty()) {
    return fallback;
  }
  const Conformer &conf = mol.getConformer();
  const double currentTarget = MolTransforms::getDihedralDeg(conf, i, j, k, l);
  const unsigned int maxAngles =
      d_maxAngles ? d_maxAngles : std::numeric_limits<unsigned int>::max();
  const double minSeparation = 10.0;

  const int eJ = mol.getAtomWithIdx(j)->getAtomicNum();
  const int eK = mol.getAtomWithIdx(k)->getAtomicNum();
  std::vector<std::size_t> candidates;
  const auto &ambiguousRules = d_table->ambiguousRules;
  const auto bit =
      d_table->unambiguousRules.find(centralElementPairKey(eJ, eK));
  if (bit == d_table->unambiguousRules.end()) {
    candidates = ambiguousRules;
  } else if (ambiguousRules.empty()) {
    candidates = bit->second;
  } else {
    candidates.reserve(bit->second.size() + ambiguousRules.size());
    std::merge(bit->second.begin(), bit->second.end(), ambiguousRules.begin(),
               ambiguousRules.end(), std::back_inserter(candidates));
  }

  // SMIRNOFF precedence: the LAST matching SMIRKS wins (candidates are in file
  // order).
  const Proper *best = nullptr;
  std::array<unsigned int, 4> bestAngleAtoms{{0, 0, 0, 0}};
  for (const std::size_t ri : candidates) {
    const Proper &p = d_table->propers[ri];
    std::array<unsigned int, 4> angleAtoms{};
    if (matchCentralBondAnchored(
            mol, *p.query, p.queryAtoms, j, k,
            /*collectAll=*/false,
            [&](const std::array<unsigned int, 4> &la) { angleAtoms = la; })) {
      best = &p;
      bestAngleAtoms = angleAtoms;
    }
  }
  if (!best) {
    return fallback;
  }
  const std::vector<double> minima =
      smirnoffTorsionMinima(best->terms, d_minimaScanStepDeg, minSeparation);
  if (minima.empty()) {
    return fallback;
  }
  // Map library-frame minima into the driver's i-j-k-l frame (the matched
  // SMIRKS may use different terminal reference atoms), exactly as the Hamburg
  // sampler does.
  const double currentAngle =
      MolTransforms::getDihedralDeg(conf, bestAngleAtoms[0], bestAngleAtoms[1],
                                    bestAngleAtoms[2], bestAngleAtoms[3]);
  const double frameOffset = currentTarget - currentAngle;
  std::vector<double> angles;
  for (const double m : minima) {
    if (angles.size() >= maxAngles) {
      break;
    }
    addUniqueAngle(angles, m + frameOffset, minSeparation);
  }
  return angles.empty() ? fallback : angles;
}

}  // namespace RDKit
