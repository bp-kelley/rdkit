//  Copyright (C) 2026 Glysade Inc and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SynthonSearch/EnumerateSynthons3D.h"

#include <set>
#include <thread>
#include <atomic>
#include <istream>
#include <ostream>

#include "SynthonZipUtils.h"
#include "Utils/ParamsIO.h"

#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/RDLog.h>

#include <limits>

namespace RDKit {

namespace {
//! Heavy atoms excluding exit dummies (atomic number 0, or the U/Np/Pu
//! placeholders the label schemes use before conversion).
unsigned int countHeavyNoDummies(const ROMol &m) {
  unsigned int n = 0;
  for (const auto atom : m.atoms()) {
    const int z = atom->getAtomicNum();
    if (z <= 1 || z == 92 || z == 93 || z == 94) {
      continue;  // H, dummy, or an exit-vector placeholder
    }
    ++n;
  }
  return n;
}
}  // namespace


const char *synthonBuildStatusMessage(SynthonBuildStatus s) {
  switch (s) {
    case SynthonBuildStatus::Ok:
      return "ok";
    case SynthonBuildStatus::BadReagentIndex:
      return "reagent index out of range for that position";
    case SynthonBuildStatus::ZipFailed:
      return "the labelled synthons could not be joined";
    case SynthonBuildStatus::NoReagentConfs:
      return "a chosen synthon has no conformers";
    case SynthonBuildStatus::ConfGenFailed:
      return "conformer generation produced nothing for the product";
  }
  return "unknown";
}

EnumerateSynthons3D::EnumerateSynthons3D(
    const EnumerationTypes::BBS &reagents, EnumerateSynthons3DParams params)
    // Delegate so a log blocker is alive while the BASE is constructed; see
    // LogBlocker.
    : EnumerateSynthons3D(LogBlocker{}, reagents, std::move(params)) {}

EnumerateSynthons3D::EnumerateSynthons3D(
    LogBlocker, const EnumerationTypes::BBS &reagents,
    EnumerateSynthons3DParams params)
    // No reaction: synthons carry exit-vector labels and are joined by molzip,
    // so the base gets an empty ChemicalReaction -- same as EnumerateSynthons.
    : EnumerateLibrary(ChemicalReaction(), reagents),
      d_params(std::move(params)) {
  d_valid = !m_bbs.empty() && SynthonZip::sniffMolzipParams(d_molzipParams, m_bbs[0]);
  if (!d_valid) {
    BOOST_LOG(rdWarningLog)
        << "EnumerateSynthons3D: synthons do not share one exit-vector label "
           "scheme; the library cannot be zipped\n";
    return;
  }
  // This library's own cache, shared by every generator working on it.
  FraglibParams flp = d_params.confgen.embedding;
  // FragmentConfGenParams::randomSeed is the public seed for the whole
  // pipeline; the nested FraglibParams seed is an implementation detail.
  flp.randomSeed = d_params.confgen.randomSeed;
  flp.ffVariant = d_params.confgen.joiner.ffVariant;
  // n.b. share the fraglib across searches for better
  //  optimization
  if (d_params.confgen.fraglib) {
    d_fraglib = d_params.confgen.fraglib;
  } else {
    d_fraglib = std::make_shared<Fraglib>(flp);
    d_params.confgen.fraglib = d_fraglib;
  }

  cacheSynthonSizes();

  if (d_params.prefillFraglib) {
    prefill();
  }
}



void EnumerateSynthons3D::cacheSynthonSizes() {
  d_fragmentHeavyCount.clear();
  for (const auto &position : getReagents()) {
    d_fragmentHeavyCount.emplace_back();
    auto &row = d_fragmentHeavyCount.back();
    row.reserve(position.size());
    for (const auto &syn : position) {
      row.push_back(syn ? countHeavyNoDummies(*syn) : 0u);
    }
  }
}

const std::vector<unsigned int> &EnumerateSynthons3D::synthonCutBonds(
    unsigned int pos, unsigned int idx) const {
  static const std::vector<unsigned int> empty;
  if (pos >= d_synthonCutBonds.size() ||
      idx >= d_synthonCutBonds[pos].size()) {
    return empty;
  }
  return d_synthonCutBonds[pos][idx];
}

void EnumerateSynthons3D::setSynthonCutBonds(unsigned int pos,
                                             unsigned int idx,
                                             std::vector<unsigned int> bonds) {
  if (d_synthonCutBonds.size() != getReagents().size()) {
    d_synthonCutBonds.resize(getReagents().size());
    for (size_t p = 0; p < getReagents().size(); ++p) {
      d_synthonCutBonds[p].resize(getReagents()[p].size());
    }
  }
  if (pos >= d_synthonCutBonds.size() ||
      idx >= d_synthonCutBonds[pos].size()) {
    return;
  }
  if (!bonds.empty()) {
    d_haveCutBonds = true;
  }
  d_synthonCutBonds[pos][idx] = std::move(bonds);
}

unsigned int EnumerateSynthons3D::fragmentHeavyCount(unsigned int pos,
                                                    unsigned int idx) const {
  if (pos >= d_fragmentHeavyCount.size() ||
      idx >= d_fragmentHeavyCount[pos].size()) {
    return 0;
  }
  return d_fragmentHeavyCount[pos][idx];
}

unsigned int EnumerateSynthons3D::productHeavyCount(
    const std::vector<unsigned int> &reagentIdx) const {
  unsigned int n = 0;
  for (unsigned int p = 0; p < reagentIdx.size(); ++p) {
    n += fragmentHeavyCount(p, reagentIdx[p]);
  }
  return n;
}

std::pair<unsigned int, unsigned int> EnumerateSynthons3D::positionSizeRange(
    unsigned int pos) const {
  if (pos >= d_fragmentHeavyCount.size()) {
    return {0, 0};
  }
  unsigned int mn = std::numeric_limits<unsigned int>::max(), mx = 0;
  for (unsigned int v : d_fragmentHeavyCount[pos]) {
    if (!v) {
      continue;  // a null synthon contributes nothing
    }
    mn = std::min(mn, v);
    mx = std::max(mx, v);
  }
  if (mn == std::numeric_limits<unsigned int>::max()) {
    mn = 0;
  }
  return {mn, mx};
}

std::pair<unsigned int, unsigned int> EnumerateSynthons3D::productSizeRange()
    const {
  unsigned int lo = 0, hi = 0;
  for (const auto &row : d_fragmentHeavyCount) {
    unsigned int mn = std::numeric_limits<unsigned int>::max(), mx = 0;
    for (unsigned int v : row) {
      if (!v) {
        continue;  // a null synthon contributes nothing
      }
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
    if (mn == std::numeric_limits<unsigned int>::max()) {
      mn = 0;
    }
    lo += mn;
    hi += mx;
  }
  return {lo, hi};
}

namespace {
//! XXX FIX ME -> this is MMFF specific and will need to be
//! Can the force field assign parameters to every atom?
//!   updated/moved in the future to the forcefield bits.
bool ffCanType(const ROMol &mol, const std::string &variant) {
  try {
    RWMol probe(mol);
    // Replace isotopes with carbons for typing
    for (auto atom : probe.atoms()) {
      if (!atom->getAtomicNum()) {
        atom->setAtomicNum(6);
        atom->setIsotope(0);
        atom->setNoImplicit(false);
        atom->setNumExplicitHs(0);
      }
    }
    RDLog::LogStateSetter blocker;  // typing failures are expected here
    MolOps::sanitizeMol(probe);
    MolOps::addHs(probe);
    MMFF::MMFFMolProperties props(probe, variant);
    return props.isValid();
  } catch (...) {
    return false;
  }
}
}  // namespace

unsigned int EnumerateSynthons3D::prefill(unsigned int numThreads) {
  if (!d_fraglib) {
    return 0;
  }
  const size_t before = d_fraglib->size();
  // Embed the final fragments of the synthon. A synthon is not
  // itself a fragment: FragmentConfGen cuts at rotatable bonds, so a synthon
  // with internal rotors may become several during embedding.

  // However: Coarse mode tries to embed as much of the synthon as it
  //  can to minimize the number of rotatable bonds during the search.
  //  The default is to make the cut bond the zip-junction bond and
  //  embed the rest using a fine embedding (low RMSD, high # attempts)
  //
  // This falls down when rings are formed by the reaction.
  //  In this case, the synthon library generator detects this
  //  and sets appropriate cut bonds that don't cut rings.

  const bool coarse = d_params.embedStyle == SynthonEmbedStyle::Coarse;
  const std::set<std::string> exitSymbols(d_molzipParams.atomSymbols.begin(),
                                          d_molzipParams.atomSymbols.end());
  const std::vector<unsigned int> noCuts;

  // Make the work list for threading
  std::vector<const ROMol *> work;
  for (const auto &position : getReagents()) {
    for (const auto &synthon : position) {
      if (synthon) {
        work.push_back(synthon.get());
      }
    }
  }
  unsigned int nThreads = numThreads ? numThreads
                                     : std::thread::hardware_concurrency();
  if (!nThreads) {
    nThreads = 1;
  }
  nThreads = std::min<unsigned int>(nThreads, static_cast<unsigned int>(
                                                  std::max<size_t>(1, work.size())));

  // atomic counter so we grab synthons in order when threading.
  std::atomic<size_t> next{0};
  auto embedSynthons = [&]() {
    FragmentConfGen gen(d_params.confgen);
    for (;;) {
      const size_t i = next++;
      if (i >= work.size()) {
        return;
      }
      const ROMol *synthon = work[i];
      const auto id = cacheFragment(*synthon);
      if (!id) {
        continue;  // will not sanitize: no key exists to record it under
      }
      try {
        if (coarse) {
          // See if we can actually embed.
          if (!ffCanType(*id, d_params.confgen.joiner.ffVariant)) {
            d_fraglib->markUnembeddable(*id);
            continue;
          }
          gen.fragmentAndEmbed(*id, nullptr, &noCuts);
        } else {
          gen.fragmentAndEmbed(*id);
        }
      } catch (...) {
        // We can't embed this one so mark it dead
        d_fraglib->markUnembeddable(*id);
      }
    }
  };

  if (nThreads <= 1) {
    embedSynthons();
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nThreads);
    for (unsigned int i = 0; i < nThreads; ++i) {
      pool.emplace_back(embedSynthons);
    }
    for (auto &thread : pool) {
      thread.join();
    }
  }
  return static_cast<unsigned int>(d_fraglib->size() - before);
}

std::unique_ptr<RWMol> EnumerateSynthons3D::cacheFragment(
    const ROMol &synthon) const {
  try {
    auto out = std::make_unique<RWMol>(synthon);
    if (d_params.embedStyle == SynthonEmbedStyle::Coarse) {
      const std::set<std::string> exitSymbols(d_molzipParams.atomSymbols.begin(),
                                              d_molzipParams.atomSymbols.end());
      for (auto atom : out->atoms()) {
        if (exitSymbols.count(atom->getSymbol())) {
          atom->setAtomicNum(0);
          atom->setIsotope(0);
          atom->setNoImplicit(true);
          atom->setNumExplicitHs(0);
        }
      }
    }
    RDLog::LogStateSetter blocker;
    MolOps::sanitizeMol(*out);
    return out;
  } catch (...) {
    return nullptr;  // cannot be keyed, so cannot be recorded either way
  }
}

bool EnumerateSynthons3D::synthonUnusable(unsigned int pos,
                                          unsigned int idx) const {
  if (!d_fraglib) {
    return false;  // nothing recorded: let the product build find out
  }
  const auto &bbs = getReagents();
  if (pos >= bbs.size() || idx >= bbs[pos].size() || !bbs[pos][idx]) {
    return false;
  }
  std::shared_ptr<UnusableCache> cache;
  {
    static std::mutex initMutex;
    std::lock_guard<std::mutex> lock(initMutex);
    if (!d_unusable) {
      d_unusable = std::make_shared<UnusableCache>();
      d_unusable->flags.resize(bbs.size());
      for (size_t p = 0; p < bbs.size(); ++p) {
        d_unusable->flags[p].assign(bbs[p].size(), -1);
      }
    }
    cache = d_unusable;
  }
  {
    std::lock_guard<std::mutex> lock(cache->mutex);
    const signed char cached = cache->flags[pos][idx];
    if (cached >= 0) {
      return cached != 0;
    }
  }

  bool dead = true;
  if (const auto frag = cacheFragment(*bbs[pos][idx])) {
    const auto n = d_fraglib->numFragmentConfs(*frag);
    dead = n.has_value() && *n == 0;
  }
  std::lock_guard<std::mutex> lock(cache->mutex);
  cache->flags[pos][idx] = dead ? 1 : 0;
  return dead;
}

ROMOL_SPTR EnumerateSynthons3D::get2D(
    const std::vector<unsigned int> &reagentIdx) const {
  const auto &bbs = getReagents();
  if (!d_valid || bbs.empty() || reagentIdx.size() != bbs.size()) {
    return {};
  }
  for (size_t i = 0; i < bbs.size(); ++i) {
    if (reagentIdx[i] >= bbs[i].size() || !bbs[i][reagentIdx[i]]) {
      return {};
    }
  }
  RDLog::LogStateSetter blocker;
  try {
    RWMol m(*bbs[0][reagentIdx[0]]);
    for (size_t i = 1; i < bbs.size(); ++i) {
      m.insertMol(*bbs[i][reagentIdx[i]]);
    }
    std::unique_ptr<ROMol> zipped(molzip(m, d_molzipParams));
    if (!zipped || !zipped->getNumAtoms()) {
      return {};
    }
    return ROMOL_SPTR(new ROMol(*zipped));
  } catch (...) {
    return {};
  }
}

SynthonProduct EnumerateSynthons3D::getProduct(
    const std::vector<unsigned int> &reagentIdx) const {
  SynthonProduct out;
  const auto &bbs = getReagents();
  if (!d_valid || bbs.empty() || reagentIdx.size() != bbs.size()) {
    out.status = SynthonBuildStatus::BadReagentIndex;
    return out;
  }
  for (size_t i = 0; i < bbs.size(); ++i) {
    if (reagentIdx[i] >= bbs[i].size() || !bbs[i][reagentIdx[i]]) {
      out.status = SynthonBuildStatus::BadReagentIndex;
      return out;
    }
  }

  // Check if the synthon is usable or not first
  for (size_t i = 0; i < reagentIdx.size(); ++i) {
    if (synthonUnusable(static_cast<unsigned int>(i), reagentIdx[i])) {
      out.status = SynthonBuildStatus::NoReagentConfs;
      return out;
    }
  }

  // Assembly can fail in some cases, so block the logs
  RDLog::LogStateSetter blocker;

  //  Note: Coarse needs to set it's own junction bonds, we can't use the defaults.
  const bool coarse = d_params.embedStyle == SynthonEmbedStyle::Coarse;
  static const std::string kSynthonPos = "_synthonPos";
  static const std::string kSynthonAtom = "_synthonAtom";
  std::unique_ptr<ROMol> zipped;
  try {
    RWMol m(*bbs[0][reagentIdx[0]]);
    if (coarse) {
      for (auto atom : m.atoms()) {
        atom->setProp<int>(kSynthonPos, 0);
        // Mark the original atoms as they will be renumbered
        atom->setProp<int>(kSynthonAtom, static_cast<int>(atom->getIdx()));
      }
    }
    for (size_t i = 1; i < bbs.size(); ++i) {
      RWMol part(*bbs[i][reagentIdx[i]]);
      if (coarse) {
        for (auto atom : part.atoms()) {
          atom->setProp<int>(kSynthonPos, static_cast<int>(i));
          atom->setProp<int>(kSynthonAtom, static_cast<int>(atom->getIdx()));
        }
      }
      m.insertMol(part);
    }
    zipped = molzip(m, d_molzipParams);
  } catch (...) {
    zipped.reset();
  }
  
  if (!zipped || !zipped->getNumAtoms()) {
    out.status = SynthonBuildStatus::ZipFailed;
    return out;
  }

  // build the conf
  std::vector<unsigned int> junctionBonds;
  if (coarse) {
    for (const auto bond : zipped->bonds()) {
      int a = -1, b = -1;
      if (bond->getBeginAtom()->getPropIfPresent(kSynthonPos, a) &&
          bond->getEndAtom()->getPropIfPresent(kSynthonPos, b) && a != b) {
        junctionBonds.push_back(bond->getIdx());
      }
    }

    // if a synthon has precomputed cut bonds honor those.
    if (d_haveCutBonds) {
      std::map<std::pair<int, int>, unsigned int> fromSynthon;
      for (const auto atom : zipped->atoms()) {
        int pos = -1, idx = -1;
        if (atom->getPropIfPresent(kSynthonPos, pos) &&
            atom->getPropIfPresent(kSynthonAtom, idx)) {
          fromSynthon[{pos, idx}] = atom->getIdx();
        }
      }
      std::set<unsigned int> cuts;
      bool nonStandardCuts = false;
      for (unsigned int p = 0; p < reagentIdx.size(); ++p) {
        const auto &want = synthonCutBonds(p, reagentIdx[p]);
        if (want.empty()) {
          continue;  // this synthon is cut at its junctions, as usual
        }
        if (p >= bbs.size() || reagentIdx[p] >= bbs[p].size() ||
            !bbs[p][reagentIdx[p]]) {
          continue;  // a null synthon: nothing to translate
        }
        nonStandardCuts = true;
        const ROMol &syn = *bbs[p][reagentIdx[p]];
        for (unsigned int b : want) {
          if (b >= syn.getNumBonds()) {
            continue;
          }
          const Bond *sb = syn.getBondWithIdx(b);
          const auto a1 = fromSynthon.find(
              {static_cast<int>(p), static_cast<int>(sb->getBeginAtomIdx())});
          const auto a2 = fromSynthon.find(
              {static_cast<int>(p), static_cast<int>(sb->getEndAtomIdx())});
          if (a1 == fromSynthon.end() || a2 == fromSynthon.end()) {
            continue;  // an exit atom consumed by the zip; nothing to cut
          }
          if (const Bond *pb =
                  zipped->getBondBetweenAtoms(a1->second, a2->second)) {
            cuts.insert(pb->getIdx());
          }
        }
      }
      
      if (nonStandardCuts) {
        // keep the junctions of any synthon that has declard cuts
        for (unsigned int b : junctionBonds) {
          const Bond *bond = zipped->getBondWithIdx(b);
          int p1 = -1, p2 = -1;
          bond->getBeginAtom()->getPropIfPresent(kSynthonPos, p1);
          bond->getEndAtom()->getPropIfPresent(kSynthonPos, p2);
          const bool lhs =
              p1 >= 0 && static_cast<size_t>(p1) < reagentIdx.size() &&
              synthonCutBonds(p1, reagentIdx[p1]).empty();
          const bool rhs =
              p2 >= 0 && static_cast<size_t>(p2) < reagentIdx.size() &&
              synthonCutBonds(p2, reagentIdx[p2]).empty();
          if (lhs && rhs) {
            cuts.insert(b);
          }
        }
        junctionBonds.assign(cuts.begin(), cuts.end());
      }
    }
  }

  FragmentConfGen gen(d_params.confgen);
  FragmentConfGenResult res =
      coarse ? gen.build(*zipped, &junctionBonds) : gen.build(*zipped);

  // If we fail the coarse embedding, maybe for steric issues,
  //  fall back to a full confgen.
  if (coarse && res.conformers.empty() &&
      res.resultType() != FragConfGenResultType::TIMED_OUT) {
    res = gen.build(*zipped);
    if (!res.conformers.empty()) {
      out.usedCoarseFallback = true;
    }
  }

  if (res.conformers.empty()) {
    // Set the output status
    out.status = res.resultType() == FragConfGenResultType::FF_FAIL
                     ? SynthonBuildStatus::NoReagentConfs
                     : SynthonBuildStatus::ConfGenFailed;
    return out;
  }

  auto product = boost::make_shared<RWMol>(*res.conformers.front());
  for (size_t i = 1; i < res.conformers.size(); ++i) {
    if (res.conformers[i]->getNumConformers()) {
      product->addConformer(new Conformer(res.conformers[i]->getConformer()),
                            true);
    }
  }
  out.mol = product;
  out.status = SynthonBuildStatus::Ok;
  return out;
}

void EnumerateSynthons3D::toStream(std::ostream &ss) const {
  //  XXX FIX ME, we should make a proper serializer not this
  //   text + binary nonsense
  EnumerateLibrary::toStream(ss);
  ss.put('\n');
  // Serialize the parameters
  {
    const std::string cfg = fragmentConfGenParamsToString(d_params.confgen);
    const auto n = static_cast<std::uint64_t>(cfg.size());
    ss.write(reinterpret_cast<const char *>(&n), sizeof(n));
    ss.write(cfg.data(), static_cast<std::streamsize>(n));
    const auto style = static_cast<std::uint8_t>(d_params.embedStyle);
    ss.write(reinterpret_cast<const char *>(&style), sizeof(style));
    const bool prefilled = d_params.prefillFraglib;
    ss.write(reinterpret_cast<const char *>(&prefilled), sizeof(prefilled));
  }
  // serialize the extra cut bonds
  {
    const bool haveCuts = d_haveCutBonds;
    ss.write(reinterpret_cast<const char *>(&haveCuts), sizeof(haveCuts));
    if (haveCuts) {
      const auto nPos = static_cast<std::uint64_t>(d_synthonCutBonds.size());
      ss.write(reinterpret_cast<const char *>(&nPos), sizeof(nPos));
      for (const auto &row : d_synthonCutBonds) {
        const auto nSyn = static_cast<std::uint64_t>(row.size());
        ss.write(reinterpret_cast<const char *>(&nSyn), sizeof(nSyn));
        for (const auto &bonds : row) {
          const auto nB = static_cast<std::uint64_t>(bonds.size());
          ss.write(reinterpret_cast<const char *>(&nB), sizeof(nB));
          if (nB) {
            ss.write(reinterpret_cast<const char *>(bonds.data()),
                     static_cast<std::streamsize>(nB * sizeof(unsigned int)));
          }
        }
      }
    }
  }
  const bool haveLib = d_params.storeFraglib && d_fraglib;
  ss.write(reinterpret_cast<const char *>(&haveLib), sizeof(haveLib));
  if (haveLib) {
    d_fraglib->serialize(ss);
  }
}

void EnumerateSynthons3D::initFromStream(std::istream &ss) {
  EnumerateLibrary::initFromStream(ss);
  d_valid = !m_bbs.empty() &&
            SynthonZip::sniffMolzipParams(d_molzipParams, m_bbs[0]);
  ss >> std::ws;  // see toStream: text archive -> raw binary boundary
  {
    std::uint64_t n = 0;
    ss.read(reinterpret_cast<char *>(&n), sizeof(n));
    std::string cfg(static_cast<size_t>(n), '\0');
    if (n) {
      ss.read(&cfg[0], static_cast<std::streamsize>(n));
    }
    const std::string err =
        fragmentConfGenParamsFromString(cfg, d_params.confgen);
    if (!err.empty()) {
      BOOST_LOG(rdWarningLog)
          << "EnumerateSynthons3D: stored parameters did not parse (" << err
          << "); falling back to defaults\n";
    }
    std::uint8_t style = 0;
    ss.read(reinterpret_cast<char *>(&style), sizeof(style));
    d_params.embedStyle = static_cast<SynthonEmbedStyle>(style);
    bool prefilled = false;
    ss.read(reinterpret_cast<char *>(&prefilled), sizeof(prefilled));
    d_params.prefillFraglib = prefilled;
  }
  // Make the fraglib now that we have the settings
  FraglibParams flp = d_params.confgen.embedding;
  flp.randomSeed = d_params.confgen.randomSeed;
  flp.ffVariant = d_params.confgen.joiner.ffVariant;
  if (d_params.confgen.fraglib) {
    d_fraglib = d_params.confgen.fraglib;
  } else {
    d_fraglib = std::make_shared<Fraglib>(flp);
  }
  d_synthonCutBonds.clear();
  d_haveCutBonds = false;
  {
    bool haveCuts = false;
    ss.read(reinterpret_cast<char *>(&haveCuts), sizeof(haveCuts));
    if (haveCuts) {
      std::uint64_t nPos = 0;
      ss.read(reinterpret_cast<char *>(&nPos), sizeof(nPos));
      d_synthonCutBonds.resize(nPos);
      for (auto &row : d_synthonCutBonds) {
        std::uint64_t nSyn = 0;
        ss.read(reinterpret_cast<char *>(&nSyn), sizeof(nSyn));
        row.resize(nSyn);
        for (auto &bonds : row) {
          std::uint64_t nB = 0;
          ss.read(reinterpret_cast<char *>(&nB), sizeof(nB));
          bonds.resize(nB);
          if (nB) {
            ss.read(reinterpret_cast<char *>(bonds.data()),
                    static_cast<std::streamsize>(nB * sizeof(unsigned int)));
          }
        }
      }
      d_haveCutBonds = true;
    }
  }
  bool haveLib = false;
  ss.read(reinterpret_cast<char *>(&haveLib), sizeof(haveLib));
  if (haveLib && ss.good()) {
    d_fraglib->initFromStream(ss);
  }
  d_params.confgen.fraglib = d_fraglib;
  cacheSynthonSizes();
}

std::vector<MOL_SPTR_VECT> EnumerateSynthons3D::get(
    const EnumerationTypes::RGROUPS &pos) const {
  std::vector<unsigned int> idx(pos.begin(), pos.end());
  const SynthonProduct p = getProduct(idx);
  if (!p) {
    return {}; // nada
  }
  std::vector<MOL_SPTR_VECT> res(1);
  res[0].push_back(p.mol);
  return res;
}

}  // namespace RDKit
