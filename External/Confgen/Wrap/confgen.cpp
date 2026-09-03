//
//  Copyright (c) 2026, Glysade Inc.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>
#include <Confgen/SynthonSearch/EnumerateSynthons3D.h>
#include <Confgen/SynthonSearch/SynthonSearch.h>
#include <boost/python/stl_iterator.hpp>
#include <cstdint>
#include <fstream>

namespace python = boost::python;

namespace RDKit {

template <class T>
std::vector<RDKit::MOL_SPTR_VECT> ConvertToVectVect(T bbs) {
  std::vector<RDKit::MOL_SPTR_VECT> vect;
  unsigned int num_bbs = python::len(bbs);
  vect.resize(num_bbs);
  for (unsigned int i = 0; i < num_bbs; ++i) {
    unsigned int len1 = python::len(bbs[i]);
    RDKit::MOL_SPTR_VECT &reacts = vect[i];
    reacts.reserve(len1);
    for (unsigned int j = 0; j < len1; ++j) {
      auto mol = python::extract<RDKit::ROMOL_SPTR>(bbs[i][j]);
      if (mol.check()) {
        reacts.push_back(mol);
      } else {
        throw_value_error("reaction called with non molecule reactant");
      }
    }
  }
  return vect;
}

template <class T>
std::vector<std::vector<std::string>> ConvertToSynthonsVect(T bbs) {
  std::vector<std::vector<std::string>> vect;
  unsigned int num_bbs = python::len(bbs);
  vect.resize(num_bbs);
  for (unsigned int i = 0; i < num_bbs; ++i) {
    unsigned int len1 = python::len(bbs[i]);
    auto &reacts = vect[i];
    reacts.reserve(len1);
    for (unsigned int j = 0; j < len1; ++j) {
      python::extract<std::string> extractor(bbs[i][j]);
      if (extractor.check()) {
        reacts.push_back(extractor());
      } else {
        throw_value_error("conversion called with non string synthon");
      }
    }
  }
  return vect;
}

template <class U, class T>
std::vector<U> ConvertToVect(T bbs) {
  std::vector<U> vect;
  unsigned int len = python::len(bbs);
  vect.reserve(len);
  for (unsigned int i = 0; i < len; ++i) {
    auto v = python::extract<U>(bbs[i]);
    if (v.check()) {
      vect.push_back(v);
    } else {
      throw_value_error("Could not convert python to vector");
    }
  }
  return vect;
}

namespace {
template <typename T>
inline std::vector<T> to_std_vector(const python::object &iterable) {
  return std::vector<T>(python::stl_input_iterator<T>(iterable),
                        python::stl_input_iterator<T>());
}
}  // namespace

void ToBBS(EnumerationStrategyBase &rgroup, ChemicalReaction &rxn,
           python::list ob) {
  rgroup.initialize(rxn, ConvertToVectVect(ob));
}

class EnumerateSynthons3DWrap : public RDKit::EnumerateSynthons3D {
 public:
  ~EnumerateSynthons3DWrap() override {}
  EnumerateSynthons3DWrap() : RDKit::EnumerateSynthons3D() {}
  //! Synthons are joined by their exit-vector LABELS -- there is no reaction
  //! and no Molzip3DParams; conformer generation is configured through
  //! EnumerateSynthons3DParams.
  explicit EnumerateSynthons3DWrap(const python::list ob)
      : RDKit::EnumerateSynthons3D(ConvertToVectVect(ob)) {}
  EnumerateSynthons3DWrap(const python::list ob,
                          const RDKit::EnumerateSynthons3DParams &params)
      : RDKit::EnumerateSynthons3D(ConvertToVectVect(ob), params) {}
  explicit EnumerateSynthons3DWrap(const python::tuple ob)
      : RDKit::EnumerateSynthons3D(ConvertToVectVect(ob)) {}
  EnumerateSynthons3DWrap(const python::tuple ob,
                          const RDKit::EnumerateSynthons3DParams &params)
      : RDKit::EnumerateSynthons3D(ConvertToVectVect(ob), params) {}
};

namespace {

//! reagent indices: python sequence -> vector
std::vector<unsigned int> ToIdx(python::object o) {
  return to_std_vector<unsigned int>(o);
}

SynthonProduct GetProductHelper(const EnumerateSynthons3DWrap &self,
                                python::object idx) {
  return self.getProduct(ToIdx(idx));
}

ROMOL_SPTR ZipProductHelper(const EnumerateSynthons3DWrap &self,
                            python::object idx) {
  return self.zipProduct(ToIdx(idx));
}

//! score() takes a non-const ROMol (alignment mutates conformers), so the
//! caller's molecule is copied rather than modified underneath them.
python::object ScoreHelper(const ShapeProductScorer &self,
                           const ROMol &product) {
  ROMol copy(product);
  const auto v = self.score(copy);
  return v ? python::object(*v) : python::object();
}

python::list ReagentsOf(const SynthonSearchResult &r) {
  python::list l;
  for (auto v : r.reagents) {
    l.append(v);
  }
  return l;
}

python::list BestOf(const SynthonSearchResult &r) {
  python::list l;
  for (const auto &h : r.best) {
    python::list idx;
    for (auto v : h.reagents) {
      idx.append(v);
    }
    l.append(python::make_tuple(idx, h.score, h.mol));
  }
  return l;
}

SynthonSearchResult RefineHelper(const EnumerateSynthons3DWrap &lib,
                                 const SynthonProductScorer &scorer,
                                 python::object seed,
                                 unsigned int maxIters, unsigned int numThreads,
                                 int numBestProducts) {
  return refineSynthons(lib, scorer, ToIdx(seed), maxIters, numThreads,
                        numBestProducts);
}

//! Load a saved library from a file -- the production path (genSynthonLib
//! writes it; a search loads it warm rather than re-embedding).
EnumerateSynthons3DWrap *LoadSynthonLibrary(const std::string &path) {
  std::ifstream in(path, std::ios_base::binary);
  if (!in) {
    throw_value_error("could not open " + path);
  }
  auto *lib = new EnumerateSynthons3DWrap();
  try {
    lib->initFromStream(in);
  } catch (...) {
    delete lib;
    throw;
  }
  return lib;
}

void SaveSynthonLibrary(const EnumerateSynthons3DWrap &self,
                        const std::string &path) {
  std::ofstream out(path, std::ios_base::binary);
  if (!out) {
    throw_value_error("could not write " + path);
  }
  self.toStream(out);
}

SynthonSearchResult ThompsonHelper(const EnumerateSynthons3DWrap &lib,
                                   const SynthonProductScorer &scorer,
                                   const ThompsonSynthonParams &params) {
  return thompsonSynthonSearch(lib, scorer, params);
}

}  // namespace

struct confgen_wrapper {
  static void wrap() {
    std::string docString;

    RegisterVectorConverter<MOL_SPTR_VECT>("VectMolVect");
    RegisterVectorConverter<std::vector<std::string>>("VectSynthonVect");

    docString =
        "EnumerateSynthons3DParams\n\
Controls 3D enumeration of a labelled-synthon library.\n\
Options:\n\
  prefillFraglib [default True]\n\
    Embed this library's synthon fragments up front rather than on first use.\n\
  storeFraglib [default True]\n\
    Persist the fragment cache inside the serialized library.\n\
";

    python::class_<RDKit::EnumerateSynthons3DParams,
                   boost::shared_ptr<RDKit::EnumerateSynthons3DParams>,
                   RDKit::EnumerateSynthons3DParams &>(
        "EnumerateSynthons3DParams", docString.c_str(),
        python::init<>(python::args("self")))
        .def_readwrite("prefillFraglib",
                       &RDKit::EnumerateSynthons3DParams::prefillFraglib)
        .def_readwrite("storeFraglib",
                       &RDKit::EnumerateSynthons3DParams::storeFraglib)
        .def_readwrite("embedStyle",
                       &RDKit::EnumerateSynthons3DParams::embedStyle,
                       "Full cuts at every rotatable bond; Coarse cuts ONLY at "
                       "the synthon junctions (far fewer rotors per product). "
                       "A library must be built for one style or the other.")
        .add_property(
            "numOutputConfs",
            +[](const RDKit::EnumerateSynthons3DParams &p) {
              return p.confgen.numOutputConfs;
            },
            +[](RDKit::EnumerateSynthons3DParams &p, int v) {
              p.confgen.numOutputConfs = v;
            },
            "Conformers generated per assembled product (-1 = AUTO).");

    python::enum_<RDKit::SynthonEmbedStyle>("SynthonEmbedStyle")
        .value("Full", RDKit::SynthonEmbedStyle::Full)
        .value("Coarse", RDKit::SynthonEmbedStyle::Coarse)
        .export_values();

    python::class_<RDKit::SynthonProduct>("SynthonProduct", python::no_init)
        .def_readonly("mol", &RDKit::SynthonProduct::mol)
        .def_readonly("usedCoarseFallback",
                      &RDKit::SynthonProduct::usedCoarseFallback)
        .add_property("ok", +[](const RDKit::SynthonProduct &p) {
          return static_cast<bool>(p);
        });

    python::class_<RDKit::SynthonProductScorer, boost::noncopyable>(
        "SynthonProductScorer",
        "Abstract base: scores an assembled product.", python::no_init);

    python::class_<RDKit::ShapeProductScorer, boost::noncopyable,
                   python::bases<RDKit::SynthonProductScorer>>(
        "ShapeProductScorer",
        "Shape+colour overlay against a fixed query.  The query shape is built "
        "once and reused, so one scorer should serve a whole search.",
        python::init<const RDKit::ROMol &, python::optional<int, bool>>(
            python::args("self", "query", "queryConfId", "allCarbonRadii")))
        .def("Score", &RDKit::ScoreHelper,
             "Best-aligning conformer wins.  Returns None if the product could "
             "not be scored.",
             python::args("self", "product"))
        .def("LastShape", &RDKit::ShapeProductScorer::lastShape,
             python::args("self"))
        .def("LastColour", &RDKit::ShapeProductScorer::lastColour,
             python::args("self"));

    python::class_<RDKit::SynthonSearchResult>("SynthonSearchResult",
                                               python::no_init)
        .add_property("reagents", &RDKit::ReagentsOf)
        .def_readonly("score", &RDKit::SynthonSearchResult::score)
        .def_readonly("evaluations", &RDKit::SynthonSearchResult::evaluations)
        .def_readonly("unscorable", &RDKit::SynthonSearchResult::unscorable)
        .def_readonly("cacheHits", &RDKit::SynthonSearchResult::cacheHits)
        .def_readonly("confgenMs", &RDKit::SynthonSearchResult::confgenMs)
        .def_readonly("scoreMs", &RDKit::SynthonSearchResult::scoreMs)
        .add_property("best", &RDKit::BestOf,
                      "[(reagentIdx, score, mol), ...] highest score first.");

    python::class_<RDKit::ThompsonSynthonParams>(
        "ThompsonSynthonParams", python::init<>(python::args("self")))
        .def_readwrite("budget", &RDKit::ThompsonSynthonParams::budget)
        .def_readwrite("priorMean", &RDKit::ThompsonSynthonParams::priorMean)
        .def_readwrite("priorVar", &RDKit::ThompsonSynthonParams::priorVar)
        .def_readwrite("noiseVar", &RDKit::ThompsonSynthonParams::noiseVar)
        .def_readwrite("batchSize", &RDKit::ThompsonSynthonParams::batchSize)
        .def_readwrite("randomSeed", &RDKit::ThompsonSynthonParams::randomSeed)
        .def_readwrite("numThreads", &RDKit::ThompsonSynthonParams::numThreads)
        .def_readwrite("cacheScores", &RDKit::ThompsonSynthonParams::cacheScores)
        .def_readwrite("numBestProducts",
                       &RDKit::ThompsonSynthonParams::numBestProducts);

    python::def("ThompsonSynthonSearch", &RDKit::ThompsonHelper,
                (python::arg("lib"), python::arg("scorer"),
                 python::arg("params") = RDKit::ThompsonSynthonParams()),
                "Thompson sampling over (position, reagent) arms.");

    python::def("RefineSynthons", &RDKit::RefineHelper,
                (python::arg("lib"), python::arg("scorer"), python::arg("seed"),
                 python::arg("maxIters") = 5, python::arg("numThreads") = 0,
                 python::arg("numBestProducts") = 10),
                "Coordinate descent from a seed combination.");

    docString = "";
    python::class_<EnumerateSynthons3DWrap, boost::noncopyable,
                   python::bases<RDKit::EnumerateLibraryBase>>(
        "EnumerateSynthons3D", docString.c_str(),
        python::init<>(python::args("self")))
        .def(python::init<python::list>(python::args("self", "reagents")))
        .def(python::init<python::list,
                          const RDKit::EnumerateSynthons3DParams &>(
            python::args("self", "reagents", "params")))
        .def(python::init<python::tuple>(python::args("self", "reagents")))
        .def(python::init<python::tuple,
                          const RDKit::EnumerateSynthons3DParams &>(
            python::args("self", "reagents", "params")))

        .def("Prefill", &RDKit::EnumerateSynthons3D::prefill,
             "Embed this library's synthon fragments into its cache. Returns "
             "the number newly embedded (0 if already warm).",
             python::args("self"))

        .def(
            "GetReagents", &RDKit::EnumerateSynthons3D::getReagents,
            "Return the reagents used in this library.",
            python::return_internal_reference<
                1, python::with_custodian_and_ward_postcall<0, 1>>(),
            python::args("self"))

        .def("Arity", &RDKit::EnumerateSynthons3DWrap::arity, python::args("self"))
        .def("NumReagents", &RDKit::EnumerateSynthons3DWrap::numReagents,
             python::args("self", "position"))
        .def("IsValid", &RDKit::EnumerateSynthons3DWrap::isValid,
             python::args("self"))
        .def("GetProduct", &RDKit::GetProductHelper,
             "Assemble one product WITH conformers.",
             python::args("self", "reagentIdx"))
        .def("ZipProduct", &RDKit::ZipProductHelper,
             "The product GRAPH only -- no conformers.  Use this to build a "
             "query independently, so it does not inherit the candidates' "
             "conformers.",
             python::args("self", "reagentIdx"))
        .def("Save", &RDKit::SaveSynthonLibrary,
             "Write this library, fragment cache included, to a file.",
             python::args("self", "path"));

    python::def("LoadSynthonLibrary", &RDKit::LoadSynthonLibrary,
                python::return_value_policy<python::manage_new_object>(),
                "Load a library written by Save() or genSynthonLib.",
                python::arg("path"));
  }
};

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdFragmentConfGen) {
  python::scope().attr("__doc__") = "Conformation Generation Playground";
  RDKit::confgen_wrapper::wrap();
}
