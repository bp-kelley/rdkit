
rdkit_library(FragmentConfGen
		    FragmentConfGen.cpp
		    Embedder/Fraglib.cpp
		    Search/RotorDriver.cpp
		    Search/InterFragScore.cpp
		    Search/TreeSearch.cpp
		    Search/RigidRotorSearch.cpp
		    Joiner/FragmentJoiner.cpp
		    Search/ThompsonSamplingSearch.cpp
		    SynthonSearch/EnumerateSynthons3D.cpp
		    SynthonSearch/SynthonSearch3D.cpp
		    Search/RotorRefine.cpp
		    Search/SystematicSearch.cpp
		    Search/RigidRotorSearchExperimental.cpp
		    Utils/RotorTopology.cpp
		    Utils/TheobaldRmsd.cpp
		    Utils/SymmetricRmsd.cpp
		    Utils/ParamsIO.cpp
		    Sampler/TorsionSampler.cpp
		    Sampler/TorsionSymmetries.cpp
              LINK_LIBRARIES
              ChemReactions DistGeomHelpers ForceFieldHelpers MolTransforms ForceField
              GaussianShape
              DataStructs Descriptors
              FileParsers SubstructMatch ChemTransforms GraphMol)
target_compile_definitions(FragmentConfGen PRIVATE RDKIT_FRAGMENT_CONFGEN_BUILD)

target_include_directories(FragmentConfGen PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

if(TARGET FragmentConfGen_static)
  target_include_directories(FragmentConfGen_static PUBLIC
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
endif()

rdkit_headers(	      FragmentConfGen.h
	      Embedder/Fraglib.h
	      Sampler/TorsionSampler.h
	      Sampler/TorsionSymmetries.h
	      SynthonSearch/EnumerateSynthons3D.h
	      SynthonSearch/SynthonSearch3D.h
	      SynthonSearch/SynthonZipUtils.h
	      Utils/RotorTopology.h
	      Utils/NonbondedLookup.h
	      Utils/TheobaldRmsd.h
	      Utils/SymmetricRmsd.h
	      Utils/DiagnosticsParams.h
	      Utils/ParamsIO.h
	      Search/SearchParams.h
	      Search/RigidRotorSearch.h
              DEST GraphMol)

rdkit_catch_test(synthon3DTestCatch test/testEnumerateSynthons3D.cpp
LINK_LIBRARIES FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign )

rdkit_catch_test(fragConfGenTestCatch test/testFragmentConfGen.cpp
LINK_LIBRARIES FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors )

# Standalone timing / measurement harness for the fragment conformer generator.
find_package(Threads REQUIRED)

# Built but NOT registered as a test (it is slow and needs benchmark data);
# run a profile by name, e.g. `fragConfGenBench breakdown`.  The `platinum`
# profile parallelises across molecules, hence Threads::Threads.
if(RDK_BUILD_CPP_TESTS)
  add_executable(fragConfGenBench Tools/fragConfGenBench.cpp)
  target_link_libraries(fragConfGenBench
    FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
    FileParsers Threads::Threads)
endif()

add_executable(genFragLib Tools/genFragLib.cpp)
target_link_libraries(genFragLib
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  FileParsers Threads::Threads)

# rotorProfile: DIAGNOSTIC.  For the widest-moving rotors of one molecule, print
# the MMFF profile over a full turn and a clash census bucketed by topological
# distance from the rotor bond (local strain vs global).
add_executable(rotorProfile Tools/rotorProfile.cpp)
target_link_libraries(rotorProfile
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  FileParsers Threads::Threads)

# genSynthonLib: build a 3D synthon library ONCE (with its fragment cache warm)
# and save it.  Production loads a prepared library rather than embedding at
# search time; doing the same here keeps prefill out of the measurement.
add_executable(genSynthonLib Tools/genSynthonLib.cpp)
target_link_libraries(genSynthonLib
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  GaussianShape FileParsers Threads::Threads)

# synthonSearchBench: RECOVERY BENCHMARK for the 3D synthon search.  Plants a
# random product from a synthon library, uses its conformers as the shape query,
# and reports whether the selected search finds it again -- with logical and
# physical work, confgen time, and shape time split out. Ships without data:
# point it at a synthon CSV or a prepared 3D library.
# synthonReachability: given (planted, found) pairs from a benchmark run, ask
# whether the missed product was REACHABLE -- walk to it through single-reagent
# substitutions and report whether the shape score rises along the way.  A test
# utility: it explains recovery failures, it does not search.
# synthonScan: is the search's answer actually the best product in there?
# Exhaustive column scans plus a uniform random sample of the product space.
add_executable(synthonScan Tools/synthonScan.cpp)
target_link_libraries(synthonScan
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  GaussianShape FileParsers Threads::Threads)

add_executable(synthonReachability Tools/synthonReachability.cpp)
target_link_libraries(synthonReachability
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  GaussianShape FileParsers Threads::Threads)

add_executable(synthonSearchBench Tools/synthonSearchBench.cpp)
target_link_libraries(synthonSearchBench
  FragmentConfGen ChemReactions DistGeomHelpers ForceFieldHelpers MolAlign Descriptors
  GaussianShape FileParsers Threads::Threads)

if(RDK_BUILD_PYTHON_WRAPPERS)
 add_subdirectory(Wrap)
endif()
