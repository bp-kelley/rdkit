if(NOT RDK_BUILD_NORMALIZED_DESCRIPTORS)
  return()
endif()

rdkit_library(NormalizedDescriptors NormalizedDescriptors.cpp SHARED
  LINK_LIBRARIES RDGeneral)
target_compile_definitions(NormalizedDescriptors PRIVATE RDKIT_NORMALIZEDDESCRIPTORS_BUILD)

rdkit_headers(NormalizedDescriptors.h DEST GraphMol/NormalizedDescriptors)

rdkit_catch_test(normalizedDescriptorsTest catch_tests.cpp
  LINK_LIBRARIES NormalizedDescriptors)

if(RDK_BUILD_BOOST_PYTHON_WRAPPERS)
  add_subdirectory(Wrap)
endif()
