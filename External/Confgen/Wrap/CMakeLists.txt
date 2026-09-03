remove_definitions(-DRDKIT_CHEMREACTIONS_BUILD)
rdkit_python_extension(rdFragmentConfGen
                       confgen.cpp
                       DEST Chem
                       LINK_LIBRARIES
		       FragmentConfGen
                       MarvinParser
                       ChemReactions )


add_pytest(pyChemFragmentConfgen
         ${CMAKE_CURRENT_SOURCE_DIR}/testConfgen.py)


