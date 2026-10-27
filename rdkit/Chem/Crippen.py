# Copyright (C) 2002-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Atom-based calculation of LogP and MR using Crippen's approach


    Reference:
      S. A. Wildman and G. M. Crippen *JCICS* _39_ 868-873 (1999)


"""

import os

from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolDescriptors

# this is the file containing the atom contributions
defaultPatternFileName = os.path.join(RDConfig.RDDataDir, 'Crippen.txt')


_GetAtomContribs = rdMolDescriptors._CalcCrippenContribs


MolLogP = lambda *x, **y: rdMolDescriptors.CalcCrippenDescriptors(*x, **y)[0]
MolLogP.version = rdMolDescriptors._CalcCrippenDescriptors_version
MolLogP.__doc__ = """ Wildman-Crippen LogP value

  Uses an atom-based scheme based on the values in the paper:
     S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

  **Arguments**

    - inMol: a molecule

    - addHs: (optional) toggles adding of Hs to the molecule for the calculation.
      If true, hydrogens will be added to the molecule and used in the calculation.

"""

MolMR = lambda *x, **y: rdMolDescriptors.CalcCrippenDescriptors(*x, **y)[1]
MolMR.version = rdMolDescriptors._CalcCrippenDescriptors_version
MolMR.__doc__ = """ Wildman-Crippen MR value

  Uses an atom-based scheme based on the values in the paper:
     S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

  **Arguments**

    - inMol: a molecule

    - addHs: (optional) toggles adding of Hs to the molecule for the calculation.
      If true, hydrogens will be added to the molecule and used in the calculation.

"""

if __name__ == '__main__':
  import sys

  if len(sys.argv):
    ms = []
    verbose = 0
    if '-v' in sys.argv:
      verbose = 1
      sys.argv.remove('-v')
    for smi in sys.argv[1:]:
      ms.append((smi, Chem.MolFromSmiles(smi)))

    for smi, m in ms:
      print('Mol: %s' % (smi))
      logp = MolLogP(m, verbose=verbose)
      print('----')
      mr = MolMR(m, verbose=verbose)
      print('Res:', logp, mr)
      newM = Chem.AddHs(m)
      logp = MolLogP(newM, addHs=0)
      mr = MolMR(newM, addHs=0)
      print('\t', logp, mr)
      print('-*-*-*-*-*-*-*-*')
