#
#  Copyright (C) 2026 RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Generates reference values from the python implementations of the
descriptors that were ported to C++ (GraphDescriptors.h, EState.h,
MiscDescriptors.h, FragmentDescriptors.h, QED.h and SpacialScore.h).

The output is used by catch_pyports.cpp. The reference file in this
directory was generated with the RDKit 2026.03.1 wheel from PyPI:

  python make_py_descriptor_reference.py input.smi py_descriptor_reference.tsv

The input is a SMILES file with the SMILES in the first column. The SMILES
are written out unchanged so that the C++ test parses exactly the same input.
"""
import gzip
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog('rdApp.*')

names = [
  'MaxAbsEStateIndex', 'MaxEStateIndex', 'MinAbsEStateIndex', 'MinEStateIndex', 'qed', 'SPS',
  'HeavyAtomMolWt', 'NumValenceElectrons', 'NumRadicalElectrons', 'MaxPartialCharge',
  'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge', 'FpDensityMorgan1',
  'FpDensityMorgan2', 'FpDensityMorgan3', 'AvgIpc', 'BalabanJ', 'BertzCT', 'Chi0', 'Chi1', 'Ipc'
]
names += [f'EState_VSA{i}' for i in range(1, 12)] + [f'VSA_EState{i}' for i in range(1, 11)]
names += [n for n, _ in Descriptors.descList if n.startswith('fr_')]


def main(inName, outName):
  fns = dict(Descriptors.descList)
  opener = gzip.open if inName.endswith('.gz') else open
  with opener(inName, 'rt') as inF, open(outName, 'w') as outF:
    outF.write('SMILES\t' + '\t'.join(names) + '\n')
    for line in inF:
      if not line.strip() or line.startswith('#'):
        continue
      smi = line.split()[0]
      m = Chem.MolFromSmiles(smi)
      if m is None or not m.GetNumAtoms():
        continue
      vals = []
      for nm in names:
        try:
          vals.append(repr(float(fns[nm](Chem.Mol(m)))))
        except Exception:
          vals.append('ERR')
      outF.write(smi + '\t' + '\t'.join(vals) + '\n')


if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2])
