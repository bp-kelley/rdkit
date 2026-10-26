# Copyright (C) 2002-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" functions to match a bunch of fragment descriptors from a file

No user-servicable parts inside.  ;-)

"""
import os

from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolDescriptors

defaultPatternFileName = os.path.join(RDConfig.RDDataDir, 'FragmentDescriptors.csv')


def _CountMatches(mol, patt, unique=True):
  return len(mol.GetSubstructMatches(patt, uniquify=unique))


fns = []


_cxxFragmentNames = set(rdMolDescriptors.GetFragmentDescriptorNames())


def _LoadPatterns(fileName=None):
  useCxx = fileName is None
  if fileName is None:
    fileName = defaultPatternFileName
  try:
    with open(fileName, 'r') as inF:
      for line in inF.readlines():
        if len(line) and line[0] != '#':
          splitL = line.split('\t')
          if len(splitL) >= 3:
            name = splitL[0]
            descr = splitL[1]
            sma = splitL[2]
            descr = descr.replace('"', '')
            patt = Chem.MolFromSmarts(sma)
            if not patt or patt.GetNumAtoms() == 0:
              raise ImportError('Smarts %s could not be parsed' % (repr(sma)))
            name = name.replace('=', '_')
            name = name.replace('-', '_')
            if useCxx and name in _cxxFragmentNames:
              # the default patterns are also available in C++
              fn = lambda mol, countUnique=True, pattern=patt, name=name: (
                rdMolDescriptors.CalcFragmentDescriptor(mol, name)
                if countUnique else _CountMatches(mol, pattern, unique=False))
            else:
              fn = lambda mol, countUnique=True, pattern=patt: _CountMatches(
                mol, pattern, unique=countUnique)
            fn.__doc__ = descr
            fns.append((name, fn))
  except IOError:
    pass


_LoadPatterns()
for name, fn in fns:
  exec('%s=fn' % (name))
fn = None
