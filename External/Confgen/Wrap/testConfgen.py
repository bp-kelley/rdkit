#  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import copy
import itertools
import os
import pickle
import sys
import time
import unittest

import numpy as np

from rdkit import Chem, Geometry, RDConfig, rdBase
from rdkit.Chem import AllChem, rdFragmentConfGen, rdChemReactions


def log(s):
  rdBase.LogErrorMsg("== " + s)


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ChemReactions', 'testData')

  def testSynthons3D(self):
    synthons = [
      [Chem.MolFromSmiles(m) for m in ['CCN[U]','C#CCN[U]']],
      [Chem.MolFromSmiles(m) for m in ['O=C([U])[C@@H]1CO1','C/C=C/CC(=O)[U]']],
    ]
    
    p = Chem.MolzipParams()
    p.setAtomSymbols(['U'])
    p.label = Chem.MolzipLabel.AtomType

    smiresults = set()
    for a in synthons[0]:
      for b in synthons[1]:
        z = Chem.molzip(a,b,p)
        smiresults.add(Chem.MolToSmiles(Chem.RemoveHs(z)))
        
    en = rdFragmentConfGen.EnumerateSynthons3D(synthons)

    results = []
    for result in en:
      for prodSet in result:
        for mol in prodSet:
          assert mol.GetNumConformers() > 0
          results.append(Chem.MolToSmiles(Chem.RemoveHs(mol)))

    print("molzip", smiresults)
    print("3d gen", results)
        
    self.assertEqual(set(results), set(smiresults))

    if rdChemReactions.EnumerateLibraryCanSerialize():
      pickle = en.Serialize()
      enumerator2 = rdFragmentConfGen.EnumerateSynthons3D()
      enumerator2.InitFromString(pickle)
      enumerator2.ResetState()

      results = []
      for result in enumerator2:
        for prodSet in result:
          for mol in prodSet:
            results.append(Chem.MolToSmiles(Chem.RemoveHs(mol)))

      self.assertEqual(set(results), set(smiresults))
    

      
      

  def test3DSearchAPI(self):
    """Build a small 3D library in memory, then drive the search from Python.

    Covers the surface needed to check this against an external implementation:
    product assembly, shape scoring, Thompson sampling and refinement.
    """
    synthons = [
        [Chem.MolFromSmiles(s) for s in ('CC(=O)[U]', 'CCC(=O)[U]', 'CCCC(=O)[U]')],
        [Chem.MolFromSmiles(s) for s in ('[U]NCC', '[U]NCCC', '[U]NC1CC1')],
    ]
    params = rdFragmentConfGen.EnumerateSynthons3DParams()
    params.numOutputConfs = 3
    lib = rdFragmentConfGen.EnumerateSynthons3D(synthons, params)
    self.assertTrue(lib.IsValid())
    self.assertEqual(lib.Arity(), 2)
    self.assertEqual([lib.NumReagents(i) for i in range(2)], [3, 3])

    # the graph-only zip must agree with the assembled product
    target = [1, 2]
    graph = lib.ZipProduct(target)
    product = lib.GetProduct(target)
    self.assertTrue(product.ok)
    self.assertGreater(product.mol.GetNumConformers(), 0)
    self.assertEqual(Chem.MolToSmiles(graph),
                     Chem.MolToSmiles(Chem.RemoveHs(product.mol)))

    # a product scored against one of its OWN poses overlays perfectly
    confId = product.mol.GetConformer(0).GetId()
    scorer = rdFragmentConfGen.ShapeProductScorer(product.mol, confId)
    self.assertAlmostEqual(scorer.Score(product.mol), 1.0, places=3)

    # ... and the search finds it again
    tp = rdFragmentConfGen.ThompsonSynthonParams()
    tp.budget = 60
    tp.numThreads = 1
    tp.numBestProducts = 3
    res = rdFragmentConfGen.ThompsonSynthonSearch(lib, scorer, tp)
    self.assertGreater(res.evaluations, 0)
    self.assertEqual(list(res.reagents), target)
    self.assertLessEqual(len(res.best), 3)

    refined = rdFragmentConfGen.RefineSynthons(lib, scorer, res.reagents, 3, 1, 3)
    self.assertEqual(list(refined.reagents), target)


if __name__ == '__main__':
  unittest.main()
