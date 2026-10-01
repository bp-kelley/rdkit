#
#  Copyright (C) 2026 The RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import csv
import math
import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import rdNormalizedDescriptors as rdnd

refFile = os.path.join(RDConfig.RDBaseDir, 'External', 'NormalizedDescriptors', 'test_data',
                       'reference_descriptors.tsv')


class TestCase(unittest.TestCase):

  def testCDFTable(self):
    table = rdnd.CDFTable(0.0, 10.0, [0.0, 2.0, 4.0, 8.0], [0.1, 0.5, 0.7, 0.9])
    self.assertAlmostEqual(table.Normalize(1.0), 0.3)
    self.assertAlmostEqual(table(3.0), 0.6)
    self.assertAlmostEqual(table.Normalize(-5.0), 0.1)
    self.assertAlmostEqual(table.Normalize(100.0), 0.9)
    self.assertEqual(table.Normalize(math.nan), 0.0)
    self.assertEqual(table.GetXs(), (0.0, 2.0, 4.0, 8.0))
    self.assertEqual(table.GetCDF(), (0.1, 0.5, 0.7, 0.9))
    self.assertEqual(table.GetMin(), 0.0)
    self.assertEqual(table.GetMax(), 10.0)
    self.assertEqual(table.GetDistribution(), "")
    self.assertEqual(rdnd.CDFTable(0, 1, [0, 1], [0, 1], distribution="norm").GetDistribution(),
                     "norm")
    with self.assertRaises(ValueError):
      rdnd.CDFTable(0.0, 1.0, [1.0, 0.0], [0.0, 1.0])

  def testCDFTableSet(self):
    tables = rdnd.CDFTableSet()
    tables.LoadFromString("# comment\ndescriptor foo norm 0 10 3\n0 0\n5 0.5\n10 1\n")
    tables.AddTable("bar", rdnd.CDFTable(-1, 1, [-1, 1], [0.25, 0.75]))
    self.assertEqual(len(tables), 2)
    self.assertEqual(tables.GetNames(), ["bar", "foo"])
    self.assertAlmostEqual(tables.Normalize("foo", 2.5), 0.25)
    self.assertAlmostEqual(tables.Normalize("bar", 0.0), 0.5)
    self.assertEqual(tables.Normalize("baz", 1.0), 0.0)
    self.assertFalse(tables.HasTable("baz"))
    with self.assertRaises(KeyError):
      tables.GetTable("baz")
    tables2 = rdnd.CDFTableSet()
    tables2.LoadFromString(tables.ToString())
    self.assertEqual(tables2.GetNames(), tables.GetNames())
    self.assertEqual(tables2.GetTable("foo").GetCDF(), tables.GetTable("foo").GetCDF())
    with self.assertRaises(ValueError):
      rdnd.CDFTableSet().LoadFromString("descriptor foo norm 0 1 2\n0 0\n")

  def testDefaultTables(self):
    tables = rdnd.GetDefaultTables()
    self.assertGreaterEqual(len(tables), 201)
    self.assertTrue(tables.GetTable("MolLogP").GetDistribution())
    # reference values from descriptastorus' RDKit2DNormalized
    self.assertAlmostEqual(rdnd.NormalizeDescriptor("MolLogP", 2.5), 0.28147, places=3)
    self.assertAlmostEqual(tables.Normalize("ExactMolWt", 350.1), 0.34379, places=3)
    self.assertEqual(rdnd.NormalizeDescriptor("NotADescriptor", 1.0), 0.0)

  def testCalcNormalizedDescriptors(self):
    names = rdnd.GetNormalizedDescriptorNames()
    self.assertEqual(len(names), 217)
    with open(refFile) as inf:
      rows = list(csv.reader(inf, delimiter='\t'))
    self.assertEqual(tuple(rows[0][1:]), names)
    tables = rdnd.GetDefaultTables()
    for row in rows[1:]:
      mol = Chem.MolFromSmiles(row[0])
      raw = rdnd.CalcDescriptorValues(mol)
      normalized = rdnd.CalcNormalizedDescriptors(mol)
      self.assertEqual(len(raw), len(names))
      self.assertEqual(rdnd.CalcNormalizedDescriptors(mol, tables), normalized)
      for name, val, ref, nval in zip(names, raw, row[1:], normalized):
        ref = float(ref)
        self.assertAlmostEqual(val, ref, delta=1e-4 * max(1.0, abs(ref)), msg=f'{row[0]} {name}')
        self.assertAlmostEqual(nval, tables.Normalize(name, ref), delta=1e-4,
                               msg=f'{row[0]} {name}')
    # a table set without a descriptor gives 0.0 for it
    empty = rdnd.CalcNormalizedDescriptors(Chem.MolFromSmiles('CCO'), rdnd.CDFTableSet())
    self.assertEqual(set(empty), {0.0})


if __name__ == '__main__':
  unittest.main()
