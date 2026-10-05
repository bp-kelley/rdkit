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

  def testHistogramTable(self):
    table = rdnd.HistogramTable([0.0, 2.0, 4.0, 8.0], [0.1, 0.5, 0.7, 0.9])
    self.assertEqual(table.Normalize(2.0), 0.5)
    self.assertEqual(table.Normalize(1.0), 0.5)
    self.assertEqual(table(3.0), 0.7)
    self.assertEqual(table.Normalize(-5.0), 0.1)
    self.assertEqual(table.Normalize(100.0), 1.0)
    self.assertEqual(table.Normalize(math.nan), 0.0)
    self.assertEqual(table.GetEdges(), (0.0, 2.0, 4.0, 8.0))
    self.assertEqual(table.GetFractions(), (0.1, 0.5, 0.7, 0.9))
    with self.assertRaises(ValueError):
      rdnd.HistogramTable([1.0, 0.0], [0.0, 1.0])

  def testHistogramTableSet(self):
    tables = rdnd.HistogramTableSet()
    tables.LoadFromString("# comment\nhistogram foo 3\n0 0.2\n5 0.5\n10 1\n")
    tables.AddTable("bar", rdnd.HistogramTable([-1, 1], [0.25, 0.75]))
    self.assertEqual(len(tables), 2)
    self.assertEqual(tables.GetNames(), ["foo", "bar"])
    self.assertEqual(tables.GetTableIndex("bar"), 1)
    self.assertEqual(tables.GetTableIndex("baz"), -1)
    self.assertEqual(tables.GetTable(1).GetFractions(), tables.GetTable("bar").GetFractions())
    self.assertEqual(tables.Normalize(0, 2.5), 0.5)
    with self.assertRaises(IndexError):
      tables.GetTable(2)
    self.assertEqual(tables.Normalize("foo", 2.5), 0.5)
    self.assertEqual(tables.Normalize("bar", 0.0), 0.75)
    self.assertEqual(tables.Normalize("baz", 1.0), 0.0)
    self.assertFalse(tables.HasTable("baz"))
    with self.assertRaises(KeyError):
      tables.GetTable("baz")
    tables2 = rdnd.HistogramTableSet()
    tables2.LoadFromString(tables.ToString())
    self.assertEqual(tables2.GetNames(), tables.GetNames())
    self.assertEqual(tables2.GetTable("foo").GetEdges(), tables.GetTable("foo").GetEdges())
    self.assertEqual(tables2.GetTable("foo").GetFractions(),
                     tables.GetTable("foo").GetFractions())
    with self.assertRaises(ValueError):
      rdnd.HistogramTableSet().LoadFromString("histogram foo 2\n0 0\n")

  def testDefaultTables(self):
    tables = rdnd.GetDefaultTables()
    self.assertGreaterEqual(len(tables), 201)
    # reference values from descriptastorus' RDKit2DHistogramNormalized
    self.assertEqual(rdnd.NormalizeDescriptor("MolLogP", 2.5), 0.29897092796495756)
    self.assertEqual(tables.Normalize("ExactMolWt", 350.1), 0.37919654375806305)
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
        self.assertEqual(nval, tables.Normalize(name, val), msg=f'{row[0]} {name}')
      self.assertEqual(tables.NormalizeDescriptors(raw), normalized)
    with self.assertRaises(ValueError):
      tables.NormalizeDescriptors([1.0])
    # a table set without a descriptor gives 0.0 for it
    empty = rdnd.CalcNormalizedDescriptors(Chem.MolFromSmiles('CCO'), rdnd.HistogramTableSet())
    self.assertEqual(set(empty), {0.0})


if __name__ == '__main__':
  unittest.main()
