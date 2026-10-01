#
#  Copyright (C) 2026 The RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import math
import unittest

from rdkit.Chem import rdNormalizedDescriptors as rdnd


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


if __name__ == '__main__':
  unittest.main()
