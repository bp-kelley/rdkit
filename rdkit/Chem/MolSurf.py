#
# Copyright (C) 2001-2026 greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Exposes functionality for MOE-like approximate molecular surface area
descriptors.

  The MOE-like VSA descriptors are also calculated here

"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ptable = Chem.GetPeriodicTable()
bondScaleFacts = {
  Chem.BondType.AROMATIC: 0.1,
  Chem.BondType.SINGLE: 0.0,
  Chem.BondType.DOUBLE: 0.2,
  Chem.BondType.TRIPLE: 0.3,
}
#bondScaleFacts = [.1, 0, .2, .3]  # aromatic,single,double,triple


def _LabuteHelper(mol, includeHs=True, force=False):
  """ *Internal Use Only*
    helper function for LabuteASA calculation
    returns an array of atomic contributions to the ASA

  **Note:** Changes here affect the version numbers of all ASA descriptors

  """
  if not force:
    try:
      res = mol._labuteContribs
    except AttributeError:
      pass
    else:
      if res:
        return res
  tpl = rdMolDescriptors._CalcLabuteASAContribs(mol, bool(includeHs))
  ats, hs = tpl
  Vi = [hs] + list(ats)
  mol._labuteContribs = Vi
  return Vi


# def SMR_VSA(mol,bins=[0.11,0.26,0.35,0.39,0.44,0.485,0.56]):
# original default bins from assuming Labute values are logs
# mrBins=[1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63]
mrBins = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0]


SMR_VSA_ = rdMolDescriptors.SMR_VSA_

#
# Original bins (from Labute paper) are:
#  [-0.4,-0.2,0,0.1,0.15,0.2,0.25,0.3,0.4]
#
logpBins = [-0.4, -0.2, 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]


SlogP_VSA_ = rdMolDescriptors.SlogP_VSA_

chgBins = [-.3, -.25, -.20, -.15, -.10, -.05, 0, .05, .10, .15, .20, .25, .30]


PEOE_VSA_ = rdMolDescriptors.PEOE_VSA_


# -------------------------------------------------
# install the various VSA descriptors in the namespace
def _InstallDescriptors():
  for i in range(len(mrBins)):
    fn = lambda x, y=i: SMR_VSA_(x, force=False)[y]
    if i > 0:
      fn.__doc__ = "MOE MR VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, mrBins[i - 1],
                                                                        mrBins[i])
    else:
      fn.__doc__ = "MOE MR VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, mrBins[i])
    name = "SMR_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: SMR_VSA_(x, force=False)[y]
  fn.__doc__ = "MOE MR VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, mrBins[i - 1])
  fn.version = "1.0.1"
  name = "SMR_VSA%d" % (i + 1)
  globals()[name] = fn

  for i in range(len(logpBins)):
    fn = lambda x, y=i: SlogP_VSA_(x, force=False)[y]
    if i > 0:
      fn.__doc__ = "MOE logP VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, logpBins[i - 1],
                                                                          logpBins[i])
    else:
      fn.__doc__ = "MOE logP VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, logpBins[i])
    name = "SlogP_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: SlogP_VSA_(x, force=False)[y]
  fn.__doc__ = "MOE logP VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, logpBins[i - 1])
  fn.version = "1.0.1"
  name = "SlogP_VSA%d" % (i + 1)
  globals()[name] = fn

  for i in range(len(chgBins)):
    fn = lambda x, y=i: PEOE_VSA_(x, force=False)[y]
    if i > 0:
      fn.__doc__ = "MOE Charge VSA Descriptor %d (% 4.2f <= x < % 4.2f)" % (i + 1, chgBins[i - 1],
                                                                            chgBins[i])
    else:
      fn.__doc__ = "MOE Charge VSA Descriptor %d (-inf < x < % 4.2f)" % (i + 1, chgBins[i])
    name = "PEOE_VSA%d" % (i + 1)
    fn.version = "1.0.1"
    globals()[name] = fn
  i += 1
  fn = lambda x, y=i: PEOE_VSA_(x, force=False)[y]
  fn.version = "1.0.1"
  fn.__doc__ = "MOE Charge VSA Descriptor %d (% 4.2f <= x < inf)" % (i + 1, chgBins[i - 1])
  name = "PEOE_VSA%d" % (i + 1)
  globals()[name] = fn
  fn = None
  # Change log for the MOE-type descriptors:
  #  version 1.0.1: optimizations, values unaffected


_InstallDescriptors()


# Change log for LabuteASA:
#  version 1.0.1: optimizations, values unaffected
LabuteASA = lambda *x, **y: rdMolDescriptors.CalcLabuteASA(*x, **y)
LabuteASA.version = rdMolDescriptors._CalcLabuteASA_version


TPSA = lambda *x, **y: rdMolDescriptors.CalcTPSA(*x, **y)
TPSA.version = rdMolDescriptors._CalcTPSA_version

if __name__ == '__main__':  # pragma: nocover
  smis = ['C', 'CC', 'CCC', 'CCCC', 'CO', 'CCO', 'COC']
  smis = ['C(=O)O', 'c1ccccc1']
  for smi in smis:
    m = Chem.MolFromSmiles(smi)
    # print(smi, LabuteASA(m))
    print('-----------\n', smi)
    # print('M:',['% 4.2f'%x for x in SMR_VSA_(m)])
    # print('L:',['% 4.2f'%x for x in SlogP_VSA_(m)])
    print('P:', ['% 4.2f' % x for x in PEOE_VSA_(m)])
    print('P:', ['% 4.2f' % x for x in PEOE_VSA_(m)])
    print()
