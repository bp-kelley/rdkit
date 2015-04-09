//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RDCHEM_INCL_
#define _RDCHEM_INCL_

#include <GraphMol/RDKitBase.h>

namespace RDKit {
  class ConformerException;

  class ReadWriteMol : public RWMol {
  public:
    ReadWriteMol(const ROMol &m,bool quickCopy=false,int confId=-1) : RWMol(m,quickCopy,confId){
    };

    void RemoveAtom(unsigned int idx){
      removeAtom(idx);
    };
    void RemoveBond(unsigned int idx1,unsigned int idx2){
      removeBond(idx1,idx2);
    };
    int AddBond(unsigned int begAtomIdx,
                 unsigned int endAtomIdx,
                 Bond::BondType order=Bond::UNSPECIFIED)
    {
      return addBond(begAtomIdx,endAtomIdx,order);
    };
    int AddAtom(Atom *atom){
      PRECONDITION(atom,"bad atom");
      return addAtom(atom,true,false);
    };
    void ReplaceAtom(unsigned int idx,Atom *atom){
      replaceAtom(idx,atom);
    };
    ROMol *GetMol() const{
      ROMol *res=new ROMol(*this);
      return res;
    }
  };
  
}
void rdExceptionTranslator(RDKit::ConformerException const& x);

#endif

