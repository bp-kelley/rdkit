DO $$ BEGIN CREATE EXTENSION rdkit; EXCEPTION WHEN duplicate_object THEN NULL; END $$;

SELECT mol_numheavyatoms(mol_from_cdxml($cdxml$<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML><page><fragment><n id="1" p="168 142" Element="6"/><n id="2" p="193.98 127" Element="6"/><b id="3" B="1" E="2" Order="1"/></fragment></page></CDXML>$cdxml$));
SELECT mol_to_smiles(mol_from_cdxml($cdxml$<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML><page><fragment><n id="1" p="168 142" Element="6"/><n id="2" p="193.98 127" Element="6"/><b id="3" B="1" E="2" Order="1"/></fragment></page></CDXML>$cdxml$));
SELECT mol_to_smarts(qmol_from_cdxml($cdxml$<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML><page><fragment><n id="1" p="168 142" Element="6"/><n id="2" p="193.98 127" Element="6"/><b id="3" B="1" E="2" Order="1"/></fragment></page></CDXML>$cdxml$));
SELECT mol_numheavyatoms(mol_from_cdxml($cdxml$<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML><page><fragment><n id="1" p="168 142" Element="6"/><n id="2" p="193.98 127" Element="6"/><b id="3" B="1" E="2" Order="1"/></fragment></page></CDXML>$cdxml$));