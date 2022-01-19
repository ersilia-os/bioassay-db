from src.pubchem import PubChemBioAssayRecord

c = PubChemBioAssayRecord(1851)
c.save_json(".")
