from src.pubchem import PubChemBioAssayRecord

c = PubChemBioAssayRecord(400)
c.save_json("./examples")
