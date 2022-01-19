from src.json2df import PubChemBioAssayJsonConverter

c = PubChemBioAssayJsonConverter("./", "PUBCHEM1851.json")
df = c.get_all_results()
c.save_df(df, ".")
