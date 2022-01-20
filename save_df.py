from src.json2df import PubChemBioAssayJsonConverter

c = PubChemBioAssayJsonConverter("./examples", "PUBCHEM400.json")
df = c.get_all_results()
c.save_df(df, "./examples")
c.get_description("./examples")
