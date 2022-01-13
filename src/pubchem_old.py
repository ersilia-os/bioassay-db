import json
import requests
from tqdm import tqdm

HEADER = "PC_AssaySubmit"
PUBCHEM_PREFIX = "PUBCHEM"

COMPOUND_QUERY_BATCH_SIZE = 100


class PubChemBioAssayRecordFromJson(object):

    def __init__(self, assay_id=None, assay_json_file=None, batch_size=None):
        if assay_json_file is not None:
            with open(assay_json_name, "r") as f:
                self.record = json.load(f)
        elif assay_id is not None:
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{0}/JSON".format(assay_id)
            print(url)
            result = requests.get(url)
            self.record = json.loads(result.text)
        else:
            raise Exception
        self.header = HEADER
        self.prefix = PUBCHEM_PREFIX
        if batch_size is None:
            self.batch_size = COMPOUND_QUERY_BATCH_SIZE
        else:
            self.batch_size = batch_size

    def get_id(self, record):
        return record[self.header]["assay"]["descr"]["aid"]["id"]

    def get_sids(self, record):
        sids = []
        for d in record[self.header]["data"]:
            sids += [d["sid"]]
        return sids

    def chunker(self, seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    def get_cids_from_sids(self, sids):
        s = ",".join([str(sid) for sid in sids])
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{0}/cids/JSON".format(s)
        r = requests.get(url)
        data = json.loads(r.text)
        if "InformationList" not in data:
            return {}
        data = data["InformationList"]["Information"]
        sid2cids = {}
        for d in data:
            if "CID" not in d:
                continue
            sid2cids[d["SID"]] = d["CID"]
        return sid2cids

    def get_smiles_from_cids(self, cids):
        s = ",".join([str(cid) for cid in list(set(cids))])
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{0}/property/CanonicalSmiles,IsomericSmiles/JSON".format(s)
        print(url)
        r = requests.get(url)
        data = json.loads(r.text)["PropertyTable"]["Properties"]
        cid2smiles = {}
        for d in data:
            if "CanonicalSMILES" in d:
                cid2smiles[d["CID"]] = d["CanonicalSMILES"]
            elif "IsomericSMILES" in d:
                cid2smiles[d["CID"]] = d["IsomericSMILES"]
            else:
                continue
        return cid2smiles

    def get_substances(self, sids):
        sid2cids = self.get_cids_from_sids(sids)
        all_cids = list(set([x for _,v in sid2cids.items() for x in v]))
        if not all_cids:
            return [(sid, None, None) for sid in sids]
        cid2smiles = self.get_smiles_from_cids(all_cids)
        compounds = []
        for sid in sids:
            if sid not in sid2cids:
                compounds += [(sid, None, None)]
            else:
                for cid in sid2cids[sid]:
                    found = False
                    if not found and cid in cid2smiles:
                        compounds += [(sid, cid, cid2smiles[cid])]
                        found = True
                if not found:
                    compounds += [(sid, None, None)]
        return compounds

    def get_compounds_from_sids(self, sids):
        done_sids = set()
        compounds = {}
        todo_sids = set(sids)
        i = 0
        while len(todo_sids) > 0:
            print("Attempt {0}".format(i))
            sids_ = list(todo_sids)
            for j, chunk in enumerate(self.chunker(sids_, self.batch_size)):
                print("Chunk {0}. Done: {1}".format(j, len(compounds)))
                for subs in tqdm(self.get_substances(chunk)):
                    compounds[subs[0]] = (subs[1], subs[2])
                    done_sids.update([subs[0]])
            todo_sids = todo_sids.difference(done_sids)
            i += 1
        compounds = [compounds[sid] for sid in sids]
        return compounds

    def get_description(self, record):
        return record[self.header]["assay"]["descr"]

    def get_data(self, record):
        return record[self.header]["data"]

    def get_data_with_compounds(self, record):
        sids = self.get_sids(record)
        compounds = self.get_compounds_from_sids(sids)
        data = self.get_data(record)
        for i in range(len(data)):
            cid = compounds[i][0]
            smiles = compounds[i][1]
            data[i]["cid"] = cid
            data[i]["smiles"] = smiles
        return data

    def get(self):
        result = {
            "assay_id": self.get_id(self.record),
            "description": self.get_description(self.record),
            "data": self.get_data_with_compounds(self.record)
        }
        return {"AssayId": "{0}{1}".format(PUBCHEM_PREFIX, result["assay_id"]), "Description": result["description"], "Data": result["data"]}
