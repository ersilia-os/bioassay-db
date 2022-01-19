import json
import requests
import os

PUBCHEM_PREFIX = "PUBCHEM"


class PubChemBioAssayRecord(object):

    def __init__(self, assay_id, batch_size=10000):
        self.batch_size=batch_size
        self.url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{}".format(assay_id)
        r = requests.get("{0}/sids/json".format(self.url))
        self.record = json.loads(r.text)

    def _chunker(self, seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    def _get_id(self, record):
        return record["InformationList"]["Information"][0]["AID"]

    def _get_sids(self, record):
        sids = record["InformationList"]["Information"][0]["SID"]
        sids.sort()
        return sids

    def _get_description(self, record):
        req = requests.get("{}/description/json".format(self.url))
        record = json.loads(req.text)
        return record["PC_AssayContainer"][0]["assay"]["descr"]

    def _get_data(self, record):
        data = []
        for chunk in self._chunker(self._get_sids(record), self.batch_size):
            s = ",".join([str(sid) for sid in chunk])
            r = requests.post("{0}/json".format(self.url), data={'sid': s})
            result = json.loads(r.text)["PC_AssaySubmit"]["data"]
            data.append(result)
        data = [d for chunk in data for d in chunk]
        return data

    def _get_cids_from_sids(self, sids):
        sid2cids = {}
        for chunk in self._chunker(sids, self.batch_size):
            s = ",".join([str(sid) for sid in chunk])
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/cids/json"
            r = requests.post(url, data={'sid': s})
            result = json.loads(r.text)
            if "InformationList" not in result:
                return {}
            result = result["InformationList"]["Information"]
            for res in result:
                if "CID" not in res:
                    continue
                sid2cids[res["SID"]] = res["CID"]
        return sid2cids


    def _get_smiles_from_cids(self, cids):
        cid2smiles = {}
        for chunk in self._chunker(cids, self.batch_size):
            s = ",".join([str(cid) for cid in list(set(cids))])
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/CanonicalSmiles,IsomericSmiles/JSON".format(s)
            r = requests.post(url, data={"cid": s})
            result = json.loads(r.text)["PropertyTable"]["Properties"]
            for res in result:
                if "CanonicalSMILES" in res:
                    cid2smiles[res["CID"]] = res["CanonicalSMILES"]
                elif "IsomericSMILES" in res:
                    cid2smiles[res["CID"]] = res["IsomericSMILES"]
                else:
                    continue
        return cid2smiles

    def _get_substances(self, sids):
        sid2cids = self._get_cids_from_sids(sids)
        all_cids = list(set([c for cid in list(sid2cids.values()) for c in cid]))
        if not all_cids:
            return [(sid, None, None) for sid in sids]
        cid2smiles = self._get_smiles_from_cids(all_cids)
        compounds = {}
        for sid in sids:
            if sid not in sid2cids:
                compounds[sid] = [None, None]
            else:
                for cid in sid2cids[sid]:
                    found = False
                    if not found and cid in cid2smiles:
                        compounds[sid] = [cid, cid2smiles[cid]]
                        found = True
                if not found:
                    compounds[sid] = [None, None]
        return compounds


    def _get_data_with_compounds(self, record):
        sids = self._get_sids(record)
        compounds = self._get_substances(sids)
        data = self._get_data(record)
        for i in range(len(data)):
            sid = data[i]["sid"]
            data[i]["cid"] = compounds[sid][0]
            data[i]["smiles"] = compounds[sid][1]
        return data


    def get(self):
        result = {
            "assay_id": self._get_id(self.record),
            "description": self._get_description(self.record),
            "data": self._get_data_with_compounds(self.record)
        }
        return {"AssayId": "{0}{1}".format(PUBCHEM_PREFIX, result["assay_id"]),
                "Description": result["description"],
                "Data": result["data"]}

    def save_json(self, path):
        data = self.get()
        assay_id = data["AssayId"]
        with open(os.path.join(path, "{}.json".format(assay_id)), 'w') as outfile:
            json.dump(data, outfile)
