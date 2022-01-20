import json
import pandas as pd
from collections import OrderedDict
import os


class PubChemBioAssayJsonConverter(object):
    def __init__(self, path, json_file):
        self.name = json_file.split(".")[0]
        with open(os.path.join(path, json_file), "r") as f:
            self.record = json.load(f)

    def _get_assay_id(self, record):
        assay_id = record["AssayId"]
        return assay_id

    def _get_description(self, record):
        description = record["Description"]
        return description

    def _get_data(self, record):
        data = record["Data"]
        return data

    def _get_tid_names(self):
        tid2name = {}
        results = self._get_description(self.record)["results"]
        for i in range(len(results)):
            tid = results[i]["tid"]
            name = results[i]["name"]
            tid2name[tid] = name
        return tid2name

    def _get_sid_cid_smiles(self):
        data = self._get_data(self.record)
        substances = {}
        for i in range(len(data)):
            substances[i] = [data[i]["sid"], data[i]["cid"], data[i]["smiles"]]
        return substances

    def _get_tid_results(self):
        data = self._get_data(self.record)
        sid2value = {}
        for i in range(len(data)):
            results = data[i]["data"]
            tid2value = OrderedDict()
            for i2 in range(len(results)):
                value = results[i2]["value"]
                for k, v in value.items():
                    val = v
                tid2value[results[i2]["tid"]] = val
            sid2value[i] = tid2value

        tid2name = self._get_tid_names()
        tid_int = [k for k in tid2name.keys()]

        for t in tid_int:
            for k, v in sid2value.items():
                tid_res = []
                for k2, v2 in v.items():
                    tid_res += [k2]
                if t not in tid_res:
                    v[t] = None

        for k, v in sid2value.items():
            for key in tid_int:
                v[key] = v.pop(key)

        return sid2value

    def _get_outcome(self):
        data = self._get_data(self.record)
        sid2outcome = {}
        for i in range(len(data)):
            sid2outcome[i] = data[i]["outcome"]
        return sid2outcome

    def _substances_to_df(self):
        substances = self._get_sid_cid_smiles()

        df = pd.DataFrame.from_dict(
            data=substances, orient="index", columns=["sid", "cid", "smiles"]
        )
        return df

    def _outcome_to_df(self):
        outcome = self._get_outcome()
        df = pd.DataFrame.from_dict(data=outcome, orient="index", columns=["outcome"])
        return df

    def _tid_to_df(self):
        sid2value = self._get_tid_results()
        tid2name = self._get_tid_names()
        tid_names = [k for k in tid2name.values()]
        df = pd.DataFrame.from_dict(data=sid2value, orient="index")
        df.columns = tid_names
        return df

    def get_description(self, path):
        descr = self._get_description(self.record)
        with open(os.path.join(path, "{}_descr.txt".format(self.name)), "w") as f:
            for k, v in descr.items():
                f.write("%s:%s\n" % (k, v))

    def get_outcome(self):
        df1 = self._substances_to_df()
        df2 = self._outcome_to_df
        df = pd.concat([df1, df2], axis=1)
        return df

    def get_all_results(self):
        df1 = self._substances_to_df()
        df2 = self._outcome_to_df()
        df3 = self._tid_to_df()
        df = pd.concat([df1, df2, df3], axis=1)
        return df

    def save_df(self, df, path):
        df.to_csv(os.path.join(path, "{}.csv".format(self.name)), index=False)
