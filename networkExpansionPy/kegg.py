## kegg.py
## Functions adapted from https://github.com/ELIFE-ASU/ecg/blob/master/ecg/kegg.py

import os
import re
import json
from Bio.KEGG import REST #, Enzyme, Compound, Map
import Bio.TogoWS as TogoWS
from tqdm import tqdm
from datetime import datetime
import shutil

def download_kegg(path=None):

    if path==None:
        curr_time = datetime.now().strftime('%Y.%m.%d-%H.%M.%S')
        path = os.path.join(os.path.abspath(os.getcwd()),"KEGG",curr_time)
    
    if not os.path.exists(path):
        os.makedirs(path)

    print("KEGG data will be downloaded to %s"%path)

    print("Downloading lists...")
    _download_lists(path)
    print("Downloading entries...")
    _download_entries(path)
    print("Detailing compounds...")
    _detail_compounds(path)
    print("Detailing reactions...")
    _detail_reactions(path)

def _download_lists(path,dbs=["reaction","compound"]):
    """
    Returns db.json of ids and names in dbs (default: rn, cpd).

    Used to grab entries.
    """

    lists_path = os.path.join(path, "lists")
    if not os.path.exists(lists_path):
        os.makedirs(lists_path)

    lists = _retrieve_lists(dbs)

    ## Write json of all entry ids and names
    for db in dbs:
        list_path = os.path.join(lists_path, db+".json")
        with open(list_path, 'w') as f:   
            json.dump(lists[db], f, indent=2)

def _retrieve_lists(dbs):
        
        lists = dict()
        for db in dbs:
        
            ## Retreive all entry ids and names
            id_name_dict = dict()
            raw_list = REST.kegg_list(db)
            id_name_list = [s.split('\t') for s in raw_list.read().splitlines()]
            for i in id_name_list:
                id_name_dict[i[0]] = i[1]

            lists[db] = list(id_name_dict.keys())
            
        return lists

def _download_entries(path,dbs=["reaction","compound"]):
        """
        Returns jsons of entries of dbs (default: rn, cpd).
        """

        ## Create dir to store entries in
        entries_path = os.path.join(path,"entries")
        if not os.path.exists(entries_path):
            os.makedirs(entries_path)

        ## Get entries on all kegg types
        for db in dbs:

            print("Downloading %s entries..."%db)

            ## Read list of all kegg ids
            list_path = os.path.join(path,"lists",db+".json")
            with open(list_path) as f:    
                list_data = json.load(f) #[0]

            ## Grab each entry in list
            entries = dict()
            for entry in tqdm(list_data):
                
                entry_id = entry.split(":")[1]
                # entry_fname = entry_id+".json"
                # entry_path = os.path.join(entries_path, entry_fname)

                while entry_id not in entries:
                    try:
                        handle = TogoWS.entry(db, entry_id, format="json") ## Will always return something, even if entry doesn't exist
                        entries[entry_id] = json.loads(handle.read())[0]
                        # with open(entry_path, 'a') as f:
                        #     f.write(handle.read())
                    except:
                        pass
            
            with open(os.path.join(entries_path,db+".json"), 'w') as f:
                json.dump(entries, f, indent=2)

def _detail_compounds(path):
    """
    Add information about elements in compounds
    """

    compound_path = os.path.join(path,'entries','compound.json')

    with open(compound_path) as f:
        compounds = json.load(f)

    for k,v in compounds.items():
        elements = re.findall(r"([A-Z][a-z]?)",v['formula'])
        v["elements"] = list(set(elements))

    ## Create detailed entry path if it doesn't exist
    entries_path_detailed = os.path.join(path,"entries_detailed")
    if not os.path.exists(entries_path_detailed):
        os.makedirs(entries_path_detailed)

    ## Write new file with added detail
    with open(os.path.join(entries_path_detailed,'compound.json'), 'w') as f:
        json.dump(compounds, f, indent=2)


def _detail_reactions(path):
        """
        Add reaction details in convenient fields for rn jsons.
        """

        reaction_path = os.path.join(path,'entries','reaction.json')
        compound_path_detailed = os.path.join(path,'entries_detailed','compound.json')

        with open(compound_path_detailed) as f:    
            compound_dict = json.load(f)#[0]

        with open(reaction_path) as f:
            reaction_dict = json.load(f)
        # for path in glob.glob(compound_path+"*.json"):
        #     with open(path) as f:
        #         compound_json = json.load(f) #[0]
        #         compound_dict[compound_json["entry_id"]] = compound_json
        
        for k,v in reaction_dict.items():
                
            equation = v["equation"]

            if re.search(r'(G\d+)',equation) == None: ## Only find entries without glycans

                for i, side in enumerate(equation.split(" <=> ")):

                    compounds = []
                    stoichiometries = []

                    ## match (n+1) C00001, (m-1) C00001 or similar
                    matches = re.findall(r'(\(\S*\) C\d+)',side)
                    # print matches
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match 23n C00001, m C00001 or similar
                    matches = re.findall(r'(\d*[n,m] C\d+)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\d*[n,m])',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match C06215(m+n), C06215(23m) or similar
                    matches = re.findall(r'(C\d+\(\S*\))',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "3 C00002" or similar (but NOT C00002 without a number)
                    matches = re.findall(r'(\d+ C\d+)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = match.split(' '+compound)[0]# re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "C00001 "at the start of the line (no coefficients)
                    matches = re.findall(r'(^C\d+) ',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = '1'
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "+ C00001 " (no coefficients)
                    matches = re.findall(r'(\+ C\d+ )',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "+ C00001" at the end of the line (no coefficients)
                    matches = re.findall(r'(\+ C\d+$)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "C00001" which is at the start and end of the line
                    matches = re.findall(r'(^C\d+$)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    if i==0:
                        v["left"] = compounds
                        v["left_stoichiometries"] = stoichiometries
                        v["left_elements"] = set()
                        ## Add element data
                        for c in compounds:
                            if c in compound_dict:
                                v["left_elements"] = v["left_elements"].union(compound_dict[c]['elements'])
                            else:
                                v["left_elements"] = v["left_elements"].union('missing_cid')
                                v["contains_missingcid"] = True
                    elif i==1:
                        v["right"] = compounds
                        v["right_stoichiometries"] = stoichiometries
                        v["right_elements"] = set()
                        ## Add element data
                        for c in compounds:
                            if c in compound_dict:
                                v["right_elements"] = v["right_elements"].union(compound_dict[c]['elements'])
                            else:
                                v["right_elements"] = v["right_elements"].union('missing_cid')
                                v["contains_missingcid"] = True
                    
                if "contains_missingcid" not in v:
                    v["contains_missingcid"] = False

                if v["left_elements"] != v["right_elements"]:
                    v["element_conservation"] = False
                    v["elements_mismatched"] = list(v["left_elements"]^v["right_elements"])
                else:
                    v["element_conservation"] = True
                    v["elements_mismatched"] = list()
                
                assert len(compounds) == len(stoichiometries)
                v["glycans"] = False

            else:

                v["glycans"] = True

        ## Create detailed entry path if it doesn't exist
        entries_path_detailed = os.path.join(path,"entries_detailed")
        if not os.path.exists(entries_path_detailed):
            os.makedirs(entries_path_detailed)

        ## Rewrite file with added detail
        with open(os.path.join(entries_path_detailed,'reaction.json'), 'w') as f:
            json.dump(reaction_dict, f, indent=2, default=_serialize_sets)

def _serialize_sets(obj):
    if isinstance(obj, set):
        return list(obj)

    return obj

def create_shutl_zips(path):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".json"):
                shutil.make_archive(os.path.join(root, file), 'zip', root, file)



if __name__ == "__main__":
    # download_kegg("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52")
    # with open("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52/entries_detailed/compound.json.zip") as f:
    #     myjson = json.load(f.decode("utf-8"))

    create_shutl_zips("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52/")

    ## Write zipped files
    # create_zips("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52/")

    ## Read zipped json
    # with zipfile.ZipFile("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52/entries_detailed/compound.json.zip","r") as z:
    #     myjson = json.loads(z.read(z.infolist()[0]).decode())

    # print(len(myjson))
    # print(myjson["C00001"])    
    # with open("/Users/harrison/ELSI/bioxp/networkExpansionPy/KEGG/2021.05.31-18.06.52/entries_detailed/compound.pkl", 'wb') as f:
    #     pickle.dump(myjson, f, pickle.HIGHEST_PROTOCOL)