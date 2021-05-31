## kegg.py
## Functions adapted from https://github.com/ELIFE-ASU/ecg/blob/master/ecg/kegg.py

import os
import re
import json
import copy
import glob
import itertools
import warnings
import argparse
from Bio.KEGG import REST #, Enzyme, Compound, Map
import Bio.TogoWS as TogoWS
from tqdm import tqdm
from datetime import datetime

def download_kegg(path=None):

    if path==None:
        curr_time = datetime.now().strftime('%Y.%m.%d-%H.%M.%S')
        path = os.path.join(os.path.abspath(os.getcwd()),"KEGG",curr_time)
    
    if not os.path.exists(path):
        os.makedirs(path)

    print("KEGG data will be downloaded to %s"%path)

    _download_lists(path)
    _download_entries(path)

def _download_lists(path,dbs=["reaction","compound"]):
    """
    Returns db.json of ids and names in dbs (default: rn, cpd).

    Used to grab entries.
    """

    lists_path = os.path.join(path, "lists")
    if not os.path.exists(lists_path):
        os.makedirs(lists_path)

    lists = __retrieve_lists(dbs)

    ## Write json of all entry ids and names
    for db in dbs:
        list_path = os.path.join(lists_path, db+".json")
        with open(list_path, 'w') as f:   
            json.dump(lists[db], f, indent=2)

def __retrieve_lists(dbs):
        
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

            ## Read list of all kegg ids
            list_path = os.path.join(path,"lists",db+".json")
            with open(list_path) as f:    
                list_data = json.load(f) #[0]

            ## Grab each entry in list
            entries = dict()
            for entry in tqdm(list_data):
                
                entry_id = entry.split(":")[1]
                entry_fname = entry_id+".json"
                # entry_path = os.path.join(entries_path, entry_fname)

                while entry_fname not in os.listdir(entries_path):
                    try:
                        handle = TogoWS.entry(db, entry_id, format="json")
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


def serialize_sets(obj):
    if isinstance(obj, set):
        return list(obj)

    return obj