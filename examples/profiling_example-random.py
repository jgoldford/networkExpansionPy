import warnings
warnings.filterwarnings("ignore")
import networkExpansionPy.folds as nf
import networkExpansionPy.lib as ne
import pandas as pd
from pathlib import PurePath, Path
import random
from copy import copy

RANDOMSEED = 3144
random.seed(RANDOMSEED)
asset_path = nf.asset_path
ALGORITHM = "random_fold_order" #"random_fold_order"
WRITE = True # write result to disk
WRITE_TMP = False # write after each iteration
CUSTOM_WRITE_PATH = None # if writing result, custom path to write to
STR_TO_APPEND_TO_FNAME = "random_fold_ordering_%s-fixed-rn-seeds"%RANDOMSEED # if writing result, str to append to filename

METABOLISM_PATH = PurePath(asset_path, "metabolic_networks","metabolism.23Aug2022.pkl") # path to metabolism object pickle
RN2RULES_PATH = PurePath(asset_path, "rn2fold","rn2rules.20230224.pkl") # path to rn2rules object pickle
SEED_CPDS_PATH = PurePath(asset_path, "compounds", "seeds.Goldford2022.csv") # path to seed compounds csv

ORDERED_OUTCOME = False # ignore random seed and always choose folds based on sort order
IGNORE_REACTION_VERSIONS = True # when maximizing for reactions, don't count versioned reactions

## Metabolism
metabolism = pd.read_pickle(METABOLISM_PATH)

## FoldRules
rn2rules = pd.read_pickle(RN2RULES_PATH)
foldrules = nf.FoldRules.from_rn2rules(rn2rules)
popular_folds = set([
    2002,
    2007,
    7560,
    543,
    210,
    325,
    205,
    282,
    246,
    109])
popular_folds = set(str(i) for i in popular_folds)
foldrules = foldrules.subset_from_folds(popular_folds)

## Modify seeds with AA and GATP_rns
aa_cids = set(["C00037",
    "C00041",
    "C00065",
    "C00188",
    "C00183",
    "C00407",
    "C00123",
    "C00148",
    "C00049",
    "C00025"])

GATP_rns = {'R00200_gATP_v1',
    'R00200_gATP_v2',
    'R00430_gGTP_v1',
    'R00430_gGTP_v2',
    'R01523_gATP_v1',
    'R04144_gATP_v1',
    'R04208_gATP',
    'R04463_gATP',
    'R04591_gATP_v1',
    'R06836_gATP',
    'R06974_gATP',
    'R06975_gATP_v1'}

## Seed
seed = nf.Params(
    rns = set(metabolism.network["rn"]) - set(rn2rules) | GATP_rns,
    cpds = set((pd.read_csv(SEED_CPDS_PATH)["ID"])) | aa_cids,
    folds = set(['spontaneous'])
)

## Inititalize fold metabolism
fm = nf.FoldMetabolism(metabolism, foldrules, seed)
## Run fold expansion
result = fm.rule_order(algorithm=ALGORITHM, write=WRITE, write_tmp=WRITE_TMP, path=CUSTOM_WRITE_PATH, str_to_append_to_fname=STR_TO_APPEND_TO_FNAME, ordered_outcome=ORDERED_OUTCOME, ignore_reaction_versions=IGNORE_REACTION_VERSIONS)