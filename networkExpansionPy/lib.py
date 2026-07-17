import json
import multiprocessing
import os
import pickle
import zipfile
from concurrent.futures import ProcessPoolExecutor
from copy import copy, deepcopy

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Try to import Rust acceleration module
try:
    import netexprs

    _HAS_RUST = True
except ImportError:
    _HAS_RUST = False

_rust_announced = False

def _announce_rust():
    global _rust_announced
    if not _rust_announced:
        print("🦀 netexprs")
        _rust_announced = True

# define asset path
asset_path, filename = os.path.split(os.path.abspath(__file__))
asset_path = asset_path + "/assets"


def load_metabolism(fname):
    return pd.read_pickle(asset_path + "/metabolic_networks/" + fname)


def load_compounds(fname):
    return pd.read_csv(asset_path + "/compounds/" + fname)


def netExp(R, P, x, b):
    k = np.sum(x)
    k0 = 0
    n_reactions = np.size(R, 1)
    y = csr_matrix(np.zeros(n_reactions))
    while k > k0:
        k0 = np.sum(x)
        y = np.dot(R.transpose(), x) == b
        y = y.astype("int")
        x_n = np.dot(P, y) + x
        x_n = x_n.astype("bool")
        x = x_n.astype("int")
        k = np.sum(x)
    return x, y


# single step of network expansion
def netExp_step(R, P, x, b):
    k = np.sum(x)
    k0 = 0
    n_reactions = np.size(R, 1)
    y = csr_matrix(np.zeros(n_reactions))

    k0 = np.sum(x)
    y = np.dot(R.transpose(), x) == b
    y = y.astype("int")
    x_n = np.dot(P, y) + x
    x_n = x_n.astype("bool")
    x = x_n.astype("int")
    k = np.sum(x)
    return x, y


# define a new network expansion, s.t. stopping criteria is now no new compounds or reactions can be added at subsequent iterations
def netExp_cr(R, P, x, b):
    k = np.sum(x)
    k0 = 0
    n_reactions = np.size(R, 1)
    y = csr_matrix(np.zeros(n_reactions))
    l = 0
    l0 = 0

    while (k > k0) | (l > l0):
        k0 = np.sum(x)
        l0 = np.sum(y)
        y = np.dot(R.transpose(), x) == b
        y = y.astype("int")
        x_n = np.dot(P, y) + x
        x_n = x_n.astype("bool")
        x = x_n.astype("int")
        k = np.sum(x)
        l = np.sum(y)
    return x, y


def netExp_trace(R, P, x, b):

    X = []
    Y = []

    X.append(x)
    k = np.sum(x)
    k0 = 0
    n_reactions = np.size(R, 1)
    y = csr_matrix(np.zeros(n_reactions))
    Y.append(y)

    while k > k0:
        k0 = np.sum(x)
        y = np.dot(R.transpose(), x) == b
        y = y.astype("int")
        x_n = np.dot(P, y) + x
        x_n = x_n.astype("bool")
        x = x_n.astype("int")
        k = np.sum(x)
        X.append(x)
        Y.append(y)
    return X, Y


def netContract(R, P, b, x_ac, y, yext, x_seeds=None):
    # code for running network contraction algorithm
    Xac = []
    Yac = []

    # x_ac = x.copy()
    Xac.append(x_ac)
    k0 = np.sum(x_ac)
    n_reactions = np.size(R, 1)
    # yactive = csr_matrix(np.multiply((y.toarray()),1-yext.toarray()))
    yactive = yext.__lt__(1).multiply(y)
    k0 = yactive.sum()
    Yac.append(yactive)

    # If no reactions are active, only seeds survive
    if k0 == 0 and x_seeds is not None:
        Xac.append(x_seeds)
        Yac.append(yactive)
        return Xac, Yac

    k = 0
    while k < k0:
        # find
        k0 = yactive.sum()
        x_ex = np.dot(P, yext).astype("bool").astype("int")
        x_ac = np.dot(P, yactive).astype("bool").astype("int")

        # Preserve seed compounds
        if x_seeds is not None:
            x_ac = (x_ac + x_seeds).astype("bool").astype("int")

        # find all extinct metabolites
        # metabolite has to be not active and extinct
        x_ex = x_ac.__lt__(1).multiply(x_ex)

        # compute extinct reactions
        yext = np.dot(R.transpose(), x_ex).astype("bool").astype("int")
        # compute active reactions

        # yactive = csr_matrix(np.multiply((yactive.toarray()),1-yext.toarray()))
        yactive = yext.__lt__(1).multiply(yactive)
        # k = np.sum(x_ac);
        k = yactive.sum()
        Xac.append(x_ac)
        Yac.append(yactive)

    return Xac, Yac


def parse_reaction_trace(reaction_trace, network):
    rxns_list = []
    for i in range(1, len(reaction_trace)):
        idx = reaction_trace[i].nonzero()[0]
        rxns = list(network.iloc[:, idx])
        rxns = pd.DataFrame(rxns, columns=["rn", "direction"])
        rxns["iter"] = i
        rxns_list.append(rxns)
    rxns_list = pd.concat(rxns_list, axis=0)
    return rxns_list


def isRxnCoenzymeCoupled(rxn, cosubstrate, coproduct):
    g = rxn[rxn.cid.isin([cosubstrate, coproduct])]
    out = False
    if len(g) > 1:
        if g.s.sum() == 0:
            out = True
    return out


def load_ecg_network(ecg):
    network_list = []
    consistent_rids = []
    for rid, v in ecg["reactions"].items():
        cids = v["left"] + v["right"]
        try:  ## This skips all reactions with n stoichiometries
            stoichs = [-int(i) for i in v["metadata"]["left_stoichiometries"]] + [
                int(i) for i in v["metadata"]["right_stoichiometries"]
            ]
            network_list += list(zip(cids, [rid for _ in range(len(stoichs))], stoichs))
        except:
            pass
        if v["metadata"]["element_conservation"] == True:
            consistent_rids.append(rid)
    return pd.DataFrame(network_list, columns=("cid", "rn", "s")), pd.DataFrame(
        consistent_rids, columns=["rn"]
    )


def load_json_network(rdict):
    network_list = []
    consistent_rids = []
    for rid, v in rdict.items():
        if v["glycans"] == False:  ## This skips all reactions with glycans
            cids = v["left"] + v["right"]
            try:  ## This skips all reactions with n stoichiometries
                stoichs = [-int(i) for i in v["left_stoichiometries"]] + [
                    int(i) for i in v["right_stoichiometries"]
                ]
                network_list += list(
                    zip(cids, [rid for _ in range(len(stoichs))], stoichs)
                )
            except:
                pass
            if v["element_conservation"] == True:
                consistent_rids.append(rid)
    return pd.DataFrame(network_list, columns=("cid", "rn", "s")), pd.DataFrame(
        consistent_rids, columns=["rn"]
    )


def _load_tuple_network(tlist):
    """
    Load a simple network, defined by 2-tuples of reactant lists and product lists

    Note: Intended for DEV only
    Note: Stoichiometry information in output is only accurate for directionality

        tlist_example = [
             (["A","B"],["C"]),
             (["C","D"],["E","F"]),
             (["E","F"],["G"]),
             (["G","H"],["I"]),
             (["A","J"],["I"])]
    """
    rows = list()
    for i, d in enumerate(tlist):
        if len(d) != 2:
            raise ValueError("reactions must be two-tuples")
        for cid in d[0]:
            rows.append({"rn": i, "cid": cid, "s": -1})
        for cid in d[1]:
            rows.append({"rn": i, "cid": cid, "s": 1})
    return pd.DataFrame(rows)


class GlobalMetabolicNetwork:
    """A metabolic network supporting expansion and contraction algorithms.

    Metabolic networks are loaded from a named database (default: "KEGG_OG")
    and typically prepared with::

        net = GlobalMetabolicNetwork()
        net.pruneInconsistentReactions()
        net.convertToIrreversible()

    After ``convertToIrreversible()``, each reaction is split into directed
    copies and indexed internally as ``(rn_id, direction)`` tuples, e.g.
    ``("R00001", "forward")``.

    **Reaction ID conventions**

    All public methods that accept reaction ID sets (``reaction_mask``,
    ``extinctReactions``, ``reactionScope``, etc.) accept either form:

    - Tuple: ``(rn_id, direction)`` — targets one directed copy
    - String: plain ``rn_id`` — targets both directed copies

    Normalization is handled by ``initialize_reaction_vector()``.

    **Internal vector conventions**

    Binary reaction vectors use opposite sign conventions depending on context:

    - Expansion *mask* vectors: ``1`` = reaction is **allowed**, ``0`` = excluded
    - Contraction *extinct* vectors: ``1`` = reaction **is extinct**, ``0`` = active

    The inversion between user-supplied ID sets (always "reactions to remove")
    and internal mask vectors (1 = keep) happens once, inside
    ``initialize_reaction_vector()``.

    **Rust acceleration**

    When ``netexprs`` is installed, compute-intensive methods automatically
    dispatch to a Rust backend with Rayon-based parallelism. All algorithms
    including ``trace`` use Rust when available. Falls back to Python
    silently if ``netexprs`` is unavailable.
    """
    

    def __init__(self, metabolism="KEGG_OG"):
        # load the data
        if metabolism == "KEGG_OG":
            network = pd.read_csv(asset_path + "/KEGG/network_full.csv")
            cpds = pd.read_csv(asset_path + "/compounds/cpds.txt", sep="\t")
            thermo = pd.read_csv(
                asset_path + "/reaction_free_energy/kegg_reactions_CC_ph7.0.csv",
                sep=",",
            )
            self._network = network
            self.thermo = thermo
            self.compounds = cpds  ## Includes many compounds without reactions

        elif metabolism == "ecg":
            with open(
                os.path.join(asset_path, "ecg", "master_from_kegg_2021-01-05.json")
            ) as f:
                ecg = json.load(f)
            network, consistent_rxns = load_ecg_network(ecg)
            self._network = network
            self.consistent_rxns = consistent_rxns
            self.compounds = pd.DataFrame(
                network["cid"].unique(), columns=["cid"]
            )  ## Only includes compounds with reactions

        elif metabolism == "KEGG":
            with zipfile.ZipFile(
                os.path.join(
                    asset_path,
                    "KEGG",
                    "2021.05.31-18.06.52",
                    "entries_detailed",
                    "reaction.json.zip",
                ),
                "r",
            ) as z:
                rdict = json.loads(z.read(z.infolist()[0]).decode())
            network, consistent_rxns = load_json_network(rdict)
            self._network = network
            self.consistent_rxns = consistent_rxns
            self.compounds = pd.DataFrame(
                network["cid"].unique(), columns=["cid"]
            )  ## Only includes compounds with reactions

        elif metabolism == "dev":
            ## Just for testing, etc.
            self._network = None

        else:
            raise (ValueError("'metabolism' must be one of 'KEGG_OG, 'ecg', 'KEGG'"))

        self.metabolism = metabolism
        self.temperature = 25
        self.seedSet = None
        self.rid_to_idx = None
        self.idx_to_rid = None
        self.cid_to_idx = None
        self.idx_to_cid = None
        self.S = None
        self._rust_arrays = None  # Lazily initialized cache for Rust arrays

    def __getstate__(self):
        state = self.__dict__.copy()
        # Ensure we save as _network
        if "network" in state and "_network" not in state:
            state["_network"] = state.pop("network")
        return state

    def __setstate__(self, state):
        # Handle old pickles that have 'network' instead of '_network'
        if "network" in state and "_network" not in state:
            state["_network"] = state.pop("network")
        # Handle old pickles that don't have '_rust_arrays'
        if "_rust_arrays" not in state:
            state["_rust_arrays"] = None
        self.__dict__.update(state)

    @property
    def network(self):
        return self._network

    @network.setter
    def network(self, value):
        self._network = value
        self._invalidate_rust_cache()

    def _ensure_dicts(self):
        """Ensure compound and reaction dictionaries are initialized."""
        if self.rid_to_idx is None or self.idx_to_rid is None:
            self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        if self.cid_to_idx is None or self.idx_to_cid is None:
            self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        if not hasattr(self, '_cid_list') or self._cid_list is None:
            self._cid_list = [self.idx_to_cid[i] for i in range(len(self.idx_to_cid))]
        if not hasattr(self, '_rid_list') or self._rid_list is None:
            self._rid_list = [self.idx_to_rid[i] for i in range(len(self.idx_to_rid))]
        # Object-dtype arrays holding the same (interned) ID objects as the
        # lists above. Fancy-indexing these with a nonzero-index array and
        # calling .tolist() copies the existing object references in C, which
        # is ~3x faster than a Python list comprehension when marshalling
        # large batches back to ID lists. The np.empty + slice-assign trick
        # avoids numpy interpreting the (rn, direction) reaction tuples as a
        # 2D array.
        if not hasattr(self, '_cid_arr') or self._cid_arr is None:
            self._cid_arr = np.empty(len(self._cid_list), dtype=object)
            self._cid_arr[:] = self._cid_list
        if not hasattr(self, '_rid_arr') or self._rid_arr is None:
            self._rid_arr = np.empty(len(self._rid_list), dtype=object)
            self._rid_arr[:] = self._rid_list
        if not hasattr(self, '_cid_set') or self._cid_set is None:
            self._cid_set = set(self.cid_to_idx.keys())

    def _ensure_rust_ready(self):
        """Lazily prepare arrays needed for Rust acceleration."""
        if self._rust_arrays is not None:
            return

        self._ensure_dicts()

        # Build R and P matrices
        R, P = self.create_RP_from_irreversible_network()
        R = csr_matrix(R)
        P = csr_matrix(P)

        # Compute b vector
        b = np.asarray(R.sum(axis=0)).ravel()

        # Convert R^T to CSR format for Rust
        R_T = R.T.tocsr()

        # Cache all arrays needed by Rust functions
        self._rust_arrays = {
            "rt_data": R_T.data.astype(np.float64),
            "rt_indices": R_T.indices.astype(np.int32),
            "rt_indptr": R_T.indptr.astype(np.int32),
            "p_data": P.data.astype(np.float64),
            "p_indices": P.indices.astype(np.int32),
            "p_indptr": P.indptr.astype(np.int32),
            "n_reactions": R.shape[1],
            "n_compounds": P.shape[0],
            "b": b.astype(np.float64),
            # Also cache sparse matrices for Python fallback
            "R": R,
            "P": P,
            "b_sparse": csr_matrix(b).T,
        }

    def _invalidate_rust_cache(self):
        """Invalidate cached Rust arrays. Call this after modifying the network."""
        self._rust_arrays = None
        self.rid_to_idx = None
        self.idx_to_rid = None
        self.cid_to_idx = None
        self.idx_to_cid = None
        self._cid_list = None
        self._rid_list = None
        self._cid_arr = None
        self._rid_arr = None
        self._cid_set = None

        global _rust_announced
        _rust_announced = False

    def _x_to_compounds(self, x_arr):
        """Convert compound boolean array to list of compound IDs."""
        if isinstance(x_arr, np.ndarray):
            cidx = np.nonzero(x_arr)[0]
        else:
            cidx = np.nonzero(x_arr.toarray().ravel())[0]
        # Object-array take + .tolist() copies the interned ID references in C,
        # far faster than a Python list comprehension for large scopes.
        return self._cid_arr[cidx].tolist()

    def _y_to_reactions(self, y_arr):
        """Convert reaction boolean array to list of reaction IDs."""
        if isinstance(y_arr, np.ndarray):
            ridx = np.nonzero(y_arr)[0]
        else:
            ridx = np.nonzero(y_arr.toarray().ravel())[0]
        return self._rid_arr[ridx].tolist()

    def _iters_to_dict(self, iter_arr, idx_to_id):
        """Convert a Rust iteration vector (int32) to {id: iteration} dict.

        Entries with value -1 (never reached) are excluded.
        """
        return {
            idx_to_id[int(i)]: int(iter_arr[i])
            for i in range(len(iter_arr))
            if iter_arr[i] >= 0
        }

    def _call_expand_batch(self, x_init_batch, masks=None):
        """Low-level batch expansion on raw numpy arrays (no ID conversion).

        Use this to chain operations without marshalling overhead. Build input
        arrays via initialize_metabolite_vector/matrix and
        initialize_reaction_vector/matrix, then pass uint8 arrays directly.

        Args:
            x_init_batch: (N × n_compounds) uint8 — seed compounds per row
            masks: optional (N × n_reactions) uint8 — 1=allowed, 0=masked

        Returns:
            (x_batch, y_batch) — (N × n_compounds) and (N × n_reactions) uint8
        """
        self._ensure_rust_ready()
        ra = self._rust_arrays
        return netexprs.expand_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], ra["n_reactions"],
            ra["p_data"], ra["p_indices"], ra["p_indptr"], ra["n_compounds"],
            x_init_batch, ra["b"], masks,
        )

    def _call_expand_trace_batch(self, x_init_batch, masks=None):
        """Call netexprs.expand_trace_batch with cached CSR arrays.

        Args:
            x_init_batch: (N × n_compounds) uint8 array
            masks: optional (N × n_reactions) uint8 array

        Returns:
            (c_iters_batch, r_iters_batch) — 2D int32 numpy arrays
        """
        self._ensure_rust_ready()
        ra = self._rust_arrays
        return netexprs.expand_trace_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], ra["n_reactions"],
            ra["p_data"], ra["p_indices"], ra["p_indptr"], ra["n_compounds"],
            x_init_batch, ra["b"], masks,
        )

    def _call_contract_batch(self, x_active_batch, y_active_batch, y_extinct_batch, x_seeds=None):
        """Low-level batch contraction on raw numpy arrays (no ID conversion).

        Use this to chain operations without marshalling overhead. Build input
        arrays via initialize_metabolite_vector/matrix and
        initialize_reaction_vector/matrix, then pass uint8 arrays directly.

        Args:
            x_active_batch: compound scope — either (n_compounds,) uint8 shared
                across all rows or (N × n_compounds) uint8 per row. Passing a 1D
                vector lets the Rust kernel broadcast it and avoids tiling.
            y_active_batch: reaction scope — either (n_reactions,) uint8 shared
                across all rows or (N × n_reactions) uint8 per row.
            y_extinct_batch: (N × n_reactions) uint8 — 1=extinct, 0=active
                (determines the batch size N)
            x_seeds: optional (n_compounds,) uint8 — seed compounds preserved
                each iteration (1D, shared across all batch rows)

        Returns:
            (x_batch, y_batch) — (N × n_compounds) and (N × n_reactions) uint8
        """
        self._ensure_rust_ready()
        ra = self._rust_arrays
        return netexprs.contract_batch(
            ra["rt_data"], ra["rt_indices"], ra["rt_indptr"], ra["n_reactions"],
            ra["p_data"], ra["p_indices"], ra["p_indptr"], ra["n_compounds"],
            x_active_batch, y_active_batch, y_extinct_batch, x_seeds,
        )

    def copy(self):
        return deepcopy(self)

    def set_ph(self, pH):
        if ~(type(pH) == str):
            pH = str(pH)
        if self.metabolism == "KEGG_OG":
            try:
                thermo = pd.read_csv(
                    asset_path
                    + "/reaction_free_energy/kegg_reactions_CC_ph"
                    + pH
                    + ".csv",
                    sep=",",
                )
                self.thermo = thermo
            except Exception as error:
                print("Failed to open pH files (please use 5.0-9.0 in 0.5 increments)")
        elif self.metabolism == "ecg":
            try:
                self.thermo = self.load_ecg_thermo(self.metabolism, pH)
            except:
                raise ValueError(
                    "Try another pH, that one appears not to be in the ecg json"
                )
        else:
            raise (
                NotImplementedError(
                    "pH not yet implemented for metabolism = %s" % self.metabolism
                )
            )

    def load_ecg_thermo(self, ph=9):
        thermo_list = []
        for rid, v in self.metabolism["reactions"].items():

            phkey = str(ph) + "pH_100mM"

            if v["metadata"]["dg"][phkey]["standard_dg_prime_value"] == None:
                dg = np.nan
            else:
                dg = v["metadata"]["dg"][phkey]["standard_dg_prime_value"]

            if v["metadata"]["dg"][phkey]["standard_dg_prime_error"] == None:
                dgerror = np.nan
            else:
                dgerror = v["metadata"]["dg"][phkey]["standard_dg_prime_error"]

            if v["metadata"]["dg"][phkey]["is_uncertain"] == None:
                note = "uncertainty is too high"
            else:
                note = np.nan

            thermo_list.append(
                (
                    rid,
                    dg,
                    dgerror,
                    v["metadata"]["dg"][phkey]["p_h"],
                    v["metadata"]["dg"][phkey]["ionic_strength"] / 1000,
                    v["metadata"]["dg"][phkey]["temperature"],
                    note,
                )
            )

        return pd.DataFrame(
            thermo_list,
            columns=(
                "!MiriamID::urn:miriam:kegg.reaction",
                "!dG0_prime (kJ/mol)",
                "!sigma[dG0] (kJ/mol)",
                "!pH",
                "!I (mM)",
                "!T (Kelvin)",
                "!Note",
            ),
        )

    def pruneInconsistentReactions(self):
        # remove reactions with qualitatively different sets of elements in reactions and products
        if self.metabolism == "KEGG_OG":
            consistent = pd.read_csv(
                asset_path + "/reaction_sets/reactions_consistent.csv"
            )
            self.network = self.network[self.network.rn.isin(consistent.rn.tolist())]
        elif self.metabolism == "KEGG" or self.metabolism == "ecg":
            self.network = self.network[
                self.network.rn.isin(self.consistent_rxns.rn.tolist())
            ]
        else:
            raise (
                NotImplementedError(
                    "Function not yet implemented for metabolism = %s" % self.metabolism
                )
            )
        self._invalidate_rust_cache()

    def pruneUnbalancedReactions(self):
        # only keep reactions that are elementally balanced
        if self.metabolism == "KEGG_OG":
            balanced = pd.read_csv(asset_path + "/reaction_sets/reactions_balanced.csv")
            self.network = self.network[self.network.rn.isin(balanced.rn.tolist())]
        else:
            raise (
                NotImplementedError(
                    "Function not yet implemented for metabolism = %s" % self.metabolism
                )
            )
        self._invalidate_rust_cache()

    def subnetwork(self, rxns):
        # only keep reactions that are in list
        self.network = self.network[self.network.rn.isin(rxns)]
        self._invalidate_rust_cache()

    def addGenericCoenzymes(self):
        replace_metabolites = {
            "C00003": "Generic_oxidant",
            "C00004": "Generic_reductant",
            "C00006": "Generic_oxidant",
            "C00005": "Generic_reductant",
            "C00016": "Generic_oxidant",
            "C01352": "Generic_reductant",
        }
        coenzyme_pairs = {}
        coenzyme_pairs["NAD"] = ["C00003", "C00004"]
        coenzyme_pairs["NADP"] = ["C00006", "C00005"]
        coenzyme_pairs["FAD"] = ["C00016", "C01352"]
        coenzyme_pairs = pd.DataFrame(coenzyme_pairs).T.reset_index()
        coenzyme_pairs.columns = ["id", "oxidant", "reductant"]
        # create reactions copies with coenzyme pairs
        new_rxns = []
        new_thermo = []
        for idx, rxn in self.network.groupby("rn"):
            z = any(
                [
                    isRxnCoenzymeCoupled(rxn, row.oxidant, row.reductant)
                    for x, row in coenzyme_pairs.iterrows()
                ]
            )
            if z:
                new_rxn = (
                    rxn.replace(replace_metabolites)
                    .groupby(["cid", "rn"])
                    .sum()
                    .reset_index()
                )
                new_rxn["rn"] = new_rxn["rn"] = idx + "_G"
                new_rxns.append(new_rxn)
                t = self.thermo[
                    self.thermo["!MiriamID::urn:miriam:kegg.reaction"] == idx
                ].replace({idx: idx + "_G"})
                new_thermo.append(t)

        new_rxns = pd.concat(new_rxns, axis=0)
        new_thermo = pd.concat(new_thermo, axis=0)

        self.network = pd.concat([self.network, new_rxns], axis=0)
        self.thermo = pd.concat([self.thermo, new_thermo], axis=0)
        self._invalidate_rust_cache()

    def convertToIrreversible(self):
        nf = self.network.copy()
        nb = self.network.copy()
        nf["direction"] = "forward"
        nb["direction"] = "reverse"
        nb["s"] = -nb["s"]
        net = pd.concat([nf, nb], axis=0)
        net = net.set_index(["cid", "rn", "direction"]).reset_index()
        self.network = net
        self._invalidate_rust_cache()

    def setMetaboliteBounds(self, ub=1e-1, lb=1e-6):

        self.network["ub"] = ub
        self.network["lb"] = lb

    def pruneThermodynamicallyInfeasibleReactions(self, keepnan=False):
        fixed_mets = ["C00001", "C00080"]

        if not hasattr(self, "thermo"):
            raise (AttributeError("Metabolism has no thermo data."))

        RT = 0.008309424 * (273.15 + self.temperature)
        rns = []
        dirs = []
        dgs = []
        for (rn, direction), dff in self.network.groupby(["rn", "direction"]):
            effective_deltaG = np.nan
            if rn in self.thermo["!MiriamID::urn:miriam:kegg.reaction"].tolist():
                deltaG = self.thermo[
                    self.thermo["!MiriamID::urn:miriam:kegg.reaction"] == rn
                ]["!dG0_prime (kJ/mol)"].values[0]
                if direction == "reverse":
                    deltaG = -1 * deltaG

                dff = dff[~dff["cid"].isin(fixed_mets)]
                subs = dff[dff["s"] < 0]
                prods = dff[dff["s"] > 0]
                k = np.dot(subs["ub"].apply(np.log), subs["s"]) + np.dot(
                    prods["lb"].apply(np.log), prods["s"]
                )

                effective_deltaG = RT * k + deltaG

            dgs.append(effective_deltaG)
            dirs.append(direction)
            rns.append(rn)

        res = pd.DataFrame({"rn": rns, "direction": dirs, "effDeltaG": dgs})
        if keepnan:
            # change effective free energy to negative number, so it passes the next filter
            res["effDeltaG"] = res["effDeltaG"].fillna(-1)
        else:
            res = res.dropna()

        # res = res[res['effDeltaG'] < 0].set_index(['rn','direction'])
        res = res[~(res["effDeltaG"] > 0)].set_index(["rn", "direction"])
        res = res.drop("effDeltaG", axis=1)
        self.network = res.join(
            self.network.set_index(["rn", "direction"])
        ).reset_index()
        self._invalidate_rust_cache()

    def pruneReactionsFromMetabolite(self, cpds):
        # find all reactions that use metabolites in cpds list
        reactions_to_drop = (
            self.network[self.network.cid.isin(cpds)].rn.unique().tolist()
        )
        reactions_to_keep = [
            x for x in self.network.rn.unique().tolist() if x not in reactions_to_drop
        ]
        # only keep reactions that do not use that metabolite
        self.subnetwork(reactions_to_keep)

    def initialize_metabolite_vector(self, seedSet):
        if seedSet is None:
            print("No seed set")
        else:
            x0 = np.zeros(len(self.cid_to_idx), dtype=int)
            indices = [self.cid_to_idx[x] for x in set(seedSet) & self._cid_set]
            if indices:
                x0[indices] = 1
            return x0

    def initialize_reaction_vector(self, reactionSet):
        """Build a binary reaction indicator vector from a set of reaction IDs.

        This is the single normalization gateway for all reaction ID sets entering
        the library. Accepts either ID representation:

        - Tuple ``(rn_id, direction)`` — sets that one directed copy to 1
        - String ``rn_id`` — sets both directed copies (forward and reverse) to 1

        Tuple keys are tried first; unmatched items are retried as plain string
        IDs expanded to both directions. Unknown IDs are silently ignored.

        Returns a vector of length ``n_reactions`` with dtype int, where 1
        indicates membership. The semantic meaning of that 1 (allowed, extinct,
        in-scope, etc.) is determined by the caller.
        """
        if reactionSet is None:
            print("No reactions in set")
            return
        x0 = np.zeros([len(self.rid_to_idx)], dtype=int)
        for rid in reactionSet:
            if rid in self.rid_to_idx:
                # exact match — tuple or any future string-keyed network
                x0[self.rid_to_idx[rid]] = 1
            else:
                # plain string → expand to both directed copies
                for direction in ("forward", "reverse"):
                    key = (rid, direction)
                    if key in self.rid_to_idx:
                        x0[self.rid_to_idx[key]] = 1
        return x0

    def initialize_reaction_matrix(self, reaction_id_sets):
        """Convert a list of reaction ID sets to a 2D binary matrix.

        Batch version of initialize_reaction_vector(). Each row corresponds to
        one set, using the same ID normalization rules (plain string expands to
        both forward and reverse directed copies).

        Args:
            reaction_id_sets: list of lists/sets of reaction IDs

        Returns:
            (N × n_reactions) uint8 numpy array
        """
        self._ensure_dicts()
        n = len(reaction_id_sets)
        mat = np.zeros((n, len(self.rid_to_idx)), dtype=np.uint8)
        rid_to_idx = self.rid_to_idx
        for i, id_set in enumerate(reaction_id_sets):
            for rid in id_set:
                if rid in rid_to_idx:
                    mat[i, rid_to_idx[rid]] = 1
                else:
                    for direction in ("forward", "reverse"):
                        key = (rid, direction)
                        if key in rid_to_idx:
                            mat[i, rid_to_idx[key]] = 1
        return mat

    def initialize_metabolite_matrix(self, compound_id_sets):
        """Convert a list of compound ID sets to a 2D binary matrix.

        Batch version of initialize_metabolite_vector().

        Args:
            compound_id_sets: list of lists/sets of compound IDs

        Returns:
            (N × n_compounds) uint8 numpy array
        """
        self._ensure_dicts()
        n = len(compound_id_sets)
        mat = np.zeros((n, len(self.cid_to_idx)), dtype=np.uint8)
        cid_to_idx = self.cid_to_idx
        cid_set = self._cid_set
        for i, id_set in enumerate(compound_id_sets):
            indices = [cid_to_idx[x] for x in set(id_set) & cid_set]
            if indices:
                mat[i, indices] = 1
        return mat

    def create_reaction_dicts(self):
        rids = set(zip(self.network["rn"], self.network["direction"]))
        rid_to_idx = dict()
        idx_to_rid = dict()
        for v, k in enumerate(rids):
            rid_to_idx[k] = v
            idx_to_rid[v] = k

        return rid_to_idx, idx_to_rid

    def create_compound_dicts(self):
        cids = set(self.network["cid"])
        cid_to_idx = dict()
        idx_to_cid = dict()
        for v, k in enumerate(cids):
            cid_to_idx[k] = v
            idx_to_cid[v] = k

        return cid_to_idx, idx_to_cid

    def create_S_from_irreversible_network(self):

        S = np.zeros([len(self.cid_to_idx), len(self.rid_to_idx)])

        for c, r, d, s in zip(
            self.network["cid"],
            self.network["rn"],
            self.network["direction"],
            self.network["s"],
        ):
            S[self.cid_to_idx[c], self.rid_to_idx[(r, d)]] = s

        return S

    def create_RP_from_irreversible_network(self):

        # only use entries in teh network dataframe with s < 0:
        reactant_network = self.network[self.network.s < 0]

        R = np.zeros([len(self.cid_to_idx), len(self.rid_to_idx)])

        for c, r, d, s in zip(
            reactant_network["cid"],
            reactant_network["rn"],
            reactant_network["direction"],
            reactant_network["s"],
        ):
            R[self.cid_to_idx[c], self.rid_to_idx[(r, d)]] = 1

        # only use entries in teh network dataframe with s > 0:
        product_network = self.network[self.network.s > 0]

        P = np.zeros([len(self.cid_to_idx), len(self.rid_to_idx)])

        for c, r, d, s in zip(
            product_network["cid"],
            product_network["rn"],
            product_network["direction"],
            product_network["s"],
        ):
            P[self.cid_to_idx[c], self.rid_to_idx[(r, d)]] = 1

        return R, P

    def create_iteration_dict(self, M, idx_to_id):
        idx_iter = dict()
        for i, row in enumerate(M):
            idxs = np.nonzero(row.toarray().T[0])[0]
            for idx in idxs:
                if idx not in idx_iter:
                    idx_iter[idx] = i

        id_iter = dict()
        for idx, i in idx_iter.items():
            id_iter[idx_to_id[idx]] = i

        return id_iter

    def expand(self, seedSet, algorithm="naive", excluded_reactions=None):
        """
        Run network expansion from a seed set of compounds.

        Args:
            seedSet: List of compound IDs to start expansion from
            algorithm: 'naive', 'cr', 'trace', or 'step'
            excluded_reactions: Optional list of reaction IDs to remove from
                the network before expanding. Accepts either tuple
                (rn_id, direction) or plain string rn_id (both directed
                copies are excluded). See initialize_reaction_vector().

        Returns:
            For 'naive', 'cr', 'step': (compounds, reactions) - lists of IDs in scope
            For 'trace': (compound_dict, reaction_dict) - dicts mapping ID to iteration
        """
        if _HAS_RUST and algorithm.lower() in ("naive", "trace"):
            _announce_rust()
            if algorithm.lower() == "trace":
                return self._expand_trace_rust(seedSet, excluded_reactions)
            else:
                return self._expand_rust(seedSet, excluded_reactions)
        else:
            return self._expand_python(seedSet, algorithm, excluded_reactions)

    def _expand_rust(self, seedSet, excluded_reactions=None):
        """Rust-accelerated network expansion (naive algorithm).

        Wraps a single seed set into a 1-row batch call to expand_batch.
        """
        self._ensure_rust_ready()
        ra = self._rust_arrays

        x0 = self.initialize_metabolite_vector(seedSet).astype(np.uint8)
        x_init_batch = x0.reshape(1, -1)

        masks = None
        if excluded_reactions is not None and len(excluded_reactions) > 0:
            exclude_vec = self.initialize_reaction_vector(excluded_reactions).astype(np.uint8)
            mask = (1 - exclude_vec).astype(np.uint8)
            masks = mask.reshape(1, -1)

        x_batch, y_batch = self._call_expand_batch(x_init_batch, masks)

        compounds = self._x_to_compounds(x_batch[0])
        reactions = self._y_to_reactions(y_batch[0])
        return compounds, reactions

    def _expand_trace_rust(self, seedSet, excluded_reactions=None):
        """Rust-accelerated trace expansion.

        Wraps a single seed set into a 1-row batch call to expand_trace_batch.
        Returns (compound_dict, reaction_dict) mapping IDs to iteration numbers.
        """
        self._ensure_rust_ready()

        x0 = self.initialize_metabolite_vector(seedSet).astype(np.uint8)
        x_init_batch = x0.reshape(1, -1)

        masks = None
        if excluded_reactions is not None and len(excluded_reactions) > 0:
            exclude_vec = self.initialize_reaction_vector(excluded_reactions).astype(np.uint8)
            mask = (1 - exclude_vec).astype(np.uint8)
            masks = mask.reshape(1, -1)

        c_iters, r_iters = self._call_expand_trace_batch(x_init_batch, masks)

        compound_dict = self._iters_to_dict(c_iters[0], self.idx_to_cid)
        reaction_dict = self._iters_to_dict(r_iters[0], self.idx_to_rid)
        return compound_dict, reaction_dict

    def _expand_python(self, seedSet, algorithm="naive", excluded_reactions=None):
        """Pure Python expansion (supports all algorithms).
 
        Internally uses a *keep* mask vector: 1 = reaction allowed, 0 = excluded.
        ``excluded_reactions`` (IDs to remove) is inverted once via
        ``initialize_reaction_vector()`` before being applied.
        """
        self._ensure_dicts()
 
        x0 = self.initialize_metabolite_vector(seedSet)
        R, P = self.create_RP_from_irreversible_network()
        b = sum(R)
 
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()
 
        if excluded_reactions is not None and len(excluded_reactions) > 0:
            exclude_vec = self.initialize_reaction_vector(excluded_reactions)
            mask = csr_matrix(np.diag(1 - exclude_vec))
            P = P * mask
            R = R * mask
 
        x0 = csr_matrix(x0)
        x0 = x0.transpose()
 
        if algorithm.lower() == "naive":
            x, y = netExp(R, P, x0, b)
        elif algorithm.lower() == "cr":
            x, y = netExp_cr(R, P, x0, b)
        elif algorithm.lower() == "trace":
            X, Y = netExp_trace(R, P, x0, b)
        elif algorithm.lower() == "step":
            x, y = netExp_step(R, P, x0, b)
        else:
            raise ValueError(
                "algorithm needs to be naive (compound stopping criteria) or cr (reaction/compound stopping criteria)"
            )
 
        if algorithm.lower() == "trace":
            compound_iteration_dict = self.create_iteration_dict(X, self.idx_to_cid)
            reaction_iteration_dict = self.create_iteration_dict(Y, self.idx_to_rid)
            return compound_iteration_dict, reaction_iteration_dict
        else:
            compounds = self._x_to_compounds(x)
            reactions = self._y_to_reactions(y)
            return compounds, reactions

    def contract(self, seedSet, reactionScope, compoundScope, extinctReactions):
        """
        Run network contraction starting from an expanded scope.

        Args:
            seedSet: List of compound IDs that are always available (preserved
                during contraction). These are compounds that entered the network
                as initial conditions and aren't produced by any reaction.
            reactionScope: List of reaction IDs in the current scope
            compoundScope: List of compound IDs in the current scope
            extinctReactions: List of reaction IDs to remove (extinct)

        Returns:
            (compounds, reactions) - lists of IDs remaining after contraction
        """
        if _HAS_RUST:
            _announce_rust()
            return self._contract_rust(reactionScope, compoundScope, extinctReactions, seedSet)
        else:
            return self._contract_python(reactionScope, compoundScope, extinctReactions, seedSet)

    def _contract_rust(self, reactionScope, compoundScope, extinctReactions, seedSet=None):
        """Rust-accelerated contraction (single, via 1-row batch)."""
        self._ensure_rust_ready()
        ra = self._rust_arrays

        x_active = self.initialize_metabolite_vector(compoundScope).astype(np.uint8).reshape(1, -1)
        y_active = self.initialize_reaction_vector(reactionScope).astype(np.uint8).reshape(1, -1)
        y_extinct = self.initialize_reaction_vector(extinctReactions).astype(np.uint8).reshape(1, -1)

        x_seeds = None
        if seedSet is not None:
            x_seeds = self.initialize_metabolite_vector(seedSet).astype(np.uint8)

        x_batch, y_batch = self._call_contract_batch(x_active, y_active, y_extinct, x_seeds)

        compounds = self._x_to_compounds(x_batch[0])
        reactions = self._y_to_reactions(y_batch[0])
        return compounds, reactions

    def _contract_python(self, reactionScope, compoundScope, extinctReactions, seedSet=None):
        """Pure Python contraction."""
        self._ensure_dicts()

        # create vectors for reactionScope, compoundScope
        xactive = self.initialize_metabolite_vector(compoundScope)
        yactive = self.initialize_reaction_vector(reactionScope)
        yextinct = self.initialize_reaction_vector(extinctReactions)

        R, P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        xactive = csr_matrix(xactive).transpose()
        yactive = csr_matrix(yactive).transpose()
        yextinct = csr_matrix(yextinct).transpose()

        x_seeds = None
        if seedSet is not None:
            x_seeds = csr_matrix(self.initialize_metabolite_vector(seedSet)).transpose()

        # run contraction algorithm
        X, Y = netContract(R, P, b, xactive, yactive, yextinct, x_seeds)
        x = X[-1]
        y = Y[-1]

        compounds = self._x_to_compounds(x)
        reactions = self._y_to_reactions(y)
        return compounds, reactions

    def run_expansions_batch(self, seedSets, maskedReactionSets=None, algorithm="naive"):
        """
        Run expansions with varying seeds and/or varying masks.

        Accepts flexible input shapes — a single seed set or mask set is
        broadcast across the other dimension automatically:

            - N seed sets, no masks → N unmasked expansions
            - N seed sets, N mask sets → N paired expansions
            - 1 seed set, N mask sets → seed broadcast across all masks
            - N seed sets, 1 mask set → mask broadcast across all seeds

        Args:
            seedSets: A single seed set (list of compound IDs) or a list of seed sets.
            maskedReactionSets: Optional. A single reaction set to remove, or a list
                of reaction sets to remove. If a list, length must match seedSets
                or be 1 (broadcast).
            algorithm: 'naive' or 'trace'

        Returns:
            For 'naive': (compoundScopes, reactionScopes) - lists of ID lists
            For 'trace': (compoundScopes, reactionScopes) - lists of {ID: iteration} dicts
        """
        # ── Normalize inputs to lists of lists ──
        # Detect single seed set: a flat list of compound IDs (strings)
        if len(seedSets) > 0 and not isinstance(seedSets[0], (list, set, frozenset)):
            seedSets = [seedSets]

        if maskedReactionSets is not None:
            if len(maskedReactionSets) > 0 and not isinstance(maskedReactionSets[0], (list, set, frozenset)):
                maskedReactionSets = [maskedReactionSets]

        # ── Broadcast to matching lengths ──
        n_seeds = len(seedSets)
        n_masks = len(maskedReactionSets) if maskedReactionSets is not None else 0

        if maskedReactionSets is not None:
            if n_seeds == 1 and n_masks > 1:
                seedSets = seedSets * n_masks
                n = n_masks
            elif n_masks == 1 and n_seeds > 1:
                maskedReactionSets = maskedReactionSets * n_seeds
                n = n_seeds
            elif n_seeds == n_masks:
                n = n_seeds
            else:
                raise ValueError(
                    f"seedSets (len {n_seeds}) and maskedReactionSets (len {n_masks}) "
                    f"must have matching lengths, or one must have length 1."
                )
        else:
            n = n_seeds

        # ── Dispatch ──
        if _HAS_RUST and algorithm.lower() in ("naive", "trace"):
            _announce_rust()
            self._ensure_rust_ready()
            ra = self._rust_arrays

            x_init_batch = np.zeros((n, ra["n_compounds"]), dtype=np.uint8)
            for i, s in enumerate(seedSets):
                x_init_batch[i] = self.initialize_metabolite_vector(s).astype(np.uint8)

            masks = None
            if maskedReactionSets is not None:
                masks = np.ones((n, ra["n_reactions"]), dtype=np.uint8)
                for i, rxns_removed in enumerate(maskedReactionSets):
                    exclude_vec = self.initialize_reaction_vector(rxns_removed).astype(np.uint8)
                    masks[i] = 1 - exclude_vec

            if algorithm.lower() == "trace":
                c_iters, r_iters = self._call_expand_trace_batch(x_init_batch, masks)
                compoundScopes = [self._iters_to_dict(c_iters[i], self.idx_to_cid) for i in range(n)]
                reactionScopes = [self._iters_to_dict(r_iters[i], self.idx_to_rid) for i in range(n)]
            else:
                x_batch, y_batch = self._call_expand_batch(x_init_batch, masks)
                compoundScopes = [self._x_to_compounds(x_batch[i]) for i in range(n)]
                reactionScopes = [self._y_to_reactions(y_batch[i]) for i in range(n)]

            return compoundScopes, reactionScopes
        else:
            compoundScopes = []
            reactionScopes = []
            for i, seedSet in enumerate(seedSets):
                excluded = maskedReactionSets[i] if maskedReactionSets else None
                cpds, rxns = self._expand_python(seedSet, algorithm, excluded)
                compoundScopes.append(cpds)
                reactionScopes.append(rxns)
            return compoundScopes, reactionScopes
    
    
    def run_expansions(self, seedSets, algorithm="naive"):
        """Run expansion for multiple seed sets.

        Delegates to run_expansions_batch with no masks.
        """
        return self.run_expansions_batch(seedSets, algorithm=algorithm)

    def run_expansions_reactionMasks(self, seedSet, maskedReactionSets):
        """Run masked expansions with a shared seed set.

        Delegates to run_expansions_batch — the single seed set is
        automatically broadcast across all masks.
        """
        return self.run_expansions_batch(seedSet, maskedReactionSets)

    def _run_expansions_python(self, seedSets, algorithm="naive"):
        """Pure Python batch expansion."""
        self._ensure_dicts()

        R, P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        compoundScopes = []
        reactionScopes = []
        for seedSet in seedSets:
            x0 = self.initialize_metabolite_vector(seedSet)
            x0 = csr_matrix(x0)
            x0 = x0.transpose()

            if algorithm.lower() == "naive":
                x, y = netExp(R, P, x0, b)
                compounds = self._x_to_compounds(x)
                reactions = self._y_to_reactions(y)
            elif algorithm.lower() == "trace":
                X, Y = netExp_trace(R, P, x0, b)
                compounds = self.create_iteration_dict(X, self.idx_to_cid)
                reactions = self.create_iteration_dict(Y, self.idx_to_rid)
            else:
                raise ValueError("please define algorithm (trace or naive)")

            compoundScopes.append(compounds)
            reactionScopes.append(reactions)
        return compoundScopes, reactionScopes

    def _run_expansions_reactionMasks_python(self, seedSet, maskedReactionSets):
        """Pure Python batch masked expansion."""
        self._ensure_dicts()

        R, P = self.create_RP_from_irreversible_network()
        b = sum(R)
        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        x0 = self.initialize_metabolite_vector(seedSet)
        x0 = csr_matrix(x0)
        x0 = x0.transpose()

        compoundScopes = []
        reactionScopes = []

        for rxns_removed in maskedReactionSets:
            # build new R and P matrices with Masks
            exclude_vec = self.initialize_reaction_vector(rxns_removed)
            reaction_mask = csr_matrix(np.diag(1 - exclude_vec))
            Pstar = P * reaction_mask
            Rstar = R * reaction_mask
            # run expansion algorithm
            x, y = netExp(Rstar, Pstar, x0, b)

            compounds = self._x_to_compounds(x)
            reactions = self._y_to_reactions(y)
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)

        return compoundScopes, reactionScopes

    def run_contractions(
        self, seedSet, reactionScope, compoundScope, extinctReactionSets
    ):
        """
        Run contraction for multiple extinction sets.

        Args:
            seedSet: List of compound IDs that are always available (preserved
                during contraction). These are compounds that entered the network
                as initial conditions and aren't produced by any reaction.
            reactionScope: List of reaction IDs in the current scope
            compoundScope: List of compound IDs in the current scope
            extinctReactionSets: List of extinction sets (each is a list of reaction IDs to remove)

        Returns:
            (compoundScopes, reactionScopes) - lists of results for each extinction set
        """
        if _HAS_RUST:
            _announce_rust()
            return self._run_contractions_rust(
                reactionScope, compoundScope, extinctReactionSets, seedSet
            )
        else:
            return self._run_contractions_python(
                reactionScope, compoundScope, extinctReactionSets, seedSet
            )

    def run_contractions_raw(
        self, seedSet, reactionScope, compoundScope, extinctReactionSets
    ):
        """Batch contraction returning raw uint8 arrays (no ID marshalling).

        Same inputs and semantics as :meth:`run_contractions`, but returns the
        raw Rust output arrays instead of converting each row back to lists of
        ID strings. This skips the dominant cost of ``run_contractions`` for a
        large batch — building N full ``list[list[str]]`` scopes under the GIL,
        which does not parallelize — for callers that only need vectorized
        numpy readouts (e.g. scope sizes, or losses relative to a wildtype
        scope, which is all the epistasis sweep computes).

        Rust-only: raises if the netexprs backend is unavailable.

        Args:
            seedSet: Compound IDs preserved during contraction (always
                available). Pass ``None`` to disable seed preservation.
            reactionScope: Reaction IDs in scope (shared across all rows).
            compoundScope: Compound IDs in scope (shared across all rows).
            extinctReactionSets: List of extinction sets (each a list of
                reaction IDs to remove). Defines the batch size N.

        Returns:
            ``(x_batch, y_batch)`` — ``(N × n_compounds)`` and
            ``(N × n_reactions)`` uint8 arrays. Row ``i`` holds the contracted
            compound / reaction indicator for ``extinctReactionSets[i]``. Column
            order matches ``idx_to_cid`` / ``idx_to_rid``; use
            ``initialize_metabolite_vector`` / ``initialize_reaction_vector`` to
            map IDs to the same column space (e.g. to diff against a wildtype
            scope).
        """
        if not _HAS_RUST:
            raise RuntimeError(
                "run_contractions_raw requires the netexprs Rust backend"
            )
        self._ensure_rust_ready()

        x_active = self.initialize_metabolite_vector(compoundScope).astype(np.uint8)
        y_active = self.initialize_reaction_vector(reactionScope).astype(np.uint8)

        # Per-row extinction matrix (single allocation, one Python loop inside).
        y_extinct_batch = self.initialize_reaction_matrix(extinctReactionSets)

        x_seeds = None
        if seedSet is not None:
            x_seeds = self.initialize_metabolite_vector(seedSet).astype(np.uint8)

        # Shared scope vectors passed 1D; Rust broadcasts across the batch.
        return self._call_contract_batch(x_active, y_active, y_extinct_batch, x_seeds)

    def _run_contractions_rust(self, reactionScope, compoundScope, extinctReactionSets, seedSet=None):
        """Rust-accelerated batch contraction returning ID-list scopes.

        Thin wrapper over :meth:`run_contractions_raw` that marshals the raw
        arrays back to per-row ID lists to preserve the ``run_contractions``
        return type. Callers that only need numpy readouts should use
        ``run_contractions_raw`` directly to skip this marshalling.
        """
        x_batch, y_batch = self.run_contractions_raw(
            seedSet, reactionScope, compoundScope, extinctReactionSets
        )

        n_batches = len(extinctReactionSets)
        compoundScopes = [self._x_to_compounds(x_batch[i]) for i in range(n_batches)]
        reactionScopes = [self._y_to_reactions(y_batch[i]) for i in range(n_batches)]

        return compoundScopes, reactionScopes

    def _run_contractions_python(
        self, reactionScope, compoundScope, extinctReactionSets, seedSet=None
    ):
        """Pure Python batch contraction."""
        self._ensure_dicts()

        # create vectors for reactionScope, compoundScope
        xactive = self.initialize_metabolite_vector(compoundScope)
        yactive = self.initialize_reaction_vector(reactionScope)

        R, P = self.create_RP_from_irreversible_network()
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        xactive = csr_matrix(xactive).transpose()
        yactive = csr_matrix(yactive).transpose()

        x_seeds = None
        if seedSet is not None:
            x_seeds = csr_matrix(self.initialize_metabolite_vector(seedSet)).transpose()

        compoundScopes = []
        reactionScopes = []

        for extinctReactions in extinctReactionSets:
            yextinct = self.initialize_reaction_vector(extinctReactions)
            yextinct = csr_matrix(yextinct).transpose()
            # run contraction algorithm
            X, Y = netContract(R, P, b, xactive, yactive, yextinct, x_seeds)
            x = X[-1]
            y = Y[-1]

            compounds = self._x_to_compounds(x)
            reactions = self._y_to_reactions(y)
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)

        return compoundScopes, reactionScopes

    def contract_and_reexpand(
        self, seedSet, reactionScope, compoundScope, extinctReactionSets
    ):
        """Contract and re-expand in one call, without intermediate ID conversion.

        This is the fast path for the common pattern of contracting a network
        scope under multiple extinction sets, then re-expanding from the
        surviving compounds. All marshalling between contraction output and
        expansion input stays in numpy — no Python ID strings are allocated
        unless you choose to convert specific rows afterward.

        **Why one extinction set works for both phases:** Contraction only
        operates on reactions present in ``reactionScope`` (y_active). Marking
        reactions outside the current scope as extinct has no effect — they're
        already inactive. So passing the full (un-intersected) extinction sets
        is equivalent to scope-intersecting them first. For re-expansion, you
        *want* the full set masked (reactions not yet in scope should stay
        masked). A single ``extinctReactionSets`` argument therefore produces
        correct behavior for both phases without manual scope-intersection.

        **Recovering ID lists from the output:** The returned arrays can be
        converted to compound/reaction ID lists for any row using
        ``_x_to_compounds`` / ``_y_to_reactions``::

            result = metabolism.contract_and_reexpand(...)
            cpds_row0 = metabolism._x_to_compounds(result['x_reexpanded'][0])
            rxns_row0 = metabolism._y_to_reactions(result['y_reexpanded'][0])

        Only pay the conversion cost for rows you actually need. For aggregate
        statistics (counts, presence matrices), operate directly on the arrays::

            counts = result['x_reexpanded'].sum(axis=1)          # (N,)
            presence = result['x_reexpanded'][1:].sum(axis=0)    # (n_compounds,)

        Args:
            seedSet: Seed compounds preserved during contraction.
            reactionScope: Reaction IDs in the current scope.
            compoundScope: Compound IDs in the current scope.
            extinctReactionSets: List of reaction ID sets to remove. Used as-is
                for both contraction (extinct vector) and re-expansion (inverted
                to a keep-mask). No scope-intersection needed.

        Returns:
            dict with keys:
                'x_contracted': (N × n_compounds) uint8 — post-contraction compounds
                'y_contracted': (N × n_reactions) uint8 — post-contraction reactions
                'x_reexpanded': (N × n_compounds) uint8 — post-re-expansion compounds
                'y_reexpanded': (N × n_reactions) uint8 — post-re-expansion reactions
        """
        if not _HAS_RUST:
            raise RuntimeError("contract_and_reexpand requires the netexprs Rust backend")

        _announce_rust()
        self._ensure_rust_ready()
        ra = self._rust_arrays
        n_batches = len(extinctReactionSets)

        # Marshal inputs once
        x_active = self.initialize_metabolite_vector(compoundScope).astype(np.uint8)
        y_active = self.initialize_reaction_vector(reactionScope).astype(np.uint8)
        y_extinct_batch = self.initialize_reaction_matrix(extinctReactionSets)

        x_seeds = None
        if seedSet is not None:
            x_seeds = self.initialize_metabolite_vector(seedSet).astype(np.uint8)

        # Contraction — shared scope vectors passed 1D; Rust broadcasts them.
        x_contracted, y_contracted = self._call_contract_batch(
            x_active,
            y_active,
            y_extinct_batch,
            x_seeds,
        )

        # Re-expansion: seeds are the contracted compound sets, masks invert extinction
        masks = (1 - y_extinct_batch).astype(np.uint8)
        x_reexpanded, y_reexpanded = self._call_expand_batch(x_contracted, masks)

        return {
            'x_contracted': x_contracted,
            'y_contracted': y_contracted,
            'x_reexpanded': x_reexpanded,
            'y_reexpanded': y_reexpanded,
        }

    def ne_output_to_graph(self, cpds, rxns):
        # build a network constructing metabolites from prior iteration to connecting subsequent iteration
        # input: cpds (rxns) is a dict with compound id (reaction id) as the key, and iteration as the value
        # output: dataframe with target and source node and iteration for

        graph = {"target": [], "source": [], "iteration": []}
        maxIter = np.array(list(cpds.values())).max()

        cpd_iter_df = pd.DataFrame(
            {"cid": [x[0] for x in cpds.items()], "iter": [x[1] for x in cpds.items()]}
        )

        for i in range(maxIter, 0, -1):
            molecules = [x[0] for x in cpds.items() if x[1] == i]
            reactions = [x[0] for x in rxns.items() if x[1] == i]
            reactions = self.network.set_index(["rn", "direction"]).loc[reactions]
            for molecule in molecules:
                # reactions
                # find any reactions that produce this metabolite
                r_sample = reactions[
                    (reactions.cid == molecule) & (reactions.s > 0)
                ].sample(1)
                # find metabolite in reaction that was produced in previous iteration
                r_sample_cpds = (
                    reactions.loc[r_sample.index]
                    .set_index("cid")
                    .join(cpd_iter_df.set_index("cid"))
                )
                r_sample_cpds = r_sample_cpds[r_sample_cpds.iter == i - 1].sample(1)
                molecules_origin = r_sample_cpds.index.tolist()[0]
                graph["target"].append(molecule)
                graph["source"].append(molecules_origin)
                graph["iteration"].append(i)

        graph = pd.DataFrame(graph)
        return graph

    def save(self, name):
        path_to_save = asset_path + "/metabolic_networks/" + name + ".pkl"
        with open(path_to_save, "wb") as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def rxns2tuple(self, rn_list):
        t = self.network[self.network.rn.isin(rn_list)][
            ["rn", "direction"]
        ].drop_duplicates()
        rn_list_tuple = list(zip(t.rn.tolist(), t.direction.tolist()))
        return rn_list_tuple

    def run_expansions_parallel(self, seedSets, algorithm="naive"):
        """
        Run expansion for multiple seed sets in parallel.

        Note: When Rust acceleration is available, this delegates to run_expansions()
        which batches all seed sets into a single Rust call parallelized via Rayon.
        The _parallel suffix is kept for API compatibility.

        Args:
            seedSets: List of seed sets (each is a list of compound IDs)
            algorithm: 'naive' or 'trace'

        Returns:
            (compoundScopes, reactionScopes) - lists of results for each seed set
        """
        if _HAS_RUST and algorithm.lower() == "naive":
            # Rust handles parallelization internally
            return self.run_expansions(seedSets, algorithm)
        else:
            # Fall back to Python multiprocessing
            return self._run_expansions_parallel_python(seedSets, algorithm)

    def _run_expansions_parallel_python(self, seedSets, algorithm="naive"):
        """Python multiprocessing-based parallel expansion."""
        self._ensure_dicts()

        R, P = self.create_RP_from_irreversible_network()
        b = np.sum(R, axis=0)

        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        compoundScopes = []
        reactionScopes = []

        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            results = executor.map(
                expansion_helper,
                [(self, seedSet, algorithm, R, P, b) for seedSet in seedSets],
            )

        for compounds, reactions in results:
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)

        return compoundScopes, reactionScopes

    def run_expansions_reactionMasks_parallel(self, seedSet, maskedReactionSets):
        """
        Run masked expansions in parallel.

        Note: When Rust acceleration is available, this delegates to run_expansions_reactionMasks()
        which uses Rust's native parallelization (Rayon). The _parallel suffix is kept
        for API compatibility.

        Args:
            seedSet: List of compound IDs to start expansion from
            maskedReactionSets: List of reaction sets to REMOVE

        Returns:
            (compoundScopes, reactionScopes) - lists of results for each mask
        """
        if _HAS_RUST:
            # Rust handles parallelization internally
            return self.run_expansions_reactionMasks(seedSet, maskedReactionSets)
        else:
            # Fall back to Python multiprocessing
            return self._run_expansions_reactionMasks_parallel_python(
                seedSet, maskedReactionSets
            )

    def _run_expansions_reactionMasks_parallel_python(
        self, seedSet, maskedReactionSets
    ):
        """Python multiprocessing-based parallel masked expansion."""
        self._ensure_dicts()

        R, P = self.create_RP_from_irreversible_network()
        b = np.sum(R, axis=0)

        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        compoundScopes = []
        reactionScopes = []

        with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            results = executor.map(
                expansion_helper_reaction_masks,
                [
                    (self, seedSet, rxns_removed, R, P, b)
                    for rxns_removed in maskedReactionSets
                ],
            )

        for compounds, reactions in results:
            compoundScopes.append(compounds)
            reactionScopes.append(reactions)

        return compoundScopes, reactionScopes


# define an example helper function for parrellel execution of netowrk expansion code
def expansion_helper(args):
    instance, seedSet, algorithm, R, P, b = args
    x0 = instance.initialize_metabolite_vector(seedSet)
    x0 = csr_matrix(x0)
    x0 = x0.transpose()

    if algorithm.lower() == "naive":
        x, y = netExp(R, P, x0, b)
        if x.toarray().sum() > 0:
            cidx = np.nonzero(x.toarray().T[0])[0]
            compounds = [instance.idx_to_cid[i] for i in cidx]
        else:
            compounds = []

        if y.toarray().sum() > 0:
            ridx = np.nonzero(y.toarray().T[0])[0]
            reactions = [instance.idx_to_rid[i] for i in ridx]
        else:
            reactions = []

    elif algorithm.lower() == "trace":
        X, Y = netExp_trace(R, P, x0, b)  # Assuming this function is already defined
        compounds = instance.create_iteration_dict(X, instance.idx_to_cid)
        reactions = instance.create_iteration_dict(Y, instance.idx_to_rid)
    else:
        raise ValueError("please define algorithm (trace or naive)")

    return compounds, reactions


# define an example helper function for parrellel execution of network expansion with reaction masks
def expansion_helper_reaction_masks(args):
    instance, seedSet, rxns_removed, R, P, b = args
    x0 = instance.initialize_metabolite_vector(seedSet)
    x0 = csr_matrix(x0)
    x0 = x0.transpose()

    exclude_vec = instance.initialize_reaction_vector(rxns_removed)
    reaction_mask = csr_matrix(np.diag(1 - exclude_vec))
    Pstar = P * reaction_mask
    Rstar = R * reaction_mask

    x, y = netExp(Rstar, Pstar, x0, b)

    if x.toarray().sum() > 0:
        cidx = np.nonzero(x.toarray().T[0])[0]
        compounds = [instance.idx_to_cid[i] for i in cidx]
    else:
        compounds = []

    if y.toarray().sum() > 0:
        ridx = np.nonzero(y.toarray().T[0])[0]
        reactions = [instance.idx_to_rid[i] for i in ridx]
    else:
        reactions = []

    return compounds, reactions