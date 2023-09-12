import pandas as pd # only used in example main
from pathlib import PurePath, Path
from copy import deepcopy
import timeit
from pprint import pprint
from datetime import datetime
import random
import pickle
import gzip
import re

asset_path = PurePath(__file__).parent / "assets"

def get_versionless_reactions(reactions):
    versionless_reactions = list()
    for i in reactions:
        versionless_reactions.append(get_versionless_reaction(i))
    return set(versionless_reactions)

def get_versionless_reaction(reaction):
    match = re.match(r'(.+)_v\d', reaction)
    if match!=None:
        return match[1]
    else:
        return reaction

class ImmutableParams:
    """
    Stores parameters that commonly need to be kept track of together.

    After initialization attributes are immutable, to avoid unintentionally changing parameters like the seed or scope.
    """

    def __init__(self, folds=None, cpds=None, cpd_iteration_dict=None, rns=None, rn_iteration_dict=None, rules=None):
        self._folds = deepcopy(folds)
        self._cpds = deepcopy(cpds)
        self._cpd_iteration_dict = deepcopy(cpd_iteration_dict)
        self._rns = deepcopy(rns)
        self._rn_iteration_dict = deepcopy(rn_iteration_dict)
        self._rules = deepcopy(rules)

    @property
    def folds(self):
        return self._folds

    @property
    def cpds(self):
        return self._cpds

    @property
    def cpd_iteration_dict(self):
        return self._cpd_iteration_dict

    @property
    def rns(self):
        return self._rns

    @property
    def rn_iteration_dict(self):
        return self._rn_iteration_dict

    @property
    def rules(self):
        return self._rules

class Params(ImmutableParams):
    """
    Mutable child class of Params.

    Can be used to store parameters requiring updating, like during the main rule order loop.
    """

    @ImmutableParams.folds.setter
    def folds(self, value):
        self._folds = deepcopy(value)

    @ImmutableParams.cpds.setter
    def cpds(self, value):
        self._cpds = deepcopy(value)

    @ImmutableParams.cpd_iteration_dict.setter
    def cpd_iteration_dict(self, value):
        self._cpd_iteration_dict = deepcopy(value)

    @ImmutableParams.rns.setter
    def rns(self, value):
        self._rns = deepcopy(value)

    @ImmutableParams.rn_iteration_dict.setter
    def rn_iteration_dict(self, value):
        self._rn_iteration_dict = deepcopy(value)

    @ImmutableParams.rules.setter
    def rules(self, value):
        self._rules = deepcopy(value)

class Metadata:
    """
    A class for storing data that's useful to know but not critical for success of the expansion.

    Attributes:
        max_effects (dict): A dictionary containing all possible foldsets which could be injected after each iteration.
        size2foldsets (dict): A dictionary containing remaining foldsets keyed by size after each iteration.
    """
    def __init__(self):
        self.size2foldsets = None
        self.max_effects = None

class Result:
    """
    A class for storing data from the run.

    Attributes:
        iteration (int): The current iteration number.
        cpds (dict): A dictionary containing compounds and what iteration they appear.
        rns (dict): A dictionary containing reactions and what iteration they appear.
        folds (dict): A dictionary containing folds and what iteration they appear.
        rules (dict): A dictionary containing rules and what iteration they appear.
        start_datetime (str): A string representing the start date and time of the run.
        start_time (float): The start time of the run in seconds.
        iteration_time (dict): A dictionary containing iteration times.
        final_path (str): A string representing the path to the final result file.
        temp_path (str): A string representing the path to the temporary result file.
        max_effects (dict): A dictionary containing all possible foldsets which could be injected after each iteration.
        size2foldsets (dict): A dictionary containing remaining foldsets keyed by size after each iteration.

    Methods:
        first_update: Updates the object attributes on the first iteration.
        update: Updates the object attributes on any other iteration.
        update_cpds: Updates the cpds attribute.
        update_rns: Updates the rns attribute.
        update_folds: Updates the folds attribute.
        update_rules: Updates the rules attribute.
        update_iteration_time: Updates the iteration_time attribute.
        update_max_effects: Updates the max_effects attribute.
        update_size2foldsets: Updates the size2foldsets attribute.
        update_iter: Updates the iteration attribute.
        get_path: Returns the path to the result file.
        temp_write: Writes temporary results to a file.
        final_write: Writes final results to a file.
    """

    def __init__(self, scope):
        self.scope = deepcopy(scope)
        self.iteration = 0
        self.iteration_cum = 0
        self.cpds_folditer = dict()
        self.cpds_subiter = dict()
        self.cpds_cumiter = dict()
        self.rns_folditer = dict()
        self.rns_subiter = dict()
        self.rns_cumiter = dict()
        self.folds_folditer = {"fold_independent":0}
        self.folds_cumiter = {"fold_independent":0}
        # self.rules = dict() # activated; not simply possible
        self.rules_folditer = dict()
        self.rules_cumiter = dict()
        self.start_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.start_time = timeit.default_timer()
        self.iteration_time = dict()
        self.final_path = None
        self.temp_path = None

        self.max_effects = dict()
        self.size2foldsets = dict()

    def first_update(self, current, metadata=None, write=False, path=None, str_to_append_to_fname=None):
        self.update_cpds(current)
        self.update_rns(current)
        self.update_iteration_time()
        if write==True:
            self.temp_write(path=path, str_to_append_to_fname=str_to_append_to_fname)
        self.update_iter_cum(current)

    def update(self, current, metadata=None, write=False, path=None, str_to_append_to_fname=None):
        self.update_iter()
        self.update_folds(current)
        self.update_rules(current)
        if metadata != None:
            self.update_max_effects(metadata)
            self.update_size2foldsets(metadata)
        self.first_update(current, write=write, path=path, str_to_append_to_fname=str_to_append_to_fname)

    def update_cpds(self, current):
        if len(set(current.cpd_iteration_dict.keys()) - set(self.cpds_folditer.keys())) != 0:
            self.cpds_subiter[self.iteration] = deepcopy(current.cpd_iteration_dict)

        for i in current.cpds:
            if i not in self.cpds_folditer:
                self.cpds_folditer[i] = self.iteration
                self.cpds_cumiter[i] = self.iteration_cum + self.cpds_subiter[self.iteration][i]

    def update_rns(self, current):
        if len(set(current.rn_iteration_dict.keys()) - set(self.rns_folditer.keys())) != 0:
            self.rns_subiter[self.iteration] = deepcopy(current.rn_iteration_dict)

        for i in current.rns:
            if i not in self.rns_folditer:
                self.rns_folditer[i] = self.iteration
                self.rns_cumiter[i] = self.iteration_cum + self.rns_subiter[self.iteration][i]

    def update_folds(self, current):
        for i in current.folds:
            if i not in self.folds_folditer:
                self.folds_folditer[i] = self.iteration
                self.folds_cumiter[i] = self.iteration_cum+1 #+ self.rns_subiter[self.iteration][i]

    def update_rules(self, current):
        for i in current.rules.ids:
            if i not in self.rules_folditer:
                # self.rules[i] = self.iteration
                self.rules_folditer[i] = self.iteration
                self.rules_cumiter[i] = self.iteration_cum+1

    def update_iteration_time(self):
        self.iteration_time[self.iteration] =  timeit.default_timer() - self.start_time

    def update_max_effects(self, metadata):
        self.max_effects[self.iteration] = metadata.max_effects

    def update_size2foldsets(self, metadata):
        self.size2foldsets[self.iteration] = metadata.size2foldsets

    def update_iter(self):
        self.iteration+=1

    def update_iter_cum(self, current):
        if not len(current.rn_iteration_dict) == 0:
            self.iteration_cum += max(current.rn_iteration_dict.values())

    def get_path(self, path=None, str_to_append_to_fname=None):
        if str_to_append_to_fname == None:
            fname = self.start_datetime+".pkl.gz"
        else:
            fname = self.start_datetime+"_"+str_to_append_to_fname+".pkl.gz"
        
        if path==None:
            path = Path.joinpath(Path.cwd(), "fold_results", fname)
        else:
            path = Path.joinpath(path, "fold_results", fname)

        return path

    def temp_write(self, path=None, str_to_append_to_fname=None):
        ## Doesn't check if overwriting

        path = self.get_path(path, str_to_append_to_fname)
        path = Path.joinpath(path.parent, path.stem+"_tmp"+path.suffix)
        path.parent.mkdir(parents=True, exist_ok=True) 
        with gzip.open(path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        if self.temp_path == None:
            self.temp_path = str(path)
            print("Temporary results written to:\n{}".format(self.temp_path))

    def final_write(self, path=None, str_to_append_to_fname=None):
        
        path = self.get_path(path, str_to_append_to_fname)
        
        i = 0
        while path.is_file():
            path = Path.joinpath(path.parent, path.stem+"_"+str(i)+path.suffix)
            i+=1

        path.parent.mkdir(parents=True, exist_ok=True) 
        with gzip.open(path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        self.final_path = str(path)
        print("Final results written to:\n{}".format(self.final_path))

class Rule:
    """
    A class representing a fold rule.

    Attributes:
        rn: the rule's associated rn
        foldset (frozenset): a set of folds that trigger the rule
        id (tuple): a unique identifier for the rule (rn, foldset)

    Methods:
        __hash__(): returns the hash value of the rule's id
        __eq__(other): returns True if the other rule has the same id
        __repr__(): returns a string representation of the rule

    Usage:
        Create a new Rule object by providing the rule id and a set of folds that trigger the rule.
        The id attribute is automatically generated and used to compare rules for equality and hashing.
    """

    def __init__(self,rn,foldset:frozenset):
        self.rn = rn
        self.rn_versionless = get_versionless_reaction(rn)
        self.foldset = foldset
        self.id = (rn, foldset)

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if isinstance(other, Rule):
            return self.id == other.id
        return False

    def __repr__(self):
        return "id:\t\t{0} \nrn:\t\t{1} \nfoldset:\t{2}".format(self.id, self.rn, self.foldset)

class FoldRules:
    """
    Object that stores a universe of possible fold rules.

    Attributes:
        rules (list): a list of Rule objects
        rns (set): the set of unique rns in FoldRules
        folds (set): the set of unique folds in FoldRules
        foldsets (set): the set of unique foldsets in FoldRules
        ids (set): the set of unique ids in FoldRules

    Methods:
        from_rn2rules(rn2rules): create a new FoldRules object from a dictionary of rns to foldsets
        subset_from_rns(rns): create a new FoldRules object containing only rules with the given rn
        subset_from_folds(folds): create a new FoldRules object containing only rules whose foldsets are entirely within the given folds
        remaining_rules(current_folds): create a new FoldRules object containing only rules remaining after all rules possible with current_folds is accounted for
        foldset2rules(): return a dictionary mapping frozensets of folds to lists of rules
        to_list(): return a list of Rule objects
        __len__(): return the number of rules in FoldRules
        __iter__(): return an iterator over the Rule objects

    Usage:
        Create a new FoldRules object by providing a list of Rule objects. You can also use the `from_rn2rules` method
        to create a new FoldRules object from a dictionary of reactions to frozensets of folds. The `subset_from_rns`,
        `subset_from_folds`, and `remaining_rules` methods allow you to create new FoldRules objects that are subsets of
        the original universe based on reactions or folds. The `foldset2rules` method returns a dictionary that maps
        frozensets of folds to rules.
    """

    def __init__(self,rules:list):
        self._rules = rules
        self._rns = None
        self._folds = None 
        self._foldsets = None
        self._ids = None

    def __repr__(self):
        return "[\n"+",\n".join([str(i) for i in self.rules])+"]"

    @classmethod
    def from_rn2rules(cls, rn2rules):
        rules = list()
        for rn, fs_set in rn2rules.items():
            for fs in fs_set:
                rules.append(Rule(rn, fs))

        return cls(rules)

    @property
    def rules(self):
        return self._rules

    @property
    def rns(self):
        if self._rns == None:
            self._rns = set([r.rn for r in self.rules])
        return self._rns

    @property
    def folds(self):
        if self._folds == None:
            self._folds = set([f for r in self.rules for f in r.foldset])
        return self._folds

    @property 
    def foldsets(self):
        if self._foldsets == None:
            self._foldsets = set([r.foldset for r in self.rules])
        return self._foldsets

    @property 
    def ids(self):
        if self._ids == None:
            self._ids = set([r.id for r in self.rules])
        return self._ids

    @property
    def versionless(self):
        return FoldRules([Rule(r.rn_versionless, r.foldset) for r in self.rules])

    def subset_from_rns(self, rns):
        return FoldRules([r for r in self.rules if r.rn in rns])

    def subset_from_folds(self, folds):
        return FoldRules([r for r in self.rules if r.foldset <=folds])

    def remaining_rules(self, current_folds):
        return FoldRules([r for r in self.rules if len(r.foldset-current_folds)>0])

    def foldset2rules(self):
        foldset2rules = {k:list() for k in self.foldsets}
        for r in self.rules:
            foldset2rules[r.foldset].append(r)
        return foldset2rules

    def to_list(self):
        return [r for r in self.rules]

    def __len__(self):
        return len(self.ids)

    def __iter__(self):
        return iter(self.rules)

class FoldMetabolism:
    """
    A class to do fold expansion from a metabolism, foldrules, and seeds.

    Attributes:        
        metabolism (GlobalMetabolicNetwork): defines the compound/reaction rules of network expansion
        foldrules (FoldRules): defines the fold rules of fold expansion
        seed (Params): a Params object that defines the initial compounds, folds, and reactions (these are fold independent reactions)
        """

    def __init__(self, metabolism, foldrules, seed, scope=None):#, preexpansion=False):        
        self._m = metabolism ## GlobalMetabolicNetwork object
        self._f = foldrules # FoldRules object
        self._seed = ImmutableParams(folds=seed.folds, rns=seed.rns, cpds=seed.cpds) ## seed.rns == fold_independent_rns
        if scope==None:
            self._scope = self.calculate_scope(seed)
        else:
            self._scope = scope

    ## Disallow changing metabolism or foldrules after initialization b/c no setter
    @property 
    def metabolism(self):
        return self._m
    @property 
    def m(self):
        return self._m
    
    @property 
    def foldrules(self):
        return self._f
    @property 
    def f(self):
        return self._f

    @property
    def seed(self):
        return self._seed

    @property 
    def scope(self):
        return self._scope
    def calculate_scope(self, seed):
        """
        Calculates the scope of the seeds with all reactions enabled by the global fold network (including fold independent reactions i.e. seed.rns)

        Arguments:
            seed (ImmutableParams): object with seed information
        
        Outputs:
            ImmutableParams object specificying the scope
        """
        print("calculating scope...")
        rn_tup_set = set(self.m.rxns2tuple((self.f.rns | self.seed.rns)))
        scope_cpds, scope_rns = self.m.expand(seed.cpds, reaction_mask=rn_tup_set)

        scope = Params()
        scope.cpds = set(scope_cpds)
        scope.rns = set([i[0] for i in scope_rns])
        scope.rules = self.f.subset_from_rns(scope.rns)
        scope.folds = scope.rules.folds 
        print("...done.")
        return ImmutableParams(folds=scope.folds, rns=scope.rns, rules=scope.rules, cpds=scope.cpds)

    def fold_expand(self, folds, current_cpds, fold_algorithm="trace"):
        """
        Expand the given set of compounds and rules using the specified fold algorithm and a subset of the rules based on the given folds.

        Arguments:
            folds (list): collection of folds for expansion
            current_cpds (list): collection of compounds
            fold_algorithm (str, optional): The name of the algorithm to use for expanding the compounds and rules. Defaults to "naive".

        Returns:
            A tuple containing two sets from the expansion: 
            1. Set of compounds
            2. Set of reactions
        """

        possible_rules = self.f.subset_from_folds(folds)
        rn_tup_set = set(self.m.rxns2tuple(possible_rules.rns | self.seed.rns))
        if fold_algorithm=="trace":
            compound_iteration_dict, reaction_iteration_dict = self.m.expand(current_cpds | self.seed.cpds, algorithm=fold_algorithm, reaction_mask=rn_tup_set)
            return compound_iteration_dict, {k[0]:v for k,v in reaction_iteration_dict.items()}#set(cx), set([i[0] for i in rx])
        elif fold_algorithm=="step":
            cx,rx = self.m.expand(current_cpds | self.seed.cpds, algorithm=fold_algorithm, reaction_mask=rn_tup_set)
            return set(cx), set([i[0] for i in rx])

    def sort_remaining_foldsets_by_size(self, current_folds):
        """
        Sorts and returns the remaining foldsets based on their size after removing the current foldsets.

        Arguments:
            current_folds (set): A set of current folds to be removed from all foldsets.

        Returns:
            A dict of foldsets keyed by their size after removing the current_folds.
        """
        remaining_rules = self.scope.rules.remaining_rules(current_folds) ## the already found rules aren't necessarily activated
        remaining_foldsets = set([i-current_folds for i in remaining_rules.foldsets]) ## excludes folds already discovered
        return self.sort_foldsets_by_size(remaining_foldsets)

    def sort_foldsets_by_size(self, foldsets):
        """
        Sorts and returns foldsets based on their size.

        Arguments:
            foldsets (set): A set of foldsets to sort.

        Returns:
            A dict of foldsets keyed by their size.
        """
        rule_sizes = sorted(set([len(i) for i in foldsets]) - set([0]))
        size2foldsets = {size:list() for size in rule_sizes}

        foldset_tuples = sorted([sorted(tuple(i)) for i in foldsets]) ## cast as tuples for predictable sorting
        for fs in foldset_tuples: 
            size2foldsets[len(fs)].append(frozenset(fs)) ## change back to frozenset

        return size2foldsets

    def loop_through_remaining_foldsets_no_look_ahead(self, size2foldsets, current, key_to_maximize, debug=False, ordered_outcome=False, ignore_reaction_versions=False):
        """
        Loop through remaining foldsets with no look-ahead, returning the foldset(s) that maximize the given key.
        
        Arguments:
            size2foldsets (dict): A dictionary with foldset sizes as keys and sets of foldsets as values.
            current (Params): A Params object representing the current state of the fold expansion.
            key_to_maximize (str): A string specifying the key to maximize ("rns", "cpds", or "rules").
            debug (bool, optional): Whether to print debug information (default: False).
            ordered_outcome (bool, optional): Whether to return the foldsets in an ordered list (default: False).
        
        Returns:
            List of the foldset(s) that maximize the given key.
        """

        ## key_to_maximize is one of "rns", "cpds", "rules"
        if key_to_maximize=="folds":
            raise(ValueError("It doesn't make sense to choose a fold which maximizes number of folds."))

        one_step_effects = Params()
        one_step_effects.cpds, one_step_effects.rns = self.fold_expand(self.scope.folds, current.cpds, fold_algorithm="step")

        possible_next_rules = self.scope.rules.remaining_rules(current.folds).subset_from_rns(one_step_effects.rns)

        max_foldsets = list()
        max_foldset2key_counts = dict()
        for size in sorted(size2foldsets.keys()):
            
            foldset2key_count = dict() ## key_to_maximize
            for foldset in size2foldsets[size]:
                _foldset_rules = possible_next_rules.subset_from_folds(current.folds | foldset) # rules enabled after trialing the addition of a foldset

                if ignore_reaction_versions:
                    _foldset_rules = _foldset_rules.versionless

                foldset2key_count[foldset] = len(getattr(_foldset_rules, key_to_maximize))

            max_v = max(foldset2key_count.values()) 
            max_foldsets = [k for k, v in foldset2key_count.items() if v==max_v and max_v>0]
            max_foldset2key_counts = {k:v for k, v in foldset2key_count.items() if v==max_v and max_v>0}
            
            print("+++++++++++++++++")
            pprint(f"foldset2key_count: {foldset2key_count}")
            print(f"max_v: {max_v}")
            pprint(f"max_foldsets:\n\t{max_foldsets}")
            if max_v>0:
                break
        
        return max_foldset2key_counts

    def loop_through_remaining_foldsets_look_ahead(self, size2foldsets, current, key_to_maximize, debug=False, ignore_reaction_versions=False):
        """
        Loop through remaining foldsets with look-ahead, returning the foldset(s) that maximize the given key.
        
        Arguments:
            size2foldsets (dict): A dictionary with foldset sizes as keys and sets of foldsets as values.
            current (Params): A Params object representing the current state of the fold expansion.
            key_to_maximize (str): A string specifying the key to maximize ("rns", "cpds", or "rules").
            debug (bool, optional): Whether to print debug information (default: False).
        
        Returns:
            A dictionary with foldsets as keys and Params objects as values. Each foldset in the dictionary maximizes the given key_to_maximize. The Params objects represent the effects of adding the corresponding foldset and include keys for "folds", "cpds", "rns", and "rules".
        """
        
        ## key_to_maximize is one of "rns", "cpds", "rules"
        if key_to_maximize=="folds":
            raise(ValueError("It doesn't make sense to choose a fold which maximizes number of folds."))

        max_effects = dict()
        max_v = 0
        for size in sorted(size2foldsets.keys()):
            
            foldsets = size2foldsets[size]
            
            for foldset in foldsets:

                effects = Params()

                effects.folds = current.folds | set(foldset)
                effects.cpd_iteration_dict, effects.rn_iteration_dict = self.fold_expand(effects.folds, current.cpds)
                effects.cpds, effects.rns = set(effects.cpd_iteration_dict.keys()), set(effects.rn_iteration_dict.keys())
                effects.rules = self.f.subset_from_folds(effects.folds).subset_from_rns(effects.rns) ## this could include many unreachable rules because we never restricted ourselves to the present folds!

                if key_to_maximize == "rns" and ignore_reaction_versions:
                    n_new_set = len(get_versionless_reactions(effects.rns) - get_versionless_reactions(current.rns))
                elif key_to_maximize == "rules" and ignore_reaction_versions:
                    n_new_set = len(effects.rules.versionless - current.rules.versionless)
                else:
                    n_new_set = len(set(getattr(effects, key_to_maximize)) - set(getattr(current, key_to_maximize)))

                if debug:
                    print("size: ", size)
                    print("foldset: ", foldset)
                    print("len_effects / len_current: ", len(getattr(effects, key_to_maximize)),  len(getattr(current, key_to_maximize)))
                    print("set_effects / set_current: ", len(set(getattr(effects, key_to_maximize))),  len(set(getattr(current, key_to_maximize))))
                    print("set_effects - set_current: ", [i.id for i in set(getattr(effects, key_to_maximize)) - set(getattr(current, key_to_maximize))])
                    print("set_current - set_effects: ", [i.id for i in set(getattr(current, key_to_maximize)) - set(getattr(effects, key_to_maximize))])
                    print("fold_current - fold_effects: ", set(current.folds) - set(effects.rules.folds))
                    print("rn_current - rn_effects: ", set(current.rns) - set(effects.rns))
                    print("cpd_current - cpd_effects: ", set(current.cpds) - set(effects.cpds))
                    print("max_v: ", max_v)
                    print("max_v foldsets: ", max_effects.keys())
                    print("n_new_set: ", n_new_set)
                    print("~"*40)

                if (n_new_set == max_v) and (max_v > 0):
                    max_effects[foldset] = effects
                elif n_new_set > max_v:
                    max_v = n_new_set
                    max_effects = dict()
                    max_effects[foldset] = effects
                else: # n_new < max_v or n_new == max_v == 0
                    pass

            ## Don't look for longer rules if shorter rules enable new reactions
            if max_v>0:
                break
        return max_effects

    def choose_next_foldset_no_look_ahead(self, current, max_foldset2key_counts, ordered_outcome=False):
        """
        Given the current foldset, choose the next foldset to expand using the no-look-ahead algorithm.

        Args:
            current (Params): The current state of the fold expansion.
            max_foldset2key_counts (list): List of frozenset objects representing the maximum effect foldsets to consider for expansion.
            ordered_outcome (bool): Whether to select the next foldset deterministically or randomly. (default:False)

        Returns:
            Tuple containing:
            - Next foldset to expand represented as a frozenset object.
            - Dictionary containing the effects of the expansion on the model, where each key is a frozenset object representing a foldset, and the corresponding value is a Params object representing the updated model state.
        """
        if len(max_foldset2key_counts)>0:
            foldset_tuples = sorted([sorted(tuple(i)) for i in max_foldset2key_counts]) ## cast as tuples for predictable sorting
            if ordered_outcome:
                next_foldset = frozenset(foldset_tuples[0])
            else:
                next_foldset = frozenset(random.choice(foldset_tuples)) ## change back to frozenset

            ## Do expansion
            effects = Params()
            effects.folds = current.folds | set(next_foldset)
            effects.cpd_iteration_dict, effects.rn_iteration_dict = self.fold_expand(effects.folds, current.cpds)
            effects.cpds, effects.rns = set(effects.cpd_iteration_dict.keys()), set(effects.rn_iteration_dict.keys())
            effects.rules = self.f.subset_from_folds(effects.folds).subset_from_rns(effects.rns) ## this could include many unreachable rules because we never restricted ourselves to the present folds!
            
            max_foldset2key_counts[next_foldset] = effects
            return next_foldset, max_foldset2key_counts ## to mimic the structure of max_effects

        else:
            print("No foldsets remaining.")
            return frozenset(), {frozenset():deepcopy(current)}

    def choose_next_foldset_look_ahead(self, current, max_effects, ordered_outcome=False):
        """
        Given the current foldset, choose the next foldset to expand using the look-ahead algorithm.

        Args:
            current (Params): The current state of the fold expansion.
            max_effects (dict): Dictionary where each key is a frozenset object representing a foldset, and the corresponding value is a Params object representing the maximum effects of expanding that foldset.
            ordered_outcome (bool): Whether to select the next foldset deterministically or randomly. (default:False)

        Returns:
            Tuple containing:
            - Next foldset to expand represented as a frozenset object.
            - Dictionary containing the effects of the expansion on the model, where each key is a frozenset object representing a foldset, and the corresponding value is a Params object representing the updated model state.
        """
        if len(max_effects) == 0:
            next_foldset = frozenset()
            max_effects[next_foldset] = deepcopy(current)
            print("No foldsets remaining.")
        else:
            foldset_tuples = sorted([sorted(tuple(i)) for i in max_effects.keys()]) ## cast as tuples for predictable sorting
            if ordered_outcome:
                next_foldset = frozenset(foldset_tuples[0])
            else:
                next_foldset = frozenset(random.choice(foldset_tuples)) ## change back to frozenset
        return next_foldset, max_effects #[next_foldset]

    def choose_next_foldset_random(self, current):
        remaining_folds = set(self.scope.folds) - set(current.folds)
        print(len(remaining_folds))

        if len(remaining_folds) == 0:
            print("No foldsets remaining.")
            return frozenset(), {frozenset():deepcopy(current)}
            
        else:
            next_foldset = frozenset([random.choice(sorted(remaining_folds))]) ## this will be a single fold; can't sample from set
            ## Do expansion
            effects = Params()
            effects.folds = current.folds | set(next_foldset)
            effects.cpd_iteration_dict, effects.rn_iteration_dict = self.fold_expand(effects.folds, current.cpds)
            effects.cpds, effects.rns = set(effects.cpd_iteration_dict.keys()), set(effects.rn_iteration_dict.keys())
            effects.rules = self.f.subset_from_folds(effects.folds).subset_from_rns(effects.rns)

            return next_foldset, {next_foldset:effects} ## to mimic the structure of max_effects

    def choose_next_foldset(self, algorithm, size2foldsets, current, debug=False, ordered_outcome=False, ignore_reaction_versions=False):
        """
        Given the current foldset, choose the next foldset to expand using the specified algorithm.

        Arguments:
            algorithm (str): The algorithm to use for choosing the next foldset.
            size2foldsets (dict): A dict where each key is a tuple of size values, and the corresponding value is a list of frozenset objects representing the foldsets with those size values.
            current (Params): The current state of the fold expansion.
            debug (bool, optional): Whether to print debug information. (defualt: False)
            ordered_outcome (bool, optional): Whether to select the next foldset deterministically or randomly. (default: False)

        Returns:
            Tuple containing:
            - Next foldset to expand represented as a frozenset.
            - Dictionary containing the effects of the expansion on the model, where each key is a frozenset object representing a foldset, and the corresponding value is a Params object representing the updated model state.
        """

        look_ahead_algorithms = {
            "look_ahead_rules":"rules",
            "look_ahead_rns":"rns",
            "look_ahead_cpds":"cpds",
        }

        no_look_ahead_algorithms = {
            "no_look_ahead_rules":"rules",
            "no_look_ahead_rns":"rns"
        }
        
        if algorithm in look_ahead_algorithms:
            max_effects = self.loop_through_remaining_foldsets_look_ahead(size2foldsets, current, look_ahead_algorithms[algorithm], debug=debug, ignore_reaction_versions=ignore_reaction_versions)
            next_foldset, max_effects = self.choose_next_foldset_look_ahead(current, max_effects, ordered_outcome)

        elif algorithm in no_look_ahead_algorithms:
            max_foldset2key_counts = self.loop_through_remaining_foldsets_no_look_ahead(size2foldsets, current, no_look_ahead_algorithms[algorithm], debug=debug, ordered_outcome=ordered_outcome, ignore_reaction_versions=ignore_reaction_versions)
            next_foldset, max_effects = self.choose_next_foldset_no_look_ahead(current, max_foldset2key_counts, ordered_outcome)

        elif algorithm=="random_fold_order":
            ## loop_through function not needed in this case
            ## size2foldsets; ordered_outcome also unused
            next_foldset, max_effects = self.choose_next_foldset_random(current)

        else:
            raise(ValueError("algorithm not found."))

        return next_foldset, max_effects

    def keep_going(self, current, algorithm):
        """
        Determines whether the fold expansion should stop based on current state of the expansion and the specified algorithm.

        Arguments:
            current (Params):  The current state of the fold expansion.
            algorithm (str): The algorithm to use for choosing the next foldset.

        Returns:
            True if the expansion should continue.
        """

        if algorithm in ["look_ahead_rules", "no_look_ahead_rules"]:
            if set(self.scope.folds).issubset(set(current.folds)):
                print("Reached scope folds.")
                return False

        elif algorithm in ["look_ahead_rns", "no_look_ahead_rns"]:
            if set(self.scope.rns).issubset(set(current.rns)):
                print("Reached scope reactions.")
                return False
        
        elif algorithm in ["look_ahead_cpds", "no_look_ahead_cpds"]:
            if set(self.scope.cpds).issubset(set(current.cpds)):
                print("Reached scope compounds.")
                return False

        elif algorithm in ["random_fold_order"]:
            if set(self.scope.folds).issubset(set(current.folds)):
                print("Reached scope folds.")
                return False

        else:
            return True

    def rule_order(self, algorithm, write=False, write_tmp=False, path=None, str_to_append_to_fname=None, debug=False, ordered_outcome=False, ignore_reaction_versions=False):
        """
        Determine the ordering of all rules/folds.

        Arguments:
            algorithm (str): The algorithm to use for determining the rule/fold order.
            write (bool, optional): If True, write the results to a file. (default: False)
            path (str, optional): The path to the directory where the results file should be written. (default: None)
            str_to_append_to_fname (str, optional): A string to append to the end of the results file name. (default: None)
            debug (bool, optional): If True, print debug information. (default: False)
            ordered_outcome (bool, optional): Whether to select the next foldset deterministically or randomly. (default: False)

        Returns:
            A `Result` object that stores the results of the rule/fold ordering algorithm.
        """
        ## Place to store results and current state of expansion
        ## ITERATION 0 (Avoid updating folds on the 0th iteration since they don't apply until iteration=1)
        result = Result(self.scope)
        current = Params(folds=self.seed.folds, cpds=self.seed.cpds, rns=set(), rules=self.scope.rules.subset_from_folds(self.seed.folds).subset_from_rns(self.seed.rns))
        current.cpd_iteration_dict = {k:0 for k in current.cpds}
        current.rn_iteration_dict = {}
        metadata = Metadata()
        result.first_update(current, write=write_tmp, path=path, str_to_append_to_fname=str_to_append_to_fname)

        ## ITERATION 1 (using only seed folds and fold independent reactions)
        current.cpd_iteration_dict, current.rn_iteration_dict = self.fold_expand(current.folds, current.cpds)
        current.cpds, current.rns = set(current.cpd_iteration_dict.keys()), set(current.rn_iteration_dict.keys())
        current.rules = self.scope.rules.subset_from_folds(current.folds).subset_from_rns(current.rns)
        result.update(current, write=write_tmp, path=path, str_to_append_to_fname=str_to_append_to_fname)

        ## Needed in case expansion not possible at all
        keep_going = self.keep_going(algorithm, current)

        ################################################
        ## ITERATION 2+
        while keep_going:
            size2foldsets = self.sort_remaining_foldsets_by_size(current.folds)
            next_foldset, max_effects = self.choose_next_foldset(algorithm, size2foldsets, current, debug=debug, ordered_outcome=ordered_outcome, ignore_reaction_versions=ignore_reaction_versions)
            effects = max_effects[next_foldset]

            if len(next_foldset)==0:
                keep_going = False 

            if keep_going:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                current.folds = (current.folds | set(next_foldset))
                current.cpds = effects.cpds
                current.cpd_iteration_dict = effects.cpd_iteration_dict
                current.rns = effects.rns
                current.rn_iteration_dict = effects.rn_iteration_dict
                current.rules = self.scope.rules.subset_from_folds(current.folds).subset_from_rns(current.rns)
                metadata.max_effects = max_effects
                metadata.size2foldsets = size2foldsets

            ## Store when cpds and rns appear in the expansion
            result.update(current, metadata=metadata, write=write_tmp, path=path, str_to_append_to_fname=str_to_append_to_fname)
            print("="*60)
            print("rule iter: {} ({:.2} sec) {}".format(result.iteration, result.iteration_time[result.iteration], next_foldset))
            print("="*60)

        if write:
            result.final_write(path=path, str_to_append_to_fname=str_to_append_to_fname)

        return result

def example_main():
    asset_path = nf.asset_path

    ALGORITHM = "no_look_ahead_rules"
    WRITE = True # write result to disk
    WRITE_TMP = True # write after each iteration
    CUSTOM_WRITE_PATH = None # if writing result, custom path to write to
    STR_TO_APPEND_TO_FNAME = "EXAMPLE" # if writing result, str to append to filename
    ORDERED_OUTCOME = False # ignore random seed and always choose folds based on sort order
    IGNORE_REACTION_VERSIONS = True # when maximizing for reactions, don't count versioned reactions
    METABOLISM_PATH = PurePath(asset_path, "metabolic_networks","metabolism.23Aug2022.pkl") # path to metabolism object pickle
    RN2RULES_PATH = PurePath(asset_path,"rn2fold","rn2rules.20230224.pkl") # path to rn2rules object pickle
    SEED_CPDS_PATH = PurePath(asset_path, "compounds", "seeds.Goldford2022.csv") # path to seed compounds csv

    ## Metabolism
    metabolism = pd.read_pickle(METABOLISM_PATH)

    ## FoldRules
    rn2rules = pd.read_pickle(RN2RULES_PATH)
    foldrules = nf.FoldRules.from_rn2rules(rn2rules)

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
    