import pandas as pd
from pathlib import PurePath, Path
from copy import copy, deepcopy
import timeit
from pprint import pprint
from collections import Counter
from datetime import datetime
import random
import pickle

asset_path = PurePath(__file__).parent / "assets"

class ImmutableParams:
    """
    Stores parameters that commonly need to be kept track of together.

    After initialization attributes are immutable, to avoid unintentionally changing parameters like the seed or scope.
    """

    def __init__(self, folds=None, cpds=None, rns=None, rules=None):
        self._folds = deepcopy(folds)
        self._cpds = deepcopy(cpds)
        self._rns = deepcopy(rns)
        self._rules = deepcopy(rules)

    @property
    def folds(self):
        return self._folds

    @property
    def cpds(self):
        return self._cpds

    @property
    def rns(self):
        return self._rns

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

    @ImmutableParams.rns.setter
    def rns(self, value):
        self._rns = deepcopy(value)

    @ImmutableParams.rules.setter
    def rules(self, value):
        self._rules = deepcopy(value)

class Result:
    """
    Store data from the run
    """

    def __init__(self):
        self.iteration = 0
        self.runtime = dict()  
        self.cpds = dict()
        self.rns = dict()
        self.folds = {"fold_independent":0}
        self.rules = dict()
        self.start_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.start_time = timeit.default_timer()
        self.iteration_time = dict()

    def first_update(self, current, write=False, path=None, str_to_append_to_fname=None):
        self.update_cpds(current)
        self.update_rns(current)
        self.update_iter()
        self.update_iteration_time()
        if write==True:
            temp_write(self, path=path, str_to_append_to_fname=str_to_append_to_fname)

    def update(self, current, write=False, path=None, str_to_append_to_fname=None):
        self.update_folds(current)
        self.update_rules(current)
        self.first_update(current, write=write, path=path, str_to_append_to_fname=str_to_append_to_fname)

    def update_cpds(self, current):
        for i in current.cpds:
            if i not in self.cpds:
                self.cpds[i] = self.iteration

    def update_rns(self, current):
        for i in current.rns:
            if i not in self.rns:
                self.rns[i] = self.iteration

    def update_folds(self, current):
        for i in current.folds:
            if i not in self.folds:
                self.folds[i] = self.iteration

    def update_rules(self, current):
        for i in current.rules.ids:
            if i not in self.rules:
                self.rules[i] = self.iteration

    def update_iteration_time(self):
        self.iteration_time[self.iteration] =  timeit.default_timer() - self.start_time

    def update_iter(self):
        self.iteration+=1

    def get_path(self, path=None, str_to_append_to_fname=None):
        if str_to_append_to_fname == None:
            fname = self.start_datetime+".pkl"
        else:
            fname = self.start_datetime+"_"+str_to_append_to_fname+".pkl"
        
        if path==None:
            path = Path.joinpath(Path.cwd(), "fold_results", fname)
        else:
            path = Path.joinpath(path, "fold_results", fname)

        return path

    def temp_write(self, path=None, str_to_append_to_fname=None):
        ## Doesn't check if overwriting

        path = self.get_path(path, str_to_append_to_fname)
        path = Path.joinpath(path.parent, path.stem+"_tmp", path.suffix)
        with open(path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def final_write(self, path=None, str_to_append_to_fname=None):
        
        path = self.get_path(path, str_to_append_to_fname)
        
        i = 0
        while path.is_file():
            path = Path.joinpath(path.parent, path.stem+"_"+str(i), path.suffix)
            i+=1

        with open(path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

class Rule:

    def __init__(self,rn,foldset):
        self.rn = rn
        self.foldset = foldset
        self.id = (rn, foldset)

    def __repr__(self):
        return "id:\t\t{0} \nrn:\t\t{1} \nfoldset:\t{2}".format(self.id, self.rn, self.foldset)

class FoldRules:
    """
    Object that stores Rules. This is the universe of possible fold rules. 
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

    def subset_from_rns(self, rns):
        return FoldRules([r for r in self.rules if r.rn in rns])
        # return FoldRules([self.rn2rule[r] for r in rns])

    def subset_from_folds(self, folds):
        return FoldRules([r for r in self.rules if r.foldset<=folds])
        # return Rules([self.fs2rule[r] for r in rns])

    def remaining_rules(self, current_folds):
        return FoldRules([r for r in self.rules if len(r.foldset-current_folds)>0])

    def to_list(self):
        return [r for r in self.rules]

    def __len__(self):
        return len(self.ids)

class FoldMetabolism:
    """
    A class to do fold expansion from a metabolism, foldrules, and seeds.
    """

    def __init__(self, metabolism, foldrules, seed):#, preexpansion=False):
        """
        Calculates expansion scope after seed compounds are defined. 

        :param metabolism: a GlobalMetabolicNetwork object that defines the compound/reaction rules of network expansion
        :param foldrules: a FoldRules object that defines the fold rules of fold expansion
        :param seed: a Params object that defines the initial compounds, folds, and reactions (these are fold independent reactions)
        :kwarg preexpansion: NOT YET IMPLEMENTED -- but the idea is to allow n steps of regular network expansion before beginning fold expansion
        """
        
        self._m = metabolism ## GlobalMetabolicNetwork object
        self._f = foldrules # FoldRules object
        self._seed = ImmutableParams(folds=seed.folds, rns=seed.rns, cpds=seed.cpds) ## seed.rns == fold_independent_rns
        self._scope = self.calculate_scope(seed)

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
        Calculate the scope of the seed with all reactions enabled by the global fold network (including fold independent reactions i.e. seed.rns)

        :param seed: an ImmutableParams object with seed.rns and seed.cpds != None
        :return: an ImmutableParams object specificying the scope
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

    def fold_expand(self, folds, current_cpds):
        """
        Returns a set of compounds and set of reactions enabled by expansion 
        using `current_cpds` and the reactions enabled by `folds` and fold independent 
        reactions (seed.rns).

        :param current_cpds: collection of compound ids
        :param folds: collection of fold ids
        :return: set of compounds and set of reactions from the expansion
        """

        active_rules = self.f.subset_from_folds(folds)
        rn_tup_set = set(self.m.rxns2tuple(active_rules.rns | self.seed.rns))
        cx,rx = self.m.expand(current_cpds, reaction_mask=rn_tup_set)
        return set(cx), set([i[0] for i in rx])

    def sort_remaining_foldsets_by_size(self, current_folds):#, remaining_rules):

        remaining_rules = self.scope.rules.remaining_rules(current_folds)
        remaining_foldsets = set([i-current_folds for i in remaining_rules.foldsets]) ## excludes folds already discovered
        rule_sizes = sorted(set([len(i) for i in remaining_foldsets]))
        ## Organizes rules by size
        size2foldsets = {size:list() for size in rule_sizes}

        remaining_foldset_tuples = sorted([sorted(tuple(i)) for i in remaining_foldsets]) ## cast as tuples for predictable sorting
        for fs in remaining_foldset_tuples: ## sorted for reproduceability
            size2foldsets[len(fs)].append(frozenset(fs))

        return size2foldsets

    def loop_through_remaining_foldsets(self, size2foldsets, current, key_to_maximize):
        ## key_to_maximize is one of "rns", "cpds", "rule"

        max_effects = dict()
        max_v = 0
        for size, foldsets in size2foldsets.items():

            for foldset in foldsets:

                effects = Params()

                effects.cpds, effects.rns = self.fold_expand((current.folds | set(foldset)), current.cpds)
                effects.rules = self.f.subset_from_rns(effects.rns)

                n_new = len(getattr(effects, key_to_maximize)) - len(getattr(current, key_to_maximize))

                if n_new == max_v:
                    max_effects[foldset] = effects
                elif n_new > max_v:
                    max_v = n_new
                    max_effects = dict()
                    max_effects[foldset] = effects
                else: # n_new < max_v
                    pass

            ## Don't look for longer rules if shorter rules enable new reactions
            if len(max_effects)>0:
                break

        return max_effects

    def select_next_foldset(self, algorithm, size2foldsets, current):
        
        if algorithm == "max_rules":
            max_effects = self.loop_through_remaining_foldsets(size2foldsets, current, "rules")
            if len(max_effects) == 0:
                next_foldset = frozenset()
                max_effects[next_foldset] = deepcopy(current)
                print("NO max_effects REMAINING")
            else:
                next_foldset = random.choice(sorted(max_effects.keys()))
            return next_foldset, max_effects[next_foldset]

    def keep_going(self, current):

        ## It should always stop if we've found all scope cpds, rns
        ## actually it might be possible to discover all reactions and compounds but not all rules
        ## but if we've discovered all folds then we've discovered all rules
        if (current.cpds == self.scope.cpds) and (current.rns == self.scope.rns) and (current.folds == self.scope.folds):
            return False

        else:
            return True

    def rule_order(self, algorithm="max_rules", write=False, path=None, str_to_append_to_fname=None):
        """
        Determine the ordering of all rules/folds.
        """
        ## Place to store results and current state of expansion
        ## ITERATION 0 (Avoid updating folds on the 0th iteration since they don't apply until iteration=1)
        result = Result()
        current = Params(folds=self.seed.folds, cpds=self.seed.cpds, rns=self.seed.rns, rules=self.scope.rules.subset_from_folds(self.seed.folds))
        result.first_update(current, write=write, path=path, str_to_append_to_fname=str_to_append_to_fname)

        ## ITERATION 1 (using only seed folds and fold independent reactions)
        current.cpds, current.rns = self.fold_expand(current.folds, current.cpds)
        result.update(current, write=write, path=path, str_to_append_to_fname=str_to_append_to_fname)

        ## Needed in case expansion not possible at all
        keep_going = self.keep_going(current)

        ################################################
        ## ITERATION 2+
        while keep_going:
            print("rule iter: {} ({:.2} sec)".format(result.iteration, result.iteration_time[result.iteration]))
            size2foldsets = self.sort_remaining_foldsets_by_size(current.folds)
            next_foldset, effects = self.select_next_foldset(algorithm, size2foldsets, current)

            keep_going = self.keep_going(current)

            if keep_going:
                ## Update folds, rules2rns available; Update rns in expansion, cpds in expansion
                current.folds = (current.folds | set(next_foldset))
                current.rules = self.scope.rules.subset_from_folds(current.folds)
                current.cpds = effects.cpds
                current.rns = effects.rns

            ## Store when cpds and rns appear in the expansion
            result.update(current, write=write, path=path, str_to_append_to_fname=str_to_append_to_fname)

        if write:
            result.final_write(path=path, str_to_append_to_fname=str_to_append_to_fname)

        return result

    # def select_next_rule_or_fold(self, current_folds, current_cpds, current_rns, remaining_folds, remaining_rules, algorithm="maxreactionsupersets"):#rselect_func=maxreactions):
    #     """
    #     Determines the next rule and its effects on cpds, rns, and rules, based on the selection function.

    #     :param current_folds: collection of current folds
    #     :param current_cpds: collection of current compounds
    #     :param current_rns: collection of current reactions
    #     :param rselect_func: function to use for identifying the next rule
    #     :return:- next_rule and ...
    #             - ... dictionary of its effects {rule2rn:... , cpds:... , rns:...} 
    #             - and for metadata purposes, also return the number of rules checked
    #             - and the number of equal rule groups this iteration
    #     """
    #     if algorithm == "randomfold":
    #         next_fold = random.choice(remaining_folds)
    #         _fdict = dict()
    #         _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(next_fold, current_folds, current_cpds)
    #         return next_fold, _fdict, 1, 0

    #     elif algorithm == "randomrule":
    #         rule_size_dict = rule_sizes(remaining_rules)
    #         next_rule = random.choice(rule_size_dict[min(rule_size_dict)])
    #         _fdict = dict()
    #         _fdict["rule2rns"], _fdict["cpds"], _fdict["rns"] = self.effect_per_rule_or_fold(next_fold, current_folds, current_cpds)
    #         return next_rule, _fdict, 1, 0
            
    #     elif algorithm == "maxreactionsupersets":
    #         r_effects, n_rules_checked, n_rules_skipped = self.loop_through_rules(current_cpds, current_rns, current_folds, remaining_rules)
    #         if len(r_effects) == 0:
    #             next_rule = frozenset()
    #             r_effects[next_rule] = {"cpds":deepcopy(current_cpds), "rns":deepcopy(current_rns)}
    #             print("NO R_EFFECTS REMAINING")
    #         else:
    #             next_rule = random.choice(sorted(r_effects.keys()))
    #         return next_rule, r_effects[next_rule], n_rules_checked, n_rules_skipped #, er_effects

    #     ## Assuming you start with all scope compounds in the seed set, pick fold order based on fold which enables the most reactions
    #     elif algorithm == "allcompoundsseedmaxreactions":
    #         # for r in remaining_rules:
    #         #     potential_fold_set = (current_folds | set(rule))
    #         #     potential_rule2rns = subset_rule2rn_from_folds(potential_fold_set, self.scope_rules2rn)

    #         future_rule2rns = {k:(v | current_rns) for k,v in remaining_rules.items() if len(v-current_rns)>0}
    #         rule_sizes = sorted(set([len(i) for i in future_rule2rns]))

    #         ## Organizes rules by size
    #         future_rule2rns_by_size = {size:dict() for size in rule_sizes}
    #         for rule, rns in future_rule2rns.items():
    #             future_rule2rns_by_size[len(rule)][rule]=rns

    #         for size in rule_sizes:
    #             future_rule2rns_by_size_sorted = sorted(future_rule2rns_by_size[size], key=lambda k: len(future_rule2rns_by_size[size][k]), reverse=True)
    #             # future_rule2rns[rule] for rule in future_rule2rns_by_size_sorted:
    #                 # rns = future_rule2rns[rule]


    #         ## Don't look for longer rules if shorter rules enable new reactions
    #         # if len(max_r_effects)>0:
    #         #     break

    #     ## picks folds maximizing the number of rules
    #     elif algorithm == "allcompoundsseedmaxrules":

    #         rules_minus_current_folds = {k-current_folds:v for k,v in remaining_rules.items()}
    #         rule_sizes = sorted(set([len(i) for i in rules_minus_current_folds]))
    #         # print(rule_sizes)
    #         ## Do we still want rule to be added even if it only maps to reactions which are already in the network?
    #         ## case 1: want to find rule regardless of it mapping to reactions already in the network

    #         ## case 2: want to find rule only if it maps to reactions not in the network
    #         # pprint(rules_minus_current_folds)

    #         ## pick fold that immeasiately enables the most rules 
    #         fold2activatedRuleCounter = Counter([f for k in rules_minus_current_folds for f in k if len(k)==min(rule_sizes)])
    #         fold2activatedRuleMaxValue = max(fold2activatedRuleCounter.values())
    #         fold2activatedRuleMaxFold = [k for k,v in fold2activatedRuleCounter.items() if v==fold2activatedRuleMaxValue]
    #         print(fold2activatedRuleMaxFold)
    #         next_rule = random.choice(sorted(r_effects.keys()))

    #         ## need to make a conversion between fold2rule and rule2fold

    #     ## Need to choose an option that picks folds which maximize the number of rules
    #     # elif algorithm == "maxrules":
    #     #     blah blah blah

    #     # elif algorithm == "maxreactions":
    #     #     r_effects, n_rules_checked, n_equal_rule_groups = self.loop_through_rules(current_folds, current_cpds, current_rns)

def example_main():
    asset_path = nf.asset_path

    ALGORITHM = "max_rules"
    WRITE = False
    PATH = None
    STR_TO_APPEND_TO_FNAME = "EXAMPLE"

    ## Metabolism
    metabolism_path = PurePath(asset_path) / "metabolic_networks" / 'metabolism.23Aug2022.pkl'
    metabolism = pd.read_pickle(metabolism_path)

    ## FoldRules
    rn2rules = pd.read_pickle(PurePath("data","rn2fold","rn2rules_v.pkl"))
    foldrules = nf.FoldRules.from_rn2rules(rn2rules)

    ## Seed
    fold_independent_rns =  set(metabolism.network["rn"]) - set(rn2rules)
    seed_cpds_path = PurePath("data", "josh", "seed_set.csv")
    seed_cpds = set((pd.read_csv(seed_cpds_path)["ID"])) #| aa_cids
    seed = nf.Params(
        rns = fold_independent_rns,
        cpds = seed_cpds,
        folds = set(['spontaneous'])
    )

    ## Inititalize fold metabolism
    fm = nf.FoldMetabolism(metabolism, foldrules, seed)
    ## Run fold expansion
    result = fm.rule_order(algorithm=ALGORITHM, write=WRITE, path=PATH, str_to_append_to_fname=STR_TO_APPEND_TO_FNAME)
    