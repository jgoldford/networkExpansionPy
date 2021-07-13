import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_


def returnStoichStr(x):
    if np.abs(x) == 1:
        return ''
    else:
        return str(abs(x))
    
def computeFreeEnergy(rxn_df,cc):
    # contert rxn dataframe to reaction string
    rxn_df['stoich_str'] = rxn_df['s'].apply(returnStoichStr)
    S = rxn_df[rxn_df['s'] < 0]
    P = rxn_df[rxn_df['s'] > 0]
    rxn_string =  " + ".join((S.stoich_str +  ' kegg:' + S.cid).tolist()) + ' = ' +   " + ".join((P.stoich_str +  ' kegg:' + P.cid).tolist())
    reaction_cc = cc.parse_reaction_formula(rxn_string)
    dG0_prime = cc.standard_dg_prime(reaction_cc)
    return dG0_prime.value.m_as("kJ/mol")

def initialze_cc(pH=7.0,pMg=3.0,temperature=298.15,ionic_strength=0.25):
	cc = ComponentContribution()
	cc.p_h = Q_(pH)
	cc.p_mg = Q_(pMg)
	cc.ionic_strength = Q_(str(ionic_strength) + "M")
	cc.temperature = Q_(str(temperature) + "K")
	return cc


def isCoenzyme(charged,uncharged,rxnDf):
	# determine if charge and uncharged 
	isCoenzymeRxn = (all([y in rxnDf.cid.tolist() for y in [runcharged,charged]])) & (rxnDf[rxnDf.cid.isin([uncharged,charged])].s.sum() == 0)
	return isCoenzymeRxn

def replace_coenzymes(metabolism,coenzymes,pH=7.0,pMg=3.0,temperature=298.15,ionic_strength=0.25):

	cc = initialze_cc(pH=pH,pMg=pMg,temperature=temperature,ionic_strength=ionic_strength)
	#metabolism is a metabolism object from the networkExpansionPy package
	#coenzymes is a table with the follwing columns: uncharged, charged, subs_uncharged, subs_charged
	new_rxns = []
	thermo_new = {'rn': [], 'dg': []}
	for rn,dff in metabolism.network.groupby('rn'):
	    for idx,row in coenzymes.iterrows():
	        # determine if 
	        isCoenzymeRxn = (all([y in dff.cid.tolist() for y in [row.uncharged, row.charged]])) & (dff[dff.cid.isin([row.uncharged, row.charged])].s.sum() == 0)
	        if isCoenzymeRxn:
	            dff = dff.replace({row.uncharged: row.subs_uncharged,row.charged: row.subs_charged })
	            new_rid = rn + '_G' + str(idx)
	            dff = dff.replace({rn:new_rid})
	            # only rxn if no coenzyme subs are not substrates or products in reaction already
	            if len(dff.cid.unique()) == len(dff):
	                new_rxns.append(dff)
	                dg = computeFreeEnergy(dff,cc)
	                thermo_new['rn'].append(new_rid)
	                thermo_new['dg'].append(dg)

	thermo_new = pd.DataFrame(thermo_new)
	new_rxns = pd.concat(new_rxns).drop('stoich_str',axis=1)
	thermo_new.columns = ['!MiriamID::urn:miriam:kegg.reaction', '!dG0_prime (kJ/mol)']
	metabolism.thermo = pd.concat([metabolism.thermo,thermo_new],axis=0)
	metabolism.network = pd.concat([metabolism.network,new_rxns],axis=0)
	return metabolism


def substituteCoenzyme(metabolism,coenzyme_substitution_df,label,pH=7.0,pMg=3.0,temperature=298.15,ionic_strength=0.25):
	cc = initialze_cc(pH=pH,pMg=pMg,temperature=temperature,ionic_strength=ionic_strength)
	#metabolism is a metabolism object from the networkExpansionPy package
	#coenzyme_substitution_df is a table with the follwing columns: cid,s,type, where type is either o or m (original vs modified)
	new_rxns = []
	thermo_new = {'rn': [], 'dg': []}
	new_rxns = []
	cj = coenzyme_substitution_df[coenzyme_substitution_df.type == 'o'];
	cj.columns = [x + '_o' for x in cj.columns]

	for rn,dff in metabolism.network.groupby('rn'):
	    dff_c = dff[dff.cid.isin(cj.cid_o.tolist())]
	    dff_i = dff[~dff.cid.isin(cj.cid_o.tolist())]
	    if len(dff_c) == len(cj):
	        coenzyme_rxn_df = dff_c.set_index('cid').join(cj.set_index('cid_o'))
	        factor = coenzyme_rxn_df['s'] / coenzyme_rxn_df['s_o']
	        # add step here to ensure factor is the same for each original coenzymes
	        factor = factor.iloc[0]
	        new_coenzymes = coenzyme_substitution_df[coenzyme_substitution_df.type == 'm']
	        co_new = coenzyme_substitution_df[coenzyme_substitution_df.type == 'm']
	        co_new['s'] = co_new['s'].apply(lambda x: x*factor)
	        #co_new['s'] = co_new['s'] * factor
	        dff_rn = pd.concat([dff_i[['cid','s']],co_new[['cid','s']]],axis=0)
	        new_rid = rn + '_' + label
	       	dff_rn['rn'] = new_rid
	        dff_rn = dff_rn[['cid','rn','s']]
	        new_rxns.append(dff_rn)
	        dg = computeFreeEnergy(dff_rn,cc)
	        thermo_new['rn'].append(new_rid)
	        thermo_new['dg'].append(dg)

	thermo_new = pd.DataFrame(thermo_new)
	new_rxns = pd.concat(new_rxns).drop('stoich_str',axis=1)
	thermo_new.columns = ['!MiriamID::urn:miriam:kegg.reaction', '!dG0_prime (kJ/mol)']
	metabolism.thermo = pd.concat([metabolism.thermo,thermo_new],axis=0)
	metabolism.network = pd.concat([metabolism.network,new_rxns],axis=0)
	return metabolism