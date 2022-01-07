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
	try:
		reaction_cc = cc.parse_reaction_formula(rxn_string)
		dG0_prime = cc.standard_dg_prime(reaction_cc)
		return dG0_prime.value.m_as("kJ/mol")
	except:
		return np.nan

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


def substituteCoenzyme(metabolism,coenzyme_substitution_df,label,pH=7.0,pMg=3.0,temperature=298.15,ionic_strength=0.25,predefined_reaction_list=None):
	cc = initialze_cc(pH=pH,pMg=pMg,temperature=temperature,ionic_strength=ionic_strength)
	#metabolism is a metabolism object from the networkExpansionPy package
	#coenzyme_substitution_df is a table with the follwing columns: cid,s,type, where type is either o or m (original vs modified)

	new_rxns = []
	thermo_new = {'rn': [], 'dg': []}
	new_rxns = []
	cj = coenzyme_substitution_df[coenzyme_substitution_df.type == 'o'];
	cj.columns = [x + '_o' for x in cj.columns]
	reaction_network = metabolism.network.copy()

	# if you specificy only a subset of reactions to consider, trim down the base network 
	if predefined_reaction_list is not None:
		reaction_network = reaction_network[reaction_network.rn.isin(predefined_reaction_list)]

	for rn,dff in reaction_network.groupby('rn'):
		# determine if reaction follows stoich rules
		dff_c = dff[dff.cid.isin(cj.cid_o)]
		if len(dff_c) == len(cj):
			dff_cids_s = dff_c.set_index('cid')[['s']]
			dff_cids_s.columns = ['s']
			cdf_s = cj.set_index('cid_o')[['s_o']]
			x = dff_cids_s.join(cdf_s)
			# compute eigenvalues of stoichiometrhy matrix. If there is more than 1 eigenvalue, then column are not indepdnent and the reaction does not folllow the stoichiometry rules
			#lambdas, V =  np.linalg.eig(x.values)
			#g = sum(lambdas**2 > 1e-30)

			g = np.linalg.matrix_rank(x.values)
			if g == 1:

				# partition data frame into coenzyme part and coenzyme indpendent part
				dff_c = dff[dff.cid.isin(cj.cid_o.tolist())]
				dff_i = dff[~dff.cid.isin(cj.cid_o.tolist())]

				coenzyme_rxn_df = dff_c.set_index('cid').join(cj.set_index('cid_o'))
				factor = coenzyme_rxn_df['s'] / coenzyme_rxn_df['s_o']
				# add step here to ensure factor is the same for each original coenzymes
				factor = factor.iloc[0]
				new_coenzymes = coenzyme_substitution_df[coenzyme_substitution_df.type == 'm']
				co_new = coenzyme_substitution_df[coenzyme_substitution_df.type == 'm']
				#co_new['s'] = co_new['s'].apply(lambda x: x*factor)
				#co_new['s'] = co_new['s'] * factor
				co_new = co_new.assign(s=co_new['s']*factor)
				dff_rn = pd.concat([dff_i[['cid','s']],co_new[['cid','s']]],axis=0)
				new_rid = rn + '_' + label
				dff_rn['rn'] = new_rid
				dff_rn = dff_rn[['cid','rn','s']]
				new_rxns.append(dff_rn)
				dg = computeFreeEnergy(dff_rn,cc)
				thermo_new['rn'].append(new_rid)
				thermo_new['dg'].append(dg)
			else:
				print('reaction {r} not substituable becase not rank coenzyme matrix is not rank 1'.format(r=rn))


	thermo_new = pd.DataFrame(thermo_new)
	new_rxns = pd.concat(new_rxns).drop('stoich_str',axis=1)
	thermo_new.columns = ['!MiriamID::urn:miriam:kegg.reaction', '!dG0_prime (kJ/mol)']
	metabolism.thermo = pd.concat([metabolism.thermo,thermo_new],axis=0)
	metabolism.network = pd.concat([metabolism.network,new_rxns],axis=0)
	return metabolism



def computeThermodynamics(metabolism,pH=7.0,pMg=3.0,temperature=298.15,ionic_strength=0.25):
	cc = initialze_cc(pH=pH,pMg=pMg,temperature=temperature,ionic_strength=ionic_strength)
	#metabolism is a metabolism object from the networkExpansionPy package
	#coenzyme_substitution_df is a table with the follwing columns: cid,s,type, where type is either o or m (original vs modified)
	thermo_new = {'rn': [], 'dg': []}
	for rn,dff in metabolism.network.groupby('rn'):
		dg = computeFreeEnergy(dff,cc)
		thermo_new['rn'].append(rn)
		thermo_new['dg'].append(dg)

	thermo_new = pd.DataFrame(thermo_new)
	thermo_new.columns = ['!MiriamID::urn:miriam:kegg.reaction', '!dG0_prime (kJ/mol)']
	return thermo_new