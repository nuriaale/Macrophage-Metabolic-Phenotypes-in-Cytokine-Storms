#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 11:05:47 2023

@author: nuria
"""

import sys
import os
import operator
ini_dir="/home/nuria/Escritorio/Macrophage-Metabolic_Phenotypes_in_Cytokine_Storms"
os.chdir(ini_dir)

import cobra
from cobra import Model, Reaction, Metabolite
import pandas as pd
import math
import copy
import numpy as np
from Block_2_functions import *
import input2 as inp2


from general_functions import *
from Block_2_functions import *



output_prefix=inp2.output_prefix
model_dict_string=inp2.model_dict
model_dict={}
for a in model_dict_string:
    model_dict[a]=cobra.io.read_sbml_model(model_dict_string[a])
conditions_to_sample=model_dict.keys()

gene_expression_file=inp2.gene_expression_file


gpr_mode=inp2.gpr_mode
or_mode=inp2.or_mode
convert_log2=inp2.convert_log2

absent_gene_expression=inp2.absent_gene_expression
absent_reaction_expression=inp2.absent_reaction_expression
correct_for_complexes=inp2.correct_for_complexes


medium_kpcs_file=inp2.medium_kpcs_file

intracelular_metabolomics_file=inp2.intracelular_metabolomics_file

max_NA_fraction_kpc=inp2.max_NA_fraction_kpc

max_NA_fraction_pellet=inp2.max_NA_fraction_pellet

low_expression_threshold=inp2.low_expression_threshold
base_penalty=inp2.base_penalty
penalty_precision=inp2.penalty_precision
fraction_of_optimum_biomass=inp2.fraction_of_optimum_biomass
force_run=inp2.force_run
remove_reactions_below=inp2.remove_reactions_below

gim3e_fraction_optimum=inp2.gim3e_fraction_optimum

n_samples=inp2.n_samples
thinning=inp2.thinning


metabolomics_file=inp2.metabolomics_file
gene_entrez_column=inp2.gene_entrez_column 

colum_replicates=inp2.colum_replicates
col_to_erase_l=inp2.col_to_erase_l

run_fva=inp2.run_fva
percentile=inp2.percentile
gene_method=inp2.gene_method
metabolite_list_fname=inp2.metabolite_list_fname
label_model=inp2.label_model
epsilon=inp2.epsilon
add_as_constraints=inp2.add_as_constraints
boundaries_precision=inp2.boundaries_precision
all_gene_exp_4_percentile=inp2.all_gene_exp_4_percentile
convert_solution_dict_to_reversible=inp2.convert_solution_dict_to_reversible
omit_0=inp2.omit_0

out_dir=inp2.out_dir
#%%
gene_expression_pd=pd.read_excel(gene_expression_file)

gene_expression_pd = gene_expression_pd.rename(columns={gene_entrez_column: 'gene'})
gene_expression_pd.dropna(subset=['gene'], inplace=True)
gene_expression_pd.set_index('gene', inplace=True)
gene_expression_pd.drop_duplicates(keep='first')#very importaaaaant

gene_expression_pd=gene_expression_pd[~gene_expression_pd[list(conditions_to_sample)[0]].isna()]

mask = gene_expression_pd.index.duplicated(keep='first')
gene_expression_pd = gene_expression_pd[~mask]
#This time the contraint file is none

if medium_kpcs_file is not None:
    kpc_pd=pd.read_excel(medium_kpcs_file,1)
    kpc_list=[kpc_pd]
    kpc_dict=general_calculations(kpc_list, colum_replicates, col_to_erase_l, conditions_to_sample, 0.50)


if metabolomics_file is not None:
    pellets_pd=pd.read_excel(metabolomics_file)

    pellets_dict=general_calculations([pellets_pd],colum_replicates, col_to_erase_l,conditions_to_sample, 0.50)
    for a in pellets_dict:
        pellets_dict[a].index=pellets_dict[a].index + '(pellet)'

    for a in kpc_dict:
        kpc_dict[a]=pd.concat([kpc_dict[a], pellets_dict[a]])
    

    
#%%
log2=convert_log2
reaction_list=[]
sample_sampling_dict={}
stat_dict={}
fva_dict={}
low_reaction_dict={}
   
os.chdir('outputs/2')
for sample in conditions_to_sample:
    condition=sample
    out_prefix=output_prefix+sample+""
    model_xlsx_filename=out_prefix+"_gimme_model_and_ranges.xlsx"
    restricted_model_fname=out_prefix+"_gimme__restricted_model.sbml"
    gene_expression_model_fname=out_prefix+"_gene_expression_model.sbml"
    model=model_dict[sample].copy()
    if metabolomics_file is not None:
        met_sink_dict, rejected_list=add_sink_reactions_to_model(model, pellets_pd, ['Unnamed: 1', 'pmol per cell per h'], biocrates_name_compartment_dict)
    
    if kpc_dict!={}:
        modify_model_for_seahorse(model,remove=False)
        x_dict,total_xi, xi_dict=mfa(model,kpc_dict,condition,metabolite_ex_dict,min_std=1e-4,precision=6,fname=output_prefix+sample+"MFA_fit.xlsx")
        for ex_rid in xi_dict:
            reaction_object=model.reactions.get_by_id(ex_rid)
            value=x_dict[ex_rid]
            ub=min(round(value+5e-7,7),reaction_object.upper_bound)
            lb=max(round(value-5e-7,7),reaction_object.lower_bound)
            reaction_object.bounds=(lb,ub)
        
    reaction_ids_to_omit=[x.id for x in model.reactions.query("EX_")]+[x.id for x in model.reactions.query("SINK_")]+[x.id for x in model.reactions.query("RGROUP_")]
    #gimme algorithm
    penalty_dict,gene_expression_model,objective,solution_dict,rev_fva, low_reactions, reaction_expression_dict=integrate_omics_gim3e_and_remove(model,gene_expression_pd,condition, fraction_of_optimum_biomass,low_expression_threshold,absent_gene_expression,absent_reaction_expression,percentile, gene_method,metabolite_list_fname,label_model,epsilon,gim3e_fraction_optimum,run_fva,add_as_constraints, boundaries_precision, all_gene_exp_4_percentile,base_penalty, convert_solution_dict_to_reversible, gpr_mode,or_mode,omit_0,log2,reaction_list=[], penalty_mode="normal", remove_reactions_below=remove_reactions_below,penalty_precision=penalty_precision,reactions_to_keep=[],correct_for_complexes=correct_for_complexes, reaction_ids_to_omit=reaction_ids_to_omit)

    fva_dict[sample]=rev_fva
    low_reaction_dict[sample]=low_reactions
    
    
    cobra.io.write_sbml_model(model,restricted_model_fname)    
    cobra.io.write_sbml_model(gene_expression_model,gene_expression_model_fname)
    #sampling of the solution flux
    aggregated_results, reaction_ids=sampling(model,n=1000,objective=None,starts=1,return_matrix=True,method="achr",thinning=100)
    
    #resultats del sampling
    reversible_matrix=aggregated_results
    reversible_reaction_ids=reaction_ids
    
    stat_dict=calculations_and_perceniles_rows(reversible_matrix,include_absolute_val_stats=True)
    sample_sampling_dict[sample]=stat_dict
    gene_expression_dict={}
   
    output=reactions_in_excel(model,fname=model_xlsx_filename,reaction_dict=reaction_expression_dict,gene_exp_pandas=gene_expression_pd,condition=condition, optimal_solution_dict=solution_dict,fva_dict=rev_fva,sampling_stat_dict=stat_dict,reaction_list=None)
   
    
#%%
import json

data_sufix="gimme.json"

keys=conditions_to_sample
sample_names="__".join(keys)
out_prefix=output_prefix+sample_names


with open(out_prefix+"_low_reaction_dict"+data_sufix,"w") as f:
    json.dump(low_reaction_dict,f) 
   
   
with open(out_prefix+"_fva_dict"+data_sufix,"w") as f:
     json.dump(fva_dict,f) 


with open(out_prefix+"_sample_sampling_dict"+data_sufix,"w") as f:
     json.dump(sample_sampling_dict,f) 
    


