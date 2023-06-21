#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:32 2023

@author: nuria
"""

import cobra
import math
import numpy as np
import sys
import os

ini_dir='/home/nuria/Escritorio/Macrophage-Metabolic_Phenotypes_in_Cytokine_Storms'
os.chdir(ini_dir)
import copy
from warnings import warn

import scipy
from scipy.sparse import dok_matrix
try:
    from sympy import Basic, Number
except:
    class Basic:
        pass

import pandas as pd
from cplex import Cplex, SparsePair,SparseTriple
from cplex.exceptions import CplexError
from six import iteritems, string_types
from six.moves import zip
import json
import cobra
from cobra import Reaction, Metabolite


from cobra.flux_analysis import pfba, flux_variability_analysis

from general_functions import metabolite_ex_dict

from Block_3_functions import run_qMTA, write_spreadsheet
import input3
#from qMTA import run_qMTA
reaction_pathway_dict_fname=input3.reaction_pathway_dict_fname

sampling_fname=input3.sampling_fname
data_dict=input3.data_dict

metabolomics_kpc_file=input3.metabolomics_kpc_file

metabolomics_file=input3.metabolomics_file

biomass_reaction=input3.biomass_reaction
drug_information_file=input3.drug_information_file
output_signficant_genes_only=input3.output_signficant_genes_only
output_omit_reactions_with_more_than_max_genes=input3.output_omit_reactions_with_more_than_max_genes
normalize_by_scale_genes=input3.normalize_by_scale_genes
normalize_by_scale_unchanged_reactions=input3.normalize_by_scale_unchanged_reactions



gene_weight=input3.gene_weight
reaction_weight=input3.reaction_weight
met_weight=input3.met_weight
kpc_weight=input3.kpc_weight

 
output_signficant_genes_only=input3.output_signficant_genes_only
output_omit_reactions_with_more_than_max_genes=input3.output_omit_reactions_with_more_than_max_genes
normalize_by_scale_genes=input3.normalize_by_scale_genes
normalize_by_scale_unchanged_reactions=input3.normalize_by_scale_unchanged_reactions
min_flux4weight=input3.min_flux4weight
coef_precision=input3.coef_precision


gene_parameters=input3.gene_parameters
normalize_by_scale_kpc=input3.normalize_by_scale_kpc
p_weight_formula_kpc=input3.p_weight_formula_kpc
normalize_p_weight_kpc=input3.normalize_p_weight_kpc
p_adj_th_kpc=input3.p_adj_th_kpc

kpc_parameters=input3.kpc_parameters

p_weight_formula_met_pellet=input3.p_weight_formula_met_pellet
p_adj_th_met_pellet=input3.p_adj_th_met_pellet
normalize_by_scale_met=input3.normalize_by_scale_met
normalize_p_weight_met=input3.normalize_p_weight_met
met_parameters=input3.met_parameters

search_for_targets=input3.search_for_targets
restrict_factor_for_targets =input3.restrict_factor_for_targets
gene_ko_parameters=input3.gene_ko_parameters
normalize_by_ref_flux_score=input3.normalize_by_ref_flux_score
min_flux_normalization_score=input3.min_flux_normalization_score
min_flux=input3.min_flux
p_adj_th_kpc_targets=input3.p_adj_th_kpc_targets
p_weight_formula_kpc_targets=input3.p_weight_formula_kpc_targets
normalize_p_weight_kpc_targets=input3.normalize_p_weight_kpc_targets

p_adj_th_pellet_targets=input3.p_adj_th_pellet_targets
p_weight_formula_met_pellet_targets=input3.p_weight_formula_met_pellet_targets
normalize_p_weight_pellet_targets=input3.normalize_p_weight_pellet_targets


max_reactionxgene_score=10

with open(reaction_pathway_dict_fname,"r") as f:
    reaction_pathway_dict=json.load(f)


with open(sampling_fname,"r") as f:
    sampling_dict=json.load(f)
    



#%%

gene_file="inputs/dge.xlsx" #is tghe dese file
deseq_file_pd=pd.read_excel(gene_file)
deseq_file_pd= deseq_file_pd.rename(columns={'NCBI.gene..formerly.Entrezgene..ID': 'gene'})
deseq_file_pd.dropna(subset=['gene'], inplace=True)
deseq_file_pd.set_index('gene', inplace=True)
padj_str='Ckpadj'
deseq_file_pd=deseq_file_pd[~deseq_file_pd[padj_str].isna()]
deseq_file_pd=deseq_file_pd.drop_duplicates(keep='first')
#very important comand to know whenever we have duplicated indexes
duplicate_indexes = deseq_file_pd.index[deseq_file_pd.index.duplicated()]

# Identify duplicated indexes
mask = deseq_file_pd.index.duplicated(keep='first')

# Filter the dataframe to keep only rows with unique indexes
deseq_file_pd = deseq_file_pd[~mask]
output_sheet={}
condition_mta_vres={}
out_dir='outputs/3'
os.chdir(out_dir)
for key in sorted(data_dict): 
    base_model=source_model=data_dict[key]["source_model"]
    target_model=data_dict[key]["target_model"]
    gene_file=data_dict[key]["differential_gene_file"]
    source_condition=data_dict[key]["source_condition"]
    target_condition=data_dict[key]["target_condition"]
    met_parameters["target_condition"]=target_condition
    met_parameters["source_condition"]=source_condition
    kpc_parameters["target_condition"]=target_condition
    kpc_parameters["source_condition"]=source_condition
    gene_parameters_mod=copy.deepcopy(gene_parameters)
    gene_parameters_mod.update(data_dict[key]) #Get pvalue and log2FC names from condition
    vref_dict={x:sampling_dict[source_condition][x]["mean"] for x in sampling_dict[source_condition]}
    #Add fake growth gene to biomass
    base_model.reactions.get_by_id(biomass_reaction).gene_reaction_rule="GROWTH"
    target_model.reactions.get_by_id(biomass_reaction).gene_reaction_rule="GROWTH"
    output_sheet_1, vres_dict, reaction_dict_dict, variation_dict_dict,up_genes, down_genes, log2fold_change_dict,   p_value_dict ,  gene_weight_dict,signficant_met_dict,signficant_kpc=run_qMTA(target_model,base_model,gene_fname=gene_file,vref_dict=vref_dict,gene_parameters=gene_parameters_mod,gene_weight=gene_weight,unchanged_reaction_weight=reaction_weight,met_weight=met_weight,reaction_pathway_dict=reaction_pathway_dict,key=key,metabolomics_file=metabolomics_file,met_parameters=met_parameters,output_signficant_genes_only=output_signficant_genes_only,output_omit_reactions_with_more_than_max_genes=output_omit_reactions_with_more_than_max_genes,normalize_by_scale_genes=normalize_by_scale_genes,min_flux4weight=min_flux4weight,coef_precision=coef_precision,normalize_by_scale_unchanged_reactions=normalize_by_scale_unchanged_reactions,kpc_file=metabolomics_kpc_file,kpc_weight=kpc_weight,kpc_parameters=kpc_parameters)
    output_sheet.update(output_sheet_1)
    condition_mta_vres[key]=vres_dict


output_prefix=''
sheets_to_merge=[key+"_a" for key in data_dict]
for n_1,sheet in enumerate(sheets_to_merge):
    if n_1==0:
       new_sheet=copy.deepcopy(output_sheet[sheet])
       new_sheet[0]=[x+"_"+sheet for x in new_sheet[0]] 
    else:
        for n,rows in enumerate(output_sheet[sheet]):
            if n==0:
                rows=[x+"_"+sheet for x in rows]
            new_sheet[n]+=rows

output_sheet["merged_sheet"]=new_sheet

key_name=".".join(data_dict.keys())

write_spreadsheet(output_prefix+key_name+"_qMTA.xlsx",output_sheet)


