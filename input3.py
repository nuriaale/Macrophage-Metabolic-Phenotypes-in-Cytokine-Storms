#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 18:20:44 2023

@author: nuria
"""

reaction_pathway_dict_fname="inputs/reaction_pathway_dict_2021_09_03.json"
sampling_fname="outputs/2/CT__CK_sample_sampling_dictgimme.json"

import cobra
data_dict={"ct_vs_ck":{"source_model":cobra.io.read_sbml_model("outputs/1/CTmodel_blocked_low_reactions.sbml"),"target_model":cobra.io.read_sbml_model("outputs/1/CKmodel_blocked_low_reactions.sbml"),"target_condition":"CK","source_condition":"CT","differential_gene_file":"inputs/dge.xlsx","log2_str":"CKlog2fc","padj_str":"Ckpadj"}}


metabolomics_kpc_file="inputs/kpc_cobas_biocrates(medium)_seahorse.xlsx"
metabolomics_file=None
biomass_reaction="BIOMASS_reaction"
drug_information_file="inputs/drug_information.json" 
output_signficant_genes_only=False
output_omit_reactions_with_more_than_max_genes=False
normalize_by_scale_genes=True#"debug" #Debug =Warning this uses a normalization by flux value
normalize_by_scale_unchanged_reactions=True #"debug"



gene_weight=0.5    #Weight of genes

reaction_weight=0.5#Weight of uchangend reactions

met_weight=0.9 #Weight of intracelular metabolomics

kpc_weight=0.9 #Kpc weight

 
output_signficant_genes_only=False
output_omit_reactions_with_more_than_max_genes=False
normalize_by_scale_genes=True#"debug" #Debug =Warning this uses a normalization by flux value
normalize_by_scale_unchanged_reactions=True #"debug"
min_flux4weight=1e-6
coef_precision=7


gene_parameters={"log2_str":"log2FoldChange","log2_factor":1,"padj_str":"padj","p_th":0.25,"log2fc_th":0,"gene_str":"NCBI.gene..formerly.Entrezgene..ID","p_weight_formula":"-1*math.log(p_value,10)+math.log(0.25001,10)"}
normalize_by_scale_kpc=False  #Because we are using SD
p_weight_formula_kpc="1/(pow(max(sd1,0.00001),2))" #Because we are using SD
normalize_p_weight_kpc=False #Because we are using SD
p_adj_th_kpc=1 #Because we are using SD
from general_functions import metabolite_ex_dict
kpc_parameters={"kpc_name_dict":metabolite_ex_dict,"p_adj_th_kpc":p_adj_th_kpc,"p_weight_formula_kpc":p_weight_formula_kpc,"factor_kpc":1,"normalize_by_scale_kpc":normalize_by_scale_kpc,"normalize_p_weight":normalize_p_weight_kpc}

p_weight_formula_met_pellet="1/(pow(max(sd1,0.00001),2))" #Because we are using SD
p_adj_th_met_pellet=1 #Because we are using SD
normalize_by_scale_met=False #Because we are using SD
normalize_p_weight_met=False #Because we are using SD
met_parameters={"target_condition":None,"source_condition":None,"p_adj_th_met":p_adj_th_met_pellet,"convert_to_log_met":False,"p_weight_formula_met":p_weight_formula_met_pellet,"normalize_by_scale_met":normalize_by_scale_met,"normalize_p_weight":normalize_p_weight_met}

search_for_targets=True
restrict_factor_for_targets =4 #Restrict models not to increase more than X times the flux of Vref, Vres when looking for targets
gene_ko_parameters={"gpr":True,"ko_factor":0.5}
normalize_by_ref_flux_score=False
min_flux_normalization_score=1e-6
min_flux=min_flux_score=0
p_adj_th_kpc_targets=0.25 #Only consider significant variations for targets
p_weight_formula_kpc_targets="-1*math.log(p_value,10)+math.log(0.25001,10)" 
normalize_p_weight_kpc_targets=True

p_adj_th_pellet_targets=0.25 #Only consider significant variations for targets
p_weight_formula_met_pellet_targets="-1*math.log(p_value,10)+math.log(0.25001,10)" #"1/(pow(max(sd1,0.00001),2))"
normalize_p_weight_pellet_targets=True


max_reactionxgene_score=10

