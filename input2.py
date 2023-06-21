#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 18:13:47 2023

@author: nuria
"""

model_dict={"CT":"outputs/1/CTmodel_blocked_low_reactions.sbml", "CK":"outputs/1/CKmodel_blocked_low_reactions.sbml"}
output_prefix=""

gene_expression_file=file_name="inputs/fpkm.xlsx"

gpr_mode="full"
or_mode="sum"
or_mode="sum" 
convert_log2=True

absent_gene_expression=5 
absent_reaction_expression=100 
correct_for_complexes=True

gene_prefix="" 
gene_sufix=""

constraint_file=None

medium_kpcs_file="inputs/kpc_cobas_biocrates(medium)_seahorse.xlsx"

intracelular_metabolomics_file=None

max_NA_fraction_kpc=0.5 

max_NA_fraction_pellet=0.5

low_expression_threshold=100
base_penalty=1 
penalty_precision=2
correct_for_complexes=True
fraction_of_optimum_biomass=0.1
force_run=True
remove_reactions_below=0

gim3e_fraction_optimum=0.99

n_samples=1000
thinning=100 

metabolomics_file=None
gene_entrez_column="entrezgene_id" 

colum_replicates=['KPC (umol/10e6 cell · h)','umol/(hor*milion)' ]

col_to_erase_l=['condition', 'Replicado', 'Unnamed: 0']

run_fva=True
percentile=True
gene_method="average"
metabolite_list_fname=None
label_model=None
epsilon=0.0001
add_as_constraints=True
boundaries_precision=0.00001
all_gene_exp_4_percentile=False
convert_solution_dict_to_reversible=True
omit_0=True

out_dir=""