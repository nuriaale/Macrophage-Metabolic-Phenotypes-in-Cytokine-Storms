#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 09:44:14 2023

@author: nuria
"""
import os
import operator
ini_dir="/home/nuria/Escritorio/Macrophage-Metabolic_Phenotypes_in_Cytokine_Storms"
os.chdir(ini_dir)

import input1 as inp1
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.manipulation.delete import  knock_out_model_genes,  remove_genes
from cobra.flux_analysis import pfba, flux_variability_analysis
import pandas as pd 
import math
from LegacySolution import LegacySolution
from Block_1_functions import *


model_name=inp1.model_name
metabolomics_file=inp1.metabolomics_file
constraint_file=inp1.constraint_file 
col_reaction=inp1.col_reaction
medium_kpcs_file=inp1.medium_kpcs_file
colum_replicates=inp1.colum_replicates
col_to_erase_l=inp1.col_to_erase_l

#Metabolomics inputs
max_NA_fraction_pellet=inp1.max_NA_fraction_pellet
only_run_MFA=inp1.only_run_MFA

#Gene expression inputs
gene_expression_file=file_name=inp1.gene_expression_file
conditions_of_interest=inp1.conditions_of_interest

#Remove gene inputs
force_remove_th=inp1.force_remove_th
objective=inp1.objective
biomass_fraction=inp1.biomass_fraction

#Differential gene expression file
deseq_file=inp1.deseq_file
comparisson_of_interest_deseq_pvalue=inp1.comparisson_of_interest_deseq_pvalue
p_value_th=inp1.p_value_th
gene_entrez_column=inp1.gene_entrez_column
reference_condition=inp1.reference_condition
#control_condition=inp1.control_contition

#Reactions to keep and kill
#reactions_to_keep=[]
reactions_to_keep=inp1.reactions_to_keep
reactions_to_kill=inp1.reactions_to_kill

#Output directories
out_dir=inp1.out_dir

#Check essential reactions using FVA
fraction_essential_genes=inp1.fraction_essential_genes

#Find targets
drug_information_file=inp1.drug_information_file

model=cobra.io.read_sbml_model(model_name)

constraint_pd=pd.read_excel(constraint_file)

boundary_precision=1e-6
model=read_flux_constraints(model,constraint_pd,boundary_precision, col_reaction)        

flag_blocked,fva=check_fva(model,reactions_to_keep,tolerance=1e-7) 

for rid in reactions_to_kill:
    if rid in model.reactions:
       model.reactions.get_by_id(rid).remove_from_model()
       flag_blocked,fva=check_fva(model,reactions_to_keep)
       if flag_blocked:
          print( rid, "blocks reactions to keep" )
          raise Exception('Reactions to keep is incompatible with reactions to kill('+rid+')')

kpc_dict={}


if medium_kpcs_file is not None:
    kpc_pd=pd.read_excel(medium_kpcs_file,1)
    kpc_list=[kpc_pd]
    kpc_dict=general_calculations(kpc_list, colum_replicates, col_to_erase_l, conditions_of_interest, 0.50)



if metabolomics_file is not None:
    pellets_pd=pd.read_excel(metabolomics_file)
    pellets_dict=general_calculations([pellets_pd],colum_replicates, col_to_erase_l,conditions_of_interest, 0.50)
    met_sink_dict, rejected_list=add_sink_reactions_to_model(model, pellets_pd, col_to_erase_l, biocrates_name_compartment_dict)

    for a in pellets_dict:
        pellets_dict[a].index=pellets_dict[a].index + '(pellet)'

    for a in kpc_dict:
        kpc_dict[a]=pd.concat([kpc_dict[a], pellets_dict[a]])
    
    metabolite_ex_dict.update(met_sink_dict)
    

gene_expression_pd=pd.read_excel(gene_expression_file)


gene_expression_pd = gene_expression_pd.rename(columns={gene_entrez_column: 'gene'})
gene_expression_pd.dropna(subset=['gene'], inplace=True)#some times we have not the id, since we are not capable to look info about that gene 
#its better to eraseit. 
gene_expression_pd.set_index('gene', inplace=True)

gene_expression_pd=gene_expression_pd[~gene_expression_pd['CT'].isna()]

mask = gene_expression_pd.index.duplicated(keep='first')

gene_expression_pd = gene_expression_pd[~mask]

       
gene_to_remove_condition={}
gene_to_remove_condition_2={}
for condition in conditions_of_interest:
     gene_to_remove_condition[condition]=set()
     gene_to_remove_condition_2[condition]={}
     for gene in gene_expression_pd.index:
         if str(int(gene)) in model.genes and gene_expression_pd.loc[gene][condition]<=force_remove_th:
                 gene_to_remove_condition[condition].add(str(int(gene)))
                 gene_to_remove_condition_2[condition][str(int(gene))]=gene_expression_pd.loc[gene][condition]

#%%

biomass_dict={}
condition_essential_dict={}
putative_genes_to_remove_for_condition_dict={}

condition_gene_to_remove_actual={}
condition_needed_gens={}

model_dict={}

for condition in conditions_of_interest:
    condition_essential_dict[condition]=[]
    model_copy=model.copy()
    modify_model_for_seahorse(model_copy,remove=False)
    out_mfa=condition+"MFA_fit.xlsx"
    x_dict,total_xi, xi_dict=mfa(model_copy,kpc_dict,condition,metabolite_ex_dict,min_std=1e-4,precision=6,fname=out_mfa) 
    if kpc_dict!={}:
        constraint_dict={}
        for ex_rid in xi_dict:
            constraint_dict[ex_rid]={"ub":x_dict[ex_rid]+1e-6,"lb":x_dict[ex_rid]-1e-6} 
            add_kpcs(model_copy,constraint_dict,always_include_0=False,add_upper_bound=True)    
    reference_biomass=model_copy.optimize() 
    biomass_dict[condition]=reference_biomass
    #gene essensiality analysis
    gene_ko_fba=simulate_gene_ko(model_copy,reference_fluxes={},objective=objective,genes_list=None,moma=False,fba_essential_genes=[],return_skipped=True,gpr=True,return_gene_dict=True,gene_expression_dict={},reaction_expression_dict={},expression_th=2,increased_reactions_flag=False, normalize_biomass=True,normalize_biomass_th=1.5,increased_flux_th=1e-5,verbose=True,skip_reactions_bellow=1e-8,number_iteration=None)
    for gene in gene_ko_fba:
       if gene_ko_fba[gene]<reference_biomass.objective_value*biomass_fraction:  
          condition_essential_dict[condition].append(gene)
          
    genes_to_remove_ordered=()
    local_dict=gene_to_remove_condition_2[condition]
    genes_to_remove_ordered_tupple = sorted(local_dict.items(), key=operator.itemgetter(1))
    for x in genes_to_remove_ordered_tupple:
        if x not in condition_essential_dict[condition]:
            genes_to_remove_ordered+=(x[0],)
        else:
            print('error',x)
    putative_genes_to_remove_for_condition_dict[condition]=genes_to_remove_ordered
    
    gene_to_remove_actual=[]
    needed_gens=[]
    min_biomass=model_copy.optimize().objective_value*biomass_fraction 
    gene_to_remove_ordered=putative_genes_to_remove_for_condition_dict[condition]
    #simulated gene ko, in a model copy
    for gene in gene_to_remove_ordered:
        gene_to_remove_actual.append(gene)
        with model_copy as model_ko:
             reactions=knock_out_model_genes(model_ko,gene_to_remove_actual)

             solution=model_ko.optimize()
             if solution.status=="optimal":
                flag_blocked,fva=check_fva(model_ko,reactions_to_keep)
                
                if flag_blocked:
                    remove_flag=True#solution is optimal but fva reported an error
                elif solution.fluxes[objective]<min_biomass:    
                    remove_flag=True#solution is optimal but fluxes dont pass a threslhold
                else:
                    remove_flag=False#the gene can't be removed, this ko has an impact in biomass production
        
             else:#solution not optimal
                flag_blocked=True 
                remove_flag=True
                
                
             if remove_flag:#the gene can be removed, this ko NO has an impact in biomass production
               gene_to_remove_actual.remove(gene)
               needed_gens.append(gene)
               
    condition_gene_to_remove_actual[condition]=gene_to_remove_actual  
    condition_needed_gens[condition]=needed_gens
    
    #actual gene ko and reaction-associated elimination
    gene_to_remove_actual=condition_gene_to_remove_actual[condition]   
    model_metabo=model.copy()
    reactions=knock_out_model_genes(model_metabo,gene_to_remove_actual)
    cobra.manipulation.remove_genes(model_metabo, gene_to_remove_actual, remove_reactions=False)
    fname_out=condition+"model_blocked_low_reactions.sbml"
    cobra.io.write_sbml_model(model_metabo,fname_out)
    model_dict[condition]=model_metabo
          
#%%
#Create an extra output file that shows which whenever a reaction is set as essential, active or inactive in the model
             
fva_dict, rows_pd =check_active_innactive_reactions(model_dict, model, objective,fname="active_innactive_reactions.xlsx", fraction=fraction_essential_genes,tolerance_feasibility=1e-6)

