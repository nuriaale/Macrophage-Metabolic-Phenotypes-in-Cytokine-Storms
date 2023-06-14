#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inputs for script 1
"""




model_name="inputs/recon3d.sbml" #base model used for analysis
metabolomics_file=None#Set to none if not present
constraint_file="inputs/constraints.xlsx" #Basic description of  medium composition, normalizyed by cell number. Restricts the maximum uptake of cells 
col_reaction='Unnamed: 0' #name of the colum that stores reaction ids
medium_kpcs_file="inputs/kpc_cobas_biocrates(medium)_seahorse.xlsx"#Measured consumption and production rates. Metabolites being consumed should be negative, metabolites being produced should be positive
colum_replicates=['KPC (umol/10e6 cell · h)']#name of the colum that indicates the condition of the replicates. 
col_to_erase_l=['condition', 'Replicado', 'Unnamed: 0']
#Metabolomics inputs
max_NA_fraction_kpc=0.5 #if max_NA_fraction_kpc is 0.5, fluxes with more than 50% NA will be removed
max_NA_fraction_pellet=0.5 #if max_NA_fraction_pellet is 0.5, metabolites with more than 50% NA will be removed

only_run_MFA=False #Run simulation to test metabolimics. WILL NOT PRUNE THE MODELS IF SET TO TRUE

#Gene expression inputs
gene_expression_file=file_name="inputs/fpkm.xlsx"  #Absolute gene expression file, in FPKM or equivalent units
conditions_of_interest=["CT",	"CK",	]#What experimental conditions we want to analyse


#Remove gene inputs
force_remove_th=1 #Genes with expression bellow 1 FPKM (absolute) will be removed
objective="BIOMASS_reaction" #Objective reactions
biomass_fraction=0.1 #Whats is the minimum fraction of bioamss allowed when removing genes 

#Differential gene expression file
deseq_file="inputs/dge.xlsx" #Differential gene expression
#P values for differential expression
comparisson_of_interest_deseq_pvalue={"CKpadj":"CK"}
p_value_th=0.05 #P value to consider differential expression
gene_entrez_column="entrezgene_id" #Name of the column with the entrez id
#reference_condition="WT" #Control condition, used for other comparissons
reference_condition=control_condition="CT" 

#Reactions to keep and kill
#Reactions that you dont want to be removed from the model
reactions_to_keep=["r1435","THFtm","GHMT2rm","MTHFDm","MTHFD2m","MTHFCm","FTHFLmi","FORtm","PDHm","ATPS4mi","CYOR_u10mi","NADH2_u10mi","CYOOm3i","CYOOm2i","ETF","ETFQO"]
#Reactions that form unfeasible loops and should be removed
reactions_to_kill=["HMR_4957","r1109", "CITL", "r2381","MALOAAtm","GALt1r","RE1342C","SBTD_D2","BALAPAT1tc2","CALAtr","PCRNte","LNELDCCRNte","ODECRNte","STCRNte","PMTCRNte","HDCECRNte", "AKGICITtm", "r2386", "r2387", "r2388", "r2389", "r2390", "r2391", "r2392", "r2393", "r2394", "r2371", "r2372", "r2375", "r2376", "r2377", "r2378", "r2379", "r2380", "MALITm", "CITtm", "r2385", "CITtcm", "MALICITtm", "AKGCITtm", "r2374", "r2382", "r0915", "HMR_4964"]


#Output directories
out_dir="OUTPUT" #Where the output files will be saved. Directory must exist

#Check essential reactions using FVA
fraction_essential_genes=0.1 #Reactions that must be active to allow growth above this fraction will be considered essential when looking for targets.

#Find targets
drug_information_file="inputs/drug_information.json" #File with annotations of known drugs and their targets. WARNING results should be manually checked as 

