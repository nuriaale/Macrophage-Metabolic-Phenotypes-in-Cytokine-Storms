#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 17:15:13 2023

@author: nuria
"""
import math
from cobra.flux_analysis import pfba, flux_variability_analysis
import pandas as pd 
from cobra import Model, Reaction, Metabolite
import copy
from cobra.manipulation.delete import  knock_out_model_genes,  remove_genes
from general_functions import *


def read_flux_constraints(model,constraint_pd,boundary_precision, col_reaction):
    precision=int(-1*(math.log10(boundary_precision)))
    for row in constraint_pd.iloc:
        if row[col_reaction] in model.reactions: #name of the column where the reaction id is kept
            reaction=model.reactions.get_by_id(row[col_reaction])
            reaction.lower_bound=round(float(row['Lower bound']),precision)
            reaction.upper_bound=round(float(row['Upper bound']),precision)
            reaction.objective_coefficient=round(float(row['Objective coefficient']),precision)
    return model

def check_fva(model,reaction_list,tolerance=1e-7):
    blocked_flag=False
    reaction_list_in_model=[]
    for rid in reaction_list:
        if rid not in model.reactions:
            print( rid,"missing in model" )
            blocked_flag=True
        else:
            reaction_list_in_model.append(rid)
    fva=flux_variability_analysis(model,fraction_of_optimum=0,reaction_list=reaction_list_in_model)
    for rid in reaction_list_in_model:
        max_abs_flux=max(abs(fva["maximum"][rid]), abs(fva["minimum"][rid])) 
        if max_abs_flux<tolerance:
           print( rid,"blocked",  fva["maximum"][rid], fva["minimum"][rid]  )
           blocked_flag=True
    return blocked_flag, fva

def simulate_multiple_gene_ko(model,reference_fluxes={},objective="BIOMASS",genes_lists=[],moma=True,fba_essential_genes=[],return_skipped=True,gpr=True,return_gene_dict=True,gene_expression_dict={},reaction_expression_dict={},expression_th=4,return_sl_reactions_dicts=False):
    constrained_model=model.copy()
    if gpr and gene_expression_dict!={} and reaction_expression_dict!={} and expression_th!=None:
       quantitative_gpr=True
    else:
       quantitative_gpr=False
    ko_dict={}
    if isinstance(reference_fluxes, pd.Series):
       reference_fluxes=reference_fluxes.to_dict()
    if reference_fluxes=={}:
       reference_fluxes=model.optimize().fluxes
       
    default_biomass=reference_fluxes[objective]
    if moma:
        constrained_model=create_euclidian_moma_model(model,reference_fluxes, wt_model=None)

    if genes_lists==[]:
        genes_lists=[[x] for x in constrained_model.genes]
    reactions_to_ko_gene_dict={}
    for gene_list in genes_lists:
        if gene_list in fba_essential_genes:
            ko_dict[str(gene_list)]=0
            print( "Gene is already essential with FBA, no need to check, skipping" )
            continue
        if quantitative_gpr:
           reactions_to_ko=[]
           new_reaction_expression_dict,value_list2,expression_dict=get_gene_exp(model,file_name="",gene_method="",gene_prefix="",gene_sufix="",expression_dict=gene_expression_dict,genes_to_ko=gene_list,omit_reflections=False)
           for reaction_id in new_reaction_expression_dict:
               old_expression_value=reaction_expression_dict[reaction_id]
               new_expression_value=new_reaction_expression_dict[reaction_id]
               if (new_expression_value+expression_th)<old_expression_value:
                   reactions_to_ko.append(model.reactions.get_by_id(reaction_id))
        
        elif gpr:
           gene_list_present=[x for x in  gene_list if x in constrained_model.genes]
           reactions_to_ko=knock_out_model_genes(constrained_model,gene_list_present)
           
        else:
            reactions_to_ko=[]
            for gene in gene_list:
                reactions_to_ko+=list(model.genes.get_by_id(gene).reactions)
        
        #print reactions_to_ko
        if reactions_to_ko==[]:
           if return_skipped:
               ko_dict[str(gene_list)]=default_biomass
               print( "Gene KO has no effect, skipping" )
           continue
        total_flux=sum([abs(reference_fluxes[x.id]) for x in reactions_to_ko])
        print( gene_list,reactions_to_ko )
        if total_flux<1e-9: #If the reactions do not carry flux you can skip it
           if return_skipped:
              print( "Reactions have no flux, skipping" )
              ko_dict[str(gene_list)]=default_biomass
           continue
        reaction_str=str(sorted([x.id for x in reactions_to_ko]))
        if reaction_str not in reactions_to_ko_gene_dict:
           reactions_to_ko_gene_dict[reaction_str]=[]
        reactions_to_ko_gene_dict[reaction_str].append(str(gene_list))
    ntotal=len(reactions_to_ko_gene_dict)
    for nr,reactions_to_ko in enumerate(reactions_to_ko_gene_dict):
        reaction_list=[constrained_model.reactions.get_by_id(x) for x in  eval(reactions_to_ko)]
        with constrained_model:  #Outised this context changes are reverted automatically
          for reaction in reaction_list:
            reaction.lower_bound=0
            reaction.upper_bound=0
          if moma:
              print( "MOMA" )
              solution=solve_moma_model(constrained_model,objective)
          else:
             solution=constrained_model.optimize()
          if solution.status!="optimal":
             print( "unfeasible" )
             result=0
          else:
             result=solution.fluxes[objective]
             
          for gene_id in reactions_to_ko_gene_dict[reactions_to_ko]:
              ko_dict[gene_id]={"result":result}#,"mod_reactions":moddified_reactions}
    
          print( str(nr)+"/"+str(ntotal),reactions_to_ko_gene_dict[reactions_to_ko],[x.id for x in reaction_list],str(round(result,4))+"/"+str(round(default_biomass,4))," gpr"+str(gpr), " qgpr "+str(quantitative_gpr) )
    if return_sl_reactions_dicts:
       reverse_reactions_to_ko_gene_dict={}
       for reactions in reactions_to_ko_gene_dict:
           for sl in reactions_to_ko_gene_dict[reactions]:
               reverse_reactions_to_ko_gene_dict[sl]=reactions
       return ko_dict,reverse_reactions_to_ko_gene_dict
    else:
       return ko_dict
   
def simulate_gene_ko(model,reference_fluxes={},objective="BIOMASS",genes_list=None,moma=True,fba_essential_genes=[],return_skipped=True,gpr=True,return_gene_dict=True,gene_expression_dict={},reaction_expression_dict={},expression_th=2,increased_reactions_flag=False, normalize_biomass=True,normalize_biomass_th=1.5,increased_flux_th=1e-5,verbose=True,skip_reactions_bellow=1e-8,number_iteration=None):
    constrained_model=model.copy()
    ko_dict={}
    gene_ko_reactions_dict={}
    increase_reaction_dict={}
    solution_x_dict={}
    total_flux_dict={}
    if gpr and gene_expression_dict!={} and reaction_expression_dict!={} and expression_th!=None:
       quantitative_gpr=True
    else:
        quantitative_gpr=False
    
    if isinstance(reference_fluxes, pd.Series):
       reference_fluxes=reference_fluxes.to_dict()
    if reference_fluxes=={}:
       
       reference_fluxes=model.optimize()                   
       reference_fluxes=reference_fluxes.fluxes            
    default_biomass=reference_fluxes[objective]
    if genes_list in (None,[]):
        genes_list=[x.id for x in constrained_model.genes]
    reactions_to_ko_gene_dict={}
    for gene_id in genes_list:
      if gene_id not in model.genes:
          print( "not found "+gene_id )
      if gene_id in model.genes:
        gene=model.genes.get_by_id(gene_id)
        if gene.id in fba_essential_genes:
            ko_dict[gene.id]=0
            print( "Gene is already essential with FBA, no need to check, skipping" )
            continue
        reaction_bound_dicts={}
        if quantitative_gpr:
           reactions_to_ko=[]
           new_reaction_expression_dict,value_list2,expression_dict=get_gene_exp(model,file_name="",gene_method="",gene_prefix="",gene_sufix="",expression_dict=gene_expression_dict,genes_to_ko=[gene.id],omit_reflections=False)
           for reaction_id in new_reaction_expression_dict:

               print( new_reaction_expression_dict )
               print( reaction_expression_dict[reaction_id] )
               old_expression_value=reaction_expression_dict[reaction_id]
               new_expression_value=new_reaction_expression_dict[reaction_id]
               if (new_expression_value+expression_th)<old_expression_value:
                   reactions_to_ko.append(model.reactions.get_by_id(reaction_id))
        elif gpr:
           reactions_to_ko=knock_out_model_genes(constrained_model,[gene.id])
        else:
            reactions_to_ko=gene.reactions
        gene_ko_reactions_dict[gene_id]=[x.id for x in reactions_to_ko]
        if reactions_to_ko==[]:
           if return_skipped:
               ko_dict[gene_id]=default_biomass
           continue
        total_flux=0
        for x in reactions_to_ko:
            if x.id in reference_fluxes:
                total_flux+=abs(reference_fluxes[x.id])
        actual_flux=abs(reference_fluxes[x.id])
        if total_flux<skip_reactions_bellow: 
           if return_skipped:
              ko_dict[gene.id]=default_biomass
           continue
        reaction_str=str(sorted([x.id for x in reactions_to_ko]))
        if reaction_str not in reactions_to_ko_gene_dict:
           reactions_to_ko_gene_dict[reaction_str]=[]
        reactions_to_ko_gene_dict[reaction_str].append(gene.id)
    ntotal=len(reactions_to_ko_gene_dict)
    for reactions_to_ko in reactions_to_ko_gene_dict:
        reaction_list=[constrained_model.reactions.get_by_id(x) for x in  eval(reactions_to_ko)]
        increased_reactions=[]
        with constrained_model: 
          for reaction in reaction_list:
            reaction.lower_bound=0
            reaction.upper_bound=0
      
          solution=constrained_model.optimize()

          if solution.status!="optimal":
             result=0
          else:
             result= solution.fluxes[objective]       
             if increased_reactions_flag:
                for reaction in model.reactions:
                   
                    if abs(solution.fluxes[reaction.id])-abs(reference_fluxes[reaction.id])>increased_flux_th:   
                          
                        if verbose:
                             
                             print( "increased_reactions. solution and reference_flux" , solution.fluxes[reaction.id], reference_fluxes[reaction.id], 'for reaction',reaction.id )   # Sergio: nou
                        increased_reactions.append(reaction.id)
                    elif normalize_biomass and abs(result)>1e-6 and reference_fluxes[reaction.id]!=0 and (reference_fluxes[objective]-result)>1e-3:
                         
                            if abs(round(solution.fluxes[reaction.id],5)/result)/(max(abs(round(reference_fluxes[reaction.id],5)),1e-6)/reference_fluxes[objective])>normalize_biomass_th:
                              if verbose:
                               
                                 print( "increased_reactions IN ELIF. solution and reference_flux", solution.fluxes[reaction.id], reference_fluxes[reaction.id] )
                              increased_reactions.append(reaction.id)

          for gene_id in reactions_to_ko_gene_dict[reactions_to_ko]:
              ko_dict[gene_id]=result
             
              solution_x_dict[gene_id]=copy.deepcopy(solution.fluxes)      
              increase_reaction_dict[gene_id]=copy.deepcopy(increased_reactions)
         
    if increased_reactions_flag:
       return ko_dict, increase_reaction_dict, gene_ko_reactions_dict,solution_x_dict
    else:
       return ko_dict
   

    
def changing_condition_name_list(pd_l, conditioname_l):
    new=[]
    for pd_x in pd_l:
        for change in conditioname_l:
            if change in pd_x.columns:
                pd_x=pd_x.rename(columns={change:'condition'})
                new.append(pd_x)
    return(new)
                
def general_calculations(pd_l, conditioname_l, col_to_erase_l, cond_to_keep, percet_thres):
    conditions_dict={}
    pd_l=changing_condition_name_list(pd_l, conditioname_l)     
    for condition in cond_to_keep:
        actual=pd.DataFrame()
        for pd_y in pd_l:
            pd_cond=pd_y[pd_y['condition']==condition]
            
            for a in col_to_erase_l:
                if a in pd_cond.columns:
                    pd_cond=pd_cond.drop(columns=a)
    
            for col in pd_cond.columns:
                non_null_count=pd_cond[col].count()
                if non_null_count/pd_cond.shape[0]<percet_thres:
                    pd_cond.drop(col, axis=1, inplace=True)
            mean=pd_cond.mean(skipna=True)
            std = pd_cond.std(ddof=0)
            se=std/pd_cond.count()
            calc = pd.concat([mean, std, se], axis=1)
            calc.columns = ['mean', 'std', 'se']
            
            if actual.empty:
                actual = calc
            else:
                actual = pd.concat([actual, calc])
        conditions_dict[condition]=actual
        
    return conditions_dict


def only_genes_in_model(model, pandas, id_col):
    genes_to_keep=[]
    for gene in pandas[id_col].dropna():
        if str(int(gene)) in model.genes:
            genes_to_keep.append(gene)
    pandas=pandas[pandas[id_col].isin(genes_to_keep)]
    return (pandas)
    

def add_kpcs(model,reaction_constraint_dict,always_include_0=False,add_upper_bound=True):
    for rid in reaction_constraint_dict:
        #print( rid )
        reaction=model.reactions.get_by_id(rid)
        reaction.lower_bound=reaction_constraint_dict[rid]["lb"]
        if  add_upper_bound:
            reaction.upper_bound=reaction_constraint_dict[rid]["ub"]
        if always_include_0:
           reaction.lower_bound=min(0,reaction.lower_bound)
           reaction.upper_bound=max(0,reaction.upper_bound)
        #print( reaction_constraint_dict[rid], reaction.bounds )      
        
def check_active_innactive_reactions(model_dict,reference_model,objective,fname, fraction=0.1,tolerance_feasibility=1e-6):
    fva_dict={}
    blocked_dict={}
    for condition in model_dict:
        condition_model=model_dict[condition].copy()
        condition_model.reactions.get_by_id(objective).objective_coefficient=1#change the objective cof of the reaction
        
        blocked_dict[condition]=[x.id for x in condition_model.reactions if (x.lower_bound==0 and x.upper_bound==0)]
        fva_condition=flux_variability_analysis(condition_model,fraction_of_optimum=fraction)
        
        fva_dict[condition]=fva_condition
    
    #Depending on the fluxes reported by the fva, the reaction is set as inactive normal or essential
   # reaction_dict={}
    reaction_dict2={}
    #reaction_list=[]
    for condition in fva_dict:
      fva_condition=fva_dict[condition]
      for rid in fva_condition.index:#go thougtyh all the reactions
          #if rid not in reaction_dict:
          if rid not in reaction_dict2:
              #reaction_list.append(rid)
              reaction_dict2[rid]=tuple()#since it is a tupple the order of the values is guaranted the same as the keys of fva dict
              
          max_flux=fva_condition["maximum"][rid]
          min_flux=fva_condition["minimum"][rid]
          if max(abs(max_flux),abs(min_flux))<tolerance_feasibility:
             if rid in blocked_dict[condition]:
                flag="innactive(direct)"
             else:
                flag="innactive(indirect)"
          
          elif min_flux>tolerance_feasibility:
               flag="essential_forward"
          elif max_flux<-tolerance_feasibility:
               flag="essential_reverse"
          else:
               flag="normal"
          #reaction_dict[rid][condition]=flag
          reaction_dict2[rid]+=(flag,)
          
    cols=["rid","name", "reaction", "GPR"]+sorted(fva_dict.keys())
    cols+=["info"]
    output_pd=pd.DataFrame(columns=cols)
    for rid in reaction_dict2:
        reaction_object=reference_model.reactions.get_by_id(rid)
        row=[rid,reaction_object.name,reaction_object.reaction,reaction_object.gene_reaction_rule]
        same=True
        comp=reaction_dict2[rid][0]
        innactive=False
        if comp in ("innactive(direct)","innactive(indirect)"):
            innactive=True
        for flag in reaction_dict2[rid]:
            row.append(flag)
            if innactive and flag in ("innactive(direct)","innactive(indirect)"):
                #here the reaction is incative and it is still
                innactive=True
            elif same and comp==flag:
                #here the reaction is still the same
                same=True
            else:#art least one time the reaction is not the same
                same=False
                inactive=False
                
        if same or inactive:
            total='It is the same'
        else:
            total='OMG it has different states'
         
        row.append(total) 
   
        output_pd=output_pd.append(pd.Series(row, index=output_pd.columns), ignore_index=True)
        
    output_pd.to_excel(fname)
    return fva_dict, output_pd

def clean_str(string_to_process,elements_to_remove=["[u'","'","[","]"]):
    for x in elements_to_remove:
        string_to_process=string_to_process.replace(x,"") 
    return(string_to_process)

