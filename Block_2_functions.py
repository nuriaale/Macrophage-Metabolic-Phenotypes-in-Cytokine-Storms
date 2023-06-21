#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 18:24:18 2023

@author: nuria
"""

from cobra.sampling import OptGPSampler
from cobra.sampling.achr import ACHRSampler
from general_functions import *
from cobra.flux_analysis import pfba, flux_variability_analysis
import numpy as np
from general_functions import *

def average_expression(model,gene_expression_pd,condition, absent_gene_expression=50,percentile=True):
    condition = str(condition)
    value_list=list(gene_expression_pd[condition])
    if percentile==True:
        absent_gene_expression=np.percentile(value_list,absent_gene_expression)
    
    reaction_expression_dict={}
    for reaction in model.reactions: 
        genes_present=[x.id for x in reaction.genes if int(x.id) in gene_expression_pd.index ]
        if len(genes_present)==0:
            continue
        if reaction.gene_reaction_rule != "" :
            value_sum=0
            for gene in reaction.genes: 
                if gene.id not in genes_present:#meaning they are not in the gene_expression excel
                    value_sum+=absent_gene_expression
                else:
                    value_sum+=gene_expression_pd[condition][int(gene.id)]
                   
            reaction_expression_dict[reaction.id]=value_sum
             
    return reaction_expression_dict, value_list

    

def create_gim3e_model(cobra_model,gene_expression_pd,condition, metabolite_list=[],label_model=None,epsilon=0.0001,gene_method="average",low_expression_threshold=25,absent_gene_expression=100,absent_reaction_expression=100,percentile=True,all_gene_exp_4_percentile=False,base_penalty=0,milp=True,or_mode="max",gpr_mode="full",omit_0=True,log2=False,penalty_mode="normal", penalty_precision=4,reaction_ids_to_omit=[],correct_for_complexes=False):
    #if gpr_mode=="average":
    reaction_expression_dict, value_list = average_expression(cobra_model,gene_expression_pd,condition, absent_gene_expression, percentile)
   
    #elseeeeeeee   
    if log2=="reverse":#this time not happening
        for x in reaction_expression_dict:
            reaction_expression_dict[x]=pow(2,reaction_expression_dict[x]) 
        value_list=[pow(2,x) for x in value_list]
    
    if log2: #it is true, this is happening
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=math.log(reaction_expression_dict[x]+0.1,2)
       value_list=[]
       for num in gene_expression_pd[condition]:
           value_list.append(math.log(num+0.1,2))
           
    #if absent_reaction_expression!=100:##por aqui no pasa
    
    if percentile==True:
        low_expression_threshold=np.percentile(value_list,low_expression_threshold)
        
    penalty_dict={}
    if base_penalty!=0:
        for reaction in cobra_model.reactions:
            penalty_dict[reaction.id]=base_penalty
            
    for reaction_id in  reaction_expression_dict:
        gene_expression_value=reaction_expression_dict[reaction_id]
        if reaction_id not in penalty_dict:
              penalty_dict[reaction_id]=0
              
        if penalty_mode=="division":#no pasas x aqui
           penalty_dict[reaction_id]=round(1000/gene_expression_value,3) 
           
        else:
            if gene_expression_value< low_expression_threshold:
              if reaction_id not in penalty_dict:
                 penalty_dict[reaction_id]=0
              gene_expression_penalty=round(-gene_expression_value+low_expression_threshold,4)
              penalty_dict[reaction_id]+=gene_expression_penalty
    
    if correct_for_complexes: #esto pasa
        for reaction in cobra_model.reactions:
            min_coef=min([abs(reaction.metabolites[x]) for x in reaction.metabolites])     
            if reaction.id in penalty_dict:
               penalty_dict[reaction.id]*=min_coef 
               
    for rid in reaction_ids_to_omit:
        del(penalty_dict[rid])
        
    convert_to_irreversible_with_indicators(cobra_model,penalty_dict.keys(),metabolite_list=[], mutually_exclusive_directionality_constraint = milp,label_model=label_model)
    #i think mlipo is false, label model is none. 
    if len(metabolite_list)>0:
       add_turnover_metabolites(cobra_model, metabolite_id_list=metabolite_list, epsilon=epsilon,label_model=label_model)
    
    objective_reaction = Reaction('gim3e_objective')
    gim3e_indicator = Metabolite('gim3e_indicator',formula='',name='',compartment='')
    objective_reaction.add_metabolites({gim3e_indicator:-1})
    cobra_model.add_reactions([objective_reaction])
    objective_reaction.objective_coefficient=-1
    total_bound=0         
    for reaction_id in penalty_dict:
           gene_expression_penalty=round(penalty_dict[reaction_id],penalty_precision)
           reaction=cobra_model.reactions.get_by_id(reaction_id)
           reaction.add_metabolites({gim3e_indicator:gene_expression_penalty})
           total_bound+=gene_expression_penalty*reaction.upper_bound
           if "reflection" in reaction.notes:
               reflection_id=reaction.notes["reflection"]
               reflection=cobra_model.reactions.get_by_id(reflection_id)
               reflection.add_metabolites({gim3e_indicator:gene_expression_penalty})
               total_bound+=gene_expression_penalty*reflection.upper_bound
    objective_reaction.lower_bound=0.0
    objective_reaction.upper_bound=total_bound
    return penalty_dict

def get_innactive_reactions(model,optimal_solution_dict,gene_expression_pd,condition, absent_gene_expression, percentile,log2=False):
   
    reaction_expression_dict, value_list  = average_expression(model,gene_expression_pd,condition=condition,absent_gene_expression=absent_gene_expression, percentile=percentile )
    if log2=="reverse":
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=pow(2,reaction_expression_dict[x]) 
    if log2==True:
       for x in reaction_expression_dict:
           reaction_expression_dict[x]=math.log(reaction_expression_dict[x]+0.1,2)
 
    core_reactions=[x for x in optimal_solution_dict.keys() if abs(optimal_solution_dict[x])>1e-8]  # Sergio: nou

    reflection_dict={}
    for x in model.reactions:
        if "reflection" in x.notes:
            reflection_dict[x.id]=x.notes["reflection"]
    
            
    low_th=np.percentile(value_list,percentile)
    low_reactions=[x for x in reaction_expression_dict if reaction_expression_dict[x]< low_th]
    low_reactions=[x for x in low_reactions if x not in core_reactions and reflection_dict.get(x) not in core_reactions]
    

    return low_reactions, core_reactions,low_th,reaction_expression_dict, value_list

def flux_minimization_fva(model,solver,reaction_list,lp_tolerance_feasibility,objective,fraction_optimum,reaction_ids_to_omit):
  irreversible_model=model.copy() 
  if objective!=None and fraction_optimum!=None:
     objetctive_reaction=irreversible_model.reactions.get_by_id(objective)
     #
     #objetctive_reaction_ub=objetctive_reaction.upper_bound#
     #
     objetctive_reaction.objective_coefficient=-1 
     
     irreversible_model.optimize()
     xval=objetctive_reaction.x
     objetctive_reaction.objective_coefficient=0 
     objetctive_reaction.upper_bound=round_sig(xval/fraction_optimum,6)
     
  elif objective!=None:
      objetctive_reaction=irreversible_model.reactions.get_by_id(objective)
      objetctive_reaction.objective_coefficient=0
      
  reaction2test=[]
  irreversible_fva={}
  model_reaction_n_dict={}
  if reaction_list!=[]:
       for reaction_id in reaction_list:
           reaction=irreversible_model.reactions.get_by_id(reaction_id)
           if reaction.bounds==(0,0):
              irreversible_fva[reaction] 
              continue 
           reaction2test.append(reaction_id)  
           if "reflection" in reaction.notes:
               reflection= reaction.notes["reflection"]
               if reflection not in reaction2test:
                  reaction2test.append(reflection)
       
  else:
     for reaction in irreversible_model.reactions:
         if  "IRRMILP_" not in reaction.id:
             reaction2test.append(reaction.id) 
  for reaction in reaction2test:
    if reaction in reaction_ids_to_omit:
          reaction2test.remove(reaction)

  for n,reaction in enumerate(irreversible_model.reactions):
      model_reaction_n_dict[reaction.id]=n
  normal_fva_list=[]
  counter=0
  print( "starting analysis" )
  for reaction_id in reaction2test:
      reaction=irreversible_model.reactions.get_by_id(reaction_id)
      if "reflection" in reaction.notes :
               counter+=1
               irreversible_model2=irreversible_model
               reaction=irreversible_model2.reactions.get_by_id(reaction_id)
               reflection_id= reaction.notes["reflection"]
               n=model_reaction_n_dict[reaction_id]
               n_reflection=model_reaction_n_dict[reflection_id]
               reflection=irreversible_model2.reactions.get_by_id(reflection_id)
               with irreversible_model:
                   #Test max
                   irreversible_model2.reactions.get_by_id(reflection_id).objective_coefficient=-100.0 #make the reflection inacive 
                   irreversible_model2.reactions.get_by_id(reaction_id).objective_coefficient=1.0
                   
                   max_sol=irreversible_model2.optimize()                                                                # Sergio: nou, optimize sin argumentos
                   
                   
                   #max_value=max_sol.x_dict[reaction_id]     # Sergio: antic
                   max_value=max_sol.fluxes[reaction_id]      # Sergio: nou
                   
                  
                   status1=irreversible_model2.optimize()    
               with irreversible_model:
                    irreversible_model3=irreversible_model
                    irreversible_model3.reactions.get_by_id(reflection_id).objective_coefficient=0 #make the reflection inacive 
                    irreversible_model3.reactions.get_by_id(reaction_id).objective_coefficient=-1.0
                   
                    min_sol=irreversible_model3.optimize()                                                                # Sergio: nou.
                    
                    #min_value=min_sol.x_dict[reaction_id]    # Sergio: antic
                    min_value=min_sol.fluxes[reaction_id]     # Sergio: nou
                    
        
               irreversible_fva[reaction_id]={"minimum":min_value,"maximum":max_value}
               #print('C,R_id', counter,reaction_id,irreversible_fva[reaction_id],max_sol.status,min_sol.status )

      else:
        normal_fva_list.append(reaction_id)
  if normal_fva_list!=[]:
     #print(  irreversible_model.optimize() )
     try:
       #print( "running FVA" )
       normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0)
       
     except:
       #print( "using cplex",  lp_tolerance_feasibility )
       normal_fva=flux_variability_analysis(irreversible_model,reaction_list=normal_fva_list,fraction_of_optimum=0.0,tolerance_feasibility=lp_tolerance_feasibility,solver="cplex") #changed
    
     #print('## normal_fva',normal_fva)
     for reaction_id in normal_fva['minimum'].keys():   # Sergio: check ?????   
     
         irreversible_fva[reaction_id]={'minimum':normal_fva['minimum'][reaction_id],'maximum':normal_fva['maximum'][reaction_id]}     # Sergio: nou, prime
  reversible_fva={}
  reverese_reactions=irreversible_model.reactions.query("_reverse") 
  #print irreversible_fva
  for reaction_id in irreversible_fva:
        #print('#$ Reaction_id',reaction_id)
        reaction=irreversible_model.reactions.get_by_id(reaction_id)
        if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
           continue
        reversible_fva[reaction.id]=copy.deepcopy(irreversible_fva[reaction.id])
        if "reflection" in reaction.notes:
            reverse_reaction= reaction.notes["reflection"]
            if irreversible_fva[reverse_reaction]["maximum"]!=0:
               reversible_fva[reaction.id]["minimum"]=-irreversible_fva[reverse_reaction]["maximum"]
            if irreversible_fva[reverse_reaction]["minimum"]!=0: 
               reversible_fva[reaction.id]["maximum"]=-irreversible_fva[reverse_reaction]["minimum"] 
  #print reversible_fva
  #objetctive_reaction.upper_bound=objetctive_reaction_ub  
  return reversible_fva,irreversible_fva  

def integrate_omics_gim3e_and_remove(metabolic_model,gene_expression_pd,condition, fraction_of_optimum=1,low_expression_threshold=25,absent_gene_expression=50,absent_reaction_expression=100,percentile=True,gene_method="average",metabolite_list_fname=None,label_model=None,epsilon=0.0001,gim3e_fraction_optimum=0.75,run_fva=True,add_as_constraints=True,boundaries_precision=0.001,all_gene_exp_4_percentile=False,base_penalty=0,convert_solution_dict_to_reversible=True,gpr_mode="full",or_mode="max",omit_0=True,log2=False,reaction_list=[],penalty_mode="normal",abs_max_objective=None,lp_tolerance_feasibility=1e-6,remove_reactions_below=33,remove_percentile=True,gene_expression_model=None,penalty_precision=4,reaction_ids_to_omit=[],reactions_to_keep=[],correct_for_complexes=False):
    low_reactions=[]
    core_reactions=[]   
    precision=int(-1*(math.log10(boundaries_precision)))
    
    if gene_expression_model==None:
      gene_expression_model=copy.deepcopy(metabolic_model)
      if metabolite_list_fname!=None and metabolite_list_fname!="": 
          metabolite_list=read_metabolomics_data(gene_expression_model,metabolite_list_fname)
      else:
           metabolite_list=[]
      for reaction in gene_expression_model.reactions:
          if reaction.objective_coefficient!=0:
              #reaction_l_act=[reaction]
              #fva1=flux_variability_analysis(gene_expression_model,fraction_of_optimum,reaction_l_act)
              fva1=flux_variability_analysis(gene_expression_model,fraction_of_optimum=fraction_of_optimum,reaction_list=[reaction])
              reaction.lower_bound=max(round_down(fva1["minimum"][reaction.id]-boundaries_precision/2.0,precision),reaction.lower_bound)   # Sergio: nou
              reaction.upper_bound=min(round_up(fva1["maximum"][reaction.id]+boundaries_precision/2.0,precision),reaction.upper_bound)     # Serigo: nou
              reaction.objective_coefficient=0
            
      milp=False
      penalty_dict=create_gim3e_model(gene_expression_model,gene_expression_pd,condition, metabolite_list,label_model,epsilon,gene_method,low_expression_threshold,absent_gene_expression,absent_reaction_expression,percentile,all_gene_exp_4_percentile,base_penalty,milp,or_mode,gpr_mode,omit_0,log2,penalty_mode,penalty_precision,reaction_ids_to_omit,correct_for_complexes)
    #else: #per el que no pasa          
    solution=gene_expression_model.optimize()
    objective=-1*float(solution.objective_value)
    solution_dict={}
    
    irreversible_optimal_solution_dict=solution.fluxes
    
    #reverese_reactions=gene_expression_model.reactions.query("_reverse")
    if convert_solution_dict_to_reversible: 
      for reaction in gene_expression_model.reactions:
         reaction_id=reaction.id
         if reaction_id not in metabolic_model.reactions:
            continue 
         """reaction=gene_expression_model.reactions.get_by_id(reaction_id)
         if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
            continue"""
         solution_dict[reaction_id]=irreversible_optimal_solution_dict[reaction_id]
         if "reflection" in reaction.notes:
             reverse_reaction= reaction.notes["reflection"]
             solution_dict[reaction_id]-=irreversible_optimal_solution_dict[reverse_reaction]
    #else#per el que no pasa
    if remove_reactions_below>0:
       low_reactions, core_reactions,low_th, reaction_expression_dict,value_list=get_innactive_reactions(gene_expression_model,irreversible_optimal_solution_dict,gene_expression_pd,condition, absent_gene_expression,percentile=remove_reactions_below,log2=log2)
       for rid in reactions_to_keep:
           if rid in low_reactions:
              low_reactions.remove(rid)
       for rid in low_reactions:
           gene_expression_model.reactions.get_by_id(rid).bounds=(0,0)
    else:
                                                                                                     
        low_reactions, core_reactions,low_th, reaction_expression_dict,value_list=get_innactive_reactions(gene_expression_model,irreversible_optimal_solution_dict,gene_expression_pd,condition, absent_gene_expression,percentile=remove_reactions_below,log2=log2)
        #reaction_expression_dict={}
        #gene_expression_dict={}
        #value_list={} 
        low_th=remove_reactions_below
    
    if run_fva==True:
      if  abs_max_objective!=None:
          gene_expression_model.reactions.get_by_id("gim3e_objective").upper_bound=abs_max_objective
          gim3e_fraction_optimum=None
      
                #flux_minimization_fva(model,solver,reaction_lis,lp_tolerance_feasibility,objective,fraction_optimum,reaction_ids_to_omit):
      solver=None
      rev_fva,irrevfva=flux_minimization_fva(gene_expression_model,solver,reaction_list,lp_tolerance_feasibility,objective="gim3e_objective",fraction_optimum=gim3e_fraction_optimum,reaction_ids_to_omit=reaction_ids_to_omit)
      if add_as_constraints==True:
        for reaction_id in rev_fva: 
          if reaction_id in metabolic_model.reactions:
             if "RGROUP_" in reaction_id:
                 continue
             reaction=metabolic_model.reactions.get_by_id(reaction_id)
             lower_bound=rev_fva[reaction_id]["minimum"]
             upper_bound=rev_fva[reaction_id]["maximum"] 
             reaction.lower_bound=max(round_down(lower_bound,precision),reaction.lower_bound)
             reaction.upper_bound=min(round_up(upper_bound,precision),reaction.upper_bound)       
             reaction.objective_coefficient=0 
    else:
      rev_fva={}           
    return penalty_dict,gene_expression_model,objective,solution_dict,rev_fva, low_reactions, reaction_expression_dict

def sample_0133(model, n, method="optgp", thinning=100, processes=1, seed=None):
    if method == "optgp":
        sampler = OptGPSampler(model, processes, thinning, seed=None)
    elif method == "achr":
        print('about to sample')
        sampler = ACHRSampler(model, thinning=thinning, seed=seed)
        print('sampling done')
    else:
        raise ValueError("method must be 'optgp' or 'achr'!")

    return pd.DataFrame(columns=[rxn.id for rxn in model.reactions],
                            data=sampler.sample(n))


def sampling(model,n=100,processes=6,objective=None,starts=1,return_matrix=False,method="optgp",thinning=100):
    print("Method:", method, "Thinning:", thinning )
    reaction_ids=[x.id for x in model.reactions]
    if objective!=None:
        print( model.reactions.get_by_id(objective).lower_bound )
    flux_dict_list=[]
    for i in range(0,starts):
       result_matrix = sample_0133(model, n,processes=processes,method=method,thinning=thinning)#.as_matrix() #Valid methods are optgp and achr. Process is only used in optgp. Thinning (“Thinning” means only recording samples every n iterations) is only used in achr
       if not return_matrix:
         for row in result_matrix:
          flux_dict={}
          for n_flux,flux in enumerate(row):
            flux_dict[reaction_ids[n_flux]]=flux
          flux_dict_list.append(flux_dict)
          if objective!=None:
             print( flux_dict[objective] )
       elif return_matrix:
            if i==0:
               aggregated_results=result_matrix
            else:
               aggregated_results=np.vstack((aggregated_results,result_matrix)) 
    if not return_matrix:
       return flux_dict_list
    else:
       return np.transpose(aggregated_results), reaction_ids
def calculations_and_perceniles_rows(matrix,include_absolute_val_stats=False,percentiles=[0.25,0.50,0.75]):
    mean=matrix.mean(axis=1)
    std=matrix.std(axis=1, ddof=0)
    stat_dict={}
    if include_absolute_val_stats:
        abs_matrix=matrix.abs()
        abs_mean=abs_matrix.mean(axis=1)
        abs_std=abs_matrix.std(axis=1, ddof=0)
    
        percentile_dict={}
        abs_percentile_dict={}
        for per in percentiles:
            percentile_dict[str(per)]= matrix.apply(lambda x: x.quantile(per), axis=1)
            abs_percentile_dict[str(per)]= abs_matrix.apply(lambda x: x.quantile(per), axis=1)

        for reac in matrix.index:
            stat_dict[reac]={"mean":mean[reac],"std":std[reac],"max":max(matrix.loc[reac]),"min":min(matrix.loc[reac])}
            stat_dict[reac]["percentiles"]={'0.25':percentile_dict['0.25'].loc[reac], '0.5':percentile_dict['0.5'].loc[reac], '0.75':percentile_dict['0.75'].loc[reac]}

            stat_dict[reac]["abs_percentiles"]={'0.25':abs_percentile_dict['0.25'].loc[reac], '0.5':abs_percentile_dict['0.5'].loc[reac], '0.75':abs_percentile_dict['0.75'].loc[reac] }
            stat_dict[reac]["abs_mean"]=abs_mean[reac]
            stat_dict[reac]["abs_std"]=abs_std[reac]
            
    else:
        percentile_dict={}
        for per in percentiles:
            percentile_dict[str(per)]= matrix.apply(lambda x: x.quantile(per), axis=1)
            
        for reac in matrix.index:
             stat_dict[reac]={"mean":mean[reac],"std":std[reac],"max":max(matrix.loc[reac]),"min":min(matrix.loc[reac])}
             stat_dict[reac]["percentiles"]={'0.25':percentile_dict['0.25'].loc[reac], '0.5':percentile_dict['0.5'].loc[reac], '0.75':percentile_dict['0.75'].loc[reac]}
    return stat_dict
    
def reactions_in_excel(model,fname,reaction_dict={},gene_exp_pandas={},condition=None, optimal_solution_dict={},fva_dict={},sampling_stat_dict={},reaction_list=None):
    output_pd=pd.DataFrame(columns=["Reaction id", "Reaction", "Reaction Specified", "Reaction name", "Optimal solution", "Minimum", "Maximum","Gene expression(fpkm)", "Gene reaction rule", "gene expression", "Innactive(direct or indirect)", "mean", "SD"])
    if reaction_list in ([],None):
       reactions=model.reactions
    else:
       reactions=[model.reactions.get_by_id(x) for x in reaction_list] 
       
    for reaction in reactions:
        row=[reaction.id, reaction.reaction]
        
        reaction_str=" "+reaction.reaction+" "
        for x in reaction.metabolites:
            if x.compartment!=None:
              compartment_str="["+x.compartment+"]"
            else:
              compartment_str=""
            if x.name==None:
               x.name=x.id   
            reaction_str=reaction_str.replace(" "+x.id+" "," "+x.name+compartment_str+" ")
        row.append(reaction_str)
        row.append(reaction.name)
        try:
           row_sol=optimal_solution_dict.get(reaction.id)
        except:
           row_sol="-"   
        row.append(row_sol)
        
        try:
           row_min=fva_dict.get[reaction.id]["minimum"]
           row_max=fva_dict.get[reaction.id]["maximum"]
        except:
           row_min="-"
           row_max="-"
        row.append(row_min)
        row.append(row_max)
        #mirara next two line
        row.append(reaction_dict.get(reaction.id))
        
        row.append(reaction.gene_reaction_rule)
        
        local_dict={}
        for gene in reaction.genes: 
            try:
                local_dict[gene.id]=gene_exp_pandas[condition][int(gene.id)]
            except:
                local_dict={}
            
        row.append(str(local_dict))
        
        try:
          row.append(str(fva_dict.get(reaction.id)["minimum"]>-1e-6 and fva_dict.get(reaction.id)["maximum"]<1e-6))
        except:
            row.append('-')
        try:
          mean=sampling_stat_dict[reaction.id]["mean"]
          std=sampling_stat_dict[reaction.id]["std"]
        except:
            mean="-"
            std="-"
            
        row.append(mean)
        row.append(std)
        
        #print(len(row))
        output_pd=output_pd.append(pd.Series(row, index=output_pd.columns), ignore_index=True)
        
        
    output_pd.to_excel(fname)
    return output_pd

def convert_to_irreversible_with_indicators(cobra_model,reaction_id_list,metabolite_list, mutually_exclusive_directionality_constraint = False,label_model=None,reactions_with_no_indicators=[]):
    #Function modified from the work by : """Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    """Will break all of the reversible reactions into two separate irreversible
     reactions with different directions.  This function call modified from
     a version in the core cobra to facilitate the MILP formulation and
     include gene_reaction_rules with the reverse reaction
   
     Arguments:
      cobra_model: A model object which will be modified in place.
      mutually_exclusive_directionality_constraint: Boolean.  If True, turnover 
       reactions are constructed to serve as MILP constraints to prevent loops.
      
     Returns:
      None, cobra_model is modified in place
    
    
    
    
     """
    reactions_to_add = []
    #from cobra.core.Reaction import Reaction
    #from cobra.core import Metabolite
    reactions_to_make_irreversible=[]
    for x in reaction_id_list:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))
    """for x in lexs:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))"""
    
    #If a label model object  is provided make sure all experimentally measured metabolites (or at least one of the metabolites in the pool) is produced
    full_metabolite_list=copy.copy(metabolite_list)
    print( full_metabolite_list)
    if label_model!=None:
       emus=[]
       for condition in label_model.experimental_dict:
           for emu in label_model.experimental_dict[condition]:
               if emu not in emus:
                  emus.append(emu)
       measured_metabolite_dict={}
       
       for emu in emus:
           iso_id=str(label_model.emu_dict[emu]["met_id"])
           #print label_model.id_isotopomer_object_dict
           #isotopomer_object=label_model.id_isotopomer_object_dict[iso_id]
           metabolites=label_model.isotopomer_id_metabolite_id_dict[iso_id]
           print( [iso_id,label_model.isotopomer_id_metabolite_id_dict[iso_id]] )
           if isinstance(metabolites,list):
              for metabolite in metabolites:
                  full_metabolite_list.append(metabolite)
           else:
              full_metabolite_list.append(metabolites)
    
    for metabolites in full_metabolite_list:
       print( '- Metabolite: ',metabolites )
       if not isinstance(metabolites,list):
          metabolites=[metabolites]
       for metabolite in metabolites:
          print('-- Metabolite:', metabolite )
          the_metabolite=cobra_model.metabolites.get_by_id(metabolite)
          for x in the_metabolite.reactions:
             if x not in reactions_to_make_irreversible:
              reactions_to_make_irreversible.append(x)    
                  
    for reaction in reactions_to_make_irreversible:
        # Potential artifact because a reaction might run backwards naturally
        # and this would result in adding an empty reaction to the
        # model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            #reverse_reaction = copy.deepcopy(reaction)
            print( "adding reverse reaction for: ", reaction )
            print( 'upper and lower bounds:',reaction.upper_bound,reaction.lower_bound)
            reverse_reaction.gene_reaction_rule=reaction.gene_reaction_rule
            reverse_reaction.id = reaction.id + "_reverse"
            reverse_reaction.lower_bound = max(0,-1*reaction.upper_bound)
            reverse_reaction.upper_bound = reaction.lower_bound * -1.
            # reaction.lower_bound = 0    # Sergio: antic. Canvi ordre.
            if reaction.upper_bound<0:
               reaction.upper_bound=0
            reaction.lower_bound = 0      # Sergio: nou. Pq ara no deixa fer lower_bound=0 si el upper_bound és negatiu
               
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {}
            current_metabolites = [x for x in reaction.metabolites]
            for the_metabolite in current_metabolites:
                reaction_dict[the_metabolite] = -1 * reaction.get_coefficient(the_metabolite.id)
            reverse_reaction.add_metabolites(reaction_dict)
            reactions_to_add.append(reverse_reaction)
            # Also: GPRs should already copy
            # reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            # reverse_reaction._genes = reaction._genes#
            
            if mutually_exclusive_directionality_constraint and reaction.id not in reactions_with_no_indicators:
                # A continuous reaction bounded by 0., 1.
                # Serves as a source for the indicator metabolites
                tmp_source = Reaction('IRRMILP_direction_constraint_source_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_source.upper_bound = 1.
                tmp_source.lower_bound = 0.
                # The reverse indicator reaction is
                # an integer-valued reaction bounded by 0,1
                # that activates flux to the reverse reaction
                # and deactivates the forward reaction only when it is
                # turned on to 1
                tmp_indicator = Reaction('IRRMILP_reverse_indicator_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_indicator.upper_bound = 1
                tmp_indicator.lower_bound = 0
                tmp_indicator.variable_kind = 'integer'                    
                flux_constraint_forward = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reaction.id)
                flux_constraint_reverse = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reverse_reaction.id)
                flux_constraint_reverse._constraint_sense = 'G'
                flux_constraint_reverse._bound = 0.
                
                tmp_source.add_metabolites({flux_constraint_forward: 1})
                
                tmp_indicator.add_metabolites({flux_constraint_forward: -1,
                                      flux_constraint_reverse: 1})
                """NEW"""                      
                reverse_reaction.add_metabolites({flux_constraint_reverse: -1e-5})
                reaction.add_metabolites({flux_constraint_forward: -1e-5})
                """ENDNEW"""
                                      
                """OLD if reaction.upper_bound != 0:
                        reaction.add_metabolites({flux_constraint_forward: -1./reaction.upper_bound})
                else:
                    # could put 1.01 X the tolerance here,
                    # This is arbitrary.  Use 0.001
                    # since 1000 is a typical upper bound
                    reaction.add_metabolites({flux_constraint_forward: -0.00001})
                if reverse_reaction.upper_bound != 0:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -1./reverse_reaction.upper_bound})
                else:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -0.00001})"""
                reactions_to_add.append(tmp_indicator)
                reactions_to_add.append(tmp_source)
    cobra_model.add_reactions(reactions_to_add)