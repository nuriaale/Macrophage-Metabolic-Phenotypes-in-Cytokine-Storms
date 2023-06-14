#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 17:49:40 2023

@author: nuria
"""

from cobra import Model, Reaction, Metabolite
from cplex import Cplex, SparsePair,SparseTriple
from six import iteritems, string_types
from sympy import Basic, Number
from LegacySolution import LegacySolution
import copy
import pandas as pd

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
            calc = pd.concat([mean, std, se], axis=1)
            calc.columns = ['mean', 'std', 'se']
            
            if actual.empty:
                actual = calc
            else:
                actual = pd.concat([actual, calc])
        conditions_dict[condition]=actual
        
    return conditions_dict



def add_sink_reactions_to_model(model, pellets_pd, exclude, biocrates_name_compartment_dict):
    rejected_list=[]
    met_sink_dict= {}
    for column_name in pellets_pd.columns:
        if column_name not in exclude:
            if column_name.lower() in biocrates_name_compartment_dict:
                metabolite_id=biocrates_name_compartment_dict[column_name.lower()]['met_id'] #Metabolite ID without compartment tag
                rid='SINK_'+metabolite_id
                met=column_name.lower()+'(pellet)'
                met_sink_dict[met]=rid
                if rid not in model.reactions:
                    for compartment in biocrates_name_compartment_dict[column_name.lower()]["compartment"]:
                        met_id=metabolite_id+"_"+compartment
                        if met_id in model.metabolites:
                            met_actual=model.metabolites.get_by_id(met_id)
                            reporter_metabolite=Metabolite("reporter_"+column_name.lower(), formula=met_actual.formula, name=met_actual.name, compartment=compartment)
                            sink_rid='SINK_'+met_id
                            sink_reaction=Reaction(sink_rid)
                            sink_reaction.name = 'sink of'+met_id
                            sink_reaction.subsystem = 'sink'
                            sink_reaction.add_metabolites({reporter_metabolite: 1.0,model.metabolites.get_by_id(met_id):-1})
                            model.add_reactions([sink_reaction])
                        
                    reaction = Reaction(rid)
                    reaction.name = 'sink of'+metabolite_id
                    reaction.subsystem = 'sink'
                    reaction.add_metabolites({reporter_metabolite: -1.0})   
                    model.add_reactions([reaction])
                else:#rid is in model reactions
                   reaction=model.reactions.get_by_id(rid)
            else:
             rejected_list.append(column_name)
    return met_sink_dict, rejected_list

def mfa(cobra_model,stat_dict,condition,metabolite_ex_dict,min_std=1e-4,precision=6,fname=None):
    coefficient_dict, quadaratic_dict,numeric_qp_dict,numeric_qp=set_cofficients_mfa(cobra_model,stat_dict,condition,metabolite_ex_dict,min_std,precision)
    mfa_problem=create_rMTA_problem_quadratic(cobra_model,quadaratic_dict, coefficient_dict,out_name="MFA.lp",verbose=False)
    mfa_problem.solve()
    solution=cplex_format_solution(mfa_problem, cobra_model)
    total_xi,xi_dict=get_xi2(stat_dict,condition,solution.x_dict,min_std=min_std,verbose=True,fname=fname)
    return solution.x_dict,  total_xi, xi_dict

def set_cofficients_mfa(cobra_model,stat_dict,condition,metabolite_ex_dict,min_std=1e-4,precision=6):
 numeric_qp_dict={}
 numeric_qp=0
 coefficient_dict={}
 quadaratic_dict={}
 #counter=0
 for met in stat_dict[condition].index:
     if met.lower() in metabolite_ex_dict:
         rid=metabolite_ex_dict[met.lower()]
         if rid in cobra_model.reactions:
             #counter+=1
             #print(rid, counter, )
             try:
                 target_vref=stat_dict[condition]['mean'][met]
                 std=max(stat_dict[condition]['std'][met],min_std) 
                 weight=pow(std,-2)
                 weight=round_sig(weight,precision)
                 if rid not in coefficient_dict or rid not in quadaratic_dict:
                    coefficient_dict[rid]=0
                    quadaratic_dict[rid]=0
                 
                 coefficient_dict[rid]+=-2*target_vref*weight#-2ab*weight
                 quadaratic_dict[rid]+=weight*2 #b2 Multipy by 2
                 numeric_qp_dict[rid]=pow(target_vref,2)*weight
                 numeric_qp+=pow(target_vref,2)*weight
                 #print(rid)
                
             except:
                 print('error') 
 return  coefficient_dict, quadaratic_dict,numeric_qp_dict,numeric_qp

parameter_defaults = {'objective_sense': 'maximize',
                      'tolerance_optimality': 1e-9,
                      'tolerance_feasibility': 1e-9,
                      'tolerance_integer': 1e-9,
                      'lp_method': 1,
                      'tolerance_barrier': 1e-9,
                      'verbose': False,
                      'qpmethod': 1}

def create_rMTA_problem_quadratic(cobra_model,quadaratic_dict, coefficient_dict,out_name="qrMTA.lp",**kwargs ):
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)
        print( the_parameters  )
    if 'relax_b' in the_parameters:
        relax_b = the_parameters.pop("relax_b")
        warn('need to reimplement relax_b')
        relax_b = False
    else:
        relax_b = False
    
    lp = Cplex()
    
    
    for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)
    
    
    objective_coefficients=[]
    for x in cobra_model.reactions:
        if x.id in coefficient_dict:
           obj_coefficent=coefficient_dict[x.id]
        else:
           obj_coefficent=0
        objective_coefficients.append(obj_coefficent)

    lower_bounds = [_float(x.lower_bound) for x in cobra_model.reactions]
    upper_bounds = [_float(x.upper_bound) for x in cobra_model.reactions]
    variable_names = cobra_model.reactions.list_attr("id")
    
    variable_kinds = [variable_kind_dict['continuous'] for x   
                      in cobra_model.reactions]
    
    lp.variables.add(obj=objective_coefficients, lb=lower_bounds, ub=upper_bounds, names=variable_names, types=variable_kinds)
    constraint_sense = []; constraint_names = []; constraint_limits = []
    for x in cobra_model.metabolites:
        if "MOMA_wt_" in x.id:
            continue
        constraint_sense.append('E')   #  E de Equal. Per a tots.
        constraint_names.append(x.id)
        constraint_limits.append(float(x._bound))
     
    the_linear_expressions = []
    for the_metabolite in cobra_model.metabolites:
        if "MOMA_wt_" in the_metabolite.id:
            continue
        variable_list = []
        coefficient_list = []
        for the_reaction in the_metabolite._reaction:
            variable_list.append(str(the_reaction.id))
            coefficient_list.append(_float(the_reaction._metabolites[the_metabolite]))
        the_linear_expressions.append(SparsePair(ind=variable_list,val=coefficient_list))
    for n,x in enumerate(cobra_model.reactions):
        if x.id in quadaratic_dict:
           coef=quadaratic_dict[x.id]
        else:
            coef=0
        lp.objective.set_quadratic_coefficients(n, n, coef)
    if relax_b:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,rhs=constraint_limits, range_values=list(range_values), senses=constraint_sense, names=constraint_names)
                                  
    else:
        lp.linear_constraints.add(lin_expr=the_linear_expressions,rhs=constraint_limits, senses=constraint_sense, names=constraint_names)
    
    lp.set_problem_type(Cplex.problem_type.QP)
    lp.objective.set_sense(1) 
    set_parameter(lp,"verbose",5)
    if out_name not in ("",None) and False:
       pass 
    return(lp)

def cplex_format_solution(lp, cobra_model, **kwargs):
    from cplex import Cplex, SparsePair
    from cplex.exceptions import CplexError
    status = lp.solution.get_status_string().lower()
    if status in ('optimal', 'time_limit', 'non-optimal',"integer optimal solution","integer optimal, tolerance"):
        objective_value = lp.solution.get_objective_value()
        x_dict = dict(zip(lp.variables.get_names(),
                          lp.solution.get_values()))

        
        x = lp.solution.get_values()
        if lp.get_problem_type() in (Cplex.problem_type.MIQP, Cplex.problem_type.MILP):
            y = y_dict = None
        else:
            y_dict = dict(zip(lp.linear_constraints.get_names(),lp.solution.get_dual_values()))
            y = lp.solution.get_dual_values()
    else:
        x = y = x_dict = y_dict = objective_value = None
    return LegacySolution(objective_value, x=x, x_dict=x_dict, status=status,y=y, y_dict=y_dict)

def round_sig(x, sig=2):
  if x==0:
    value=0
  else:
     value=round(x, sig-int(math.floor(math.log10(abs(x))))-1)
  return value

def get_xi2(stat_dict,condition,x_dict,min_std=1e-6,verbose=True,fname=None):
  total_xi=0
  xi_dict={}
  if stat_dict is not None and stat_dict!={}:
      output_pd=pd.DataFrame(columns=["met","rid", "simulated flux","target_flux","sd", "Xi2"])
      for met in stat_dict[condition].index:
          if met.lower() in metabolite_ex_dict:
              rid=metabolite_ex_dict[met.lower()]
              if rid not in x_dict:
                 out_str=rid+"("+met+")"+" not in model"
            
                 row=[out_str, '--', '--', '--', '--', '--']
            
                 output_pd=output_pd.append(pd.Series(row, index=output_pd.columns), ignore_index=True)
                 continue 
              vres=x_dict[rid]
              vref=stat_dict[condition]['mean'][met]
              sd=max(stat_dict[condition]['std'][met],min_std)
              xi2= pow((vref-vres)/sd,2)
              xi_dict[rid]=xi2
              total_xi+=xi2
              if fname!=None:
                  row=[met,rid,round(vres,6), round(vref,6),round(sd,6), round(xi2,3)]
            #rows.append() 
                  output_pd=output_pd.append(pd.Series(row, index=output_pd.columns), ignore_index=True)
          else: 
              out_str=met+" not in metabolite_ex_dict"
              row=[out_str, '--', '--', '--', '--', '--']

              output_pd=output_pd.append(pd.Series(row, index=output_pd.columns), ignore_index=True)
      if fname!=None:
          output_pd.to_excel(fname)
      
  return total_xi,xi_dict  

def modify_model_for_seahorse(model,remove=False):

   if "CYOOm_token_atpsyn" not in model.metabolites and "ATPS4mi_token" not in model.metabolites and "NADH2_u10mi_token" not in model.metabolites: 
      CYOOm_token_atpsyn=Metabolite("CYOOm_token_atpsyn"); CYOOm_token_atpsyn.name="CYOOm_token_atpsyn"
      CYOOm_token_atpsyn.compartment="c" 
      ATPS4mi_token=Metabolite("ATPS4mi_token"); ATPS4mi_token.name="ATPS4mi_token"
      ATPS4mi_token.compartment="c"
      NADH2_u10mi_token=Metabolite("NADH2_u10mi_token"); NADH2_u10mi_token.name="NADH2_u10mi_token"
      NADH2_u10mi_token.compartment="c"
      fadhdh_like_token=Metabolite("fadhdh_like_token"); fadhdh_like_token.name="fadhdh_like_token"
      fadhdh_like_token.compartment="c"
      sulfox_like_token=Metabolite("sulfox_like_token"); sulfox_like_token.name="sulfox_like_token"
      sulfox_like_token.compartment="c"
      
      model.add_metabolites([CYOOm_token_atpsyn,ATPS4mi_token])

      #Nadhdh
      model.reactions.NADH2_u10mi.add_metabolites({NADH2_u10mi_token:1}) #	5.0 h_m + nadh_m + q10_m --> 4.0 h_i + nad_m + q10h2_m + NADH2_u10mi_token
      RGROUP_non_atpsyn_NADH2_u10mi=reaction_from_string(model,"RGROUP_non_atpsyn_NADH2_u10mi","NADH2_u10mi_token ->",bounds=None,gpr="")
      RGROUP_atpsyn_NADH2_u10mi=reaction_from_string(model,"RGROUP_atpsyn_NADH2_u10mi","NADH2_u10mi_token -> 0.5 CYOOm_token_atpsyn + 2.5 ATPS4mi_token ",bounds=None,gpr="")
      
      FADH_Like_list=[]
      for reaction in model.metabolites.get_by_id("q10_m").reactions:
          if reaction.metabolites[model.metabolites.get_by_id("q10_m")]==-1:
             if reaction.metabolites.get(model.metabolites.get_by_id("q10h2_m"))==1:
                 if model.metabolites.get_by_id("h_i") not in reaction.metabolites:

                    FADH_Like_list.append(reaction)
      
      for reaction in FADH_Like_list:
          reaction.add_metabolites({fadhdh_like_token:1}) #  etfrd_m + q10_m --> etfox_m + q10h2_m + ETFQO_token
      
      RGROUP_non_atpsyn_fadhdh_like=reaction_from_string(model,"RGROUP_non_atpsyn_fadhdh_like","fadhdh_like_token ->",bounds=None,gpr="") 
      RGROUP_atpsyn_fadhdh_like=reaction_from_string(model,"RGROUP_atpsyn_fadhdh_like","fadhdh_like_token -> 0.5 CYOOm_token_atpsyn + 1.5 ATPS4mi_token ",bounds=None,gpr="")
      ###sulfox Like
      SULFOX_Like_list=[]
      for reaction in model.metabolites.get_by_id("ficytC_m").reactions:
          if reaction.metabolites[model.metabolites.get_by_id("ficytC_m")]==-2:
             if reaction.metabolites[model.metabolites.get_by_id("focytC_m")]==2:
                 if model.metabolites.get_by_id("h_i") not in reaction.metabolites:
                    #print( reaction.id, reaction.reaction  )
                    SULFOX_Like_list.append(reaction)
      
      for reaction in SULFOX_Like_list:
          reaction.add_metabolites({sulfox_like_token:1}) #  etfrd_m + q10_m --> etfox_m + q10h2_m + ETFQO_token
        
      RGROUP_non_atpsyn_sulfox_like=reaction_from_string(model,"RGROUP_non_atpsyn_sulfox_like","sulfox_like_token ->",bounds=None,gpr="") 
      RGROUP_atpsyn_sulfox_like=reaction_from_string(model,"RGROUP_atpsyn_sulfox_like","sulfox_like_token -> 0.5 CYOOm_token_atpsyn + 0.5 ATPS4mi_token ",bounds=None,gpr="")
      #ATP synthase
      model.reactions.ATPS4mi.add_metabolites({ATPS4mi_token:-1}) 
        
      #O2 linked to ATPsynthase
      RGROUP_CYOOm_atpsyn=reaction_from_string(model,"RGROUP_CYOOm_atpsyn","CYOOm_token_atpsyn -> ",bounds=None,gpr="")
      RGROUP_CYOOm3i_CYOOm2i=define_reaction_group(model,{"CYOOm3i":1,"CYOOm2i":1},group_reaction_id="RGROUP_CYOOm3i_CYOOm2i",lower_bound=None,upper_bound=None,objective_coefficient=0)
      #identfy proton transport reactions: 
      proton_transport_dict={}
      for reaction in model.metabolites.get_by_id("h_e").reactions:
          if model.metabolites.get_by_id("h_c") in reaction.metabolites: 
             factor=reaction.metabolites[model.metabolites.get_by_id("h_e")]

             proton_transport_dict[reaction.id]=factor
      
      RGROUP_h_transport=define_reaction_group(model,proton_transport_dict,group_reaction_id="RGROUP_h_transport",lower_bound=None,upper_bound=None,objective_coefficient=0)
   if remove:
       #Mets
       CYOOm_token_atpsyn=model.metabolites.CYOOm_token_atpsyn
       ATPS4mi_token=model.metabolites.ATPS4mi_token
       NADH2_u10mi_token=model.metabolites.NADH2_u10mi_token
       fadhdh_like_token=model.metabolites.fadhdh_like_token
       sulfox_like_token=model.metabolites.sulfox_like_token
       mRGROUP_CYOOm3i_CYOOm2i=model.metabolites.mRGROUP_CYOOm3i_CYOOm2i
       mRGROUP_h_transport=model.metabolites.mRGROUP_h_transport
       
       mets2remove=[CYOOm_token_atpsyn,ATPS4mi_token,NADH2_u10mi_token,fadhdh_like_token,sulfox_like_token,mRGROUP_CYOOm3i_CYOOm2i,mRGROUP_h_transport]

       for met in mets2remove:
           met.remove_from_model(destructive=False)
       #Reactions
       reaction_id=["RGROUP_non_atpsyn_NADH2_u10mi","RGROUP_atpsyn_NADH2_u10mi","RGROUP_non_atpsyn_fadhdh_like","RGROUP_atpsyn_fadhdh_like","RGROUP_non_atpsyn_sulfox_like","RGROUP_atpsyn_sulfox_like","RGROUP_CYOOm_atpsyn","RGROUP_h_transport","RGROUP_CYOOm3i_CYOOm2i"]

       for rid in  reaction_id:
           reaction=model.reactions.get_by_id(rid)
           reaction.remove_from_model()  
           
def reaction_from_string(model,rid,r_string,bounds=None,gpr=""):
     if rid in model.reactions:
        new_reaction=model.reactions.get_by_id(rid)
     else:    
       new_reaction=Reaction(str(rid))
       reaction_list_provisional=[new_reaction]
       model.add_reactions(reaction_list_provisional)
       new_reaction.build_reaction_from_string(r_string)
     if bounds!=None:
        new_reaction.bounds=bounds
     if gpr!=None:
        new_reaction.gene_reaction_rule=gpr
     return new_reaction
 
def define_reaction_group(model,reaction_dict,group_reaction_id=None,lower_bound=None,upper_bound=None,objective_coefficient=0):
    new_reaction_id="RGROUP"
    if group_reaction_id!=None:
       if "RGROUP" in group_reaction_id:
           new_reaction_id=group_reaction_id
       else:
           new_reaction_id+="_"+group_reaction_id
    else:
       for reaction_id in reaction_dict:
           new_reaction_id+="_"+reaction_id
    if new_reaction_id in model.reactions:
       model.reactions.get_by_id(new_reaction_id).remove_from_model()
    new_reaction_name="Reaction Group:"
    for reaction_id in reaction_dict:
        if  reaction_dict[reaction_id]>0:
            new_reaction_name+="+"+reaction_id
        else:
            new_reaction_name+="-"+reaction_id
    metabolite = Metabolite("m"+new_reaction_id,formula='',name="mGROUP"+new_reaction_id,compartment='gr')
    group_reaction = Reaction(new_reaction_id)
    group_reaction.name = new_reaction_name
    group_reaction.subsystem = 'Reaction group'
    if upper_bound!=None:
       group_reaction.upper_bound=upper_bound
    group_reaction.add_metabolites({metabolite:-1})
    if objective_coefficient==None:
        group_reaction.objective_coefficient=0
    reaction_list_provisional=[group_reaction]
    model.add_reactions(reaction_list_provisional)
    #model.add_reaction(group_reaction)
    group_reaction.objective_coefficient=objective_coefficient
    theoretical_lower_bound=0
    theoretical_upper_bound=0
    for reaction_id in reaction_dict:
        coef=reaction_dict[reaction_id]
        reaction=model.reactions.get_by_id(reaction_id)
        reaction.add_metabolites({metabolite:coef})
        if coef>=0:
           theoretical_upper_bound+=reaction.upper_bound
           theoretical_lower_bound+=reaction.lower_bound
        else:
           theoretical_upper_bound-=reaction.lower_bound
           theoretical_lower_bound-=reaction.upper_bound
    if lower_bound==None:
        group_reaction.lower_bound=min(round_down(theoretical_lower_bound,2),0)
    else:
        group_reaction.lower_bound=lower_bound
    if upper_bound==None:
        group_reaction.upper_bound=max(round_up(theoretical_upper_bound,2),1000)
    else:
        group_reaction.upper_bound=upper_bound
    return group_reaction

import math

def round_up(number,positions):
    exponent=pow(10,positions)
    new_number=math.ceil(number*exponent)/exponent
    return new_number


def round_down(number,positions):
    if number==0.0:
       return 0
    exponent=pow(10,positions)
    return math.floor(number*exponent-0.0001)/exponent


biocrates_name_compartment_dict={u'ala': {'compartment': ['c', 'x', 'm'], 'met_id': u'ala__L'},
 'alpha-aaa': {'compartment': ['c', 'm'], 'met_id': 'L2aadp'},
 u'arg': {'compartment': ['c', 'm'], 'met_id': u'arg__L'},
 u'asn': {'compartment': ['c', 'm'], 'met_id': u'asn__L'},
 u'asp': {'compartment': ['c', 'm'], 'met_id': u'asp__L'},
 'c0': {'compartment': ['c', 'x', 'm'], 'met_id': 'crn'},
 'c10': {'compartment': ['c', 'x'], 'met_id': 'c10dc'},
 'c12': {'compartment': ['c', 'x'], 'met_id': 'c12dc'},
 'c14': {'compartment': ['c', 'm'], 'met_id': 'ttdcrn'},
 'c16': {'compartment': ['c', 'x', 'm'], 'met_id': 'pmtcrn'},
 'c16:1': {'compartment': ['c', 'm'], 'met_id': 'hdcecrn'},
 'c18': {'compartment': ['c', 'm'], 'met_id': 'stcrn'},
 'c18:1': {'compartment': ['c', 'm'], 'met_id': 'odecrn'},
 'c2': {'compartment': ['c', 'x', 'm'], 'met_id': 'acrn'},
 'c3': {'compartment': ['c', 'x', 'm'], 'met_id': 'pcrn'},
 'c4': {'compartment': ['c', 'x', 'm'], 'met_id': 'c4crn'},
 'c6': {'compartment': ['c', 'x'], 'met_id': 'c6crn'},
 'c8': {'compartment': ['c', 'x', 'm'], 'met_id': 'c8crn'},
 u'carn': {'compartment': ['c'], 'met_id': u'carn'},
 'carnosine': {'compartment': ['c'], 'met_id': 'carn'},
 u'cit': {'compartment': ['c', 'm'], 'met_id': u'citr__L'},
 u'creatinine': {'compartment': ['c', 'm'], 'met_id': 'creat'},
 'dopa': {'compartment': ['c'], 'met_id': '34dhphe'},
 'dopamine': {'compartment': ['c'], 'met_id': 'dopa'},
 u'gln': {'compartment': ['c', 'm'], 'met_id': u'gln__L'},
 u'glu': {'compartment': ['c', 'm'], 'met_id': u'glu__L'},
 u'gly': {'compartment': ['c', 'x', 'm'], 'met_id': u'gly'},
 u'his': {'compartment': ['c', 'm'], 'met_id': u'his__L'},
 'histamine': {'compartment': ['c'], 'met_id': 'hista'},
 u'ile': {'compartment': ['c', 'm'], 'met_id': u'ile__L'},
 'kynurenine': {'compartment': ['c'], 'met_id': 'Lkynr'},
 u'leu': {'compartment': ['c', 'm'], 'met_id': u'leu__L'},
 u'lys': {'compartment': ['c', 'x', 'm'], 'met_id': u'lys__L'},
 u'met': {'compartment': ['c', 'm'], 'met_id': u'met__L'},
 u'orn': {'compartment': ['c', 'm'], 'met_id': u'orn'},
 u'phe': {'compartment': ['c', 'm'], 'met_id': u'phe__L'},
 u'pro': {'compartment': ['c', 'm'], 'met_id': u'pro__L'},
 'putrescine': {'compartment': ['c', 'm'], 'met_id': 'ptrc'},
 u'ser': {'compartment': ['c', 'x', 'm'], 'met_id': u'ser__L'},
 'serotonin': {'compartment': ['c'], 'met_id': 'srtn'},
 u'spermidine': {'compartment': ['c'], 'met_id': u'spmd'},
 u'spermine': {'compartment': ['c'], 'met_id': u'sprm'},
 u'srtn': {'compartment': ['c'], 'met_id': u'srtn'},
 't4-oh-pro': {'compartment': ['c', 'm'], 'met_id': '4hpro_LT'},
 'taurine': {'compartment': ['c', 'x'], 'met_id': 'taur'},
 u'thr': {'compartment': ['c', 'm'], 'met_id': u'thr__L'},
 u'trp': {'compartment': ['c'], 'met_id': u'trp__L'},
 u'tyr': {'compartment': ['c', 'm'], 'met_id': u'tyr__L'},
 u'val': {'compartment': ['c', 'm'], 'met_id': u'val__L'}}



metabolite_ex_dict={u'ala': u'EX_ala__L_e',
 'alpha-aaa': 'EX_L2aadp_e',
 u'arg': u'EX_arg__L_e',
 u'asn': u'EX_asn__L_e',
 u'asp': u'EX_asp__L_e',
 'c0': 'EX_crn_e',
 u'carn': u'EX_carn_e',
 'carnosine': 'EX_carn_e',
 u'cit': u'EX_citr__L_e',
 u'creatinine': 'EX_crtn_e',
 u'gln': u'EX_gln__L_e',
 u'glu': u'EX_glu__L_e',
 u'gly': u'EX_gly_e',
 u'his': u'EX_his__L_e',
 'histamine': 'EX_hista_e',
 u'ile': u'EX_ile__L_e',
 'kynurenine': 'EX_Lkynr_e',
 u'leu': u'EX_leu__L_e',
 u'lys': u'EX_lys__L_e',
 u'met': u'EX_met__L_e',
 u'orn': u'EX_orn_e',
 u'phe': u'EX_phe__L_e',
 u'pro': u'EX_pro__L_e',
 'putrescine': 'EX_ptrc_e',
 u'ser': u'EX_ser__L_e',
 'serotonin': 'EX_srtn_e',
 u'spermidine': u'EX_spmd_e',
 u'spermine': u'EX_sprm_e',
 u'srtn': u'EX_srtn_e',
 't4-oh-pro': 'EX_4hpro_LT_e',
 'taurine': 'EX_taur_e',
 u'thr': u'EX_thr__L_e',
 u'trp': u'EX_trp__L_e',
 u'tyr': u'EX_tyr__L_e',
 u'val': u'EX_val__L_e',
 "glc": "EX_glc__D_e",
 "glucose":"EX_glc__D_e",
 "lac":"EX_lac__L_e",
 "lactate":"EX_lac__L_e",
 "per":"RGROUP_h_transport",
 "ocr":"EX_o2_e",
 "basal respiration":"RGROUP_CYOOm3i_CYOOm2i",
 "atp production":"RGROUP_CYOOm_atpsyn"
 }

def set_parameter(lp, parameter_name, parameter_value):
    if parameter_name == 'objective_sense':
        parameter_value = getattr(lp.objective.sense, parameter_value)
        lp.objective.set_sense(parameter_value)
        return
    elif parameter_name == 'the_problem':
        warn('option the_problem removed')
        return
    elif parameter_name == 'verbose':
        return
      
    try:
        cplex_name = parameter_mappings.get(parameter_name, parameter_name)
        cplex_value = parameter_value
        #print( cplex_name, cplex_value )
        param = lp.parameters
        for i in cplex_name.split("."):
            param = getattr(param, i)
        if isinstance(cplex_value, string_types) and \
                hasattr(param.values, cplex_value):
            cplex_value = getattr(param.values, cplex_value)
        param.set(cplex_value)
    except (CplexError, AttributeError) as e:
        raise ValueError("Failed to set %s to %s: %s" %
                         (parameter_name, str(parameter_value), repr(e)))
        
        
parameter_mappings = {'lp_method': 'lpmethod',
                      'lp_parallel': 'threads',
                      'qpmethod':'qpmethod',
                      'threads': 'threads',
                      'objective_sense': 'objective_sense',
                      'time_limit': 'timelimit',
                      'iteration_limit': 'simplex.limits.iterations',
                      'tolerance_barrier': 'barrier.convergetol',
                      'tolerance_feasibility': 'simplex.tolerances.feasibility',
                      'tolerance_markowitz': 'simplex.tolerances.markowitz',
                      'tolerance_optimality': 'simplex.tolerances.optimality',
                      'tolerance_integer': 'mip.tolerances.integrality',
                      'MIP_gap_abs': 'mip.tolerances.absmipgap',
                      'MIP_gap': 'mip.tolerances.mipgap'}

variable_kind_dict ={'continuous': "C", 'integer': "I"}

def _float(value):
    if isinstance(value, Basic) and not isinstance(value, Number):
        return 0.
    else:
        return float(value)