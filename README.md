# Macrophage-Metabolic-Phenotypes-in-Cytokine-Storms

Computational workflow developed to generate macrophages GSMM before and after treatment with cytokines 

## Description of main files
- 1_model_preparation.py
  The main purpose of this file is to build Condition-Specific models
  
- 2_sampling.py
  This main purpose of this file is to characterize the cell's metabolism by integrating transcriptomics data
  
-3_qMTA.py
  his main purpose of this file is to detect metabolic differences between two cell cultures

  Block_1_functions.py, Block_2_functions.py and Block_3functions.py contains required modules for the main scripts.
  input1.py, input2.py, input3.py describe the relatives paths of the inputs needed by the main scripts

## Requirements

cobrapy module
CPLEX
math package
pandas module
copy package
numpy package
scipy package
statsmodels module
openpyxl package
string package
six package
sympy module
