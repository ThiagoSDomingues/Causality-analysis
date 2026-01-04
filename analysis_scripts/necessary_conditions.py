### Author: OptimusThi
import numpy as np


""" 
Necessary causality conditions 

[1] 

"""

# Nonlinear necessary conditions (6)

def n_1(C_eta, R_bulk, lambda_1, e, p):

    condition = 2/C_eta + (6/5)*R_bulk - (5/7)*np.abs(lambda_1)/(e + p)
    return condition >= 0 
    
def n_2(C_eta, R_bulk, lambda_3, e, p):
    
    condition = 1 - 1/C_eta + (2/5)*R_bulk - (5/14)*(lambda_3/(e + p))
    return condition >= 0    

def n_3(C_eta, R_bulk, lambda_3, e, p):
    
    condition = 1/C_eta + (3/5)*R_bulk  - (5/14)*(lambda_3)/(e + p)
    return condition >= 0    

def n_4(C_eta, R_bulk, lambda_a, lambda_d, e, p):
    
    condition = 1 - 1/C_eta + (2/5)*R_bulk + (9/14)*(lambda_a)/(e + p) - (5/14)*(lambda_d)/(e + p)   
    return condition >= 0

def n_5(C_eta, C_zeta, cs_2, R_bulk, lambda_1, e, p):
    
    condition = cs_2 + (4/3)*(1/C_eta) + 1/C_zeta + (22/15 + cs_2)*R_bulk + (38/21 + (8/5)*(1/3 - cs_2) + cs_2)*(lambda_1)/(e + p)
    return condition >= 0

def n_6(C_eta, C_zeta, cs_2, R_bulk, lambda_3, e, p):
    
    condition = 1 - (cs_2 + (4/3)*(1/C_eta) + 1/C_zeta) + (-7/15 - cs_2)*(R_bulk) + (-17/21 - (8/5)*(1/3 - cs_2) - cs_2)*(lambda_3)/(e + p)  
    return condition >= 0    
