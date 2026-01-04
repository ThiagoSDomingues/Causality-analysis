### Author: OptimusThi
# _*_ coding: utf-8 _*_

import numpy as np

"""
Sufficient conditions to perform a causality analysis in heavy-ion collisions simulations based on DNMR equations of motion 

References: 

[1] Cheng Chiu, and Chun Shen 
    Exploring theoretical uncertainties in the hydrodynamic description of relativistic
    heavy-ion collisions, Phys. Rev. C 103, 064901 (2021)     
"""

def s_1(C_eta, R_bulk, lambda_1, lambda_3, e, p):

    condition = 1/C_eta - np.abs(lambda_1)/(e+p) + (2/5)*(R_bulk) - (5/7)*(lambda_3/(e+p))   
    return condition >= 0 
    
def s_2(C_eta,  R_bulk, lambda_1, e, p):

    condition = 1/C_eta + (3/5)*R_bulk - (5/7)*np.abs(lambda_1)/(e+p) 
    return condition >= 0
       
# s_3 and s_4 are always satisfied !   

def s_5(C_eta, C_zeta, cs2, R_bulk, lambda_1, lambda_3, e, p):
    
    condition = (1 + R_bulk)*(1 - cs2) - (4/(3*C_eta) + 1/C_zeta + (22/15)*R_bulk + (38/21 + (8/5)*(1/3 - cs2) + cs2)*(lambda_3)/(e+p) + np.abs(lambda_1)/(e+p) + (17/14)*((8/5)*(1/3 - cs2) + cs2 - 5/42)*((lambda_3/(e+p) + np.abs(lambda_1)/(e+p))**2)/(1-1/C_eta + (2/5)*R_bulk - np.abs(lambda_1)/(e+p) - (5/7)*lambda_3/(e+p)))  
    return condition >= 0 
    
def s_6(C_eta, C_zeta, cs2, R_bulk, lambda_1, e, p):
    
    condition = 1/(3*C_eta) + 1/C_zeta + cs2 + (13/15 + cs2)*(R_bulk) + (-23/21 + (8/5)*(1/3 - cs2) - cs2)*np.abs(lambda_1)/(e+p)     
    return condition >= 0     
    
def s_7(C_eta, cs2, R_bulk, lambda_1, lambda_3, e, p):

    condition = ( 1/C_eta + (3/5)*(R_bulk) - (5/7)*np.abs(lambda_1)/(e+p) )**2 - (17/14)*((8/5)*(1/3 - cs2) + cs2 - 5/42)*(lambda_3/(e+p) + np.abs(lambda_1)/(e+p))**2     
    return condition >= 0 
    
def s_8(C_eta, C_zeta, cs2, R_bulk, lambda_1, lambda_2, lambda_3, e, p):

    condition = (4/3)*(1/C_eta) + 1/C_zeta + cs2 + (22/15 + cs2)*R_bulk - (22/21 - (8/5)*(1/3 - cs2) + cs2)*np.abs(lambda_1)/(e+p) - (1 + R_bulk + lambda_2/(e+p))*(1 + R_bulk + lambda_3/(e+p))/(3*(1 + R_bulk - lambda_1/(e+p))**2) 
    return condition >= 0                 


