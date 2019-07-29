# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 18:30:30 2019

@author: yogit
"""
import numpy as np
from sympy import *
from z3 import *
import math
import time
def ConstraintGenerator(B,EqX0,EqX1,po,DataSamples,x,yy,u,w,nvar): #Substituting Data Samples in B(x) and E(B(f(x,u)))
    cc0=[]
    cc1=[]
    cc2=[]
    cc3=[]    
    Bo=B.copy()
    for ap in range(len(Bo)):
            Bo[ap]=Bo[ap].subs((yy[q],x[q]) for q in range(nvar))  # Replacing the temporary variables (yy) with true variables (x)
    
    
    for ii in range(len(DataSamples[0,:])): #Substitute Data Sample value in initial condition
        X0=np.zeros([len(EqX0),len(EqX0[0])])      # stores the evaluated intial region for each data sample
        for jj in range(len(EqX0)):
            for ll in range(len(EqX0[0])):
                X0[jj,ll]=EqX0[jj][ll].subs((x[tt],DataSamples[tt,ii]) for tt in range(nvar))   # the evaluation of initial region for each data sample and storing in X0
        
        X1=np.zeros([len(EqX1),len(EqX1[0])]) #Substitute Data sample value in unsafe condition
        for jj in range(len(EqX1)):
            for ll in range(len(EqX1[0])):
                 X1[jj,ll]=EqX1[jj][ll].subs((x[tt],DataSamples[tt,ii]) for tt in range(nvar))  # the evaluation of unsafe region for each data sample and storing in X1
                 
        if any(all(X0[mm][nn]<=1e-5 for nn in range(X0.shape[1])) for mm in range(X0.shape[0])):   #checking if the evaluated data sample lies in intial region or not
            BB=Matrix(Bo).T #for those data samples in initial,
            for aqq in range(len(BB)):
                for jjj in range(nvar):
                    BB[aqq]=BB[aqq].subs(x[jjj],DataSamples[jjj,ii]) #substitute Data sample in B(x)
                c1=list(BB)
            cc1.append(c1)    # cc1 stores the monomials of Barrier Polynomial substituted with data samples which lie in region X0
            
        elif any(all(X1[mm][nn]<=1e-5 for nn in range(X1.shape[1])) for mm in range(X1.shape[0])):    ##checking if the evaluated data sample lies in unsafe region or not
            BB=Matrix(Bo).T
            for aqq in range(len(BB)): #Same for unsafe
                for jjj in range(nvar):
                    BB[aqq]=BB[aqq].subs(x[jjj],DataSamples[jjj,ii])
                c2=list(BB)
            cc2.append(c2)     ## cc2 stores the monomials of Barrier Polynomial substituted with data samples which lie in region X1
     
    BBB=Matrix(Bo).T     #required because cannot subtract an entire list just by its name      
    
    for iii in range(len(DataSamples[0,:])): #For all data samples, the E(B(f(x,u)))-B(x) value
        for ss in range(len(po[:,0])):    #for each input mode
            VV=po[ss,:]-BBB
            for jj in range(nvar):
                VV=VV.subs(x[jj],DataSamples[jj,iii])
            c3=list(VV)
            cc3.append(c3)
    
    for iiii in range(len(DataSamples[0,:])): #B(x) for all data samples
        BBB=Matrix(Bo).T           #to reduce one for loop converted the list to matrix 
        for jj in range(nvar):
            BBB=BBB.subs(x[jj],DataSamples[jj,iiii])
        c0=list(BBB)
        cc0.append(c0)
        
    v=np.float64(cc0)        #substituted all Data samples in B(x)
    x0=np.float64(cc1)       #substituted Data samples lying in region X0 in B(x)
    x1=np.float64(cc2)       #substituted Data samples lying in region X1 in B(x)
    pop=np.float64(cc3)      #substituted all Data samples E(B(f(x)))-f(x)
    
    
    return v,x0,x1,pop #Return those values in an array. pop will have size nm*(no of barrier function terms)