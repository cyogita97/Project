# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 19:04:57 2019

@author: yogit
"""
import numpy as np
from sympy import *
from z3 import *
import math
import time
start_time = time.time()
def CoefficientGenerator(B,v,x0,x1,pop,ca,ga,x,yy,u,w,nml,nsamp):
    pp=len(B)
    p=[]
    for j1 in range(pp):
        p.append(Real('p{}'.format(j1))) #Define coefficients for barrier polynomial p0*1+p1*x1+... (depending on the order of terms)
        
    s = Solver() #Initialise z3 solver
    
    
    aa=[]
    for i in range(0,len(v)):
        aaa=0
        for kk in range(pp): 
            aaa = aaa+v[i,kk]*p[kk]                  #aaa represnts the inequality with coefficients for each data sample in entire region     
        aa.append(aaa)
    A=z3.And([aa[m]>=0 for m in range(len(aa))])       # the inequalities should be satisfied for all data samples
    s.add(A)    #Condition 1 B(x)>0
   
    # =============================================================================
        
    a=[] #Condition 2: B(x)<=ga for initial condition
    for i in range(0,len(x0)):
        aaa=0
        for kk in range(pp):
            aaa = aaa+x0[i,kk]*p[kk]             #aaa represents the expression for each data sample in x0
        a.append(aaa)
    Z=z3.And([a[m]<=ga for m in range(len(a))])     # the inequalities should be satisfied for all X0 data samples
    s.add(Z)
    # =============================================================================
    b=[] #Condition 3 B(x)>=1 for unsafe condition
    for i in range(0,len(x1)):
        aaa=0
        for kk in range(pp):
            aaa = aaa+x1[i,kk]*p[kk]           #aaa represents the expression B(x) for each data sample in x1
        b.append(aaa)
    F=z3.And([b[m]>=1 for m in range(len(b))])    # the inequalities should be satisfied for all X1 data samples
    s.add(F)
    
    npo=nml
    for i in range(0,nsamp):
        d1=[]
        for k in range(npo):
            d=0
            for kk in range(pp):
                d=d+pop[i*npo+k,kk]*p[kk]     #d represents the expression E(B(f(x)))-B(x) for each input mode and each data sample in  x
            d1.append(d)
        e=z3.Or([d1[m]<=ca for m in range(len(d1))])   # the inequalities should be satisfied for atleast one of the input modes
        s.add(e)
    
            
    if s.check()==sat:
     #   print("Generating Barrier Polynomial","--- %s seconds ---" % (time.time() - start_time))   
        return p,s.check(),s.model()
    else:
        return p,s.check(),'model not available'