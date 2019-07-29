# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 11:17:48 2019

@author: yogit
"""
import numpy as np
from sympy import *
from z3 import *
import math
import time
start_time = time.time()
from scipy.optimize import minimize
import scipy as sp
from tkinter import *
import tkinter.filedialog
import tkinter.messagebox
import tkinter.simpledialog
from platypus.problems import Problem
from platypus.algorithms import NSGAII
from platypus.problems import DTLZ2
import platypus.types as tt
import matplotlib.pyplot as plt
import Polygenerator
import ConstrGenerator
import CoeffGenerator
import CountExample
def discrete(nsamples,constant,gamma,EqX0,EqX1,funcn,B,nm,nu,nvar,nw,inp,bounds,opti,analy,x,yy,w,u,dyn):
 #   print("this is nm",nm)
    nsamp=nsamples
    DataSamples=np.zeros([nvar,nsamp])
    function=funcn
    for sd in range(nvar):
        c=(bounds[sd][1]-bounds[sd][0])*np.random.rand(1,nsamp)
    if dyn==1:  # if dynamics of system is in discrete form
       (po,nml,nul,inpt)=Polygenerator.PolynomialGenDisc(B,function,analy,nm,nu,inp,x,yy,w,u,nvar,nw) #Generating Expectation Polynomial
    else:   # if dynamics of system is in continuous form
       (po,nml,nul,inpt)=Polygenerator.PolynomialGenCont(B,function,analy,nm,nu,inp,x,yy,w,u,nvar,nw) 
    ca=constant #Fix the value of c
    ga=gamma#Fix the value of gamma
        
    cap=ca
    gam=ga    
    
    #We need to change gamma and c value alternatively. First start with ca, if the model exists for a particular value 
    #ca without counterexamples, subsequently change 'turn' to false and reduce ga and continue this alternatively
    #until an unsatisfiability
    #for a particular value of ca or ga is reached. Then fix with corresponding value (of ca or ga) and keep reducing the 
    #other parameter (ca or ga) until unsatisfiability is reached or a minimum value of the parameter is reached.
    turn=True
    it=sat
    while ca>=1e-4 and ga>=1e-4:
        while turn and it==sat:
        #    print(ca,ga)
            (v,x0,x1,pop)=ConstrGenerator.ConstraintGenerator(B,EqX0,EqX1,po,DataSamples,x,yy,u,w,nvar) #Generating constraints for constructing barrier polynomial
       #     print("Generating constraints", "--- %s seconds ---" % (time.time() - start_time)) 
            (p,s,m)=CoeffGenerator.CoefficientGenerator(B,v,x0,x1,pop,ca,ga,x,yy,u,w,nml,nsamp) #Generating barrier polynomial coefficients
            #print(m)
            if s==sat:
                if opti==1:
                   (model,cex)=CountExample.CounterExampleSLSQP(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga)
                else:
                   (model,cex)=CountExample.CounterExamplePlatypus(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga,nml)
      #          print(cex)
     #           print("Generating counterexamples", "--- %s seconds ---" % (time.time() - start_time))   
                if cex.size!=0:
                    DataSamples=np.hstack([DataSamples,cex])
                    #print(DataSamples)
                    nsamp=len(DataSamples[:][0])
                else:
                    cap=ca
                    mod=model
                    ga=ga/2
                    turn=False
                    break     
            if s==unsat:
                model=mod
                ga=ga/2
                turn=False
                it=unsat
                ca=cap           
                break      
        
        while turn and it==unsat:
       #     print(ca,ga)
            (v,x0,x1,pop)=ConstrGenerator.ConstraintGenerator(B,EqX0,EqX1,po,DataSamples,x,yy,u,w,nvar) #Generating constraints for constructing barrier polynomial
      #      print("Generating constraints", "--- %s seconds ---" % (time.time() - start_time)) 
            (p,s,m)=CoeffGenerator.CoefficientGenerator(B,v,x0,x1,pop,ca,ga,x,yy,u,w,nml,nsamp) #Generating barrier polynomial coefficients
            #print(m)
            if s==sat:
                if opti==1:
                   (model,cex)=CountExample.CounterExampleSLSQP(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga)
                else:
                   (model,cex)=CountExample.CounterExamplePlatypus(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga,nml)
      #          print(cex)
      #          print("Generating counterexamples", "--- %s seconds ---" % (time.time() - start_time))   
                cap=ca
                if cex.size!=0:
                    DataSamples=np.hstack([DataSamples,cex])
                    #print(DataSamples)
                    nsamp=len(DataSamples[:][0])
                else:
                    ca=ca/2
                    break     
            if s==unsat:
                ca=cap        
                break
            
            
        while not turn and it==sat:
     #       print(ca,ga)
            (v,x0,x1,pop)=ConstrGenerator.ConstraintGenerator(B,EqX0,EqX1,po,DataSamples,x,yy,u,w,nvar) #Generating constraints for constructing barrier polynomial
        #    print("Generating constraints", "--- %s seconds ---" % (time.time() - start_time)) 
            (p,s,m)=CoeffGenerator.CoefficientGenerator(B,v,x0,x1,pop,ca,ga,x,yy,u,w,nml,nsamp) #Generating barrier polynomial coefficients
            #print(m)
            if s==sat:
                if opti==1:
                   (model,cex)=CountExample.CounterExampleSLSQP(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga)
                else:
                   (model,cex)=CountExample.CounterExamplePlatypus(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga,nml) 
      #          print(cex)
             #   print("Generating counterexamples", "--- %s seconds ---" % (time.time() - start_time))   
                if cex.size!=0:
                    DataSamples=np.hstack([DataSamples,cex])
                    #print(DataSamples)
                    nsamp=len(DataSamples[:][0])
                else:
                    gam=ga
                    mod=model
                    ca=ca/2
                    turn=True
                    break  
            if s==unsat:
                turn=True
                it=unsat
                ca=ca/2
                ga=gam           
                break     
        
        while not turn and it==unsat:
    #        print(ca,ga)
            (v,x0,x1,pop)=ConstrGenerator.ConstraintGenerator(B,EqX0,EqX1,po,DataSamples,x,yy,u,w,nvar) #Generating constraints for constructing barrier polynomial
    #        print("Generating constraints", "--- %s seconds ---" % (time.time() - start_time)) 
            (p,s,m)=CoeffGenerator.CoefficientGenerator(B,v,x0,x1,pop,ca,ga,x,yy,u,w,nml,nsamp) #Generating barrier polynomial coefficients
            #print(m)
            if s==sat:
                if opti==1:
                   (model,cex)=CountExample.CounterExampleSLSQP(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga)
                else:
                   (model,cex)=CountExample.CounterExamplePlatypus(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga,nml)
     #           print(cex)
      #          print("Generating counterexamples", "--- %s seconds ---" % (time.time() - start_time))   
                if cex.size!=0:
                    DataSamples=np.hstack([DataSamples,cex])
                    #print(DataSamples)
                    nsamp=len(DataSamples[:][0])
                else:
                    gam=ga
                    ga=ga/2
                    break  
            if s==unsat:
                ga=gam          
                break
        
        if s==unsat:
            break
    print("value of c is ",cap)
    print("value of gamma is ", gam)
    print("Final Model is ",model)
    
    #This is for printing the input value according to the value of the polynomial (4th condition in z3 solver)
        
    for ap in range(len(B)):
        B[ap]=B[ap].subs((yy[q],x[q]) for q in range(nvar))
    barrier=sum(np.multiply(B,model))
    
    poNew1=[]
    for sv in range(len(po[:,0])):
        po11=sum(list(np.multiply(list(po[sv,:]),model)))
        poNew1.append(po11)
        
    rows=nul+1
    Data=[[0] * rows for i in range(len(poNew1))]
    for sc in range(len(poNew1)):
        Data[sc][0]=str(poNew1[sc]-barrier)+"<"+str(cap)
        for sz in range(1,nul+1):
          Data[sc][sz]='input'+str(sz-1)+':'+str(inpt[sz-1][sc])
    hj=[0 for i in range(len(poNew1))]
    for i in range(len(poNew1)):
        hj[i]=poNew1[i]-barrier
  #  for kl in range(len(Data)):
   #     for jk in range(len(Data[0])):
   #         print('{:<20s}'.format(Data[kl][jk]))
    return cap,gam,model,Data,hj,barrier