# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 19:34:26 2019

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
from platypus.problems import Problem
from platypus.algorithms import NSGAII
from platypus.problems import DTLZ2
import platypus.types as tt
def CounterExampleSLSQP(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga):   #finding counter examples through SLSQP adding them to data samples and again finding barrier polynomial
    Bo=B.copy()
    for ap in range(len(Bo)):
            Bo[ap]=Bo[ap].subs((yy[q],x[q]) for q in range(nvar))
    gh=[] #Converting the model values (coefficients p) to float
    for st in range(len(p)):
        gh.append(m.eval(p[st]))
    for tu in range(len(gh)):
        gh[tu]=gh[tu].as_decimal(20) #upto precision of 20
        gh[tu]=gh[tu].replace('?','')
        gh[tu]=float(gh[tu])
    for zx in range(len(gh)):
        if gh[zx]==0:
            gh[zx]=1e-19
    B1=sum(list(np.multiply(Bo,gh))) #Substitute the coeff in barrier polynomial 
    
    poNew=[] #Find expectation by substitution
    for sv in range(len(po[:,0])):
        po1=sum(list(np.multiply(list(po[sv,:]),gh)))
        poNew.append(po1)
    
    xi=[]
    for sd in range(nvar):
        xi.append(bounds[sd][0]+(bounds[sd][1]-bounds[sd][0])*np.random.rand())

#Use SciPy optimize to minimise B(x), B(x)-ga, 1-B(x) with the required bounds and constraints,
        #check the min value, if -ve, then it is a counter example
              
    
    def Objective1(ax):
        B2=B1.copy()
        B2=B2.subs((x[i],ax[i]) for i in range(len(ax)))
        return B2-0.001
    
    
    
    sol1=minimize(Objective1,xi,method="SLSQP",bounds=bounds)
    xb=list(sol1.x)
    func1min=sol1.fun
    
    xi=[]
    for sd in range(nvar):
        xi.append(bounds[sd][0]+(bounds[sd][1]-bounds[sd][0])*np.random.rand())   
    
    def Objective2(ax):
        B2=B1.copy()
        B2=B2.subs((x[i],ax[i]) for i in range(len(ax)))
        return ga-B2-0.001
    
    xc=[]
    func2min=[]
    for bn in range(len(EqX0)):
        def Constraint1(ax):
            EqX0c=EqX0.copy()
            cc=[]
            for pp in range(len(EqX0c[bn])):
                cc.append(-EqX0c[bn][pp].subs((x[i],ax[i]) for i in range(len(ax))))
            return cc
        con1= {'type': 'ineq','fun': Constraint1}
        sol2=minimize(Objective2,xi,method="SLSQP",constraints=con1)
        xc.append((list(sol2.x)))
        func2min.append(sol2.fun)      
        
    xi=[]
    for sd in range(nvar):
        xi.append(bounds[sd][0]+(bounds[sd][1]-bounds[sd][0])*np.random.rand())
    
    def Objective3(ax):
        B2=B1.copy()
        B2=B2.subs((x[i],ax[i]) for i in range(len(ax)))
        return B2-1-0.001 
    
    xd=[]
    func3min=[]
    for bm in range(len(EqX1)):
        def Constraint2(ax):
            EqX1c=EqX1.copy()
            cd=[]
            for pp in range(len(EqX1c[bm])):
                cd.append(-EqX1c[bm][pp].subs((x[i],ax[i]) for i in range(len(ax))))
            return cd
        con2= {'type': 'ineq','fun': Constraint2}
        sol3=minimize(Objective3,xi,method="SLSQP",constraints=con2)
        xd.append(list((sol3.x))) 
        func3min.append(sol3.fun)
    
    xe=[]
    func4min=[]   

    xi=[]
    for sd in range(nvar):
        xi.append(bounds[sd][0]+(bounds[sd][1]-bounds[sd][0])*np.random.rand())
    
    def Objective4(ax):
        B2=B1.copy()
        B2=B2.subs((x[i],ax[i]) for i in range(len(ax)))

        poNew1=poNew.copy()
        for bv in range(len(poNew1)):
            poNew1[bv]=poNew1[bv].subs((x[i],ax[i]) for i in range(len(ax)))
        poo=0
        for bx in range(len(poNew1)):
            poo=poo+(-(poNew1[bx]-B2-ca)) #This is a work around for condition
            #checking the sign of all the three, if all are negative, then its counter example
            #otherwise it is not (at least one OR condition is satisfying)
        return poo
    
    sol4=minimize(Objective4,xi,method="SLSQP",bounds=bounds)
    if sol4.success==True:
        xe=list(sol4.x)
        func4min=sol4.fun
   
    
    counterexamples=[]
    if func1min<=-0.01:
        counterexamples.append(xb)
    if all(func2min[bb]<=-0.01 for bb in range(len(func2min))):
        for bv in range(len(xc)):
            counterexamples.append(xc[bv])
    if all(func3min[bb]<=-0.01 for bb in range(len(func3min))):
        for bv in range(len(xd)):
            counterexamples.append(xd[bv])
    if len(xe)!=0:
        if all((poNew[bb]-B1-ca).subs((x[i],xe[i]) for i in range(len(xe)))>=0 for bb in range(len(poNew))):
            counterexamples.append(xe)
    
    cex=np.transpose(np.array(counterexamples))
    return gh, cex

def CounterExamplePlatypus(B,EqX0,EqX1,m,p,po,x,yy,u,w,bounds,nvar,ca,ga,nml):   #finding counter examples through Platypus adding them to data samples and again finding barrier polynomial
    counterexamples=[]
    boundss=[]
    for i in range(nvar):   
        boundss.append(tt.Real(bounds[i][0],bounds[i][1]))
   
    gh=[] #Converting the model values (coefficients p) to float
    for st in range(len(p)):
        gh.append(m.eval(p[st]))
    for tu in range(len(gh)):
        gh[tu]=gh[tu].as_decimal(20) #upto precision of 20
        gh[tu]=gh[tu].replace('?','')
        gh[tu]=float(gh[tu])
    for zx in range(len(gh)):
        if gh[zx]==0:
            gh[zx]=1e-19
    #Substitute the coeff in barrier polynomial 
 #   print("this is gh",gh)
    poNew=[] #Find expectation by substitution
    for sv in range(len(po[:,0])):
        po1=sum(list(np.multiply(sympify(gh),list(po[sv,:]))))
        poNew.append(po1)
        
    def min_bx(ax):
        Bo=B.copy()
        for i in range(len(B)):
            Bo[i]=Bo[i].subs((yy[q],ax[q]) for q in range(nvar))
        B1=sum(list(np.multiply(Bo,gh)))
   #     print("this is t",t)
        return [B1-0.001] #the two brackets are for multiple objective and constraints respectively.
  
    problem = Problem(nvar,1)
    problem.types[:] = boundss
    problem.function = min_bx
    aX=[]
    algorithm = NSGAII(problem)
    algorithm.run(1000)
    gX=[]
  #  print("this is algo",algorithm.result)
    for solution in algorithm.result:
        if solution.objectives[0]<-0.01 :
           aX.append(solution.objectives[0]) 
           gX.append(solution.variables)
    
    if aX:
     #  print("this is min(aX)",min(aX))
       u=aX.index(min(aX))
       counterexamples.append(gX[u])
    #   print("gX",gX[u])

    for vv in range(len(EqX0)):
        def min_bxga(ax):
            Bo=B.copy()
            for i in range(len(B)):
                Bo[i]=Bo[i].subs((yy[q],ax[q]) for q in range(nvar))
            B1=sum(list(np.multiply(Bo,gh)))
            EqX0c=EqX0.copy()
            cc=[]
            for i in range(len(EqX0c[vv]) ):
                cc.append(EqX0c[vv][i].subs((x[q],ax[q]) for q in range(nvar)))
            return [-B1+ga-0.001],cc #the two brackets are for multiple objective and constraints respectively.
       # print("this is EqqX0[vv]",len(EqX0[vv]))
        problem = Problem(nvar,1,len(EqX0[vv]))
        problem.types[:] = boundss
        problem.constraints[:] = "<=0"
        problem.function = min_bxga
        algorithm = NSGAII(problem)
        algorithm.run(1000)
        a0=[]
        g0=[]
        for solution in algorithm.result:
            if solution.objectives[0]<-0.01 :
               a0.append(solution.objectives[0]) 
               g0.append(solution.variables)
        if a0:
           print("this is g0",g0[a0.index(min(a0))])
           counterexamples.append(g0[a0.index(min(a0))])
    
        
    for vv in range(len(EqX1)):
        def min_bx1(ax):
            Bo=B.copy()
            for i in range(len(B)):
                Bo[i]=Bo[i].subs((yy[q],ax[q]) for q in range(nvar))
            B1=sum(list(np.multiply(Bo,gh)))
            EqX1c=EqX1.copy()
            cc=[]
            for i in range(len(EqX1c[vv])):
                cc.append(EqX1c[vv][i].subs((x[q],ax[q]) for q in range(nvar)))
            return [B1-1-0.001],cc #the two brackets are for multiple objective and constraints respectively.
        problem = Problem(nvar,1,len(EqX1[vv]))
        print("this is EqX1[vv]",len(EqX1[vv]))
        problem.types[:] = boundss
        problem.constraints[:] = "<=0"
        problem.function = min_bx1
        algorithm = NSGAII(problem)
        algorithm.run(1000)
        a1=[]
        g1=[]
        for solution in algorithm.result:
            if solution.objectives[0]<-0.01 :
               a1.append(solution.objectives[0]) 
               g1.append(solution.variables)
        if a1:
           print("this is g1",g1[a1.index(min(a1))])
           counterexamples.append(g1[a1.index(min(a1))])
    def exppo(ax):
        poNew1=poNew.copy()
        for bv in range(len(poNew1)):
            poNew1[bv]=poNew1[bv].subs((x[i],ax[i]) for i in range(len(ax)))
        Bo=B.copy()
        for i in range(len(B)):
            Bo[i]=Bo[i].subs((yy[q],ax[q]) for q in range(nvar))
        B1=sum(list(np.multiply(Bo,gh)))
        jk=[]
       
        for i in range(len(poNew1)):
            jk.append(poNew1[i]-B1-ca-0.001)
  #      print("this is jk",jk)
        return jk
    problem = Problem(nvar,nml)
    problem.types[:] = boundss
    problem.function = exppo
    algorithm = NSGAII(problem)
    problem.directions[:] = Problem.MAXIMIZE
    algorithm.run(1000)
    aP=[]
    gP=[]
    for solution in algorithm.result:
        if all(solution.objectives[i]>0.01 for i in range(nml)) :
           aP.append(sum(solution.objectives[i] for i in range(nml))) 
      #     print(a4)
           gP.append(solution.variables)
    #       print(g4)
    if aP:
        print("this is gP",gP[aP.index(max(aP))])
        counterexamples.append(gP[aP.index(max(aP))])

    cex=np.transpose(np.array((counterexamples)))
    print("this is gh",gh)
    return gh, cex    