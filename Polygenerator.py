# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 17:33:29 2019

@author: yogit
"""
import numpy as np
from sympy import *
from z3 import *
import math
import time
def PolynomialGenDisc(B,function,analy,nm,nu,inp,x,yy,w,u,nvar,nw): #To generate the polynomial E(B(f(x,w,u)|(x,u))<=B(x)+c
    if analy==2:
        nml=1
        nul=0
        inpt=np.zeros([nml,nml])
    else:
        nml=nm
        nul=nu
        inpt=inp
    po=zeros(nml,len(B))
    for ab in range(nml): 
        BB=B.copy()     
        f=function.copy() 
        print("this f",f)
        print("this f[0]",f[0])
        print(f)
        for y in range(len(f)):                                                                                                                                                                                     
            for z in range(len(inpt[:,ab])):
                print(f)
                f[y]=f[y].subs(u[z],inpt[z][ab]) #for each input, substitute u in f(x,u)
                
        
        for ap in range(len(B)):
            BB[ap]=BB[ap].subs((yy[q],f[q]) for q in range(nvar)) #Now substitute x=f(x,w) in B(x)
        
        for aq in range(len(BB)): #To calculate the expectation of B(f(x,w)) E=c*m
            if BB[aq]==1: #if the term is 1, both c and m value is 1
                c=[1]
                m=[1]
            else:
                temp=Poly(expand(BB[aq])) #For each term in B(f(x,w)), seperate out the constants and the variables
                
                c=temp.coeffs() #c stores the constant
                m=[prod(x1**k for x1, k in zip(temp.gens, mon)) for mon in temp.monoms()]
            print(m) #m is a list that stores variables
            
            for k in range(len(m)): 
                s1=(str(m[k]).replace("**","^")).split("*")
                chk=0
                for k2 in range(len(s1)):
                    s2=s1[k2].split("^")
                    for k3 in range(1,nw+1):
                        if s2[0]=="w"+str(k3) and len(s2)==1: #If there is an odd w term, then the expec of that term is 0 
                            m[k]=0
                            chk=1
                            break
                        if s2[0]=="w"+str(k3) and int(s2[1])%2==1:                                                                                                                      
                            m[k]=0
                            chk=1
                            break
                        if s2[0]=="w"+str(k3) and int(s2[1])%2==0: #If there is an even w term, use this formula for expectation
                            gd=int(s2[1])
                            dfk=gd/2
                            df=math.factorial(gd)/((2**dfk)*math.factorial(dfk))
                            s1[k2]=str(df)
                if chk==0:
                  m[k]=prod(sympify(s1)) #now your expectation is independent of w (contains x)
            poX=poX=sum(np.multiply(c,m)) #multiply
            po[ab,aq]=poX
    return po,nml,nul,inpt #this polynomial is returned for every input. So 3 inputs -> 3 polynomials (this po is E(B(f(x,u))))

def PolynomialGenCont(B,function,analy,nm,nu,inp,x,yy,w,u,nvar,nw): #To generate the polynomial E(B(f(x,w,u)|(x,u))<=B(x)+c
    if analy==2:
        nml=1
        nul=0
        inpt=np.zeros([nml,nml])
    else:
        nml=nm
        nul=nu
        inpt=inp
    po=zeros(nml,len(B))
    for ab in range(nml): 
        BB=B.copy()     
        f=function.copy()  
        for y in range(len(f)):                                                                                                                                                                                     
            for z in range(len(inpt[:,ab])):
                print(f)
                f[y]=f[y].subs(u[z],inpt[z][ab]) #for each input, substitute u in f(x,u)
                
        
        for ap in range(len(B)):
            BB[ap]=BB[ap].subs((yy[q],x[q]) for q in range(nvar)) #Now substitute x=f(x,w) in B(x)
        
        deriv_1=[[0 for i in range(len(BB))] for j in range(nvar)] 
        for q in range(nvar):                         # derivative of Barrier polynomial
            for i in range(len(BB)):
                deriv_1[q][i]=diff(BB[i],x[q])
        print(deriv_1)
        deriv_1=Matrix(deriv_1)
        
        a=[0 for j in range(nvar)] 
        for j in range(nvar):
            a[j]=list(deriv_1[j,:]*f[j])
        print("this is a",a)
        
        sum_deriv_1=[sum(x) for x in zip(*a)]            #derivatives summed up
        print("this is sum_deriv_1",sum_deriv_1)
        
        deriv_2=[[0 for i in range(nvar)] for j in range(nvar)] 
        for k in range(nvar):                           #finding double derivative of barrier polynomial
            for s in range(nvar):
                d=[]
                for i in range(len(BB)):
                    d.append(diff(BB[i],x[k],x[s]))
                deriv_2[k][s]=d    
        print("THIS IS deriv_2",deriv_2)
        
        G=np.asarray(g)
        G_trans=G.T
        gtmuld2=np.dot(G_trans,deriv_2)
        print("this is gtmuld2",gtmuld2)
        G=G_trans.T
        print(G)
        gtmul2g=[[0 for i in range(G.shape[1])] for j in range(G.shape[1])]  #finding G_trans*deriv_2*G
        print(gtmul2g)
        for t in range(G.shape[1]):
            for j in range(G.shape[1]):
                for k in range(nvar):
                    gtmul2g[t][j]=gtmul2g[t][j]+np.multiply(gtmuld2[t][k],G[k][j])
        hh=Matrix(gtmul2g)
        oo=[]
        for i in range(hh.shape[0]):
            oo.append(hh[i*(hh.shape[0]+1)])
        z=[sum(x) for x in zip(*oo)]
        half_trace=list(z)
        print("this is trace after",ab,half_trace)
        sum_final=sum_deriv_1+half_trace     #finding final sum
        print("this is sum before",sum_final)
        for aq in range(len(BB)):
            po[ab,aq]=sum_final[aq]
    return po,nml,nul,inpt