# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 18:56:45 2019

@author: yogit
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:42:46 2019

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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib
import Discret         # module storing the function for generating barrier certificate for system having discrete dynamics
from random import random
def disab():   # This function enables and disables the text box depending on the type of analysis
    if ana.get()==1:  # for synthesis input text box is enabled
       inpu.config(state=NORMAL)
    if ana.get()==2:  #for verification input text box is disabled
       inpu.config(state=DISABLED)    
    if p.get()==1:    # for user defined function, the text box for polynomial form is disabled
       polyno.config(state=DISABLED) 
       defin.config(state=NORMAL)
    if p.get()==2:    # for polynomial form, the text box for user defined function is disabled
       defin.config(state=DISABLED)
       polyno.config(state=NORMAL)
def button_pressed():
    global nm 
    global nu 
    global nvar 
    global nw 
    global n 
    global inp 
    global funcn 
    global bounds 
    global EqX0
    global EqX1 
    global B 
    global x 
    global yy 
    global w 
    global u 
    global constant 
    global gamma 
    global Td 
    global nsamples
    global inp
    global jjj
    global kkk
    global lll
    global mmm
    global pol
    global bar
    global inicondn
    x=[]
    yy=[]
    w=[]
    u=[]
    funcn=[]
    inicondn=[]
    bounds=eval(statesp.get())  #the bounds of the state variables taken from the user are assigned to bounds
    if ana.get()==1:#if the user selects synthesis then take input set from user and find no of input variables and input mode from the data entered
       inp=np.array(eval(inpu.get()))
       nu=inp.shape[0]
       nm=inp.shape[1]
    else:   #if user selects verification take u=[1]
       nm=1
       nu=0
       u=[1]
    nvar=len(bounds) # the no of state variables can be found from the length of bounds entered by the user
    nw=int(dist.get()) # take number of disturbance variables from user as input
    for i in range(1,nvar+1):
        ex='x'+str(i)
        x.append(symbols(ex))         #Introducing symbols for state variables
    
    for i in range(1,nvar+1):
        ey='y'+str(i)
        yy.append(symbols(ey))        #Temporary state variables for operation within functions
            
    for i in range(1,nw+1):
        ww='w'+str(i)
        w.append(symbols(ww))        #Symbol for disturbance variables
    
    for i in range(1,nu+1):
        uu='u'+str(i)
        u.append(symbols(uu))
    funcn=list(eval(dyn.get()))      #Symbol for input variables
    EqX0=eval(initi.get())           #Initial region
    EqX1=eval(reach.get())           #Reachability (unsafe) region
    nsamples=int(sam.get())          #No of data samples to be generated randomly in the bounds
    constant=eval(consta.get())      #Initial value of C
    gamma=eval(gha.get())            #Initial value of gamma
    Td=eval(hori.get())              #Temporal horizon Td
    N=500                                    
    inicondn=eval(initial.get())     #initial state for plotting
    
    if dy.get()==1:   # the user selects the discrete dynamics
       if p.get()==1:      # user selects user defined function
          B=eval(defin.get())
       elif p.get()==2:    # user selects polynomial form
            n=int(polyno.get())
            B=list(itermonomials(yy,n))
       
       (jjj,kkk,lll,mmm,pol,bar)=Discret.discrete(nsamples,constant,gamma,EqX0,EqX1,funcn,B,nm,nu,nvar,nw,inp,bounds,opt.get(),ana.get(),x,yy,w,u,dy.get())
       print("this is model",lll)
       print("this is data",mmm)
       print("this is poly",pol)
       print("this is barrier",bar)
       pr=kkk+jjj*Td  # the upper bound of the probability of reaching the unsafe region
       if pr>1:
          pr=1
       if ana.get()==2:  # if the user selects verification then no controller is required hence print not available
           mmm=' Not available'
       s='Value of c is {}\nValue of gamma is {}\nThe barrier function is {}\nThe controller designed is \n{}\n The upper bound of the probability of reaching the unsafe region is {}  '.format(jjj,kkk,bar,mmm,pr)
       tex.insert(END,s)
       tex.see(END)
       f = open("Results.txt", "w")      # create file Results
       f.close 
       f= open("Results.txt","w")        # open the created file and write the results in it          
       f.write("This is final c %s\n" %jjj)
       f.write("This is final gamma %s\n" %kkk)
       f.write("This is the final model %s\n" %lll)
       if ana.get()==1:   # if the user selects synthesis then only write the controller designed into the file
          f.write("The controller design %s\n" %mmm)
       f.write("The upper bound of probability of going in region X1 is %s" %pr)
       f.close
       f=open("Results.txt","r")
       print(f.read())
master =Tk()   # creates the master window for the GUI
dy = IntVar()  # Bullet option for dynamics 
ana = IntVar() # Bullet option for Analysis 
opt=IntVar()   # Bullet option for optimization
smt= IntVar()  # Bullet option for SMT solver
p = IntVar()   # Bullet option for type of barrier certificate
class CustomDialog(tkinter.simpledialog.Dialog):   # Message Box Example in front of each text box is created using this function

    def __init__(self, parent, title=None, text=None):
        self.data = text
        tkinter.simpledialog.Dialog.__init__(self, parent, title=title)

    def body(self, parent):
        self.geometry("810x200")
        self.text = Text(self, width=40, height=4)
        self.text.pack(fill="both", expand=True)

        self.text.insert("1.0", self.data)

        return self.text

def plot():    # function for plotting the results
    matplotlib.use('TkAgg')
    # Initialize an instance of Tk
    # Initialize matplotlib figure for graphing purposes
    fig = plt.figure(1)
    root = Tk()
    # Special type of "canvas" to allow for matplotlib graphing
    canvas = FigureCanvasTkAgg(fig, master=root)
    plot_widget = canvas.get_tk_widget()
    cc=jjj
    f=funcn.copy()
    N=5
    root.title("Plot of state v/s time")
    N_lim=10
    we=inicondn
    h=[]
    h.append(inicondn)
    for i in range(N_lim):
       poll=pol.copy()
       for kl in range(len(poll)):
           poll[kl]=poll[kl].subs((x[ad],we[ad]) for ad in range(nvar))
       f=funcn.copy() 
       for i in range(len(funcn)):
           f[i]=f[i].subs((w[b],random()) for b in range(nw))
       for kl in range(len(f)):
           f[kl]=f[kl].subs((x[ad],we[ad]) for ad in range(nvar))
       for r in range(nm):
          if eval(str(poll[r]))<cc:
              for rt in range(len(f)):
                  f[rt]=f[rt].subs((u[z],inp[z][r]) for z in range(len(inp[:,r])))
              we=eval(str(f))
              print(we)
              h.append(we)
              break
    t=[i for i in range(len(h))]       
    plt.plot(t,h)
    plot_widget.grid(row=14, column=3)
    plt.xlabel('Time (sec)')
    plt.ylabel('Temperature (in Celsius)')
    plt.show()
def help_dyn():   # these functions are called by custom dialog
    fromonk_text ="Example : The dynamics of the system can be entered as follows: If the system has dynamics of form :\n x(k + 1) = x(k) + 1*(8e-3*(15 − x(k)) + 3.6e-3*(55 − x(k))u(k)) + 0.1w(k)\n Enter Dynamics as :\n[x[0]+1*(8e-3*(15-x[0])+3.6e-3*(55-x[0])*u[0])+0.1*w[0]]"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_statesp():
    fromonk_text ="Example : The states space/bound of each state of the system are to be represented as tuples. If the state space bounds for x0 =[0,45] and x1=[23,34]\n Enter [(0,45),(23,34)]"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_inpu():
    fromonk_text ="Example : The input variables and input to the system can be entered as follows : If the system has 3 input modes and 2 input variables u0=[1,2,3] u1=[3,4,5]\n Enter [[1,2,3],[3,4,5]]"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_initi():
    fromonk_text ="Example : The initial set of the specification can be given as follows :If the initial state of x0 lies in [21,22]\n Enter text as follows :[[x[0]-22,21-x[0]]] Here the rows represent intersections and columns represents unions. All inequalities are to be entered in form of<=0"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_reach():
    fromonk_text ="Example : The reachable set of the specification can be given as follows :If the reachable state of x0 lies in [21,22]\n Enter text as follows :[[x[0]-22,21-x[0]]] Here the rows represent intersections and columns represents unions. All inequalities are to be entered in form of<=0"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_defin():
    fromonk_text ="Example : The user defined polynomial can be entered in the following form : [cosx(0),sinx(1)**2 + cosx(0),x(2)**2]\n If entering a constant number enter it as sympify(constant)\n Example : enter 2 as sympify(2)"
    CustomDialog(master, title="Example", text=fromonk_text)
def help_inistate():
    fromonk_text ="Example : The initial state can be defined here : \n If the system has 2 variables the enter initial state as [21,33]"
    CustomDialog(master, title="Example", text=fromonk_text)
    
    
Label(master, text="Select type of system :").grid(row=0,column=0,sticky=W,padx=20,pady=5)
Radiobutton(master, text="Discrete", variable=dy, value=1).grid(row=0,column=1,sticky=W,padx=10,pady=5)
Radiobutton(master, text="Continuous", variable=dy, value=2).grid(row=0,column=2,sticky=W,padx=10,pady=5)
Label(master, text="Select type of analysis :").grid(row=1,column=0,sticky=W,padx=20,pady=5)
Radiobutton(master, text="Synthesis", variable=ana,command=disab, value=1).grid(row=1,column=1,sticky=W,padx=10,pady=5)
Radiobutton(master, text="Verification", variable=ana, command=disab, value=2).grid(row=1,column=2,sticky=W,padx=10,pady=5)
Label(master, text="Select type of SMT Solver :").grid(row=2,column=0,sticky=W,padx=20,pady=5)
Radiobutton(master, text="Z3", variable=smt,command=disab, value=1).grid(row=2,column=1,sticky=W,padx=10,pady=5)
Radiobutton(master, text="MathSAT", variable=smt, command=disab, value=2).grid(row=2,column=2,sticky=W,padx=10,pady=5)
Label(master, text="Select tool to compute counter example :").grid(row=3,column=0,sticky=W,padx=20,pady=5)
Radiobutton(master, text="SLSQP", variable=opt, value=1).grid(row=3,column=1,sticky=W,padx=10,pady=5)
Radiobutton(master, text="Platypus", variable=opt, value=2).grid(row=3,column=2,sticky=W,padx=10,pady=5)
Radiobutton(master, text="dReal", variable=opt, value=3).grid(row=4,column=1,sticky=W,padx=10,pady=5)
master.title("Synthesis and Verification via control barrier function")
Label(master, text="Dynamics of the system").grid(row=10,sticky=W, padx=20,pady =5)
Button(master, text ="Example", command = help_dyn).grid(row=10,column=2,sticky=W, padx=20,pady =5)
Label(master, text="State space/bounds of the system").grid(row=11,sticky=W,padx=20, pady = 5)
Button(master, text ="Example", command = help_statesp).grid(row=11,column=2,sticky=W,padx=20, pady = 5)
Label(master, text="Input set").grid(row=12,sticky=W,padx=20, pady = 5)
Button(master, text ="Example", command = help_inpu).grid(row=12,column=2,sticky=W,padx=20, pady = 5)
Label(master, text="Number of disturbance variables").grid(row=13,sticky=W,padx=20, pady = 5)
Label(master, text="Initial set").grid(row=14,sticky=W,padx=20, pady = 5)
Button(master, text ="Example", command = help_initi).grid(row=14,column=2,sticky=W,padx=20, pady = 5)
Label(master, text="Reachable set").grid(row=15,sticky=W,padx=20, pady = 5)
Button(master, text ="Example", command = help_reach).grid(row=15,column=2,sticky=W,padx=20, pady = 5)
Label(master, text="Number of data samples").grid(row=16,sticky=W,padx=20, pady = 5)
Label(master, text="Select the type of barrier certificate :").grid(row=17,column=0,sticky=W,padx=20,pady=5)
Radiobutton(master, text="User-defined function", variable=p,command=disab, value=1).grid(row=18,column=0,sticky=W,padx=20,pady=5)
Button(master, text ="Example", command = help_defin).grid(row=18,column=2,sticky=W,padx=20,pady=5)
Radiobutton(master, text="Polynomial form (Enter degree of polynomial)", variable=p,command=disab, value=2).grid(row=19,column=0,sticky=W,padx=20,pady=5)
Label(master, text="Initial value of constant c :").grid(row=20,column=0,sticky=W,padx=20,pady=5)
Label(master, text="Initial value of gamma :").grid(row=21,column=0,sticky=W,padx=20,pady=5)
Label(master, text="Temporal horizon T:").grid(row=22,column=0,sticky=W,padx=20,pady=5)
Label(master, text="Initial state for plotting (Optional):").grid(row=23,column=0,sticky=W,padx=20,pady=5)
Button(master, text ="Example", command = help_inistate).grid(row=23,column=2,sticky=W,padx=20, pady = 5)

initial=(Entry(master))
statesp = (Entry(master))
inpu = (Entry(master))
dist = (Entry(master))
dyn = (Entry(master))
initi = (Entry(master))
reach = (Entry(master))
sam = (Entry(master))
defin = (Entry(master))
polyno = (Entry(master))
consta = (Entry(master))
gha = (Entry(master))
hori = (Entry(master))
tex = Text(master)
tex.grid(row=0,column=3,rowspan=28,sticky=N+W+E+S,padx=20,pady=20)  # this text box displays the final result with the controller and the barrier certificate
dyn.grid(row=10, column=1)
statesp.grid(row=11, column=1)
inpu.grid(row=12, column=1)
dist.grid(row=13, column=1)
initi.grid(row=14, column=1)
reach.grid(row=15, column=1)
sam.grid(row=16, column=1)
defin.grid(row=18, column=1)
polyno.grid(row=19, column=1)
consta.grid(row=20, column=1)
gha.grid(row=21, column=1)
hori.grid(row=22, column=1)
initial.grid(row=23, column=1)

Button(master, text='                         Compute result                                  ',command=button_pressed).grid(row=24, column=0, columnspan=1, pady=4)
Button(master, text='       Cancel        ',command=master.quit).grid(row=24, column=2,columnspan=1,sticky=E, pady=4)
Button(master, text='             Plot Result                      ',command=plot).grid(row=24, column=1,columnspan=1,sticky=E, pady=4)
mainloop() # this function displays the GUI window