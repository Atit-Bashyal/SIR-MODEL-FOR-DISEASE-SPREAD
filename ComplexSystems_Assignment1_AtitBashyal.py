#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:49:18 2019

@author: atit
"""
'''Import Needed Libraries'''

from scipy . integrate import odeint
import numpy as np
import matplotlib . pyplot as plt

'''Initialize Parameters'''

c = 1.0 # contact rate per susceptible individual
a = 0.025 # Probability of disease transmition
rho = 0.03 # recovery rate per infected individual
d= 0.1 # Death rate of healthy(both Susceptibe and recovered)
sigma = 0.3 # rate of recovered individuals that become susceptible
theta = 0.15 #Immigration rate of people adding to susceptible population
delta = 0.25 # death rate of infected individuals

'''Build  function for integration'''

def deriv_sir(var,t):
    # index variables from the variable list (var)
    # variables of the model are the three types of popultion represented in the model
    s = var[0]
    i = var[1]
    r= var[2]
    # write the differential equations for each variable
    dsdt = (theta) - (d*s) - (a*c*s*i) + (sigma*r)
    didt = (a*c*s*i) - (delta*i) - (rho*i)
    drdt = (rho*i) - (sigma*r) - (d*r)
   
    return [dsdt,didt,drdt]

'''Declare time points to calcualte numerical solutions'''

t = np.linspace(0,30.0) # 30 days

'''Declare initial condition list''' 

var0 = [85,10,5] # index same as variable list

'''compute the numerical solutions'''
  
soln = odeint(deriv_sir,var0,t) 
# soln will be a matrix with each column containing solution for each variable
#index of column will be same as variable list used in integration function

'''extract solution list for each variable''' 

# list of numerical solution of suseptiable population at each time point
sus = soln[:,0]

# list of numerical solution of infected population at each time point
inf = soln[:,1]

# list of numbrical solution of recovering population at each time point
rec = soln[:,2]

# include a total population list at each computational time point to look into steady states
tot = sus+inf+rec

'''plotting the solutions against time'''

plt.plot(t,sus,'b',label= 'suseptiable population') # blue for Suseptiable population
plt.plot(t,inf,'r',label='infected population') # red for Infected population
plt.plot(t,rec,'g',label='revocering population') # green for recovering
plt.plot(t,tot,'-.k',label='total population') # black for total
plt.legend(loc="best") # set location of legend to best possible visibility
plt.title('SIR Model') # title the plot
#label the axis
plt.ylabel('Number of people') 
plt.xlabel('time(days)')
# save the plot as pdf file
plt.savefig('SIR model_Assignment1_AtitBashyal')
#display the plot if wanted!!
plt.show()


## with given rates totalpopulation die within 30 days
## changing recovery rate to 0.3 does not let the total population die within 30 days
## imigration rate increasing by a factor of 10 also leads to the all population not dying
## increasing death rate of infected population also prevents collaspe of the population size


        







   







