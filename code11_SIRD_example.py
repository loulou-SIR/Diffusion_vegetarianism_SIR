import networkx as nx
import matplotlib.pyplot as plt
import EoN
import xlrd
import random
from matplotlib import pyplot as plt
import numpy as np
import time
import csv

'''This program gives a concrete example of how to implement the custom SIRD model for a long-term scenario
Because we forecast the growth of vegetarianism on very long times, we changed the unit of the timescale to year
WARNING: Before running this code, code 5 and 8 must be run'''

##Execute custom stochastic SIRD algorithm

#Network parameters
N1=10000
K1=4
P1=0.3
#R0 value
r0=1.4
#Rates in /day
beta1_day=0.009
#Total rate gamma=gamma1+gamma1
gamma_day=beta1_day/r0
#15% of transitioning quit vegetarianism whereas 85% remain life-long vegetarians
gamma1_day=0.15*gamma_day
gamma2_day=0.85*gamma_day

#Rates in /year
beta1=beta1_day*365
gamma1=gamma1_day*365
gamma2=gamma2_day*365

#Simulation number
sol_number=50

#Execute custom stochastic SIRD
t1,S1,I1,R1,D1= average_stochastic_custom(sol_number, beta1,gamma1,gamma2)



##Plot custom stochastic SIRD and deterministic SIRD on same figure

fig1 = plt.figure(facecolor='w')
fig1.suptitle('Custom SIRD algorithm N=10000 R0='+str(r0))
ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#plot custom stochastic SIRD
ax1.plot(t1, S1/N1, 'b', linestyle='-', alpha=0.8, lw=2, label='Omnivorous')
ax1.plot(t1, I1/N1, 'r', linestyle='-', alpha=0.8, lw=2, label='Transitioning')
ax1.plot(t1, R1/N1, 'y', linestyle='-',alpha=0.8, lw=2, label='Ex-vegetarian')
ax1.plot(t1, D1/N1, 'g', linestyle='-',alpha=0.8, lw=2, label='Vegetarian')


ax1.set_xlabel('Time /years') #changed the timescale to years
ax1.set_ylabel('Number (10000s)')
ax1.set_ylim(0,1.2)
ax1.yaxis.set_tick_params(length=0)
ax1.xaxis.set_tick_params(length=0)
ax1.grid(b=True, which='major', c='w', lw=2, ls='-')
legend1 = ax1.legend()
legend1.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax1.spines[spine].set_visible(False)

#show both figures
plt.show()