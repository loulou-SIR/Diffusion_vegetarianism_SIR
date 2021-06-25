import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN
'''Objective : compare outbreak size depending on R0 for analytical and stochastic solution
WARNING: code is very long to run
'''

##Outbreak size for deterministic SIR

# Mean recovery rate, gamma, is fixed in advane
gamma=0.6
step=50
outbreak_list=[]
for n in [1000,10000,100000,1000000]:
    N = n
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 100, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    #outbreak are stored in a list
    outbreak_size=0
    r_0_list=np.linspace(0.2,2.7,step)
    outbreak=[]
    for r_0 in r_0_list:#we change r_0 parameter and calculate outbreak size
        beta=r_0*gamma
        t = np.linspace(0, 160, 160)#change later

        # The SIR model differential equations.
        def deriv(y, t, N, beta, gamma):
            S, I, R = y
            dSdt = -beta * S * I / N
            dIdt = beta * S * I / N - gamma * I
            dRdt = gamma * I
            return dSdt, dIdt, dRdt

        # Initial conditions vector
        y0 = S0, I0, R0
        # Integrate the SIR equations over the time grid, t.
        ret = odeint(deriv, y0, t, args=(N, beta, gamma))
        S, I, R = ret.T

        outbreak_size=(R[len(R)-1]-I0)/S0
        outbreak.append(outbreak_size)
    outbreak_list.append(outbreak)

##Outbreak size for stochastic SIR

# Mean recovery rate, gamma, is fixed in advane
gamma=0.6
step=50
outbreak_list1=[]
for n in [1000,10000,100000,1000000]:
    print('step n= '+str(n))
    N = n
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 100, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    #outbreak are stored in a list
    outbreak_size=0
    r_0_list=np.linspace(0.2,2.7,step)
    outbreak=[]
    for r_0 in r_0_list:#we change r_0 parameter and calculate outbreak size
        print('simulation n ='+str(n)+' R0 ='+str(r_0))
        beta=r_0*gamma

        max_t=0
        max_len=0
        sol_number=10
        simulations=[]
        for i in range(sol_number):
            #creating graph
            g=nx.watts_strogatz_graph(n=N, k=4, p=0.6)
            #initializing random weights
            E = g.number_of_edges()
            w = [random.random() for i in range(E)]
            s = max(w)
            w = [ i/s for i in w ] #normalizing
            len(w)
            k = 0
            for i, j in g.edges():
                g[i][j]['weight'] = w[k]
                k+=1
            #Use of the fast SIR function of EoN module
            t1,S1,I1,R1= EoN.fast_SIR(g, tau = beta, gamma=gamma, rho = I0/N, transmission_weight="weight", return_full_data=False)
            simulations.append([t1,S1,I1,R1])
            if max(t1)>max_t:
                max_t=max(t1)
            if len(t1)>max_len:
                max_len=len(t1)

        #mean of the 50 simulations
        t_mean,S_mean,I_mean,R_mean, counter=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
        for i in range(sol_number): #we go through the 50 simulations
            sim=simulations[i]
            t,S,I,R=sim[0],sim[1],sim[2],sim[3]

            new_t=np.linspace(min(t),max(t),5000)#creating a new time vector that will be fixed in number

            new_S=np.interp(new_t,t,S)#interpolating a new S vector from new_t, t, S
            new_I=np.interp(new_t,t,I)
            new_R=np.interp(new_t,t,R)

            for j in range(5000):#calculating the sum
                t_mean[j]+=new_t[j]
                S_mean[j]+=new_S[j]
                I_mean[j]+=new_I[j]
                R_mean[j]+=new_R[j]
                counter[j]+=1
        for k in range(5000): #calculating the mean value (replace by /5000)
            t_mean[k]/=counter[k]
            S_mean[k]/=counter[k]
            I_mean[k]/=counter[k]
            R_mean[k]/=counter[k]

        outbreak_size=(R_mean[len(R_mean)-1]-I0)/S0
        outbreak.append(outbreak_size)
    outbreak_list1.append(outbreak)


##Plot outbreak size
# Plot deterministic outbreak size
fig = plt.figure(facecolor='w')
fig.suptitle('Outbreak size depending on R0 for deterministic solution')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(r_0_list,outbreak_list[0],'b', alpha=0.5, lw=2, label='N=10'+'\u00B3')
ax.plot(r_0_list,outbreak_list[1],'r', alpha=0.5, lw=2, label='N=10'+'\u2074')
ax.plot(r_0_list,outbreak_list[2],'y', alpha=0.5, lw=2, label='N=10'+'\u2075')
ax.plot(r_0_list,outbreak_list[3],'k', alpha=0.5, lw=2, label='N=10'+'\u2076')

ax.set_xlabel('R0 = beta/gamma')
ax.set_ylabel('Outbreak size')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)

# Plot stochastic outbreak size
fig1 = plt.figure(facecolor='w')
fig1.suptitle('Outbreak size depending on R0 for stochastic solution')
ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax1.plot(r_0_list,outbreak_list1[0],'b', alpha=0.5, lw=2, label='N=10'+'\u00B3')
ax1.plot(r_0_list,outbreak_list1[1],'r', alpha=0.5, lw=2, label='N=10'+'\u2074')
ax1.plot(r_0_list,outbreak_list1[2],'y', alpha=0.5, lw=2, label='N=10'+'\u2075')
ax1.plot(r_0_list,outbreak_list1[3],'k', alpha=0.5, lw=2, label='N=10'+'\u2076')

ax1.set_xlabel('R0 = beta/gamma')
ax1.set_ylabel('Outbreak size')
ax1.set_ylim(0,1.2)
ax1.yaxis.set_tick_params(length=0)
ax1.xaxis.set_tick_params(length=0)
ax1.grid(b=True, which='major', c='w', lw=2, ls='-')
legend1 = ax1.legend()
legend1.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax1.spines[spine].set_visible(False)

#Show both solutions
plt.show()
