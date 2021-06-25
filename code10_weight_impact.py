import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

'''In this program we plot the impact of the weight on the number of transitioning individuals
'''
##SIR parameters
# Total population, N.
N = 10000
# Initial number of Vegetarian and Ex vegetarian individuals, I0 and R0.
I0, R0 = 100, 0
# Everyone else, S0, is Omnivorous to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma
beta, gamma = 1, 0.6


##Weight impact on small world network

#List of average degree values
weigth_list=[0.1,0.25, 0.5,0.75,1]
#list of simulations with different weights
w_simulations=[]
#average degree of graph
k_value=4
#probability of reconnecting edges (small world)
p=0.3
for w_value in weigth_list:
    max_t=0
    max_len=0
    sol_number=25
    simulations=[]
    for i in range(sol_number):
        #creating graph
        g=nx.watts_strogatz_graph(n=N, k=k_value, p=p)
        #initializing fixed weights
        E = g.number_of_edges()
        w = [w_value for i in range(E)]
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



    #mean of the 10 simulations
    t_mean,S_mean,I_mean,R_mean, counter=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
    for i in range(sol_number): #we go through the 10 simulations
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
    w_simulations.append([t_mean,S_mean,I_mean,R_mean])


## Plot weigth impact for I list
fig = plt.figure(facecolor='w')
fig.suptitle('Influence of weights on small world network for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma)+"k = "+str(k_value)+", p = "+str(p))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#Plot the mean of every stochastic simulation
#linestyle_list=['--',':','-','-.']
color_list=['b','g','m','y','c','r','k','p','d']
for i in range(len(weigth_list)):
    w=weigth_list[i]
    sim=w_simulations[i]
    color=color_list[i]
    #style=linestyle_list[i]
    t_mean,S_mean,I_mean,R_mean=sim[0],sim[1],sim[2],sim[3]
    ax.plot(t_mean, I_mean/N, color, linestyle='-', alpha=1, lw=2, label='Transitioning for w = '+str(w))

#Parameter
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (10000s)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)


##Show
plt.show()