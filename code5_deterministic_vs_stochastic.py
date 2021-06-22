import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

##SIR parameters
# Total population, N.
N = 10000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 100, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma
beta, gamma = 1, 0.6
r_0 = beta/gamma
print("r_0 = "+str(r_0))


##Deterministic SIR model
# A grid of time points (in days)
t = np.linspace(0, int(max_t), int(max_t)+1)

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



##Stochastic model on weighted network averaged 50 times

max_t=0
max_len=0
sol_number=50
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



#mean of the 10 simulations
t_mean,S_mean,I_mean,R_mean, counter=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
for i in range(sol_number): #we go through the 10 simulations
    sim=simulations[i]
    t,S,I,R=sim[0],sim[1],sim[2],sim[3]
    '''
    if len(t2)>len(t_mean): # we keep the largest time list
        t_mean=t2
    '''
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



##Plot the solutions

## Plot deterministic solution
fig = plt.figure(facecolor='w')
fig.suptitle('Analytic solution for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/N, 'b', alpha=0.5, lw=2, label='Omnivorous')
ax.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Vegetarian')
ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Ex-vegetarian')
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


##Plot stochastic solution
fig1 = plt.figure(facecolor='w')
fig1.suptitle('Stochastic solution for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma)+" averaged "+str(sol_number)+" times")
ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#thin lines
for j in range(len(simulations)-1):
    sim=simulations[j]
    t1,S1,I1,R1=sim[0],sim[1],sim[2],sim[3]
    ax1.plot(t1, S1/N, 'b', alpha=0.5, lw=0.5)
    ax1.plot(t1, I1/N, 'r', alpha=0.5, lw=0.5)
    ax1.plot(t1, R1/N, 'g', alpha=0.5, lw=0.5)

t2,S2,I2,R2=simulations[-1][0],simulations[-1][1],simulations[-1][2],simulations[-1][3]
ax1.plot(t1, S1/N, 'b', alpha=0.5, lw=0.5, label='Omnivorous')
ax1.plot(t1, I1/N, 'r', alpha=0.5, lw=0.5, label='Vegetarian')
ax1.plot(t1, R1/N, 'g', alpha=0.5, lw=0.5, label='Ex-vegetarian')

#mean of every simulation
ax1.plot(t_mean, S_mean/N, 'k', linestyle='--', alpha=1, lw=2)
ax1.plot(t_mean, I_mean/N, 'k', linestyle='--', alpha=1, lw=2)
ax1.plot(t_mean, R_mean/N, 'k', linestyle='--',alpha=1, lw=2)

ax1.set_xlabel('Time /days')
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
