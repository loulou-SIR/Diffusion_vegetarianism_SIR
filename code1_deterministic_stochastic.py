import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

'''In this code, we implement the deterministic SIR, deterministic SIRD and stochastic SIR models
'''


##Parameters for SIR
# Total population, N.
N = 10000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0= 0.02*N,0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, recovery rate gamma 1 and death rate gamma 2
beta, gamma = 1,0.6


##Parameters for SIRD
# Total population, N1.
N_ = 10000
# Initial number of infected and recovered and dead individuals, I1 and R1 and D1
I0_, R0_, D0_= 0.02*N_,0,0
# Everyone else, S1, is susceptible to infection initially.
S0_ = N_ - I0_ - R0_ - D0_
# Contact rate, beta, recovery rate gamma 1 and death rate gamma 2
beta1, gamma1,gamma2 = 1,0.1,0.5

##Time vector
# A grid of time points (in days)
days=30
t = np.linspace(0, days, days)

##SIR deterministic model

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


##SIRD deterministic model

# The SIRD model differential equations.
def deriv(y, t, N, beta, gamma1,gamma2):
    S, I, R, D= y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - (gamma1 + gamma2) * I
    dRdt = gamma1 * I
    dDt = gamma2 * I
    return dSdt, dIdt, dRdt, dDt

# Initial conditions vector
y0_ = S0_, I0_, R0_, D0_
# Integrate the SIRD equations over the time grid, t.
ret1 = odeint(deriv, y0_, t, args=(N_, beta1, gamma1, gamma2))
S1, I1, R1, D1 = ret1.T

##SIR stochastic model

max_t=0
max_len=0
sol_number=50 #number of simulations
simulations=[]
for i in range(sol_number):
    #creating graph
    g=nx.watts_strogatz_graph(n=N, k=4, p=0.05)

    #Use of the fast SIR function of EoN module
    t2,S2,I2,R2= EoN.fast_SIR(g, tau = beta, gamma=gamma,  rho=I0/N, return_full_data=False)
    simulations.append([t2,S2,I2,R2])
    if max(t2)>max_t:
        max_t=max(t2)
    if len(t2)>max_len:
        max_len=len(t2)


#mean of the 50 simulations
t_mean,S_mean,I_mean,R_mean, counter=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
for i in range(sol_number): #we go through the 10 simulations
    sim=simulations[i]
    t2,S2,I2,R2=sim[0],sim[1],sim[2],sim[3]

    new_t=np.linspace(min(t2),max(t2),5000)#creating a new time vector that will be fixed in number

    new_S=np.interp(new_t,t2,S2)#interpolating a new S vector from new_t, t, S
    new_I=np.interp(new_t,t2,I2)
    new_R=np.interp(new_t,t2,R2)

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


##Plot the 3 solutions obtained

# Plot deterministic SIR
fig1 = plt.figure(facecolor='w')
fig1.suptitle('Deterministic SIR for N=10000 I0= '+str(I0)+' beta= '+str(beta)+' gamma = '+str(gamma))

ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax1.plot(t, S/N, 'b', alpha=0.5, lw=2, label='Omnivorous')
ax1.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Vegetarians')
ax1.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Ex-vegetarians')

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


# Plot deterministic SIRD
fig = plt.figure(facecolor='w')
fig.suptitle('Deterministic SIRD for N=10000 I0= '+str(I0_)+' beta= '+str(beta1)+' gamma1 = '+str(gamma1)+' gamma2 = '+str(gamma2))

ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S1/N_, 'b', alpha=0.5, lw=2, label='Omnivorous')
ax.plot(t, I1/N_, 'r', alpha=0.5, lw=2, label='Transitioning vegetarians')
ax.plot(t, R1/N_, 'y', alpha=0.5, lw=2, label='Ex-vegetarians')
ax.plot(t, D1/N_, 'g', alpha=0.5, lw=2, label='Vegetarians')

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
fig2 = plt.figure(facecolor='w')
fig2.suptitle('Stochastic SIR for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma)+" averaged "+str(sol_number)+" times")
ax2 = fig2.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#thin lines
for j in range(len(simulations)-1):
    sim=simulations[j]
    t2,S2,I2,R2=sim[0],sim[1],sim[2],sim[3]
    ax2.plot(t2, S2/N, 'b', alpha=0.5, lw=0.5)
    ax2.plot(t2, I2/N, 'r', alpha=0.5, lw=0.5)
    ax2.plot(t2, R2/N, 'g', alpha=0.5, lw=0.5)

t2,S2,I2,R2=simulations[-1][0],simulations[-1][1],simulations[-1][2],simulations[-1][3]
ax2.plot(t2, S2/N, 'b', alpha=0.5, lw=0.5, label='Omnivorous')
ax2.plot(t2, I2/N, 'r', alpha=0.5, lw=0.5, label='Vegetarian')
ax2.plot(t2, R2/N, 'g', alpha=0.5, lw=0.5, label='Ex-vegetarian')

#mean of every simulation
ax2.plot(t_mean, S_mean/N, 'k', linestyle='--', alpha=1, lw=2)
ax2.plot(t_mean, I_mean/N, 'k', linestyle='--', alpha=1, lw=2)
ax2.plot(t_mean, R_mean/N, 'k', linestyle='--',alpha=1, lw=2)

ax2.set_xlabel('Time /days')
ax2.set_ylabel('Number (10000s)')
ax2.set_ylim(0,1.2)
ax2.yaxis.set_tick_params(length=0)
ax2.xaxis.set_tick_params(length=0)
ax2.grid(b=True, which='major', c='w', lw=2, ls='-')
legend2 = ax2.legend()
legend2.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax2.spines[spine].set_visible(False)


##Show figures
plt.show()

