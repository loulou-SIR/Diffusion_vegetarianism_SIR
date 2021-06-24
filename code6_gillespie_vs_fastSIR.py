import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

'''in this file we want to compare the results given by Gillespie and Fast SIR algorithms for the same networks in entry
Results are plotted and the relative difference in % between the final states of S,I and R is printed in console
'''


##SIR parameters
# Total population, N.
N = 10000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 100, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma
beta, gamma = 1, 0.6


##Running algorithms

max_t=0
max_len=0
sol_number=50
simulations,simulations2=[],[]
for i in range(sol_number):
    #creating graph
    g=nx.watts_strogatz_graph(n=N, k=4, p=0.6)

    #Use of the fast SIR function of EoN module
    t1,S1,I1,R1= EoN.fast_SIR(g, tau = beta, gamma=gamma, rho = I0/N, return_full_data=False)
    simulations.append([t1,S1,I1,R1])
    if max(t1)>max_t:
        max_t=max(t1)
    if len(t1)>max_len:
        max_len=len(t1)

    #Use of the Gillespie function of EoN module
    t2,S2,I2,R2= EoN.Gillespie_SIR(g, tau = beta, gamma=gamma, rho = I0/N, return_full_data=False)
    simulations2.append([t2,S2,I2,R2])
    if max(t2)>max_t:
        max_t=max(t2)
    if len(t2)>max_len:
        max_len=len(t2)


#mean of the 50 simulations for fast SIR
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

#mean of the 50 simulations for Gillespie
t_mean2,S_mean2,I_mean2,R_mean2, counter2=np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
for i in range(sol_number): #we go through the 50 simulations
    sim=simulations2[i]
    t,S,I,R=sim[0],sim[1],sim[2],sim[3]

    new_t=np.linspace(min(t),max(t),5000)#creating a new time vector that will be fixed in number

    new_S=np.interp(new_t,t,S)#interpolating a new S vector from new_t, t, S
    new_I=np.interp(new_t,t,I)
    new_R=np.interp(new_t,t,R)

    for j in range(5000):#calculating the sum
        t_mean2[j]+=new_t[j]
        S_mean2[j]+=new_S[j]
        I_mean2[j]+=new_I[j]
        R_mean2[j]+=new_R[j]
        counter2[j]+=1
for k in range(5000): #calculating the mean value (replace by /5000)
    t_mean2[k]/=counter2[k]
    S_mean2[k]/=counter2[k]
    I_mean2[k]/=counter2[k]
    R_mean2[k]/=counter2[k]

#difference between the final states of Gillespie and fast SIR
diff_S=(S_mean[len(S_mean)-1]-S_mean2[len(S_mean2)-1])*100/N
diff_I=(I_mean[len(I_mean)-1]-I_mean2[len(I_mean2)-1])*100/N
diff_R=(R_mean[len(R_mean)-1]-R_mean2[len(R_mean2)-1])*100/N
print("diff S final state = "+str(abs(diff_S))+"%  diff I final state = "+str(abs(diff_I))+"%  diff R final state= "+str(abs(diff_R))+"%")


##Plot the solutions


#Plot average outcome of fast SIR algorithm
fig1 = plt.figure(facecolor='w')
fig1.suptitle('Fast SIR solution for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma)+" averaged "+str(sol_number)+" times")
ax1 = fig1.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#mean of every simulation
ax1.plot(t_mean, S_mean/N, 'b', linestyle='-', alpha=1, lw=2, label='Fast SIR : Omnivorous')
ax1.plot(t_mean, I_mean/N, 'r', linestyle='-', alpha=1, lw=2, label='Fast SIR : Vegetarian')
ax1.plot(t_mean, R_mean/N, 'g', linestyle='-',alpha=1, lw=2, label='Fast SIR : Ex-vegetarian')

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

#Plot average outcome of Gillespie algorithm
fig2 = plt.figure(facecolor='w')
fig2.suptitle('Gillespie solution for N= '+str(N)+", beta = "+str(beta)+", gamma = "+str(gamma)+" averaged "+str(sol_number)+" times")
ax2 = fig2.add_subplot(111, facecolor='#dddddd', axisbelow=True)

#mean of every simulation
ax2.plot(t_mean2, S_mean2/N, 'b', linestyle='--', alpha=1, lw=2, label='Gillespie : Omnivorous')
ax2.plot(t_mean2, I_mean2/N, 'r', linestyle='--', alpha=1, lw=2, label='Gillespie : Vegetarian')
ax2.plot(t_mean2, R_mean2/N, 'g', linestyle='--',alpha=1, lw=2, label='Gillespie : Ex-vegetarian')

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

#show both figures
plt.show()

