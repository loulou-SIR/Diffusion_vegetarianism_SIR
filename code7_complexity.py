import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN
import time

'''IMPORTANT: this code is very long to run
'''
## Complexity calculation for both algorithms

times_gillespie=[]
times_fastSIR=[]
step=1000
min_population=1000
max_population=100000
avg=50
population_size=[l for l in range(min_population,max_population,step)]
id_log_population=[id_log(x) for x in population_size]

for l in range(min_population ,max_population,step):
    print("evolution : "+str(l)+"/"+str(max_population))

    #creating graph
    g=nx.watts_strogatz_graph(n=l, k=4, p=0.6)

    #initializing random weights
    E = g.number_of_edges()
    w = [random.random() for i in range(E)]
    s = max(w)
    w = [ i/s for i in w ] #normalizing
    k = 0
    for i, j in g.edges():
        g[i][j]['weight'] = w[k]
        k+=1
    #SIR parameters
    I0, R0 = 0.05*N, 0
    S0 = N - I0 - R0
    beta, gamma = 1,0.6

    #Gillespie algo done 50 times
    time_elapsed_mean=0
    for i in range(avg):
        time_start= time.perf_counter() #start time of simulation
        t,S,I,R= EoN.Gillespie_SIR(g, tau = beta, gamma=gamma, rho = I0/N, transmission_weight="weight", return_full_data=False)
        time_elapsed = (time.perf_counter() - time_start) #duration of simulation
        time_elapsed_mean+=time_elapsed
    times_gillespie.append(time_elapsed_mean/avg)

    #fast SIR algo done 50 times
    time_elapsed_mean1=0
    for i in range(avg):
        time_start1= time.perf_counter() #start time of simulation
        t,S,I,R= EoN.fast_SIR(g, tau = beta, gamma=gamma, rho = I0/N, transmission_weight="weight", return_full_data=False)
        time_elapsed1 = (time.perf_counter() - time_start1) #duration of simulation
        time_elapsed_mean1+=time_elapsed1
    times_fastSIR.append(time_elapsed_mean1/avg)

##Plot

fig = plt.figure(facecolor='w')
fig.suptitle("Simulation times for N from "+str(min_population)+" nodes to "+str(max_population)+" nodes with step "+str(step))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(population_size, times_gillespie, 'b', alpha=0.5, lw=2, label='Gillespie')
ax.plot(population_size, times_fastSIR, 'g', alpha=0.5, lw=2, label='Fast SIR')
ax.set_xlabel('network size')
ax.set_ylabel('simulation time (s)')
ax.set_ylim(0,max(max(times_gillespie),max(times_fastSIR)))
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

