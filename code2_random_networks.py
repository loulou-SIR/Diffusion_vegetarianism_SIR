import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

'''In this code, we generate and plot 4 random networks

'''

##Network parameters
#Mean degree
K=4
#number of nodes
N=50
#probability p for Erdos-Renyi graph
prob=K/N
#probability of rewiring for small world graph
prob_rewire=0.3
#parameter of power law for scale-free network
m=2


##Random k-regular graph

g=nx.random_regular_graph(d=K,n=N)

pos = nx.spring_layout(g)
plt.figure(figsize=(20,20))
#colors = [f[u][v]['color'] for u,v in edges]
nx.draw_networkx(g, pos, with_labels = False, width=1,node_color='b', node_size=250, font_size=10)
plt.axis('off')
plt.title('Random 4-regular network')


##Erdos Renyi network

f=nx.erdos_renyi_graph(n=N,p=prob)

pos = nx.spring_layout(f)
plt.figure(figsize=(20,20))
#colors = [f[u][v]['color'] for u,v in edges]
nx.draw_networkx(f, pos, with_labels = False, width=1,node_color='g', node_size=250, font_size=10)
d=dict(f.degree)
mean_degree=sum(d.values())/float(len(f))
plt.axis('off')
plt.title('Erdos Renyi network with average degree np = '+str(mean_degree))

##Small world network (Watts Strogatz)

e=nx.watts_strogatz_graph(n=N,k=K,p=prob_rewire)

pos = nx.circular_layout(e)
plt.figure(figsize=(20,20))
#colors = [f[u][v]['color'] for u,v in edges]
nx.draw_networkx(e, pos, with_labels = False, width=1,node_color='y', node_size=250, font_size=10)
d=dict(e.degree)
mean_degree=sum(d.values())/float(len(e))
plt.axis('off')
plt.title('Small world network with average degree k = '+str(mean_degree))

##Scale free network (also called power law)

h=nx.barabasi_albert_graph(n=N,m=m)

pos = nx.spring_layout(h)
plt.figure(figsize=(20,20))
#colors = [f[u][v]['color'] for u,v in edges]
nx.draw_networkx(h, pos, with_labels = False, width=1,node_color='m', node_size=250, font_size=10)
d=dict(h.degree)
mean_degree=sum(d.values())/float(len(h))
plt.axis('off')
plt.title('Scale-free network with average degree k = '+str(mean_degree))


##Show graphs

plt.show()
