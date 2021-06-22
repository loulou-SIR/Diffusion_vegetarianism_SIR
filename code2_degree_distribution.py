import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import random
import EoN

##Network parameters
#Mean degree
K=15
#number of nodes
N=10000
#probability p for Erdos-Renyi graph
prob=K/N
#probability of rewiring for small world graph
prob_rewire=0.3
#parameter of power law for scale-free network
m=3


##Random regular graph

g=nx.random_regular_graph(d=K,n=N)

degree_freq = nx.degree_histogram(g)
degrees = range(len(degree_freq))
d=dict(g.degree)
mean_degree=sum(d.values())/float(len(g))
plt.figure(figsize=(20, 20))
plt.loglog(degrees[m:], degree_freq[m:],'go-')
plt.xlabel('Degree     Average degree = '+str(mean_degree))
plt.ylabel('Frequency')

plt.title('Degree distribution of random regular graph')

##Erdos Renyi network

f=nx.erdos_renyi_graph(n=N,p=prob)

degree_freq = nx.degree_histogram(f)
degrees = range(len(degree_freq))
d=dict(f.degree)
mean_degree=sum(d.values())/float(len(f))
plt.figure(figsize=(20, 20))
plt.plot(degrees[m:], degree_freq[m:],'go-')
plt.xlabel('Degree     Average degree = '+str(mean_degree))
plt.ylabel('Frequency')
plt.title('Degree distribution of Erdos Renyi network')

##Small world network (Watts Strogatz)

e=nx.watts_strogatz_graph(n=N,k=K,p=prob_rewire)

degree_freq2 = nx.degree_histogram(e)
degrees2 = range(len(degree_freq2))
d=dict(e.degree)
mean_degree=sum(d.values())/float(len(e))
plt.figure(figsize=(20, 20))
plt.plot(degrees2[m:], degree_freq2[m:],'yo-')
plt.xlabel('Degree      Average degree = '+str(mean_degree))
plt.ylabel('Frequency')
plt.title('Degree distribution of small world network')

##Scale free network (also called power law)

h=nx.barabasi_albert_graph(n=N,m=m)

degree_freq3 = nx.degree_histogram(h)
degrees3 = range(len(degree_freq3))
d=dict(h.degree)
mean_degree=sum(d.values())/float(len(h))
plt.figure(figsize=(20, 20))
plt.loglog(degrees3[m:], degree_freq3[m:],'mo-')
plt.xlabel('Degree     Average degree = '+str(mean_degree))
plt.ylabel('Frequency')
plt.title('Degree distribution of scale-free network')

##Show graphs

plt.show()

