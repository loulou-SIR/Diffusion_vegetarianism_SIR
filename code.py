import networkx as nx
import matplotlib.pyplot as plt
import EoN
import xlrd

#to do : faire en sorte que les weights dépendent des attributs de chaque noeud (graphe dirigé ?)

## Diccionary
d = {}
wb = xlrd.open_workbook('/Users/louisdretzolis/Desktop/Proyecto Python/data_graph.xlsx')
sh = wb.sheet_by_index(0)
for i in range(1,15):
    id = sh.cell(i,0).value
    name = sh.cell(i,1).value
    age = sh.cell(i,2).value
    genre = sh.cell(i,3).value #genre 1 : woman, genre 2 : man
    friend1 = sh.cell(i,4).value
    weight1 = sh.cell(i,5).value
    friend2 = sh.cell(i,6).value
    weight2 = sh.cell(i,7).value
    d[id] = [name, age, genre,friend1,weight1,friend2,weight2]
##Graph
f=nx.Graph()

for i in d:
    print(i)
    print(d[i])
    f.add_node(i)
    if d[i][3]!='':
        f.add_edge(i,d[i][3],weight=d[i][4])
    if d[i][5]!='':
        f.add_edge(i,d[i][5],weight=d[i][6])

#storing weights in a list
edges = f.edges()
weights = [f[u][v]['weight'] for u,v in edges]

##Model SIR
gamma = 0.2
beta = 1.2
r_0 = beta/gamma
print(r_0)
N = len(f) # population size
I0 = 1   # intial n° of infected individuals
R0 = 0
S0 = N - I0 -R0


nx_kwargs = {"with_labels":True, "pos": nx.spring_layout(f), "width": weights, "alpha": 0.7} #optional arguments to be passed on to the networkx plotting command.
print("doing Gillespie simulation")

sim = EoN.Gillespie_SIR(f, tau = beta, gamma=gamma, rho = I0/N, transmission_weight="weight", return_full_data=True)
print("done with simulation, now plotting")

for i in range(0,7,1):
    sim.display(time = i,  **nx_kwargs)
    plt.axis('off')
    plt.title("Iteration {}".format(i))
    plt.draw()
    plt.show()
