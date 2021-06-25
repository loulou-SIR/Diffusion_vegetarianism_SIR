import networkx as nx
import matplotlib.pyplot as plt
import EoN
import xlrd
import random as rd
from matplotlib import pyplot as plt
import numpy as np
import csv


'''This program
IMPORTANT:
    1) This program imports data from the stata_data.xlsx file, so line #20 should be changed according to the path of the file on your computer
    2) The last section allows the user to export weights list in .csv but this is optional; here again, path should be modified
    3)This program needs to be run before the custom fast SIR algorithm
'''

##import data of prevalence and odds ratiofor 24 subcategories

wb = xlrd.open_workbook('/Users/louisdretzolis/Desktop/Proyecto Python/stata_data.xlsx') #ADD YOUR OWN PATH
sh = wb.sheet_by_index(2) #open third sheet

#create 3 dictionaries:
# - one for prevalence of subcategory among whole population
# - one for % of vegetarianism among subcategory
# - one for odds ratio of vegetarianism of subcategory

dict_prevalence,dict_vegetarian,dict_odds={},{},{}

for i in range(1,25):
    if sh.cell(i,0).value==i:
        dict_prevalence[i]=sh.cell(i,5).value
        dict_vegetarian[i]=sh.cell(i,7).value
        dict_odds[i]=sh.cell(i,8).value


##Create dictionary of cumulated prevalences for 24 subcategories

dict_cumulated_prevalence={1:dict_prevalence[1]}#we initialize the dictionary with the first prevalence

for i in range(2,25):
    dict_cumulated_prevalence[i]=dict_cumulated_prevalence[i-1]+dict_prevalence[i]


###Sub-function key prob used by the graph generator
#we want a function that takes in entry a dictionary of keys with their probabilty and return the key
def key_prob(d):
    rand=rd.random()
    for k in d.keys():
        if rand<d[k]:
            return k


##Generator of custom graph (function used by custom fast SIR algorithm)

def generate_custom(n1,k1,p1):
    cust=nx.watts_strogatz_graph(n=n1,k=k1,p=p1)

    for i in range(n1):#go through all the nodes of graph
        node_subcat=key_prob(dict_cumulated_prevalence)
        #calculate vegetarian probability by considering subcategory
        vege_prob=dict_vegetarian[node_subcat]
        nx.set_node_attributes(cust,values = {i:node_subcat}, name='subcategory')
        if rd.random()<vege_prob:
            nx.set_node_attributes(cust,values = {i:1}, name='vege')
        else:
            nx.set_node_attributes(cust,values = {i:0}, name='vege')

    weight_graph(cust)

    return cust

##Function that calcultates weights of edges depending on attributes


#Global variables used by weight function
max_weight= max(dict_odds.values()) #weight in most favorable scenario
norm=1/max_weight #normalization factor

def weight_norm(node_subcategory):#gives normalized weight in function of subcategory
    weight_raw= dict_odds[node_subcategory]#raw weight
    weight_norm=weight_raw*norm #normalized weight
    return weight_norm

def weight_graph(g):
    #adds normalized weights depending on subcategory to all nodes
    for j in range(g.number_of_nodes()):
        node_subcat=g.node[j]['subcategory']
        g.node[j]['weight']=weight_norm(node_subcat)
    #adds weights equal to 0 to all edges
    for i, j in g.edges():
        g[i][j]['weight'] = 0



'''
##Export the normalized weights to .csv

W2=[]
for i in range(1,25):
    W2.append(weight_norm(i))

with open('/Users/louisdretzolis/Desktop/Proyecto Python/weighted.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(W2)
'''


