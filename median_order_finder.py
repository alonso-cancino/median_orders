import numpy as np
import networkx as nx
import pandas as pd
from pyscipopt import Model,quicksum
from random import shuffle
from networkx.algorithms.tournament import random_tournament as rt
import time
import math
from tqdm import trange

def adjacency_matrix(t,order):
    n = len(order)
    adj_t = np.zeros((n,n))
    for e in t.edges:
        adj_t[order[e[0]],order[e[1]]]= 1 
    return adj_t

def random_tournament(n):
    r_t = rt(n)
    adj_t = adjacency_matrix(r_t,list(range(n)))
    return r_t, adj_t

def ratio(adj_T):
    r = 0
    n = adj_T.shape[0]
    for i in range(n):
        for j in range(i,n):
            r += adj_T[i,j]
    return r

def find_median_order_rt(n,verbose=True):
    t,adj_t = random_tournament(n)
    
    start_setup = time.time()
    model = Model()
    p,w,r = {},{},{}

    for k in range(n):
        for l in range(n):
            p[k,l] = model.addVar(vtype='B')
            for i in range(n):
                for j in range(i,n):
                    r[i,k,l,j] = model.addVar(vtype='B')
                    w[i,k,l,j] = adj_t[k][l]

    for i in range(n):
        # Forcing p to be a permutation 
        model.addCons(quicksum(p[s,i] for s in range(n))==1)
        model.addCons(quicksum(p[i,s] for s in range(n))==1)
        for k in range(n):
            for j in range(i,n):
                for l in range(n):
                    # Setting r[i,k,l,j] = min(p^t[i,k],p[l,j])
                    model.addCons(r[i,k,l,j] <= p[k,i])
                    model.addCons(r[i,k,l,j] <= p[l,j])

    model.setObjective(quicksum(r[i,k,l,j]*w[i,k,l,j] for i in range(n) for j in range(i,n) for k in range(n) for l in range(n)), "maximize")
    end_setup = time.time()
    
    if verbose == True:
        print(f'Variables setup, took {"{:.2f}".format(end_setup-start_setup)}s')
        print(f'Input tournament had adjacency matrix: \n\n{adj_t}')
    model.data = p,r
    model.optimize()
    sol = model.getBestSol()
    Q = np.array([math.floor(model.getVal(model.data[0][key])) for key in model.data[0].keys()]).reshape([n,n])
    end_optimization = time.time()
    
    if verbose == True:
        print(f'\nOptimization ended with status {model.getStatus()} in {"{:.2f}".format(end_optimization-end_setup)}s, with {model.getObjVal()} increasing edges and optimal solution:')
        print('\n',Q)
    order = [int(x) for x in list(np.matmul(Q.T,np.array(range(n))))]
    new_adj_t = np.matmul(Q.T,np.matmul(adj_t,Q))
    
    if verbose == True:
        print(f'\nwhich induces the ordering:\n\n {order}')
        print(f'\nand induces the following adjacency matrix: \n\n {new_adj_t}')

    return adj_t,new_adj_t,order

find_median_order_rt(8)
