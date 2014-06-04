import cvxopt   as cv
import numpy    as np
import scipy    as sp
import chompack as cp
import networkx as nx
import time
import random as rdm
from cvxopt import solvers
from Decomposition import *
from interfacePP import *
from pypower.api import case9, case14, case30, case39, case57,case300, ppoption, runpf, printpf

import matplotlib.pyplot as plt

cv.printing.options['width'] = -1
cv.printing.options['dformat'] = '%3.5f'
print "Test run"



N = 10
elapsedX9 = cv.matrix(0,(N,1),tc = 'd')
els9 = cv.matrix(0,(N,1),tc = 'd')
elsX9 = cv.matrix(0,(N,1),tc = 'd')
time9 = cv.matrix(0,(N,4),tc = 'd')

elapsedX14 = cv.matrix(0,(N,1),tc = 'd')
els14 = cv.matrix(0,(N,1),tc = 'd')
elsX14 = cv.matrix(0,(N,1),tc = 'd')
time14 = cv.matrix(0,(N,4),tc = 'd')

elapsedX30 = cv.matrix(0,(N,1),tc = 'd')
els30 = cv.matrix(0,(N,1),tc = 'd')
elsX30 = cv.matrix(0,(N,1),tc = 'd')
time30 = cv.matrix(0,(N,4),tc = 'd')

elapsedX39 = cv.matrix(0,(N,1),tc = 'd')
els39 = cv.matrix(0,(N,1),tc = 'd')
elsX39 = cv.matrix(0,(N,1),tc = 'd')
time39 = cv.matrix(0,(N,4),tc = 'd')

elapsedX300 = cv.matrix(0,(N,1),tc = 'd')
els300 = cv.matrix(0,(N,1),tc = 'd')
elsX300 = cv.matrix(0,(N,1),tc = 'd')
time300 = cv.matrix(0,(N,4),tc = 'd')

#for i in range(N):
    #mySys = case300()
    #start  = time.time()
    #A, B, b, C, c = SystemInitPyPower(mySys)
    #time300[i,1] = time.time() - start
    #sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    #time300[i,0] = time.time() - start
    #time300[i,2] = t_decomp
    #time300[i,3] = t_solver

for i in range(N):
    mySys = case9()
    start  = time.time()
    A, B, b, C, c = SystemInitPyPower(mySys)
    time9[i,1] = time.time() - start
    sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    time9[i,0] = time.time() - start
    time9[i,2] = t_decomp
    time9[i,3] = t_solver

    mySys = case9()
    A, B, b, C, c = SystemInitPyPower(mySys)
    start = time.time()
    sol, e = SDPsolverX(A,C=C,c=c, B=B, b=b)
    elapsedX9[i] = (time.time()- start)
    elsX9[i] = e


for i in range(N):
    mySys = case14()
    start  = time.time()
    A, B, b, C, c = SystemInitPyPower(mySys)
    time14[i,1] = time.time() - start
    sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    time14[i,0] = time.time() - start
    time14[i,2] = t_decomp
    time14[i,3] = t_solver

    mySys = case14()
    A, B, b, C, c = SystemInitPyPower(mySys)
    start = time.time()
    sol, e = SDPsolverX(A,C=C,c=c, B=B, b=b)
    elapsedX14[i] = (time.time()- start)
    elsX14[i] = e
    

for i in range(N):
    
    mySys = case30()
    start  = time.time()
    A, B, b, C, c = SystemInitPyPower(mySys)
    time30[i,1] = time.time() - start
    sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    time30[i,0] = time.time() - start
    time30[i,2] = t_decomp
    time30[i,3] = t_solver

    mySys = case30()
    A, B, b, C, c = SystemInitPyPower(mySys)
    start = time.time()
    sol, e = SDPsolverX(A,C=C,c=c, B=B, b=b)
    elapsedX30[i] = (time.time()- start)
    elsX30[i] = e


for i in range(N):
    mySys = case39()
    start  = time.time()
    A, B, b, C, c = SystemInitPyPower(mySys)
    time39[i,1] = time.time() - start
    sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    time39[i,0] = time.time() - start
    time39[i,2] = t_decomp
    time39[i,3] = t_solver

    mySys = case39()
    A, B, b, C, c = SystemInitPyPower(mySys)
    start = time.time()
    sol, e = SDPsolverX(A,C=C,c=c, B=B, b=b)
    elapsedX39[i] = (time.time()- start)
    elsX39[i] = e


for i in range(N):
    mySys = case300()
    start  = time.time()
    A, B, b, C, c = SystemInitPyPower(mySys)
    time300[i,1] = time.time() - start
    sol, I, t_solver, t_decomp = SDPsolver(A, B=B, b=b, C=C, c=c)
    time300[i,0] = time.time() - start
    time300[i,2] = t_decomp
    time300[i,3] = t_solver
    mySys = case300()
    A, B, b, C, c = SystemInitPyPower(mySys)
    start = time.time()
    sol, e = SDPsolverX(A,C=C,c=c, B=B, b=b)
    elapsedX300[i] = (time.time()- start)
    elsX300[i] = e

import plotting
