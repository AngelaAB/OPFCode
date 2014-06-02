import cvxopt   as cv
import numpy    as np
import scipy    as sp
import chompack as cp
from copy import copy
import time

from mylib import *
        
def sparsedecomp(A,LB,LC):
    
    # check if A is square matrix  
    if A.size[0] != A.size[1]:
        print "Matrix dimensions does not match. A matrix must be square."
        return
 
    # check if other matrices have same dimension as A
    if LB != []:
        for i in LB:
            if i.size != A.size:
                print "Dimension of equality matrices (LB) do not match. Must have same dimension as A matrix"
                return

    # check if other matrices have same dimension as A
    if LC != []:
        for i in LC:
            if i.size != A.size:
                print "Dimension of inequality matrices (LC) do not match. Must have same dimension as A matrix"
                return
    


    # generate maximum sparsity pattern
    M = copy(A)
    for i in LC:
        M += abs(i)
    for i in LB:
        M += abs(i)
#print M
        
    MC = cp.symbolic(cv.sparse(M), p=cv.amd.order)
    
    # All maximal cliques of objective 
    # matrix A are extracted and returned as a list  
    dA = clique_extract(MC,A)
    
    # All maximal cliques of all (in-)equality constraint matrices are extracted and returned as a list
    dLB = []
    if LB != []:
        for i in LB:
            dLB.append(clique_extract(MC,i))


    # All maximal cliques of all (in-)equality constraint matrices are extracted and returned as a list
    dLC = []
    if LC != []:
        for i in LC:
            dLC.append(clique_extract(MC,i))
        
    # permutation matrix
    P = cv.spmatrix(1,MC.p,range(MC.sparsity_pattern().size[0]))
    
    E, H = eq_const(MC)
      
    # dA  = list containing all cliques of A [CHECK]
    # dLB = 2D-List containing all cliques of every equality constraint [CHECK]
    # list of equality contraints
    # permutation list/matrix [CHECK]

    return dA, dLB, dLC, P, E, H


def clique_extract(AC,A):

    # creates permutation matrix, to permute original matrix A
    P = cv.spmatrix(1,AC.p,range(AC.sparsity_pattern().size[0])) # permutation matrix
    AP = P.trans()*A*P
    
    # all maximal cliques of A are extracted 
    dA = []
    for i in AC.cliques():
        dA.append(AP[i,i])
    
    return dA   # list of cliques is returned


def reformSDP(A):
    
    RA=[]
    for i in A:
        for j in range(i.size[0]):
            RA.extend(cv.matrix(i[j::,j]))
    RA = cv.matrix(RA)
    
    return RA 
    
def SDPsolverX(A, B = [], b = [], C = [], c = []):
    '''
    Solve an SDP of standard form
    
    :param A: Matrix cost function.
    :param B: Matrix equality constraint.
    :param b: Vector equality cosntrint.
    :param C: Matrix inequality constraint.
    :param c: Vector inequality constraint.
    :type A: cv.matrix.
    :type B: [cv.matrix].
    :type b: [cv.matrix].
    :type C: [cv.matrix].
    :type c: [cv.matrix].
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
    '''
    from cvxopt import solvers
    
    s = A.size[0]
    A = cv.sparse(A+A.trans()-cv.spmatrix(cv.matrix(A[::s+1]),range(s),range(s)))


    myc = reformSDP([A])
    if B != []:
        for i,cl in enumerate(B):
            s = int(cl.size[0])
            B[i] = cl+cl.trans()-cv.spmatrix(cv.matrix(cl[::s+1]),range(s),range(s))

        myA = copy(myc)
        for cl in B:
            a = reformSDP([cl])
            myA = cv.matrix([[myA],[a]])
        myA = cv.sparse(myA[::,1::].trans())
        myb = cv.matrix(b,tc = 'd')
    else:
        myA = []
        myb = []
    
    if C != []:    
        for i,cl in enumerate(C):
            s = cl.size[0]
            C[i] = cl+cl.trans()-cv.spmatrix(cv.matrix(cl[::s+1]),range(s),range(s))
        myGl = copy(myc)
        for cl in C:
            g = reformSDP([cl])
            myGl = cv.matrix([[myGl],[g]])
        myGl = cv.sparse(myGl[::,1::].trans())
        myg = copy(c)
    else:
        myGl = []
        myg = []
    
    myG, myh = psd_cond([A])
    
    if myA and myGl:
        print 'A ok'
        print 'myGl ok'
        start = time.time()
        sol = solvers.sdp(c = myc, A = myA, b  = myb, Gs = myG,  hs = myh, Gl = myGl, hl = myg)
        elapsed = time.time() - start
    elif not myA and myGl:
        print 'A not ok'
        print 'myGl ok'
        start = time.time()
        sol = solvers.sdp(c = myc, Gs= myG, hs = myh, Gl = myGl, hl = myg)
        elapsed = time.time() - start
    elif myA and not myGl:
        print 'A ok'
        print 'myGl not ok'
        start = time.time()
        sol = solvers.sdp(c = myc,A=myA, b=myb, Gs=myG, hs=myh)
        elapsed = time.time() - start
    else:
        print 'A not ok'
        print 'myGl not ok'
        start = time.time()
        sol = solvers.sdp(c = myc,Gs=myG, hs=myh)
        elapsed = time.time() - start
        
    return sol, elapsed


def SDPsolver(A, B = [], b = [], C = [], c = []):
    #print 'A is: ', A
    #print 'S is: ', S
    #print 'B is: ', B
    #print 'b is: ', b
    #print 'C is: ', C
    #print 'c is: ', c
    
    from cvxopt import solvers
    start = time.time()
    # Decomposition of matrices A,B,C according to sparsity pattern
    # dA, dB, dC are lists (of lists) of submatrices of A,B,C
    # P is sparsity pattern
    for i,cl in enumerate(B):
        s = int(cl.size[0])
        B[i] = cv.sparse(cl+cl.trans()-cv.spmatrix(cv.matrix(cl[::s+1]),range(s),range(s)))
    for i,cl in enumerate(C):
        s = cl.size[0]
        C[i] = cv.sparse(cl+cl.trans()-cv.spmatrix(cv.matrix(cl[::s+1]),range(s),range(s)))
        
    s = A.size[0]
    A = cv.sparse(A+A.trans()-cv.spmatrix(cv.matrix(A[::s+1]),range(s),range(s)))
    
    
    
    dA, dB, dC, P, E, H = sparsedecomp(A,B,C)
    I = decomp_index(dA)
    # Transforms list of matrices to list of vectors with matrices in column-first order
    myc = reformSDP(dA)
    
    if dB != []:

        myA = copy(myc)
        for cl in dB:
            a = reformSDP(cl)
            myA = cv.matrix([[myA],[a]])
        myA = cv.sparse(myA[::,1::].trans())
        myb = cv.matrix(b, tc='d')
        if E!= []:
            for i in E:
                myA = cv.matrix([myA,i])
            myb = cv.matrix([cv.matrix(b),H])
    else:
        if E!= []:
            myA = copy(myc).trans()
            for i in E:
                myA = cv.matrix([myA,i])
            myA = cv.sparse(myA[1::,::])
            myb = H
        else:
            myA = []
            myb = []

    if dC != []:    
        
        myGl = copy(myc)
        for cl in dC:
            g = reformSDP(cl)
            myGl = cv.matrix([[myGl],[g]])
        myGl = cv.sparse(myGl[::,1::].trans())
        myg = copy(c)
    else:
        myGl = []
        myg = []

    myG, myh = psd_cond(dA)
    
    e1 = time.time()- start
 
    if myA and myGl:
        print 'A ok'
        print 'myGl ok'
        start = time.time()
        sol = solvers.sdp(c = myc, A = myA, b  = myb, Gs = myG,  hs = myh, Gl = myGl, hl = myg)
        elapsed = time.time()-start
    elif not myA and myGl:
        print 'A not ok'
        print 'myGl ok'
        start = time.time()
        sol = solvers.sdp(c = myc, Gs= myG, hs = myh, Gl = myGl, hl = myg)
        elapsed = time.time()-start
    elif myA and not myGl:
        print 'A ok'
        print 'myGl not ok'
        start = time.time()
        sol = solvers.sdp(c = myc,A=myA, b=myb, Gs=myG, hs=myh)
        elapsed = time.time()-start
    else:
        print 'A not ok'
        print 'myGl not ok'
        start = time.time()
        sol = solvers.sdp(c = myc,Gs=myG, hs=myh)
        elapsed = time.time()-start

    return sol, I, elapsed, e1

def psd_cond(A):
    # M: total number of optimization parameters
    # A: cliques structure
    
    N = [0]
    for i,c in enumerate(A):
        N.append(N[i]+(c.size[0]**2+c.size[0])/2)
    
    G = []
    h = []

    for i,c in enumerate(A):
        u = N[-1] - N[i+1]
        l = N[i]
        
        Y = cv.matrix(0,(1,c.size[0]**2))
        k = 0
        for y in range(c.size[0]):
            for x in range(c.size[0]):
                if x >= y:
                    Y[k] = ind(x,y,c.size[0])
                else:
                    Y[k] = ind(y,x,c.size[0])
                k += 1
        
        Z1 = cv.spmatrix([],[],[],(c.size[0]**2,l))
        Z2 = cv.spmatrix(-1,range(c.size[0]**2),Y)
        Z3 = cv.spmatrix([],[],[],(c.size[0]**2,u))
        Z = cv.sparse([[Z1],[Z2],[Z3]])
        h.append(cv.spmatrix([],[],[],(c.size[0],c.size[0])))            
        G.append(Z)
        
    return G, h


def eq_const(AC):


    # This function generates the matrix of equality constraints
    
    # index vector of offset for cliques
    # N[-1]: total length of optimization parameter x
    N = [0]
    for i,c in enumerate(AC.cliques()):
        N.append(N[i]+(len(c)**2+len(c))/2)
    
    # make tuple out of separators with true indices according to cliques
    E = []
    for cc, pc in enumerate(AC.parent()):
        if pc != -1:
            C = AC.cliques()[cc]
            P = AC.cliques()[pc]
            S = AC.separators()[cc]
            
            T = []
            for i in S:
                T.append([C.index(i),P.index(i)])
            
            e = []                    
            for i,x in enumerate(T):
                for y in T[:i+1:]:
                    e = cv.spmatrix([],[],[],(1,N[-1]))
                    e[0,N[cc] + ind(x[0],y[0],len(C))] =  1
                    e[0,N[pc] + ind(x[1],y[1],len(P))] = -1
                    E.append(e)
    
    if E != []:
        Z = cv.spmatrix([],[],[],(len(E),1))
    else:
        Z = []
        
    return E, Z

def decomp_index(dA):
    I = [i.size[0]*(i.size[0]+1)/2 for i in dA]
    J = cv.matrix(0,(len(I)+1,1), tc='i')
    
    for i in range(len(I)-1):
        J[i+1] = J[i] + I[i]
    J[-1] = I[-1] + J[-2]
    
    return J
