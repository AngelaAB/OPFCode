import cvxopt   as cv
import numpy    as np
import scipy    as sp
import chompack as cp
from copy import copy


def rank(A, eps=1e-12):
    u, s, vh = np.linalg.svd(A)
    return len([x for x in s if abs(x) > eps])
    
def ind(x,y,l):
    return sum(range(l+1))-sum(range(l+1-y)) +x-y

def symsparsmatrix(size, dens):
    den = size**2*dens*0.5
    Dx = np.random.random_integers(0,size-1,den)    # generates random x coordinates
    Dy = np.random.random_integers(0,size-1,den)    # generates random y coordinates
    De = np.random.random_integers(0,100,den)       # generates random matrix entries for coordinate (x,y)

    D = cv.spmatrix(De,Dx,Dy,(size,size)) 
    D += D.trans() + cv.spmatrix(np.random.random_integers(0,50,size),range(size),range(size))

    return D    # matrix D is a symmetric sparse matrix of dimension (size,size) with density approx. dens


def vec2mat(x):
    L = x.size[0]
    N = int(max(-0.5-sp.sqrt(1+8*L)/2,-0.5+sp.sqrt(1+8*L)/2))
    M = cv.matrix(0,(N,N),tc='d')
    
    s = cv.matrix(0,(N,1))
    for i in range(N-1):
        s[i+1] = s[i] + N - i
    for i in range(N):
        for j in range(i+1):
            e = x[s[j]+i-j,0]
            M[i,j] = e
            M[j,i] = M[i,j]
    return M


def multvec2mat(X,I):
    s = I.size
    if s[0] < s[1]:
        I = I.trans()
    if I.size[1] != 1:
        print 'Dimension of index vector not correct.\n'
    M = list()
    for i in range(I.size[0]-1):
        x = X[I[i]:I[i+1]:]
        m = vec2mat(x)
        M.append(m)

    return M
