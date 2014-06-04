import cvxopt   as cv
import numpy    as np
import scipy    as sp
import chompack as cp

def getY(SYS):
    # generates the network admittance matrix Y from SYS
    # assumption: bus numbers start at 1,...,n
    N = len(SYS['bus'][::,0])
    Y = cv.spmatrix([],[],[],(N,N), tc='z')
    for e,i in enumerate(SYS['branch'][::,0]):
        j = int(SYS['branch'][e,1])   #bus of 'to' - line
        ind_i = int(np.where(SYS['bus'][::,0]==i)[0][0] -1)
        ind_j = int(np.where(SYS['bus'][::,0]==j)[0][0] -1)
        Y[ind_i,ind_j] = -1/(SYS['branch'][e,2] +1j*SYS['branch'][e,3])-1j/2*SYS['branch'][e,4]
        Y[ind_j,ind_i] = Y[ind_i,ind_j]
        
    for e in range(N):
        Y[e,e] = -sum(Y[e,::])+SYS['bus'][e,4]+1j*SYS['bus'][e,5]
    return Y

    
def getyk(Y,k):
    # generates matrix Y_k, admittance matrix with only k-th row.
    # is used to generate (2.2a)/(2.2b) in Molzahn dissertation
    N = Y.size[0]
    ek = cv.spmatrix([1],[k],[0],(N,1))
    return ek*ek.trans()*Y
    
def getYklm(yk):
    # generates YY_k according to Molzahn dissertation (2.2a) and (2.2b)
    Yk     =  cv.sparse(0.5*cv.matrix([[(yk+yk.trans()).real(), (yk-yk.trans()).imag()],[(yk.trans()-yk).imag(), (yk+yk.trans()).real()]]))
    Yk_bar =  cv.sparse(-0.5*cv.matrix([[(yk+yk.trans()).imag(), (yk.trans()-yk).real()],[(yk-yk.trans()).real(), (yk+yk.trans()).imag()]]))
    return Yk, Yk_bar

def getJk(k,N):
    # generates matrix M_k according to Molzahn dissertation (2.2c)
    ek = cv.spmatrix([1],[k],[0],(N,1))
    D = ek*ek.trans()
    return cv.spdiag([D,D])
    
def getLlm(l,m,N):
    # generates matrix L_lm according to myThesis
    el = cv.spmatrix([1],[l-1],[0],(N,1))
    em = cv.spmatrix([1],[m-1],[0],(N,1))
    D = el*el.trans() + em*em.trans() - el*em.trans() - em*el.trans()
    return cv.spdiag([D,D])

def getNlm(l,m,Y):
    # generates matrix N and N_bar for line capacity constraint with real Power limit
    N = Y.size[0]
    el = cv.spmatrix([1],[l-1],[0],(N,1))
    em = cv.spmatrix([1],[m-1],[0],(N,1))
    G = Y.real()[l-1,m-1]
    B = Y.imag()[l-1,m-1]
    D1 = em*em.trans()-0.5*em*el.trans()-0.5*el*em.trans()
    D2 = 0.5*em*el.trans()-0.5*el*em.trans()
    Z = cv.spmatrix([],[],[],(D2.size))
    T1 = cv.spdiag([D1,D1])
    T2 = cv.sparse([[Z,D2],[-D2,Z]])
    N = G*T1
    Nb = B*T2
    return N, Nb
    
def getTlm(l,m,N,theta):
    # generates matrix T for line capacity constraint using maximum angle difference
    el = cv.spmatrix([1],[l-1],[0],(N,1))
    em = cv.spmatrix([1],[m-1],[0],(N,1))
    T = cv.spdiag([-0.5*theta*cv.matrix([el*em.trans()+em*el.trans()]),0.5*cv.matrix([el*em.trans()+em*el.trans()])])
    return cv.sparse(T)
    

def SystemInitPyPower(mySys):
    # Initializes the equality/inequality constraints and the objective function if one of PyPower cases is used
    # e.g. mySys = case9()
    #
    
    Y = getY(mySys)
    N = Y.size[0]
    C = []
    c = []
    B = []
    b = []
    for i in range(N):
        # Maximum/minimum bus voltage limit on bus i
        Vbase = mySys['bus'][i,9]
        M = getJk(i,N)
        C.append(M)
        C.append(-M)
        Vi_max = mySys['bus'][i,11]*Vbase # maximum voltage amplitude
        Vi_min = mySys['bus'][i,12]*Vbase# minimum voltage amplitude
        c.append( Vi_max**2)
        c.append(-Vi_min**2)
        
        
        # Only for buses with attached generators:
        
        bus = mySys['bus'][i,0] # bus number
        #print 'bus is: ',bus
        
        if bus in list(mySys['gen'][::,0]):
            
            ind_g = list(mySys['gen'][::,0]).index(bus)
            #print 'generator is at bus: ', (ind_g+1)
            
            yk = getyk(Y,i)
            Y_k, Y_kb = getYklm(yk)
            
            #Maximum/minimum real power netto power generation
            C.append( Y_k)
            C.append(-Y_k)
            PGk_max = mySys['gen'][ind_g,8]
            PGk_min = mySys['gen'][ind_g,9]
            PDk = mySys['bus'][i,2]
            c.append( PGk_max - PDk)
            c.append(-PGk_min + PDk)
            
            #Maximum/minimum reactive power netto generation
            C.append( Y_kb)
            C.append(-Y_kb)
            QGk_max = mySys['gen'][ind_g,3]
            QGk_min = mySys['gen'][ind_g,4]
            QDk = mySys['bus'][i,3]
            c.append( QGk_max - QDk)
            c.append(-QGk_min + QDk)
        
        # Only buses without generators
        # Power demand at each bus
        else:
            yk = getyk(Y,i)
            Y_k, Y_kb = getYklm(yk)
            B.append(Y_k)
            B.append(Y_kb)
            PD_i = mySys['bus'][i,2]
            QD_i = mySys['bus'][i,3]
            b.append(-PD_i)
            b.append(-QD_i)
    
    #for i,l in enumerate(mySys['branch'][::,0]):
        ## Line capacity according to Madani2013 (6d)
        #l = int(l)
        #m = int(mySys['branch'][i,1])
        #Llm = getLlm(l,m,N)
        #alpha = mySys['branch'][i,9]
        #dV = 2*(1-sp.cos(alpha))
        #C.append(Llm)
        #c.append(dV**2)
        
    A = cv.spdiag(cv.matrix(1,(2*N,1)))
    c = cv.matrix(np.asarray(c))

    return A, B, b, C, c


def SystemInitCSV(source):
    
    Datapath = './Data' + source
    Branch = np.genfromtxt(Datapath + '/Branch.csv', delimiter=',', skip_header=1, autostrip=1, usecols= (0,2,3,4,5,6,7,8,9,10))
    Bus    = np.genfromtxt(Datapath + '/Bus.csv', delimiter=',', skip_header=1, autostrip=1, usecols= (0,2,3,4,5,6,7,8,9,10,11,12,13))
    Gen    = np.genfromtxt(Datapath + '/Generator.csv', delimiter=',', skip_header=1, autostrip=1, usecols= (0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
    Load   = np.genfromtxt(Datapath + '/Load.csv', delimiter=',', skip_header=1, autostrip=1, usecols= (0,2,3,4,5,6,7,8,9,10,11))
    
    
    
    N_Br = len(Branch)
    N_Bu = len(Bus)
    N_G = len(Gen)
    N_L = len(Load)
    
    mySys = {'branch' : cv.matrix(0,(N_Br,13), tc = 'd'),
             'bus' : cv.matrix(0,(N_Bu,13), tc = 'd'),
             'gen' : cv.matrix(0,(N_G,21), tc = 'd')}
    
    for i in range(N_Br):
        mySys['branch'][i,0] = int(Branch[i][1])
        mySys['branch'][i,1] = int(Branch[i][2])
        mySys['branch'][i,2] = float(Branch[i][3])
        mySys['branch'][i,3] = float(Branch[i][4])
        mySys['branch'][i,4] = float(Branch[i][5])
    
    for i in range(N_Bu):
        mySys['bus'][i,0] = int(Bus[i][0])
        mySys['bus'][i,4] = float(Bus[i][2])
        mySys['bus'][i,5] = float(Bus[i][3])
        mySys['bus'][i,11] = float(Bus[i][9])
        mySys['bus'][i,12] = float(Bus[i][10])
        mySys['bus'][i,9] = float(Bus[i][7])
    
    for i in range(N_G):
        mySys['gen'][i,0] = int(Gen[i][1])
        mySys['gen'][i,3] = float(Gen[i][11])
        mySys['gen'][i,4] = float(Gen[i][12])
        mySys['gen'][i,8] = float(Gen[i][2])
        mySys['gen'][i,9] = float(Gen[i][3])
        
    for i in range(N_L):
        b = int(Load[i][1])
        mySys['bus'][b-1,2] = Load[i][2]
        mySys['bus'][b-1,3] = Load[i][3]
        
    return mySys


########################################################################
# Maybe not needed functions
def getylm(Y,SYS,l,m):
    N = Y.size[0]
    ind = myindex(SYS['branch'],[l,m])
    l -= 1
    m -= 1
    Ylm = cv.spmatrix([],[],[],(N,N),tc = 'z')
    Ylm[l,l] = 0.5*1j*SYS['branch'][ind,4] + Y[l,m]
    Ylm[l,m] = -Y[l,m]
    return Ylm
    
def myindex(L,X):
    I = []
    for e in range(len(L)):
        if (L[e,0] == X[0]) and (L[e,1] == X[1]):
            I.append(e)
            
    return I
