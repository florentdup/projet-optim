import numpy as np
import matplotlib.pyplot as plt

alpha=1.
beta=10.

lambdapdc=4000.
lambdac=1.

Dm=0.3
DM=1.5
Pm=45.
PM=68.
Wm=0.#verifier unitées!
WM=100.

#Pour les canalisations

L=np.array([0.5,0.5,1,0.25]) 
D0=np.array([0.395,0.395,0.5,0.5])
Q=np.array([0.5,0.5,-((2-0.5)-1),2-0.5])



noeuds=[(0,1),(2,3),(3,4),(0,4)]

Qsc=np.array([0.5])
noeudssc=[(2,3)]
Winit=np.array([0.])


Nsc=1
Ncana=4
Nnoeuds=5

nbvar=Ncana+Nsc+Nnoeuds
nbcomp=4*Ncana+6*Nsc


def f(Dd,W):
    return alpha*np.dot(L,Dd)+beta*np.sum(W)

def c(x):
    #x=Dd,W,P
    X=np.split(x,[Ncana, Ncana+Nsc,Ncana+Nsc+Nnoeuds])  
    Dd=X[0]
    W=X[1]
    P=X[2]
    C=np.zeros((4*Ncana+6*Nsc))
    pi=P**2

    for i in range(Ncana):
        if Dd[i]<0:
            Dd[i]=0
        C[i]=(pi[noeuds[i][1]]-pi[noeuds[i][0]])-lambdapdc*L[i]*Q[i]**2/(pow(D0[i],0.4)+pow(Dd[i],0.4))

    for i in range(Ncana):
        if Dd[i]<0:
            Dd[i]=0
        C[i+Ncana]=(pi[noeuds[i][0]]-pi[noeuds[i][1]])-lambdapdc*L[i]*Q[i]**2/(pow(D0[i],0.4)+pow(Dd[i],0.4))

    for i in range(Nsc):
        C[i+2*Ncana]=lambdac*Qsc[i]*np.log(pi[noeudssc[i][1]]/pi[noeudssc[i][0]])-Winit[i]-W[i]
    
    for i in range(Nsc):
        C[i+2*Ncana+Nsc]=pi[noeudssc[i][0]]-pi[noeudssc[i][1]]
    
    for i in range(Nsc):
        C[i+2*Ncana+2*Nsc]=P[noeudssc[i][1]]-PM
        C[i+2*Ncana+3*Nsc]=Pm-P[noeudssc[i][1]]

    for i in range(Ncana):
        C[i+2*Ncana+4*Nsc]=Dd[i]-DM
        C[i+3*Ncana+4*Nsc]=Dm-Dd[i]
        

    for i in range(Nsc):
        if W[i]<0:
            W[i]=0.
        C[i+4*Ncana+4*Nsc]=W[i]-WM
        C[i+4*Ncana+5*Nsc]=Wm-W[i]

    return C




def dcdx(X):
    J=np.zeros((nbvar,nbcomp))
    for i in range(nbvar):
        dx=np.ndarray.flatten(np.eye(1,nbvar,i)*step_calculJ)
        dcdi=(c(X+dx)-c(X))/step_calculJ
        J[i]=dcdi
    return J

dfdx=np.concatenate([alpha*L,np.array([beta]),np.zeros((Nnoeuds))])


    
#def L(Dd,W,P,lambdaa):
#    return f(Dd,W)-np.dot(lambdaa,c(Dd,W,P))



def ArrowHurwicz(Dd0,W0,P0,lambdaa0,epsilon=1e-3,alpha_arrow=1e-3,maxiter=10):
    N=0
    Dd=np.copy(Dd0)
    W=np.copy(W0)
    P=np.copy(P0)
    x=np.concatenate([Dd,W,P])
    lambdaa=np.copy(lambdaa0)
    continuer=True

    while N<maxiter and continuer:
        N+=1
        increment_x=-epsilon*(dfdx+np.dot(dcdx(x),lambdaa))
        x+=increment_x
        lambdaa=Proj(lambdaa+alpha_arrow*c(x))
        

        continuer=(np.linalg.norm(increment_x)>seuil1)

        

        
    return x,N,lambdaa




def Proj(x):
    x[x < 0] = 0.
    return x


Dd0=np.zeros((Ncana))+Dm
W0=np.zeros((Nsc))+Wm
P0=np.zeros((Nnoeuds))+Pm
lambda0=np.zeros((nbcomp))+1.


x=np.concatenate([Dd0,W0,P0])

step_calculJ=1e-4
seuil1=1e-3


x,N,lambdaa=ArrowHurwicz(Dd0,W0,P0,lambda0,1e-2,1e-2,5000) 

print("Nombre d'itérations :", N)
x=np.split(x,[Ncana, Ncana+Nsc,Ncana+Nsc+Nnoeuds])  
print("Dd: \n",x[0])
print("W: \n",x[1])
print("P: \n",x[2])
print("cout: \n",f(x[0],x[1]))




