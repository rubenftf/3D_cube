import os
import matplotlib.pylab as plt
import six
from matplotlib import colors
import numpy as np
import random

path_to_Dist=os.getcwd()
colors_ = list(six.iteritems(colors.cnames))
random.shuffle(colors_)

def graphs(E):

    step=[]
    for file in os.listdir(path_to_Dist):
        if file.endswith(".txt"):
            if "./distributions/field_distributions/calculated/{}_".format(E) in file:
                step.append(int(file.split("_")[2]))
    step=np.sort(step)
    files=[]
    for i in step:
        files.append('./distributions/field_distributions/calculated/{}_{}_.txt'.format(E,i))

    X=[]
    for i in files:
        x,y=np.loadtxt(i, unpack=True)
        X.append((x,y))

    t,P=np.loadtxt("./distributions/field_distributions/data/cube.p04", unpack=True)

    stepj=[]
    maxE=[]
    fig, ax = plt.subplots(1,4, figsize=(24, 8))
    for j,i in enumerate(X):
        ax[0].plot(i[0], i[1], color="{}".format(colors_[j][0]), 
        label="{}({:1.2f}) at {:1.4f}".format(E,max(i[1]),step[j]))
        stepj.append(step[j])
        maxE.append(max(i[1]))
        ax[0].legend(loc=1);    
    ax[0].set_xlabel('{}'.format(E))
    ax[0].set_ylabel('frequency')
    ax[0].set_title('Distribution {}'.format(E))

    x_x=[]
    for j,e in enumerate(X):
        for i in range(len(e[1])):
            if e[1][i]==max(e[1]):
                nm=i
        Eav=float(max(e[1])/2)
        val=100
        for i in range(nm):
            if np.absolute(Eav-e[1][i])<val:
                na1=i
                val=np.absolute(Eav-e[1][i])
        val=100
        for i in range(nm,len(e[1])):
            if np.absolute(Eav-e[1][i])<val:
                na2=i
                val=np.absolute(Eav-e[1][i])
        x_x.append(e[0][na2]-e[0][na1])
    ax[1].plot(stepj,x_x)
    ax[1].set_xlabel('P/Ps')
    ax[1].set_xticks(np.linspace(-0.8,0.8,3))   
    ax[1].grid(True)    
    ax[1].set_ylabel('half-width')
    ax[1].set_title('Half-width')

    ax[2].plot(t,P, label="max={:1.2f}".format(max(P)))
    ax[2].legend(loc=1);
    ax[2].set_xlim(0,max(t)+1)
    ax[2].set_ylim(1,1)
    ax[2].set_xlabel('t/t0')
    ax[2].set_ylabel('P/Ps')
    ax[2].grid(True)
    ax[2].set_title('Polarization-time')
    P_rate=[]
    t_rate=[]
    for i in range(1,len(P)):
        P_rate.append(P[i]-P[i-1])
        t_rate.append(t[i]-t[i-1])
    P_rate=np.asarray(P_rate)
    t_rate=np.asarray(t_rate)
    np.delete(t,0)
    ax[3].plot(t,P_rate/t_rate)
    ax[3].grid(True)
    ax[3].set_title('Derivative')    

    fig.tight_layout()
    fig.show()

def dist_graph(E):

    path_to_Dist="./distributions/field_distributions/calculated"
    step=[]
    for file in os.listdir(path_to_Dist):
        if file.endswith(".txt"):
            if "{}_".format(E) in file:
                step.append(float(file.split("_")[1]))
    step=np.sort(step)
    print (step)
    files=[]
    for i in step:
        files.append('./distributions/field_distributions/calculated/{}_{}_.txt'.format(E,i))
    X=[]
    for i in files:
        x,y=np.loadtxt(i, unpack=True)
        X.append((x,y))
    t,P=np.loadtxt("/nfshome/khachaturyan/Publication/3D_cube/cube_simulation/distributions/field_distributions/data/cube_up.p01", unpack=True)

    maxE=[]
    fig, ax = plt.subplots(figsize=(12, 10))
    for j,i in enumerate(X):
        ax.plot(i[0], i[1], label=r'{}$P_s$'.format(round(step[j],10)),linewidth=3)
        maxE.append(max(i[1]))
        print (max(i[1]))
        ax.legend(loc=1,fontsize=20);
    ax.set_xlabel(r'${}/E_0$'.format(E),fontsize=35)
    ax.set_ylabel('frequency',fontsize=35)
    plt.tick_params(labelsize=20)
    plt.xlim(0.4,1.4)
    plt.ylim(0, 4)
    plt.savefig("./distributions/field_distributions/calculated/{}_1.jpg".format(E),dpi=600)

dist_graph("Em")