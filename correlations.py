from __future__ import division
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime

class correlation:

    def __init__(self, numer):

        self.symmetry="Ortho"
        self.dimention_x, self.dimention_y, self.dimention_z=5, 5, 5
        self.anglestep_p=1080
        self.anglestep_t=1080
        self.sizestep=200
        self.numer=numer
        with open("//nfshome/khachaturyan/Publication/3D_cube/cube_simulation/cube/E_{}.txt".format(self.numer),"r") as f:
            content=f.readlines()
            for line in content:
                if "Time" in line:
                    self.t_n=float(line.split(' ')[-1])
        #t,p=np.loadtxt('/nfshome/khachaturyan/Publication/3D_cube/cube_simulation/calculated/{0}_cor/{0}.p01'.format(self.symmetry), unpack=True)
        #self.p_n=round(p[np.where(t==self.t_n)][0],3)
        self.p_n=0

    def readFlexpdeoutput(self,flexdeoutputfile):
        Ez=[]
        Ex=[]
        rz=[]
        rx=[]
        with open(flexdeoutputfile,'r') as f:
            s=0
            for line in f:
                s+=1
                a=line.strip()
                if s>8:
                    Ex.append(float(a.split(" ")[0]))
                    Ez.append(float(a.split("  ")[1].split("   ")[0]))
                    rx.append(1e6*float(a.split("   ")[1]))
                    rz.append(1e6*float(a.split("   ")[-1]))
            pointvalueE=zip(Ex,Ez,Ez)
            print (np.shape(pointvalueE))
            pointvalue=zip(rx,rz,rz)
        return np.reshape(pointvalueE,(self.dimention_x, self.dimention_y, self.dimention_z,8,3)),\
            np.reshape(pointvalue,(self.dimention_x, self.dimention_y, self.dimention_z,8,3))

    def corcoefficientcalc(self,E1,E2,rho):

        cor=[[] for _ in range (len(rho))]
        for i in range(len(rho)):
            cor[i]=np.corrcoef(E1[i],E2[i])[1,0]#np.cov(E1[i],E2[i])[1,0]
        return np.asarray(cor)

    def corcoefficient(self,name,rho,E):

        data = np.zeros(rho.size, dtype=[('var{}'.format(i),float ) for i in range(1,self.anglestep_p+1)])
        data['var1']=rho*1e6
        for j in range(2,self.anglestep_p+1):
            data['var{}'.format(j)]=E[j-1]
        f='{}'.format("%10.3f "*(self.anglestep_p))
        np.savetxt('{}_cor.txt'.format(name), data, fmt=f)

    def cor_graph(self,name,Zv):

        z=np.asarray(Zv)

        b=self.Zrho
        a=self.Ztheta+np.pi/2#self.Zphi
        bn=self.Zrho_new
        an=self.Ztheta_new+np.pi/2#self.Zphi_new

        Bng=[]
        for i in z:
            f=interp1d(b, i, kind='cubic')
            Bng.append(f(bn))

        z=np.array(Bng).T
        x,y=np.meshgrid(a,bn)
        m=len(x)
        val_nxest = int(max(3+np.sqrt(m/2),2*3+3))
        tck = interpolate.bisplrep(y, x, z, s=3000, kx=3, ky=3, nxest=val_nxest)
        xn,yn=np.meshgrid(an,bn)
        zn = interpolate.bisplev(yn[:,0], xn[0,:], tck)

        for i in range(len(zn)):
            for j in range(len(zn[i])):
                if zn[i][j]>1.0:
                    zn[i][j]=1
                elif zn[i][j]<-1:
                    zn[i][j]=-1

        zn=zn.T
        for a in range(3):
            for i in range(-int(len(zn)/2),int(len(zn)/2)):
                zn[i]=0.5*(zn[int(len(zn)/2)-np.sign(i)*i]+zn[i])
            for i in range(-int(len(zn)/2),int(len(zn)/2)):
                zn[i]=0.5*(zn[-i]+zn[i])
            a+=1
        zn=zn.T

        zn_min=1
        zn_max=-1
        for i in z:
            for j in i:
                if j<zn_min:
                    zn_min=j
                if j>zn_max:
                    zn_max=j

        matplotlib.use('Agg')
        plt.figure()
        xn,yn=np.meshgrid(an,bn)
        plt.subplots(subplot_kw=dict(projection='polar'))
        matplotlib.rcParams.update({'font.size': 15, 'font.family': 'serif'})

        plt.pcolor(xn, yn, zn, vmin=0, vmax=1)
        plt.ioff()
        plt.savefig('cor_pic/{}.jpg'.format(name), dpi=600, bbox_inches="tight", format="jpg")

    def correlations(self):
        
        path_to="/nfshome/khachaturyan/Publication/3D_cube/cube_simulation/cube/"
        #E,p0=self.readFlexpdeoutput(path_to+'{}_cor/E_P_cor_{}.out'.format(self.symmetry,self.numer))
        E,p0=self.readFlexpdeoutput(path_to+'E_{}.txt'.format(self.numer))
     
        FxxR=[]
        FzzR=[]

        rho=np.linspace(0,1,self.sizestep)*0.4*self.dimention_x
        phi=np.linspace(0,1.01,self.anglestep_p)*2*np.pi*0
        theta=np.linspace(0,1.01,self.anglestep_t)*2*np.pi
        rho_new=np.linspace(0,1,2*self.sizestep)*0.4*self.dimention_x
        phi_new=np.linspace(0,1.01,2*self.anglestep_p)*2*np.pi*0
        theta_new=np.linspace(0,1.01,2*self.anglestep_t)*2*np.pi

        self.Zrho=rho
        self.Zphi=phi
        self.Ztheta=theta

        self.Zrho_new=rho_new
        self.Zphi_new=phi_new
        self.Ztheta_new=theta_new

        Ex1=[[] for _ in range(self.sizestep)]
        Ez1=[[] for _ in range(self.sizestep)]
        for m in range(self.sizestep):
            for z in range(self.dimention_x):
                for y in range(self.dimention_y):
                    for x in range(self.dimention_z):
                        for i in range(8):
                            Ex1[m].append(E[z][y][x][i][0])
                            Ez1[m].append(E[z][y][x][i][2])

        for n in range(self.anglestep_t):

            Ex2=[[] for _ in range(self.sizestep)]
            Ez2=[[] for _ in range(self.sizestep)]

            for m in range(self.sizestep):

                r=rho[m]*np.array([np.sin(theta[n])*np.cos(0),np.sin(theta[n])*np.sin(0),np.cos(theta[n])])

                for z in range(self.dimention_x):
                    for y in range(self.dimention_y):
                        for x in range(self.dimention_z):
                            for i in range(8):

                                p=p0[z][y][x][i]+r

                                if p[0]<0:
                                    p[0]+=self.dimention_x
                                if p[0]>self.dimention_x:
                                    p[0]-=self.dimention_x
                                if p[1]<0:
                                    p[1]+=self.dimention_y
                                if p[1]>self.dimention_y:
                                    p[1]-=self.dimention_y
                                if p[2]<0:
                                    p[2]+=self.dimention_z
                                if p[2]>self.dimention_z:
                                    p[2]-=self.dimention_z

                                nighbor=int((p[2]-int(p[2]))/0.5)*4+int((p[1]-int(p[1]))/0.5)*2+int((p[0]-int(p[0]))/0.5)

                                Ex2[m].append(E[int(p[2])][int(p[1])][int(p[0])][nighbor][0])
                                Ez2[m].append(E[int(p[2])][int(p[1])][int(p[0])][nighbor][2])

            FxxR.append(self.corcoefficientcalc(Ex2,Ex1,rho))
            FzzR.append(self.corcoefficientcalc(Ez2,Ez1,rho))

        self.cor_graph("E_x_{}_{}".format(self.p_n,self.t_n),FxxR)
        self.cor_graph("E_z_{}_{}".format(self.p_n,self.t_n),FzzR)

start_at=datetime.datetime.now()
print ("began at: ",start_at)
Cor=correlation(1)
Cor.correlations()
print ("spent time: {}".format(datetime.datetime.now()-start_at))