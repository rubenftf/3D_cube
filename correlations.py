import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime

class correlation:

    def __init__(self, numer):
        self.symmetry="Rhombo"                                                 #choose a symmetry Tetra, Rhombo, Ortho
        self.dimention_x, self.dimention_y, self.dimention_z=5, 5, 5
        self.anglestep_p=1080
        self.anglestep_t=1080
        self.sizestep=200
        self.numer=numer
        self.path_to="/nfshome/khachaturyan/Publication/3D_cube/cube_simulation/calculated/"
        self.data_file="E_P_cor_{}.out".format(self.symmetry,self.numer)
        with open(self.path_to+self.data_file,"r") as f:
            content=f.readlines()
            for line in content:
                if "Time" in line:
                    self.t_n=float(line.split(' ')[-1])

    def correlations(self):                                                    #correlation calculation for z and x components, as y assumed to be equal to x

        P,E,p0=self.readFlexpdeoutput(self.path_to+self.data_file)

        OxxR=[]
        OzzR=[]      
        FxxR=[]
        FzzR=[]

        rho=np.linspace(0,1,self.sizestep)*0.4*self.dimention_x                #mesh creation (correaltions is calculated among these points)
        phi=np.linspace(0,1.01,self.anglestep_p)*2*np.pi*0
        theta=np.linspace(0,1.01,self.anglestep_t)*2*np.pi
        rho_new=np.linspace(0,1,2*self.sizestep)*0.4*self.dimention_x          #mesh densification for smoth plotm used in the cor_graph function 
        phi_new=np.linspace(0,1.01,2*self.anglestep_p)*2*np.pi*0
        theta_new=np.linspace(0,1.01,2*self.anglestep_t)*2*np.pi

        self.Zrho=rho                                                          #also used in the cor_graph function 
        self.Zphi=phi
        self.Ztheta=theta
        self.Zrho_new=rho_new
        self.Zphi_new=phi_new
        self.Ztheta_new=theta_new

        Px1=[[] for _ in range(self.sizestep)]                                 #generate a statich map where at each point P and E values are assigned to
        Pz1=[[] for _ in range(self.sizestep)]                                 #P1 and E2 arrays
        Ex1=[[] for _ in range(self.sizestep)]
        Ez1=[[] for _ in range(self.sizestep)]
        for m in range(self.sizestep):
            for z in range(self.dimention_x):
                for y in range(self.dimention_y):
                    for x in range(self.dimention_z):
                        for i in range(8):
                            Px1[m].append(P[z][y][x][i][0])
                            Pz1[m].append(P[z][y][x][i][2])
                            Ex1[m].append(E[z][y][x][i][0])
                            Ez1[m].append(E[z][y][x][i][2])

        for n in range(self.anglestep_t):                                      #for different angles
            Px2=[[] for _ in range(self.sizestep)]                             #generate similar to P1 and E1 arrays, but with values corresponding to
            Pz2=[[] for _ in range(self.sizestep)]                             #a point shifted on (r,theta,phi) radious form
            Ex2=[[] for _ in range(self.sizestep)]                             # the static map posicion (line 157 in the code)
            Ez2=[[] for _ in range(self.sizestep)]                             #then spatial correlations are looking between this maps and the static one

            for m in range(self.sizestep):

                r=rho[m]*np.array([np.sin(theta[n])*np.cos(0),np.sin(theta[n])*np.sin(0),np.cos(theta[n])])

                for z in range(self.dimention_x):
                    for y in range(self.dimention_x):
                        for x in range(self.dimention_x):
                            for i in range(8):

                                p=p0[z][y][x][i]+r

                                if p[0]<0:                                     #periodic boundary canditions
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

                                Px2[m].append(P[int(p[2])][int(p[1])][int(p[0])][nighbor][0])
                                Pz2[m].append(P[int(p[2])][int(p[1])][int(p[0])][nighbor][2])
                                Ex2[m].append(E[int(p[2])][int(p[1])][int(p[0])][nighbor][0])
                                Ez2[m].append(E[int(p[2])][int(p[1])][int(p[0])][nighbor][2])

            OxxR.append(self.corcoefficientcalc(Px2,Px1,rho))
            OzzR.append(self.corcoefficientcalc(Pz2,Pz1,rho))
            FxxR.append(self.corcoefficientcalc(Ex2,Ex1,rho))
            FzzR.append(self.corcoefficientcalc(Ez2,Ez1,rho))

        self.cor_graph("P_x_{}_{}".format(self.p_n,self.t_n),OxxR)
        self.cor_graph("P_z_{}_{}".format(self.p_n,self.t_n),OzzR)
        self.cor_graph("E_x_{}_{}".format(self.p_n,self.t_n),FxxR)
        self.cor_graph("E_z_{}_{}".format(self.p_n,self.t_n),FzzR)

    def readFlexpdeoutput(self,flexdeoutputfile):                              #fills the content by values Ex,Ey,Ez,Px,Py,Pz
        with open(flexdeoutputfile,'r') as f:                                  #and correspondic coordinates x,y,z
            content=f.readlines()                                              #generate arrays of pointvalueE,pointvalueP,pointvalue
            Ez=[]                                                              #each cube supposed to contain 8 points (x,y,z)
            Ey=[]
            Ex=[]
            Pz=[]
            Py=[]
            Px=[]
            pointcoordinate_z=[]
            pointcoordinate_y=[]
            pointcoordinate_x=[]
            for line in content:
                if len(line) > 2:
                    if "VAL(xcomp(E" in line:
                        Ex.append(float(line.split()[1]))
                    elif "VAL(ycomp(E" in line:
                        Ey.append(float(line.split()[1]))
                    elif "VAL(zcomp(E" in line:
                        Ez.append(float(line.split()[1]))
                    elif "VAL(px" in line:
                        Px.append(float(line.split()[1]))
                    elif "VAL(py" in line:
                        Py.append(float(line.split()[1]))
                    elif "VAL(pz" in line:
                        Pz.append(float(line.split()[1]))
                    elif "VAL(pz" in line:
                        pointcoordinate_z.append(float(line.split(",")[3].split("*")[0]))
                        pointcoordinate_y.append(float(line.split(",")[2].split("*")[0]))
                        pointcoordinate_x.append(float(line.split(",")[1].split("*")[0]))
            pointvalueE=zip(Ex,Ez,Ez)
            pointvalueP=zip(Px,Pz,Pz)
            pointvalue=zip(pointcoordinate_x,pointcoordinate_y,pointcoordinate_z)
        return np.reshape(pointvalueP,(self.dimention_x, self.dimention_y, self.dimention_z,8,3)),\
            np.reshape(pointvalueE,(self.dimention_x, self.dimention_y, self.dimention_z,8,3)),\
            np.reshape(pointvalue,(self.dimention_x, self.dimention_y, self.dimention_z,8,3))

    def cor_graph(self,name,Zv):                                               #plots correaltion graph
        z=np.asarray(Zv)

        b=self.Zrho
        a=self.Ztheta+np.pi/2                                                  #+np.pi/2 is need to orient the graph in the field direction
        bn=self.Zrho_new                                                       #so fild is directed up (correspond to a pi/2 angle)
        an=self.Ztheta_new+np.pi/2

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
        plt.savefig('{}.jpg'.format(name), dpi=600, bbox_inches="tight", format="jpg")

    def corcoefficientcalc(self,E1,E2,rho):                                    #calculates correlations alog the selected radious
        cor=[[] for _ in range (len(rho))]
        for i in range(len(rho)):
            cor[i]=np.corrcoef(E1[i],E2[i])[1,0]
            #np.cov(E1[i],E2[i])[1,0]                                          #optionaly the covarience can be calculated instead of correlation coefficient
        return np.asarray(cor)

    def corcoefficient(self,name,rho,E):                                       #creates a file with correlations graph data, if want to plote in other program
        data = np.zeros(rho.size, dtype=[('var{}'.format(i),float ) for i in range(1,self.anglestep_p+1)])
        data['var1']=rho*1e6
        for j in range(2,self.anglestep_p+1):
            data['var{}'.format(j)]=E[j-1]
        f='{}'.format("%10.3f "*(self.anglestep_p))
        np.savetxt('{}_cor.txt'.format(name), data, fmt=f)

start_at=datetime.datetime.now()
print ("began at: ",start_at)
number=33                                                                      #select the file number
Cor=correlation(number)
Cor.correlations()
print ("spent time: {}".format(datetime.datetime.now()-start_at))
