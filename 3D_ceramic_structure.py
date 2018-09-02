import numpy as np
import random
import scipy.integrate as integrate
import lxml.etree as lxml

class material_properties():

    def __init__(self,cube_data):
        parser=lxml.XMLParser()
        tree = lxml.parse(cube_data, parser)
        root=tree.getroot()
        
        self.cube_dimentions=np.asarray([int(i) for i in root[1][0].get('D').split(",")])
        self.symmetry=str(root[0][0].get('symmetry'))
        if self.symmetry=="Tetra":
            self.tau=0.75e-11*2/1.3
            self.eps=np.array([[1721,0,0],[0,1721,0],[0,0,382]])
            self.Ea=33e6
            self.Ps=0.3641
            self.v=1.8
        elif self.symmetry=="Rhombo":
            self.tau=0.17e-11
            self.eps=np.array([[529,0,0],[0,529,0],[0,0,295]])
            self.Ea=25e6
            self.Ps=0.455
            self.v=1.4
        elif self.symmetry=="Ortho":
            self.tau=0.17e-12
            self.eps=np.array([[160,0,0],[0,100,0],[0,0,55]])
            self.Ea=8e6
            self.Ps=0.42
            self.v=1
        self.v=1.8 if self.symmetry=="Tetra" else 1.4
        self.Ea=33e6 if self.symmetry=="Tetra" else 25e6
        self.Ps=0.3641 if self.symmetry=="Tetra" else 0.455

    def assignpoldirection(self):
        if self.symmetry=="Tetra":
            #tetragonal
            theta_r1=np.pi/4
            theta_max=np.arcsin(np.sqrt(2.0/3.0))
            def f1(x):
                return 1.5/np.pi
            def f2(x):
                return 6*(np.pi/4-np.arccos(1/np.tan(x)))/np.pi**2
            #print 2*np.pi*(integrate.quad(lambda x: f1(x)*np.sin(x),0,np.pi/4)[0]+integrate.quad(lambda x: f2(x)*np.sin(x),np.pi/4,np.arcsin(np.sqrt(2.0/3.0)))[0])
            j=0
            i=[]
            f_theta_int=[]
            theta=[]
            while j<=theta_max:
                if j<=np.pi/4:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,j)[0])
                    theta.append(f1(j))
                    i.append(j)
                else:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f2(x)*np.sin(x),np.pi/4,j)[0]+integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,np.pi/4)[0])
                    theta.append(f2(j))
                    i.append(j)
                j+=np.arcsin(np.sqrt(2.0/3.0))/10000.0
        elif self.symmetry=="Rhombo":
            #rhombohedra
            theta_r1=np.arctan(1/np.sqrt(2))
            theta_max=np.arctan(np.sqrt(2))
            def f1(x):
                return 2/np.pi
            def f2(x):
                return 6*(np.pi/3-np.arccos(1/(np.sqrt(2)*np.tan(x))))/np.pi**2
            #print (2*np.pi*(integrate.quad(lambda x: f1(x)*np.sin(x),0,theta_r1)[0]+integrate.quad(lambda x: f2(x)*np.sin(x),theta_r1,theta_max)[0]))
            j=0
            i=[]
            f_theta_int=[]
            theta=[]
            while j<=theta_max:
                if j<=theta_r1:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,j)[0])
                    theta.append(f1(j))
                    i.append(j)
                else:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f2(x)*np.sin(x),theta_r1,j)[0]+integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,theta_r1)[0])
                    theta.append(f2(j))
                    i.append(j)
                j+=theta_max/10000.0
        elif self.symmetry=="Ortho":
            #orthorhombic
            theta_max=np.pi/4
            theta_0=np.arctan(1/np.sqrt(2))
            phi_0=np.arctan(np.sqrt(2))
            def f1(x):
                return 3/np.pi
            def f2(x):
                return 6*(np.pi/2-2*np.arccos(1/(np.sqrt(3)*np.tan(x))))/np.pi**2
            def f3(x):
                return 6*(phi_0-np.arccos(1/(np.sqrt(3)*np.tan(x))))/np.pi**2
            print (integrate.quad(lambda x: 2*np.pi*f3(x)*np.sin(x),theta_0,np.pi/4)[0]+\
                        integrate.quad(lambda x: 2*np.pi*f2(x)*np.sin(x),np.pi/6,theta_0)[0]+\
                        integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,np.pi/6)[0])
            j=0
            i=[]
            f_theta_int=[]
            theta=[]
            while j<=theta_max:
                if j<=np.pi/6:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,j)[0])
                    theta.append(f1(j))
                    i.append(j)
                elif j>np.pi/6 and j<=theta_0:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f2(x)*np.sin(x),np.pi/6,j)[0]+\
                    integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,np.pi/6)[0])
                    theta.append(f2(j))
                    i.append(j)
                elif j>theta_0 and j<=np.pi/4:
                    f_theta_int.append(integrate.quad(lambda x: 2*np.pi*f3(x)*np.sin(x),theta_0,j)[0]+\
                    integrate.quad(lambda x: 2*np.pi*f2(x)*np.sin(x),np.pi/6,theta_0)[0]+\
                    integrate.quad(lambda x: 2*np.pi*f1(x)*np.sin(x),0,np.pi/6)[0])
                    theta.append(f2(j))
                    i.append(j)
                j+=theta_max/10000.0
        ang_sel_points=[]
        random_point=np.random.random_sample((1,int(np.prod(self.cube_dimentions))))[0]
        for i in range(np.prod(self.cube_dimentions)):
            for j in range(len(f_theta_int)-1):
                if random_point[i]>=f_theta_int[j] and random_point[i]<f_theta_int[j+1]:
                    ang_sel_points.append(j)
                    break
        ones=np.ones(np.prod(self.cube_dimentions))
        ones[:np.prod(self.cube_dimentions)/2]=-1
        random.shuffle(ones)
        self.one=ones
        self.theta=np.array(ang_sel_points)*theta_max*ones/10000.0
        self.phi=np.random.random_sample((1,np.prod(self.cube_dimentions)))[0]*np.pi*2
        self.psi=np.random.random_sample((1,np.prod(self.cube_dimentions)))[0]*np.pi*2
        p_direction=np.array([np.sin(self.theta)*np.cos(self.phi),np.sin(self.theta)*np.sin(self.phi),np.cos(self.theta)]).T
        p_direction=np.reshape(p_direction,(self.cube_dimentions[2],self.cube_dimentions[1],self.cube_dimentions[0],3))
        eps0=[]
        for i,j,k in zip(self.phi,self.theta,self.psi):
            eps0.append(self.eps_t(i,j,k))
        epsilon=np.reshape(np.asarray(eps0),(self.cube_dimentions[2],self.cube_dimentions[1],self.cube_dimentions[0],3,3))
        print (np.mean(np.cos(self.theta)))
        return p_direction, epsilon

    def eps_t(self,phi,theta,psi):

        e=np.matrix(self.eps)
        m_phi=np.array([[np.cos(phi),np.sin(phi),0],[-np.sin(phi),np.cos(phi),0],[0,0,1]])
        e_phi=np.dot(np.dot(m_phi,e),m_phi.T)
        m_theta=np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])
        e_theta_phi= np.dot(np.dot(m_theta,e_phi),m_theta.T)
        m_psi=np.array([[np.cos(psi),np.sin(psi),0],[-np.sin(psi),np.cos(psi),0],[0,0,1]])
        return np.dot(np.dot(m_psi,e_theta_phi),m_psi.T)
