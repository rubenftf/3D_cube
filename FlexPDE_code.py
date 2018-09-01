import numpy as np
import os
import subprocess
import lxml.etree as lxml
import xml.etree.ElementTree as ET
import datetime

class flex_code():

    def __init__(self, material_properties):
        self.material_properties=material_properties
        self.cube_directory=os.getcwd()+'/cube'
        self.errlim0=10
        self.Ea=self.material_properties.Ea
        self.tau0=self.material_properties.tau
        self.Ps=self.material_properties.Ps
        self.factor=self.material_properties.symmetry

    def read_start(self,start_data_file):
        self.restart=False
        angles_dist=self.material_properties.assignpoldirection()
        self.p_direcations=angles_dist[0]
        self.epsilon=angles_dist[1]
        self.cubes_num=self.material_properties.cube_dimentions
        self.p_initial=-np.ones(np.prod(self.cubes_num)).reshape(self.cubes_num[2],self.cubes_num[1],self.cubes_num[0])
        self.voltage=self.cubes_num[2]*self.material_properties.v
        self.theta=self.material_properties.theta
        self.phi=self.material_properties.phi
        self.t_initial=0

    def read_restart(self,restart_data_file):
        self.restart=True
        tree = ET.parse(restart_data_file)
        root = tree.getroot()
        self.cubes_num=[int(root[0][0].get("Dx")),int(root[0][0].get("Dy")),int(root[0][0].get("Dz"))]
        self.voltage=float(root[0][0].get("V"))
        p_direcations=[]
        theta=[]
        phi=[]
        psi=[]
        del root[0][0]
        grain=root[0]
        for i in grain:
            p_direcationsx=float(i.get('p_direcationsx'))
            p_direcationsy=float(i.get('p_direcationsy'))
            p_direcationsz=float(i.get('p_direcationsz'))
            p_direcations.append([p_direcationsx,p_direcationsy,p_direcationsz])
            theta.append(float(i.get('theta')))
            psi.append(float(i.get('psi')))
            phi.append(float(i.get('phi')))
        p_direcations=np.asarray(p_direcations).reshape(self.cubes_num[2],self.cubes_num[1],self.cubes_num[0],3)
        theta=np.asarray(theta)
        psi=np.asarray(psi)
        phi=np.asarray(phi)
        eps0=[]
        for i,j,k in zip(phi,theta,psi):
            eps0.append(self.material_properties.eps_t(i,j,k))
        self.epsilon=np.asarray(eps0).reshape(self.cubes_num[2],self.cubes_num[1],self.cubes_num[0],3,3)
        p_direction=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]).T
        self.p_direcations=np.reshape(p_direction,(self.cubes_num[2],self.cubes_num[1],self.cubes_num[0],3))
        
    def write_to_restart(self,restart_data_file):
        head=lxml.Element("restart_data")
        geometry=lxml.SubElement(head,"geometry")
        lxml.SubElement(geometry,"dimentions", V="{}".format(self.voltage),Dx="{}".format(self.cubes_num[0]), Dy="{}".format(self.cubes_num[1]), Dz="{}".format(self.cubes_num[2]))
        for z in range (0,self.cubes_num[2]):
            for y in range (0,self.cubes_num[1]):
                for x in range (0,self.cubes_num[0]):
                    lxml.SubElement(geometry, "grain",p_direcationsx="{}".format(self.p_direcations[z][y][x][0]),p_direcationsy="{}".format(self.p_direcations[z][y][x][1]),p_direcationsz="{}".format(self.p_direcations[z][y][x][2]),\
                    theta="{}".format(self.material_properties.theta[x+self.cubes_num[0]*y+self.cubes_num[0]*self.cubes_num[1]*z]), \
                    phi="{}".format(self.material_properties.phi[x+self.cubes_num[0]*y+self.cubes_num[0]*self.cubes_num[1]*z]), \
                    psi="{}".format(self.material_properties.psi[x+self.cubes_num[0]*y+self.cubes_num[0]*self.cubes_num[1]*z]), id="{}.{}.{}".format(x+1,y+1,z+1))

        material=lxml.SubElement(head,"material")
        lxml.SubElement(material,"saturation", Ps="{}".format(self.Ps))
        lxml.SubElement(material,"activation_energy", Ea="{}".format(self.Ea))
        lxml.SubElement(material,"charact_time", tau0="{}".format(self.tau0))
        lxml.SubElement(material,"eps_main_values", eps="{}".format(self.material_properties.eps))

        with open(restart_data_file, 'w') as f:
                f.write(lxml.tostring(head, pretty_print=True))

    def point_position(self,i,j):
        a0=np.array([[(i-1)+0.5],[(j-1)+0.5]])
        a1=np.array([[(i-1)+0.25],[(j-1)+0.25]])
        a2=np.array([[(i-1)+0.25],[(j-1)+0.75]])
        a3=np.array([[(i-1)+0.75],[(j-1)+0.75]])
        a4=np.array([[(i-1)+0.75],[(j-1)+0.25]])
        return np.array([a0,a1,a2,a3,a4])

    def write_pde_start(self):
        V=self.voltage
        epsilon=self.epsilon
        cubes_num=self.cubes_num
        cube_directory=self.cube_directory
        with open ('{}/{}_s.pde'.format(cube_directory,self.material_properties.symmetry),'w') as f:
            f.write('TITLE "{}_{} from {}" \n'.format(self.material_properties.symmetry,cubes_num[0],datetime.datetime.now()))
            f.write("SELECT\nerrlim=0.01\n")
            f.write("COORDINATES \n")
            f.write("Cartesian3 \n")
            f.write("VARIABLES \n")
            f.write("uv (threshold 10)\n")
            f.write("DEFINITIONS\n")
            f.write("    voltage= {}\n    l=1e-6\n     v=l^3\n".format(-V))
            f.write("    E=-grad(uv)\n   Ex=0    Ey=0    Ez=0\n   e0=8.85418781762e-12\n\
    exx=1    exy=0    exz=0    eyx=0    eyy=1    eyz=0    ezx=0    ezy=0    ezz=1\n".format(self.Ps,self.cubes_num[2]))
            f.write("EQUATIONS \n")
            f.write("uv: e0*(dx(exx*dx(uv)+exy*dy(uv)+exz*dz(uv))+dy(eyx*dx(uv)+eyy*dy(uv)+eyz*dz(uv))+dz(ezx*dx(uv)+ezy*dy(uv)+ezz*dz(uv)))=0 \n")
            f.write("EXTRUSION\n")
            for i in range (0,cubes_num[2]):
                f.write("   SURFACE '{}' z={}*l\n".format(i+1,i))
                f.write("       LAYER '{}-{}'\n".format(i+1,i+2))
            f.write("SURFACE '{}' z={}*l\n".format(cubes_num[2]+1,cubes_num[2]))
            f.write("BOUNDARIES \n")
            f.write("   SURFACE '1' value(uv)=0\n")
            f.write("   SURFACE '{}' value(uv)=voltage\n".format(cubes_num[2]+1))
            f.write("       Region '0.0.0'\n")
            f.write("       start(0,0) \n")
            f.write("       line to (0.2*l,0)\n")
            f.write("       line to ({}*l,0)\n".format(cubes_num[0]-0.2))
            f.write("       line to ({}*l,0) periodic (x-{}*l,y)\n".format(cubes_num[0],cubes_num[0]))
            f.write("       line to ({}*l,{}*l)\n".format(cubes_num[0],cubes_num[1]))
            f.write("       line to ({}*l,{}*l) periodic (x,y-{}*l)\n".format(cubes_num[0]-0.2,cubes_num[1],cubes_num[1]))
            f.write("       line to (0.2*l,{}*l)\n".format(cubes_num[1]))
            f.write("       line to (0,{}*l)\n".format(cubes_num[1]))
            f.write("       line to close\n")
            for x in range (1,cubes_num[0]+1):
                for y in range (1,cubes_num[1]+1):
                    f.write("REGION '{}.{}' \n".format(x,y))
                    for z in range (1,cubes_num[2]+1):
                        f.write("       LAYER '{}-{}'\n".format(z,z+1))
                        f.write("           exx={}\n           exy={}\n           exz={}\n".format(epsilon[z-1][y-1][x-1][0][0],epsilon[z-1][y-1][x-1][0][1],epsilon[z-1][y-1][x-1][0][2]))
                        f.write("           eyx={}\n           eyy={}\n           eyz={}\n".format(epsilon[z-1][y-1][x-1][1][0],epsilon[z-1][y-1][x-1][1][1],epsilon[z-1][y-1][x-1][1][2]))
                        f.write("           ezx={}\n           ezy={}\n           ezz={}\n".format(epsilon[z-1][y-1][x-1][2][0],epsilon[z-1][y-1][x-1][2][1],epsilon[z-1][y-1][x-1][2][2]))
                        f.write("           Ex=VOL_INTEGRAL(xcomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           Ey=VOL_INTEGRAL(ycomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           Ez=VOL_INTEGRAL(zcomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        if z!=cubes_num[2]:
                            f.write("   SURFACE '{}' natural(uv)=0\n".format(z+1))
                    f.write("   START ({}*l,{}*l) natural(uv)=0\n\
                    LINE TO ({}*l,{}*l) natural(uv)=0\n\
                    LINE TO ({}*l,{}*l) natural(uv)=0\n\
                    LINE TO ({}*l,{}*l) natural(uv)=0\n\
                    LINE TO CLOSE\n".format(x-1,y-1,x,y-1,x,y,x-1,y))
            f.write("PLOTS\n")
            f.write("transfer(uv) file='uv.dat'\n")
            """
            f.write("vector(E) on y={}.5*l \n".format(int(self.cubes_num[1]/2)))
            f.write("table(magnitude(E)) format '#1' points=100 file='Em.txt'\n")
            f.write("SUMMARY EXPORT FILE 'E_ap.out'\n")
            f.write("REPORT('avenergy')\n")
            for z in range (1,cubes_num[2]+1):
                for y in range (1,cubes_num[1]+1):
                    for x in range (1,cubes_num[0]+1):
                        f.write('REPORT(VAL(Ex,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
                        f.write('REPORT(VAL(Ey,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
                        f.write('REPORT(VAL(Ez,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
                        f.write('REPORT(VAL(ezx,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
                        f.write('REPORT(VAL(ezy,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
                        f.write('REPORT(VAL(ezz,{}*l,{}*l,{}*l))\n'.format(x-0.5,y-0.5,z-0.5))
            """
            f.write("END")
            f.close()
            print ("FlexPDE file is ready")

    def write_pde_restart(self):
        P=self.p_direcations
        epsilon=self.epsilon
        cubes_num=self.cubes_num
        cube_directory=self.cube_directory
        with open ('{}/{}.pde'.format(cube_directory,self.material_properties.symmetry),'w') as f:
            f.write('TITLE "{}_{} from {}" \n'.format(self.material_properties.symmetry,cubes_num[0],datetime.datetime.now()))
            f.write("SELECT\nerrlim=0.01\n")
            f.write("COORDINATES \n")
            f.write("Cartesian3 \n")
            f.write("VARIABLES \n")
            f.write("up (threshold 10)\n uv (threshold 10)\n ps (threshold 0.01)\n")
            f.write("DEFINITIONS\n")
            f.write("    l=1e-6\n")
            f.write("    ps_max={}\n   tau0={}\n    v=l^3\n    P=vector(0,0,0)\n\
    ps_av=-0.01\n   px=xcomp(P)    py=ycomp(P)    pz=zcomp(P)\n\
    Ea=33e6\n    E_cr=0.4e6\n  E_v=-grad(uv)\n    E_p=-grad(up)\n\
    k=(E_cr+(magnitude(E_p)-E_cr)*ustep(sign(magnitude(E_p)-E_cr)))/E_cr\n\
    E=E_v+E_p/k\n\
    Ex=xcomp(E)    Ey=ycomp(E)    Ez=zcomp(E)\n\
    e0=8.85418781762e-12\n\
    exx=1    exy=0    exz=0    eyx=0    eyy=1    eyz=0    ezx=0    ezy=0    ezz=1\n".format(self.Ps,self.tau0))
            f.write("    ps_c= if (abs(ps_av)>0.99) then 0.99*ps_av/abs(ps_av) else ps_av\n")
            f.write("    E_pn=(Ex*px+Ey*py+Ez*pz)/(ps_c*ps_max)\n")
            f.write("    tau= if (t>1e-7)\n    \
                then\n    \
                    if (Ea/abs(E_pn)>100)\n    \
                        then (tau0*exp(100))\n    \
                          else (tau0*exp(Ea/abs(E_pn)))\n    \
                else (tau0*exp(100))\n")
            f.write("n=3\n")
            f.write("transfermesh( 'to_up.dat', uv_in,up_in,ps_in)\n")
            #f.write("transfermesh( 'uv.dat', uv_in)\n")
            f.write("INITIAL VALUES\n ps=ps_in\n uv=-uv_in\n up=up_in\n")
            f.write("EQUATIONS \n")
            f.write("up: e0*(dx(exx*dx(up)+exy*dy(up)+exz*dz(up))+dy(eyx*dx(up)+eyy*dy(up)+eyz*dz(up))+dz(ezx*dx(up)+ezy*dy(up)+ezz*dz(up)))=0 \n")
            f.write("uv: dt(uv)=0 \n")
            f.write("ps: dt(ps)=(sign(E_pn)-ps)*n*(t/tau)^(n-1)/tau\n")
            f.write("EXTRUSION\n")
            for i in range (0,cubes_num[2]):
                f.write("   SURFACE '{}' z={}*l\n".format(i+1,i))
                f.write("       LAYER '{}-{}'\n".format(i+1,i+2))
            f.write("SURFACE '{}' z={}*l\n".format(cubes_num[2]+1,cubes_num[2]))
            f.write("BOUNDARIES \n")
            f.write("   SURFACE '1' value(up)=0\n")
            f.write("   SURFACE '{}' value(up)=0\n".format(cubes_num[2]+1))
            f.write("       Region '0.0.0'\n")
            f.write("       start(0,0) \n")
            f.write("       line to (0.2*l,0)\n")
            f.write("       line to ({}*l,0)\n".format(cubes_num[0]-0.2))
            f.write("       line to ({}*l,0) periodic (x-{}*l,y)\n".format(cubes_num[0],cubes_num[0]))
            f.write("       line to ({}*l,{}*l)\n".format(cubes_num[0],cubes_num[1]))
            f.write("       line to ({}*l,{}*l) periodic (x,y-{}*l)\n".format(cubes_num[0]-0.2,cubes_num[1],cubes_num[1]))
            f.write("       line to (0.2*l,{}*l)\n".format(cubes_num[1]))
            f.write("       line to (0,{}*l)\n".format(cubes_num[1]))
            f.write("       line to close\n")
            for x in range (1,cubes_num[0]+1):
                for y in range (1,cubes_num[1]+1):
                    f.write("REGION '{}.{}' \n".format(x,y))
                    for z in range (1,cubes_num[2]+1):
                        f.write("       LAYER '{}-{}'\n".format(z,z+1))
                        f.write("           exx={}\n           exy={}\n           exz={}\n".format(epsilon[z-1][y-1][x-1][0][0],epsilon[z-1][y-1][x-1][0][1],epsilon[z-1][y-1][x-1][0][2]))
                        f.write("           eyx={}\n           eyy={}\n           eyz={}\n".format(epsilon[z-1][y-1][x-1][1][0],epsilon[z-1][y-1][x-1][1][1],epsilon[z-1][y-1][x-1][1][2]))
                        f.write("           ezx={}\n           ezy={}\n           ezz={}\n".format(epsilon[z-1][y-1][x-1][2][0],epsilon[z-1][y-1][x-1][2][1],epsilon[z-1][y-1][x-1][2][2]))
                        f.write("           Ex=VOL_INTEGRAL(xcomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           Ey=VOL_INTEGRAL(ycomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           Ez=VOL_INTEGRAL(zcomp(E),'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           ps_av=VOL_INTEGRAL(ps,'{}.{}','{}-{}')/v\n".format(x,y,z,z+1))
                        f.write("           P=vector({0},{1},{2})*ps_c*ps_max\n"\
                        .format(P[z-1][y-1][x-1][0],P[z-1][y-1][x-1][1],P[z-1][y-1][x-1][2]))
                        if z!=cubes_num[2]:
                            f.write("   SURFACE '{}' natural(up)=normal(+P)\n".format(z+1))
                    f.write("   START ({}*l,{}*l) natural(up)=normal(+P)\n\
                    LINE TO ({}*l,{}*l) natural(up)=normal(+P)\n\
                    LINE TO ({}*l,{}*l) natural(up)=normal(+P)\n\
                    LINE TO ({}*l,{}*l) natural(up)=normal(+P)\n\
                    LINE TO CLOSE\n".format(x-1,y-1,x,y-1,x,y,x-1,y))
            #in case of sigma

            for y in range (1,cubes_num[1]+1):
                for x in range (1,cubes_num[0]+1):
                    f.write("FEATURE 'y_{}.{}.{}' \n".format(x-1,x,y-1))
                    f.write("start({}*l,{}*l) line to ({}*l,{}*l) \n".format(x-1,y-1,x,y-1))
                    f.write("FEATURE 'x_{}.{}.{}' \n".format(x-1,y-1,y))
                    f.write("start({}*l,{}*l) line to ({}*l,{}*l) \n".format(x-1,y-1,x-1,y))

            f.write("TIME 0 TO 1000\n")
            f.write("PLOTS\n")
            if self.material_properties.symmetry=="Tetra":
                f.write("for t=0 by 1e-7 to 1e-6 by 0.5e-6 to 1e-5 by 1e-5 to 1e-4 by 1e-4 to 1e-3 by 0.5e-3 to 1e-2 by 1e-2 to 1e-1 by 1e-1 to 1 by 1 to 10 by 10 to 100 by 100 to 1000\n")
            elif self.material_properties.symmetry=="Rhombo":
                f.write("for t=0 by 1e-7 to 1e-6 by 1e-6 to 1e-5 by 0.5e-5 to 1e-4 by 0.1e-4 to 1e-3 by 0.1e-3 to 1e-2 by 0.5e-2 to 1e-1 by 1e-1 to 1 by 1 to 10 by 10 to 100\n")
            else :
                f.write("for t=0 by 1e-7 to 1e-6 by 1e-6 to 1e-5 by 0.5e-5 to 1e-4 by 0.1e-4 to 1e-3 by 0.1e-3 to 1e-2 by 0.5e-2 to 1e-1 by 1e-1 to 1 by 1 to 10 by 10 to 100\n")
            f.write("HISTORY(VOL_INTEGRAL(zcomp(P))/({}*ps_max*v)) export format '#t#r#i'\n".format(np.prod(self.cubes_num)))
            f.write("vector(E) on y={}.5*l \n".format(int(self.cubes_num[2]/2)))
            f.write("vector(P) on y={}.5*l \n".format(int(self.cubes_num[2]/2)))
            if self.material_properties.symmetry=="Tetra":
                f.write("for t=0, 1e-7, 1e-4, 4e-4, 1e-3 by 1e-3 to 1e-2 by 5e-3 to 3e-2 by 1e-2 to 1e-1 by 2e-1 to 1, 5, 100, 1000\n")
            elif self.material_properties.symmetry=="Rhombo":
                f.write("for t=0, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4 by 1e-4 to 1e-3 by 1e-3 to 2e-3 , 2.4e-3, 3e-3 by 1e-3 to 1e-2 by 5e-2 to 1e-1 by 5e-1 to 1, 5, 10, 100\n")
            else:
                f.write("for t=0 by 1e-7 to 1e-6 by 1e-6 to 1e-5 by 0.5e-5 to 1e-4 by 0.1e-4 to 1e-3 by 0.1e-3 to 1e-2 by 0.5e-2 to 1e-1 by 1e-1 to 1 by 1 to 10 by 10 to 100\n")
            f.write("table(magnitude(E)) format '#1' points=100 file='Em.txt'\n")
            #in case of correlation
            """
            f.write("for t = 0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, time_end\n")
            f.write("SUMMARY EXPORT FILE 'E_cor.out'\n")
            f.write("REPORT('E_p_xz')\n")
            for z in range (1,cubes_num[2]+1):
                for y in range (1,cubes_num[1]+1):
                    for x in range (1,cubes_num[0]+1):
                        for k in [1,3]:
                            for j in [1,3]:
                                for i in [1,3]:
                                    f.write('REPORT(VAL(xcomp(E),{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    f.write('REPORT(VAL(ycomp(E),{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    f.write('REPORT(VAL(zcomp(E),{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    f.write('REPORT(VAL(ezx,{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    f.write('REPORT(VAL(ezy,{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    f.write('REPORT(VAL(ezz,{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    #f.write('REPORT(VAL(px,{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    #f.write('REPORT(VAL(py,{}*l,{}*l,{}*l))\n'.formatformat(x-1+i/4.,y-1+j/4.,z-1+k/4.))
                                    #f.write('REPORT(VAL(pz,{}*l,{}*l,{}*l))\n'.format(x-1+i/4.,y-1+j/4.,z-1+k/4.))
            """
            #in case of sigma

            f.write("SUMMARY EXPORT FILE 'sigma_xyz.out'\n")
            f.write("REPORT('xyz')\n")
            for z in range (2,cubes_num[2]+1):
                for y in range (1,cubes_num[1]+1):
                    for x in range (1,cubes_num[0]+1):
                        f.write("REPORT(SURF_INTEGRAL ( normal(+P), '{}', '{}.{}', '{}-{}')+".format(z,x,y,z-1,z))                   #surfaces
                        f.write("SURF_INTEGRAL ( normal(+P), '{}',    '{}.{}', '{}-{}'))\n".format(z,x,y,z,z+1))
            for z in range (1,cubes_num[2]):
                for x in range(cubes_num[1]):
                    f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}')+".format(x,x+1,0,z,z+1,x+1,0+1))
                    f.write("SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(x,x+1,0,z,z+1,x+1,cubes_num[1]))
                for y in range(1,cubes_num[1]):
                    for x in range(1,cubes_num[0]):
                        f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}')+".format(x-1,x,y,z,z+1,x,y))
                        f.write("SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(x-1,x,y,z,z+1,x,y+1))
                for y in range(1,cubes_num[1]):
                    f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}')+".format(cubes_num[0]-1,cubes_num[0],y,z,z+1,cubes_num[0],y))
                    f.write("SURF_INTEGRAL ( normal(+P), 'y_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(cubes_num[0]-1,cubes_num[0],y,z,z+1,cubes_num[0],y+1))
                for y in range(cubes_num[1]):
                    f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}')+".format(0,y,y+1,z,z+1,0+1,y+1))
                    f.write("SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(0,y,y+1,z,z+1,cubes_num[0],y+1))
                for x in range(1,cubes_num[1]):
                    for y in range(1,cubes_num[0]):
                        f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}')+".format(x,y-1,y,z,z+1,x,y))
                        f.write("SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(x,y-1,y,z,z+1,x+1,y))
                for x in range(1,cubes_num[0]):
                    f.write("REPORT(SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}')+".format(x,cubes_num[1]-1,cubes_num[1],z,z+1,x,cubes_num[1]))
                    f.write("SURF_INTEGRAL ( normal(+P), 'x_{}.{}.{}', '{}-{}', '{}.{}'))\n".format(x,cubes_num[1]-1,cubes_num[1],z,z+1,x+1,cubes_num[1]))
            for z in range (1,cubes_num[2]+1):
                for y in range (1,cubes_num[1]+1):
                    for x in range (1,cubes_num[0]+1):
                        f.write("REPORT(VOL_INTEGRAL(ps/v,'{}.{}','{}-{}'))\n".format(x,y,z,z+1))

            #f.write("transfer(uv,up,ps) file='to_up.dat'\n")
            f.write("END")
            f.close()
            print ("FlexPDE file is ready")

    def flex_run(self):
        if self.restart==False:
            self.write_pde_start()
            subprocess.check_output([os.path.join("/nfshome/khachaturyan/fpde640linux86_64","flexpde6n"),"-Q","{}/{}_s.pde".format(self.cube_directory,self.material_properties.symmetry)])
        else:
            self.write_pde_restart()
            subprocess.check_output([os.path.join("/nfshome/khachaturyan/fpde640linux86_64","flexpde6n"),"-Q","{}/{}.pde".format(self.cube_directory,self.material_properties.symmetry)])