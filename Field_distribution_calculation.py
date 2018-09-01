import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
from materials import material_properties
mp=material_properties("start.xml")
import flex_code
fc=flex_code.flex_code(mp)

class Dist:

    def __init__(self,fc,mp):

        self.fc=fc
        self.mp=mp
        self.scale=1e-6
        print (self.mp.symmetry)
        self.width_x=self.mp.cube_dimentions[0]
        self.width_y=self.mp.cube_dimentions[1]
        self.width_z=self.mp.cube_dimentions[2]

    def distribution(self,x,num,p):
        if x=='P' or x=='Py' or x=='Px':
            resolution=1000
            applied=1
            with open ('./distributions/field_distributions/data/{}_{}.txt'.format(x,num), 'r') as f:
                a=f.readlines()
                f.close()
            with open ('./distributions/field_distributions/data/{}_{}.txt'.format(x,num), 'r') as f:
                i=0
                while i<len(a):
                    if '0.000' not in a[i]:
                        f.write(a[i])
                        i+=1
        else:
            resolution=500
            applied=mp.v/self.scale
        field_values=[]
        with open ('./distributions/field_distributions/data/{}_{}.txt'.format(x,num), 'r') as f:
            s=0
            for line in f:
                s+=1
                a=line.strip()
                if s>7:
                    field_values.append(float(a)/applied)
        field_values=np.array(field_values)
        step=float((field_values.max()-field_values.min())/resolution)
        number_in_step=[]
        for i in range(resolution+1):
            number_in_step.append(len(np.where(field_values<field_values.min()+(i+1)*step)[0])-len(np.where(field_values<field_values.min()+i*step)[0]))
        normalization=sum(number_in_step)*step
        x_axis=np.array([field_values.min()+step/2+i*step for i in range(resolution+1)])
        y_axis=np.array(number_in_step)/normalization
        with open('./distributions/field_distributions/calculated/{}_{}_.txt'.format(x,p), "w") as f:
            for i,j in zip(x_axis,y_axis):
                f.write("{} {}\n".format(i,j))
            f.close()
        return x_axis,y_axis

dist_calc=Dist(fc,mp)

name="Em"

step=[]
time_steps=[]
path_to_data=os.getcwd()+"/distributions/field_distributions/data"
for file in os.listdir('{}'.format(path_to_data)):
    if file.endswith(".txt") and name in file:
        with open (path_to_data+"/"+file, "r") as f:
            a=f.readlines()
            time_steps.append(float(a[4].split(" ")[-1]))
        if int(file.split("_")[1].split(".")[0]) not in step:
            step.append(int(file.split("_")[1].split(".")[0]))

step=sorted(np.array(step))
time_steps=sorted(np.array(time_steps))
p_values=[]
t,p=np.loadtxt('{}/cube_up.p01'.format(path_to_data), unpack=True)

for s in time_steps:
    p_values.append(round(p[np.where(t==s)][0],3))

for i,j in zip(step,p_values):
    print (i, j)
    dist_calc.distribution(name,i,j)