import numpy as np
import rebound as rb
import time
import tracemalloc
# the gap is between 0.8 au and 5.7 au
Mearth = 5.9722e24 # kg
Msun = 1.989e30 # kg
Mjup = 1.898e27
Mconv2 = Mjup/Msun
Mconv = Mearth/Msun
AUd = 1.496e8 #km
Rsun = 696340 #km
Rearth = 6371 #km
Rconv = Rearth/AUd
Rconv2 = Rsun/AUd

########
#2 different sets of initial parameters
Mpaper2 = np.array([7.99,255.4,47.8,51.2,991.6])*Mconv
apaper2 = np.array([0.0154,0.1134,0.7708,0.2373,5.957])
epaper2 = np.array([0.05,0,0.08,0.03,0.13])
Rpaper2 = np.array([2.5452621503266757,13.254298756623639,\
                    8.086855328828321,8.249880523982966,13.24438689860083])*Rconv

Mpaper1 = np.array([8.37,264.75,57.209,54.38,1169.6])*Mconv
apaper1 = np.array([0.0356,0.1148,0.781,0.2403,5.74])
epaper1 = np.array([0.170,0.010,0.32,0.005,0.020])

#######
#create the different systems
sim0 = rb.Simulation()
sim0.units = ('AU','yr','Msun')
sim0.move_to_com()
sim0.integrator = "mercurius"

sim1 = rb.Simulation()
sim1.units = ('AU','yr','Msun')
sim1.move_to_com()
sim1.integrator = "mercurius"

sim2 = rb.Simulation()
sim2.units = ('AU','yr','Msun')
sim2.move_to_com()
sim2.integrator = "mercurius"

sim3 = rb.Simulation()
sim3.units = ('AU','yr','Msun')
sim3.move_to_com()
sim3.integrator = "mercurius"

sim4 = rb.Simulation()
sim4.units = ('AU','yr','Msun')
sim4.move_to_com()
sim4.integrator = "mercurius"

sim5 = rb.Simulation()
sim5.units = ('AU','yr','Msun')
sim5.move_to_com()
sim5.integrator = "mercurius"

########
#add objects to the different systems
sim0.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim0.add(m=Mpaper1[0],e=epaper1[0],a=apaper1[0],r=Rpaper2[0],inc=0,hash="planet-b") # 
sim0.add(m=Mpaper1[1],e=epaper1[1],a=apaper1[1],r=Rpaper2[1],inc=0,hash="planet-b") # 
sim0.add(m=Mpaper1[2],e=epaper1[2],a=apaper1[2],r=Rpaper2[2],inc=0,hash="planet-f") # 
sim0.add(m=Mpaper1[3],e=epaper1[3],a=apaper1[3],r=Rpaper2[3],inc=0,hash="planet-c") # 
sim0.add(m=Mpaper1[4],e=epaper1[4],a=apaper1[4],r=Rpaper2[4],inc=0,hash="planet-d") #  

sim1.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim1.add(m=Mpaper1[1],e=epaper1[1],a=apaper1[1],r=Rpaper2[1],inc=0,hash="planet-b") # 
sim1.add(m=Mpaper1[2],e=epaper1[2],a=apaper1[2],r=Rpaper2[2],inc=0,hash="planet-f") # 
sim1.add(m=Mpaper1[3],e=epaper1[3],a=apaper1[3],r=Rpaper2[3],inc=0,hash="planet-c") # 
sim1.add(m=Mpaper1[4],e=epaper1[4],a=apaper1[4],r=Rpaper2[4],inc=0,hash="planet-d") # 

sim2.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim2.add(m=Mpaper1[2],e=epaper1[2],a=apaper1[2],r=Rpaper2[2],inc=0,hash="planet-f") # 
sim2.add(m=Mpaper1[3],e=epaper1[3],a=apaper1[3],r=Rpaper2[3],inc=0,hash="planet-c") # 
sim2.add(m=Mpaper1[4],e=epaper1[4],a=apaper1[4],r=Rpaper2[4],inc=0,hash="planet-d") #  

sim3.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim3.add(m=Mpaper2[0],e=epaper2[0],a=apaper2[0],r=Rpaper2[0],inc=0,hash="planet-b") # 
sim3.add(m=Mpaper2[1],e=epaper2[1],a=apaper2[1],r=Rpaper2[1],inc=0,hash="planet-b") # 
sim3.add(m=Mpaper2[2],e=epaper2[2],a=apaper2[2],inc=0,r=Rpaper2[2],hash="planet-f") # 
sim3.add(m=Mpaper2[3],e=epaper2[3],a=apaper2[3],inc=0,r=Rpaper2[3],hash="planet-c") # 
sim3.add(m=Mpaper2[4],e=epaper2[4],a=apaper2[4],inc=0,r=Rpaper2[4],hash="planet-d") #  


sim4.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim4.add(m=Mpaper2[1],e=epaper2[1],a=apaper2[1],r=Rpaper2[1],inc=0,hash="planet-b") # 
sim4.add(m=Mpaper2[2],e=epaper2[2],a=apaper2[2],r=Rpaper2[2],inc=0,hash="planet-f") # 
sim4.add(m=Mpaper2[3],e=epaper2[3],a=apaper2[3],r=Rpaper2[3],inc=0,hash="planet-c") # 
sim4.add(m=Mpaper2[4],e=epaper2[4],a=apaper2[4],r=Rpaper2[4],inc=0,hash="planet-d") #  

sim5.add(m=0.96+Mpaper1[0],r=0.943*Rconv2,hash="mainStar")
sim5.add(m=Mpaper2[2],e=epaper2[2],a=apaper2[2],r=Rpaper2[2],inc=0,hash="planet-f") # 
sim5.add(m=Mpaper2[3],e=epaper2[3],a=apaper2[3],r=Rpaper2[3],inc=0,hash="planet-c") # 
sim5.add(m=Mpaper2[4],e=epaper2[4],a=apaper2[4],r=Rpaper2[4],inc=0,hash="planet-d") #  


simList = [sim0,sim1,sim2,sim3,sim4,sim5]

for i,s in enumerate(simList):
    s.move_to_com()
    s.integrator = "mercurius"
    s.save("cancriStartTest3.bin")
    s.dt = 0.05*s.calculate_orbits()[1].P

########
#integrate all configurations by doing a forloop through all of them
startTime = time.time()

endTime = int(1e5)
Nplanets = [5,4,3,5,4,3]
Nsteps = int(1e3)
times = np.linspace(0,endTime,Nsteps)

for n,s in enumerate(simList):
    major = np.zeros((Nplanets[n],Nsteps))
    ecc = np.zeros((Nplanets[n],Nsteps))
    incl = np.zeros((Nplanets[n],Nsteps))
    E = np.zeros(Nsteps)
    for i,t in enumerate(times):
        s.integrate(t)
        os = s.calculate_orbits()
        for j in range(Nplanets[n]):
            major[j][i] = os[j].a 
            ecc[j][i] = os[j].e
            incl[j][i] = os[j].inc
        E[i] = s.calculate_energy()
    orbitals = (major,ecc,incl)
    with open("orbitals_data_cancriTestingR"+ str(n) + ".dat","w") as fd_out:
        for i in range(len(times)):
            for j in range(len(orbitals)):
                for k in range(Nplanets[n]-1):
                    fd_out.write(str(orbitals[j][k][i])+"\t")
                fd_out.write(str(orbitals[j][Nplanets[n]-1][i])+"\n")
    fd_out.close()

    with open("energy_data_cancriTestingR"+ str(n) +".dat","w") as fd_out:
        for i in range(len(E)):
            fd_out.write(str(E[i])+"\n")
    fd_out.close()
    s.save("cancriEndTest3.bin")
#check what the time length is of the simulation
endTime = time.time()
print("time:" , endTime-startTime)

