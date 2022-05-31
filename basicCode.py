import numpy as np
import rebound as rb
import time
import random

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
#initial parameters for the planets

Mpaper1 = np.array([7.99,255.4,47.8,51.2,991.6])*Mconv
Rpaper2 = np.array([2.5452621503266757,13.254298756623639,\
                    8.086855328828321,8.249880523982966,13.24438689860083])*Rconv
apaper1 = np.array([0.0154,0.1134,0.7708,0.2373,5.957])
epaper1 = np.array([0.05,0,0.08,0.03,0.13])

########
#function for determining the distance for an object
def norm( vector ):
    return vector.x**2 + vector.y**2 + vector.z**2

#find the closes planet to a test particle
def findClosestPlanet( planets, asteroid ):
    """
    Return hash of planet (within array planets) that's closest to asteroid.
    Aux fxn used within removeCollidedAsteroids.
    """
    distances = np.array( [ norm(p-asteroid) for p in planets ] )
    idx=np.where( distances == np.min(distances) )
    idx=idx[0]
    assert len(idx) == 1 ## There can be only one closest planet to an impacting asteroid!
    return planets[idx[0]].hash 

#remove a test particle after it collided with a planet
def removeCollidedAsteroids( sim ):
    """ 
    To be called after a rebound simulation encountered a rebound.Collision exception.
    Loops over test particles with .lastcollision > 0 (i.e., they have collided with something), 
    determine the impact target, remove them from sim.
    Return value: (impacts, asteroidCopies) where:
      impacts is a dictionary: keys are asteroid hashes, values planet hash values.
      asteroidCopies are copies of the test particles before removal.  (are those needed at all?)
    Failure mode: if planets are very close-by and time steps are very coarse, impactors may get assigned 
       to the wrong planet (assignment based on position at end of time step, not at impact).
    """
    planets = sim.particles[:sim.N_active] # also include star(s)
    testParticles = sim.particles[sim.N_active:]
    # which asteroids impacted?
    asteroidTimes = np.array( [a.lastcollision for a in testParticles] )
    idx = np.where( asteroidTimes > 0 )
    idx=idx[0]
    hashList = []
    if len( idx ) == 0: ## if no asteroid is involved in the collision
        return hashList
    for i in reversed(idx):    
        ## reversed: loop backwards, particle removal could mess up indexing otherwise
        a = testParticles[i]
        pHash = findClosestPlanet(planets, a) 
        sim.remove( hash=a.hash )
        p = sim.particles[pHash]
        hashList.append((sim.t,a.hash.value,pHash.value,a.x,a.y,a.z,a.vx,a.vy,a.vz,p.x,p.y,p.z,p.vx,p.vy,p.vz))
    return hashList

# remove a particle that escaped the system
def removeEscapedAsteroids(sim,boxsize):
    planets = sim.particles[:sim.N_active] # also include star(s)
    testParticles = sim.particles[sim.N_active:]
    # which asteroids escaped?
    norms = np.array( [ norm(p) for p in testParticles] )
    idx = np.where( norms > boxsize**2 )
    idx=idx[0]
    #assert len( idx ) > 0 ### This won't work if no asteroid is involved in a collision (don't support planet-planet collisions)
    hashList = []
    for i in reversed(idx):
        ## reversed: loop backwards, particle removal could mess up indexing otherwise
        a = testParticles[i]
        hashList.append((sim.t,a.hash.value,a.x,a.y,a.z,a.vx,a.vy,a.vz))
        sim.remove( hash=a.hash )
    return hashList
    
    
    #remove a test particle if it collided with a planet but do it immediately instead
    # of after the timestep
def my_merge(sim_pointer, collided_particles_index):
    removeValue = 2
    sim = sim_pointer.contents # retreive the standard simulation object
    ps = sim.particles # easy access to list of particles

    i = collided_particles_index.p1   # Note that p1 < p2 is not guaranteed.    
    j = collided_particles_index.p2 
    if int(j)<int(i):
        i = collided_particles_index.p2
        j = collided_particles_index.p1
        removeValue = 1
    os = sim.calculate_orbits()
    pHash = sim.particles[int(i)].hash
    a = sim.particles[int(j)]
    p = sim.particles[pHash]
    
    print(i,j,pHash,a.hash)
    # This part is exciting! We can execute additional code during collisions now!
    collidedList.append((sim.t,a.hash.value,pHash.value,a.x,a.y,a.z,a.vx,a.vy,a.vz,p.x,p.y,p.z,p.vx,p.vy,p.vz,\
                        os[int(j)-1].a,os[int(j)-1].e,os[int(j)-1].inc))
    return removeValue

###########
#setup the simulation
sim = rb.Simulation()
sim.units = ('AU','yr','Msun')
########
#add objects to the simulation
sim.add(m=0.96+Mpaper1[0]+Mpaper1[1],r=0.943*Rconv2,hash=0)

sim.add(m=Mpaper1[2],e=epaper1[2],a=apaper1[2],r=Rpaper2[2],inc=0,hash=1) # 
sim.add(m=Mpaper1[3],e=epaper1[3],a=apaper1[3],r=Rpaper2[3],inc=0,hash=2) # 
sim.add(m=Mpaper1[4],e=epaper1[4],a=apaper1[4],r=Rpaper2[4],inc=0,hash=3) # 

Nplanets=3
N_testparticle = 100
file = "asteroidValues.dat"
data = np.genfromtxt(file)

iteration=0

a = data[N_testparticle*(iteration-1):N_testparticle*iteration:,0]
e = data[N_testparticle*(iteration-1):N_testparticle*iteration:,1]
inc = data[N_testparticle*(iteration-1):N_testparticle*iteration:,2]

numbers = range(4,N_testparticle+4)
for i in range(N_testparticle):
    rand = np.random.random()*2*np.pi
    sim.add(a=a[i], e=e[i], inc=inc[i], Omega=0, omega=rand, f=rand,hash=numbers[i])
#for q in range(sim.N):
#    print(sim.particles[q].hash.value)
sim.N_active = 4

collidedList = []
escapedList = []

sim.move_to_com()
sim.integrator = "mercurius"
boxsize= 931
sim.exit_max_distance = boxsize
sim.testparticle_type = 0
sim.collision = "direct"
sim.collision_resolve = my_merge

year=1
sim.dt = year/365.25*2
############

#do a foreloop from time 0 to the end of the simulation. at certain timesteps make snapshots
#these snapshots show what the simulation looks like at that moment so we can see how the simulation
# evolves over time. 
times = np.array([0,year*100,year*1000,year*10000,year*100000,year*300000,year*500000,year*700000,year*900000,year*1000000])
#times = np.array([0,year*10000])
Nsteps = len(times)
major = np.zeros((Nplanets,Nsteps))
ecc = np.zeros((Nplanets,Nsteps))
E = np.zeros(Nsteps)
for i,t in enumerate(times):
    done = False
    while not done:
        try:
            sim.integrate(t,exact_finish_time=0)
            done = True
        except rb.Escape as error:
            hashList2 = removeEscapedAsteroids(sim,boxsize)
            for q in hashList2:
                escapedList.append(q)
            testParticles = sim.particles[sim.N_active:]
            if len(testParticles) == 0:
                done=True # stop sim once we're out of test particles
    fileName = str(iteration) + "snapshot" + str(i)+".bin"
    sim.save(fileName)


###########C
#add all the values that were saved during the simulation to lists and then to files for analysation later

escapedValues = []
collidedValues = []
remainingValues = []
remainingHashes = []

for i in sim.particles[sim.N_active:]:
    remainingHashes.append(i.hash.value)

fileName = str(iteration) + "snapshot0.bin"
sim = rb.Simulation(fileName)
os = sim.calculate_orbits()

for i in remainingHashes:
    remainingValues.append((os[i-1].a,os[i-1].e,os[i-1].inc))

for i in escapedList:
    escapedValues.append((os[i[1]-1].a,os[i[1]-1].e,os[i[1]-1].inc,i[2],i[3],i[4],i[5],i[6],i[7],i[0]))
        

collidedHashes = []
collidedValuesF = []
for i in collidedList:
    collidedHashes.append((i[1],i[2]))
    if (i[1]>i[2]):
        collidedValues.append((os[i[1]-1].a,os[i[1]-1].e,os[i[1]-1].inc,i[15],i[16],i[17],i[3],i[4],i[5],i[6],\
                               i[7],i[8],i[9],i[10],i[11],i[12],i[13],i[14],i[0],i[2]))
    else:
        collidedValues.append((os[i[2]-1].a,os[i[2]-1].e,os[i[2]-1].inc,i[15],i[16],i[17],i[3],i[4],i[5],i[6],\
                               i[7],i[8],i[9],i[10],i[11],i[12],i[13],i[14],i[0],i[1]))
    if i[0]>year*1000:
        if (i[1]>i[2]):
            collidedValuesF.append((os[i[1]-1].a,os[i[1]-1].e,os[i[1]-1].inc,i[15],i[16],i[17],i[3],i[4],i[5],i[6],\
                                   i[7],i[8],i[9],i[10],i[11],i[12],i[13],i[14],i[0],i[2]))
        else:
            collidedValuesF.append((os[i[2]-1].a,os[i[2]-1].e,os[i[2]-1].inc,i[15],i[16],i[17],i[3],i[4],i[5],i[6],\
                                   i[7],i[8],i[9],i[10],i[11],i[12],i[13],i[14],i[0],i[1]))          

        
#### here we write everything to files

for countShot in range(len(times)):
    remainingBefore = []
    snapshotName = str(iteration) + "snapshot"+str(countShot) +".bin"
    sim = rb.Simulation(snapshotName)
    os = sim.calculate_orbits()
    for i in sim.particles:
        print(i.hash)
    for counting,i in enumerate(sim.particles[sim.N_active:]):
        hashing = i.hash.value
        #print(counting+3)
        remainingBefore.append((os[counting+3].a,os[counting+3].e,os[counting+3].inc,i.x,i.y,i.z,i.vx,i.vy,i.vz))
    fileName = "remaining_particles" + str(countShot) + ".dat"
    with open(fileName,"a") as fd_out:
        for i in range(len(remainingBefore)):
            fd_out.write(str(remainingBefore[i][0])+"\t"+str(remainingBefore[i][1])+"\t"\
                         +str(remainingBefore[i][2])+"\t"+str(remainingBefore[i][3])+"\t"\
                         +str(remainingBefore[i][4])+"\t"+str(remainingBefore[i][5])+"\t"\
                         +str(remainingBefore[i][6])+"\t"+str(remainingBefore[i][7])+"\t"\
                         +str(remainingBefore[i][8])+"\n")
    fd_out.close()
  
  
fileName = "escaped_particles.dat"
with open(fileName,"a") as fd_out:
    for i in range(len(escapedValues)):
        fd_out.write(str(escapedValues[i][0])+"\t"+str(escapedValues[i][1])+"\t"+\
                     str(escapedValues[i][2])+"\t"+str(escapedValues[i][3])+"\t"+\
                     str(escapedValues[i][4])+"\t"+str(escapedValues[i][5])+"\t"+\
                     str(escapedValues[i][6])+"\t"+str(escapedValues[i][7])+"\t"+\
                     str(escapedValues[i][8])+"\t"+str(escapedValues[i][9])+"\n")
fd_out.close()
fileName = "collided_particles.dat"
with open(fileName,"a") as fd_out:
    for i in range(len(collidedValues)):
        fd_out.write(str(collidedValues[i][0])+"\t"+str(collidedValues[i][1])+"\t"+\
                     str(collidedValues[i][2])+"\t"+str(collidedValues[i][3])+"\t"+\
                     str(collidedValues[i][4])+"\t"+str(collidedValues[i][5])+"\t"+\
                     str(collidedValues[i][6])+"\t"+str(collidedValues[i][7])+"\t"+\
                     str(collidedValues[i][8])+"\t"+str(collidedValues[i][9])+"\t"+\
                     str(collidedValues[i][10])+"\t"+str(collidedValues[i][11])+"\t"+\
                     str(collidedValues[i][12])+"\t"+str(collidedValues[i][13])+"\t"+\
                     str(collidedValues[i][14])+"\t"+str(collidedValues[i][15])+"\t"+\
                     str(collidedValues[i][16])+"\t"+str(collidedValues[i][17])+"\t"+\
                     str(collidedValues[i][18])+"\t"+str(collidedValues[i][19])+"\n")
fd_out.close()

fileName = "filtered_collided_particles.dat"
with open(fileName,"a") as fd_out:
    for i in range(len(collidedValuesF)):
        fd_out.write(str(collidedValuesF[i][0])+"\t"+str(collidedValuesF[i][1])+"\t"+\
                     str(collidedValuesF[i][2])+"\t"+str(collidedValuesF[i][3])+"\t"+\
                     str(collidedValuesF[i][4])+"\t"+str(collidedValuesF[i][5])+"\t"+\
                     str(collidedValuesF[i][6])+"\t"+str(collidedValuesF[i][7])+"\t"+\
                     str(collidedValuesF[i][8])+"\t"+str(collidedValuesF[i][9])+"\t"+\
                     str(collidedValuesF[i][10])+"\t"+str(collidedValuesF[i][11])+"\t"+\
                     str(collidedValuesF[i][12])+"\t"+str(collidedValuesF[i][13])+"\t"+\
                     str(collidedValuesF[i][14])+"\t"+str(collidedValuesF[i][15])+"\t"+\
                     str(collidedValuesF[i][16])+"\t"+str(collidedValuesF[i][17])+"\t"+\
                     str(collidedValuesF[i][18])+"\t"+str(collidedValuesF[i][19])+"\n")
fd_out.close()

fileName = "collided_hashes.dat"
with open(fileName,"a") as fd_out:
    for i in range(len(collidedHashes)):
        fd_out.write(str(collidedHashes[i][0])+"\t"+str(collidedHashes[i][1])+"\n")
fd_out.close()

with open("remaining_particles_After.dat","a") as fd_out:
    for i in range(len(remainingValues)):
        fd_out.write(str(remainingValues[i][0])+"\t"+str(remainingValues[i][1])+"\t"+str(remainingValues[i][2])+"\n")
fd_out.close()
