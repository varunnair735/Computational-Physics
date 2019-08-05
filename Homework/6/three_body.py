import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

def f(s,t):
    r1 = s[0,:]
    r2 = s[1,:]
    r3 = s[2,:]
    d1 = s[3,:]
    d2 = s[4,:]
    d3 = s[5,:]


    fr1 = d1
    fr2 = d2
    fr3 = d3
    fd1 = m2*(r2-r1)/((LA.norm(r2-r1))**3)\
        + m3*(r3-r1)/((LA.norm(r3-r1))**3)

    fd2 = m1*(r1-r2)/((LA.norm(r1-r2))**3)\
        + m3*(r3-r2)/((LA.norm(r3-r2))**3)

    fd3 = m1*(r1-r3)/((LA.norm(r1-r3))**3)\
        + m2*(r2-r3)/((LA.norm(r2-r3))**3)

    return np.array([fr1,fr2,fr3,fd1,fd2,fd3],float)

t0 = 0.0
t1 = 2.0
dt = 1e-6 #target accuracy per unit time
N = 1e4
h = (t1-t0)/N #starting step size

tpoints = np.arange(t0,t1,h)
x1points = []
x2points = []
x3points = []
y1points = []
y2points = []
y3points = []

#initial conditions
m1 = 150
m2 = 200
m3 = 250
s = np.array([[3,1],
              [-1,-2],
              [-1,1],
              [0,0],
              [0,0],
              [0,0]],float)

x1points.append(s[0,0])
x2points.append(s[1,0])
x3points.append(s[2,0])
y1points.append(s[0,1])
y2points.append(s[1,1])
y3points.append(s[2,1])

for t in tpoints:
    rho = 0.5 #dummy value

    #steps need to be redone for rho too small
    while rho < 1:
        #estimates x(t+2h) in two steps (x_1)
        s1 = np.copy(s)

        k1 = h * f(s1, t)
        k2 = h * f(s1+0.5*k1, t+0.5*h)
        k3 = h * f(s1+0.5*k2, t+0.5*h)
        k4 = h * f(s1+k3, t+h)
        s1 += (k1 + 2*k2 + 2*k3 + k4) / 6

        k1 = h * f(s1, t)
        k2 = h * f(s1+0.5*k1, t+0.5*h)
        k3 = h * f(s1+0.5*k2, t+0.5*h)
        k4 = h * f(s1+k3, t+h)
        s1 += (k1 + 2*k2 + 2*k3 + k4) / 6

        #estimates x(t+2h) in one step (x_2)
        s2 = np.copy(s)

        k1 = 2*h * f(s, t)
        k2 = 2*h * f(s+0.5*k1, t+h)
        k3 = 2*h * f(s+0.5*k2, t+h)
        k4 = 2*h * f(s+k3, t+2*h)
        s2 += (k1 + 2*k2 + 2*k3 + k4) / 6

        #calculates joint error in two dimensions
        epsilon_x = (s1[:,0]-s2[:,0]) / 30
        epsilon_y = (s1[:,1]-s2[:,1]) / 30
        error = np.sqrt(np.square(epsilon_x)\
              + np.square(epsilon_y))

        #bounds step size h
        for i in range(len(error)):
            if error[i] < 1e-15:
                error[i] = h*dt/2**4

        #definition of multipart error from book
        rho = h*dt/error

        rho = np.min(rho)
        if rho > 2**4:
            rho = 2**4
        #adjusts step size
        h *= rho**(1/4)

    s = np.copy(s1)
    x1points.append(s[0,0])
    x2points.append(s[1,0])
    x3points.append(s[2,0])
    y1points.append(s[0,1])
    y2points.append(s[1,1])
    y3points.append(s[2,1])




fig6, ax6 = plt.subplots(1, 1, figsize = (12, 8))

ax6.plot(x1points,y1points,'m', label='Star 1')
ax6.plot(x2points,y2points,'k', label='Star 2')
ax6.plot(x3points,y3points,'c', label='Star 3')
ax6.set_xlabel("$x$")
ax6.set_ylabel("$y$")
ax6.set_title(\
        "Trajectories of Stars in the Three Body Problem")
ax6.legend()
plt.show()
