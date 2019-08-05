sigma = 10
r = 28
b = 8/3

def f(s,t):
    x = s[0]
    y = s[1]
    z = s[2]
    fx = sigma * (y - x)
    fy = r*x - y - x*z
    fz = x*y - b*z
    return np.array([fx,fy,fz],float)

t0 = 0.0
t1 = 50.0
N = 1e4
h = (t1-t0)/N

tpoints = np.arange(t0,t1,h)
xpoints = []
ypoints = []
zpoints = []


#initial conditions
s = np.array([0,1,0],float)

xpoints.append(s[0])
ypoints.append(s[1])
zpoints.append(s[2])


for t in tpoints:
    rho = 0.5 #dummy value
    #steps need to be redone for rho too small
    #print(rho)
    while rho < 16:
        #does first two step estimation (x1)
        k1 = h * f(s, t)
        k2 = h * f(s+0.5*k1, t+0.5*h)
        k3 = h * f(s+0.5*k2, t+0.5*h)
        k4 = h * f(s+k3, t+h)
        s1 = (k1 + 2*k2 + 2*k3 + k4) / 6

        k1 = h * f(s1, t)
        k2 = h * f(s1+0.5*k1, t+0.5*h)
        k3 = h * f(s1+0.5*k2, t+0.5*h)
        k4 = h * f(s1+k3, t+h)
        s1 += (k1 + 2*k2 + 2*k3 + k4) / 6

        #does one step estimation (x2)
        k1 = 2*h * f(s, t)
        k2 = 2*h * f(s+0.5*k1, t+h)
        k3 = 2*h * f(s+0.5*k2, t+h)
        k4 = 2*h * f(s+k3, t+2*h)
        s2 = (k1 + 2*k2 + 2*k3 + k4) / 6

        #calculating joint error in x and y dimensions
        epsilon_x = (s1[:3]-s2[:3]) / 30
        #epsilon_y = (s1[:3]-s2[:3]) / 30
        epsilon_y = np.array([0,0,0])
        error = np.sqrt(np.square(epsilon_x) + np.square(epsilon_y))

        #bounds step size h
        for i in range(len(error)):
            if error[i] < 1e-10:
                error[i] = h*dt/16

        rho = h*dt/error
        for i in range(len(rho)):
            if rho[i] > 16.0:
                rho[i] = 16.0

        rho = np.mean(rho)
        h *= rho**(1/4)

    s += s1
    xpoints.append(s[0])
    ypoints.append(s[1])
    zpoints.append(s[2])


fig10, ax10 = plt.subplots(2, 1, figsize = (12, 8))

ax10[0].plot(tpoints,ypoints[:-1],'m')
ax10[1].plot(xpoints,zpoints,'k')

plt.show()
