import numpy as np
import matplotlib.pyplot as plt

tend = 1
N = 100
L = 1

def solver(tend, N, L=1, obstacle=None,rho=1,nu=0.1,P=10):
    """Will solve a 2-d channel flow system according
        to Navier-Stokes equations for an NxN grid in
        a channel of length L, height L/2, iterating
        for <tend> seconds. There are two potential
        options for obstacle which will solve the
        system for a ball or wing placed in the channel.
        The density/viscosity of the fluid, as well as
        the pressure driving the flow can be adjusted."""


    #solving parameters
    dt = 1e-3      # Time step
    epsilon = dt/1000

    a = 2*L / N #grid spacing

    #creates ball obstacle
    if obstacle == 'ball':
        x0,y0 = N/2,N/2     #origin
        r = 0.1*N      #radius
        xgrid,ygrid = np.ogrid[-x0:N-x0+1,-y0:N-y0+1]
        mask = xgrid**2 + ygrid**2 <= r**2

    #creates wing obstacle
    elif obstacle == 'wing':
        mask = np.zeros([N+1,N+1],dtype=bool)
        mask[int(N/2):int(N/2)+3,int(N/6):int(5*N/6)] = True
        mask[int(N/2)+3:int(N/2)+4,int(N/3):int(5*N/6)] = True
        mask[int(N/2)+4:int(N/2)+5,int(N/2):int(2*N/3)] = True

        """x0,y0 = N/2,N/2      #y coordinate of wing center
        h = 8         #thickness of wing
        xgrid,ygrid = np.ogrid[-x0:N-x0,-y0:N-y0]
        mask = ygrid <= h"""

    t = 0.0
    tend += epsilon

    #initial conditions for velocity and pressure
    u = np.zeros([N+1,N+1],float) #x velocity
    up = np.zeros([N+1,N+1],float)

    v = np.zeros([N+1,N+1],float) #y velocity
    vp = np.zeros([N+1,N+1],float)

    p = P*np.ones([N+1,N+1],float) #pressure (scalar)
    pp = P*np.ones([N+1,N+1],float)

    while t < tend:

        """boundary conditions at walls"""
        u[0,:] = 0.0 #at y = 0
        v[0,:] = 0.0
        u[-1,:] = 0.0 #at y = L/2
        v[-1,:] = 0.0

        """boundary conditions for object"""
        if obstacle != None:
            #velocities in object's interior = 0
            u[mask] = 0.0
            v[mask] = 0.0

        """SOLVES DIFFERENTIAL EQUATIONS"""

        """Pressure"""
        for k in range(N+1):
            #boundary condition for x = 0
            pp[1:-1, 0] = ((p[1:-1,1] + p[1:-1,-1]
                          + p[2:,0] + p[0:-2,0]) / 4
                          - a**2 / 2)

            #pressure on the interior of the region
            pp[1:-1, 1:-1] = ((p[1:-1,2:] + p[1:-1,0:-2]
                             + p[2:,1:-1] + p[0:-2,1:-1]) / 4
                             - a**2 / 2)

            #boundary condition for x = L
            pp[1:-1, -1] = ((p[1:-1,0] + p[1:-1,-2]
                           + p[2:,-1] + p[0:-2,-1]) / 4
                           - a**2 / 2)

            #boundary conditions on upper/lower wall
            #derivative of pressure = 0 at top/bottom walls
            p[-1, :] = p[-2, :]
            p[0, :] = p[1, :]

            if obstacle == 'wing':
                #derivative of pressure at object boundary = 0
                #along bottom edge
                p[int(N/2),int(N/6):int(5*N/6)] =\
                 p[int(N/2)-1,int(N/6):int(5*N/6)]

                #along top edges
                p[int(N/2)+2,int(N/6):int(N/3)] =\
                 p[int(N/2)+3,int(N/6):int(N/3)]
                p[int(N/2)+3,int(N/3):int(N/2)] =\
                 p[int(N/2)+4,int(N/3):int(N/2)]
                p[int(N/2)+4,int(N/2):int(2*N/3)] =\
                 p[int(N/2)+5,int(N/2):int(2*N/3)]
                p[int(N/2)+3,int(2*N/3):int(5*N/6)] =\
                 p[int(N/2)+4,int(2*N/3):int(5*N/6)]


            #switches arrays storing each
            p, pp = pp, p

        """Velocity"""
        up[1:-1, 1:-1] = (u[1:-1, 1:-1]
                            - u[1:-1, 1:-1] * dt / a
                            * (u[1:-1, 1:-1] - u[1:-1, 0:-2])
                            - v[1:-1, 1:-1] * dt / a
                            * (u[1:-1, 1:-1] - u[0:-2, 1:-1])
                            - dt / (2 * rho * a)
                            * (p[1:-1, 2:] - p[1:-1, 0:-2])
                            + nu * dt / a**2
                            * (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2]
                            + u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1])
                            + P * dt)

        vp[1:-1, 1:-1] = (v[1:-1, 1:-1]
                            - u[1:-1, 1:-1] * dt / a
                            * (v[1:-1, 1:-1] - v[1:-1, 0:-2])
                            - v[1:-1, 1:-1] * dt / a
                            * (v[1:-1, 1:-1] - v[0:-2, 1:-1])
                            - dt / (2 * rho * a)
                            * (p[2:, 1:-1] - p[0:-2, 1:-1])
                            + nu * dt / a**2
                            * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2]
                            + v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1]))

        """Velocity boundary conditions"""
        #for x = 0
        up[1:-1, 0] = (u[1:-1, 0] - u[1:-1, 0] * dt / a
                    * (u[1:-1, 0] - u[1:-1, -1])
                    - v[1:-1, 0] * dt / a
                    * (u[1:-1, 0] - u[0:-2, 0])
                    - dt / (2 * rho * a)
                    * (p[1:-1, 1] - p[1:-1, -1])
                    + nu * dt / a**2
                    * (u[1:-1, 1] - 2 * u[1:-1, 0] + u[1:-1, -1]
                    + u[2:, 0] - 2 * u[1:-1, 0] + u[0:-2, 0])
                    + P * dt)

        #for x = 0
        vp[1:-1, 0] = (v[1:-1, 0] - u[1:-1, 0] * dt / a
                    * (v[1:-1, 0] - v[1:-1, -1])
                    - v[1:-1, 0] * dt / a
                    * (v[1:-1, 0] - v[0:-2, 0])
                    - dt / (2 * rho * a)
                    * (p[2:, 0] - p[0:-2, 0])
                    + nu * dt / a**2
                    * (v[1:-1, 1] - 2 * v[1:-1, 0] + v[1:-1, -1]
                    + v[2:, 0] - 2 * v[1:-1, 0] + v[0:-2, 0]))

        #for x = L
        up[1:-1, -1] = (u[1:-1, -1] - u[1:-1, -1] * dt / a
                     * (u[1:-1, -1] - u[1:-1, -2])
                     - v[1:-1, -1] * dt / a
                     * (u[1:-1, -1] - u[0:-2, -1])
                     - dt / (2 * rho * a)
                     * (p[1:-1, 0] - p[1:-1, -2])
                     + nu * dt / a**2
                     * (u[1:-1, 0] - 2 * u[1:-1,-1] + u[1:-1, -2]
                     + u[2:, -1] - 2 * u[1:-1, -1] + u[0:-2, -1])
                     + P * dt)

        #for x = L
        vp[1:-1, -1] = (v[1:-1, -1] - u[1:-1, -1] * dt / a
                     * (v[1:-1, -1] - v[1:-1, -2])
                     - v[1:-1, -1] * dt / a
                     * (v[1:-1, -1] - v[0:-2, -1])
                     - dt / (2 * rho * a)
                     * (p[2:, -1] - p[0:-2, -1])
                     + nu * dt / a**2
                     * (v[1:-1, 0] - 2 * v[1:-1, -1] + v[1:-1, -2]
                     + v[2:, -1] - 2 * v[1:-1, -1] + v[0:-2, -1]))


        #switching arrays for velocity and pressure
        u, up = up, u
        v, vp = vp, v

        t += dt #advance time step

    return u, v, p

def basic_bc(u, v):
    """Sets the basic boundary conditions for the
        channel flow. u and v are 0 on upper and
        lower boundaries"""

    u[0,:] = 0.0 #at y=0
    v[0,:] = 0.0
    u[-1,:] = 0.0 #at y=L/2
    v[-1,:] = 0.0

    return u,v

def theoretical_u1(y,p,L=1,R=L/4,rho=1,nu=0.1):
    """Finds x velocity of fluid at given
        height and pressure in channel"""

    dp = np.ones(N+1)
    eta = nu*rho
    prefac = 0.5*np.square(R)/(rho*eta)
    dp[1:N] = (p[int(L/4),2:N+1] - p[int(L/4),0:N-1]) / (2*(L/N))
    height = 1 - (y/R)**2

    return prefac*dp*height

def solved_u(y,L=1,nu=0.1):
    """Finds x velocity of fluid at given
        height and pressure in channel"""

    u = -10 / nu
    u *= (np.square(y)/2 - L*y/4)
    return u
