import numpy as np
import matplotlib.pyplot as plt

def channel_flow(tend, N, L, obstacle=None):
    """Will solve a 2-d channel flow system according
        to Navier-Stokes equations for an NxN grid in
        a channel of length L, height L/2, iterating
        for <tend> seconds"""

    # physical parameters of fluid
    rho = 1        #density
    nu = 0.1       #viscosity
    P = 10         #constant term
    #solving parameters
    dt = 1e-3      # Time-step
    epsilon = dt/1000

    step = 2 / (N)# - 1)

    #initiates quiver plot
    #x = np.linspace(0, L, N)
    #y = np.linspace(0, L/2, N)
    #X, Y = np.meshgrid(x, y)

    #creates ball obstacle
    if obstacle == 'ball':
        x0,y0 = N/2,N/2     #origin
        r = 5      #radius
        xgrid,ygrid = np.ogrid[-x0:N-x0,-y0:N-y0]
        mask = xgrid**2 + ygrid**2 <= r**2
    elif obstacle == 'wing':
        x0,y0 = N/2,N/2      #y coordinate of wing center
        h = 8         #thickness of wing
        xgrid,ygrid = np.ogrid[-x0:N-x0,-y0:N-y0]
        mask = ygrid <= h

    t = 0.0
    #t1 = tend / 4
    #t2 = 2*tend/4
    #t3 = 3*tend/4
    #t4 = tend
    tend += epsilon

    #initial conditions for velocity and pressure
    u = np.zeros([N+1,N+1],float) #x velocity
    up = np.zeros([N+1,N+1],float)

    v = np.zeros([N+1,N+1],float) #y velocity
    vp = np.zeros([N+1,N+1],float)

    p = 5*np.ones([N+1,N+1],float) #pressure (scalar)
    pp = 5*np.ones([N+1,N+1],float)

    #b = np.zeros((ny, nx))


    while t < tend:

        """boundary conditions at walls"""
        u[0,:] = 0.0 #at y = 0
        v[0,:] = 0.0
        u[-1,:] = 0.0 #at y = L/2
        v[-1,:] = 0.0

        """boundary conditions at obstacle surface"""
        if obstacle != None:
            #velocities in object's interior = 0
            u[mask] = 0.0
            v[mask] = 0.0
            p[mask] = 0.0

            if obstacle == 'ball':
                #pressure in object's interior = 0
                pass
                #derivative of pressure at object boundary = 0

            elif obstacle == 'wing':
                #pressure in object's interior = 0
                pass
                #derivative of pressure at object boundary = 0



        """SOLVES DIFFERENTIAL EQUATIONS"""

        """Pressure boundary conditions""" #i dont need to multiply by step^2 on top and divide on the bottom because i used same step sizes
        for k in range(N+1):
            pp[1:-1, 1:-1] = (((p[1:-1,2:] + p[1:-1,0:-2])
                             + (p[2:,1:-1] + p[0:-2,1:-1]))
                             / 4
                             - step**2 * step**2 / (2 * (step**2 + step**2))
                             )
                             #* b[1:-1, 1:-1])

            #boundary condition for x = L
            pp[1:-1, -1] = (
                            ((p[1:-1,0] + p[1:-1,-2])
                           + (p[2:,-1] + p[0:-2,-1]))
                           / 4
                           - step**2 * step**2 / (2 * (step**2 + step**2))
                           )
                           #* b[1:-1, -1])

            #boundary condition for x = 0
            pp[1:-1, 0] = (((p[1:-1,1] + p[1:-1,-1])
                          + (p[2:,0] + p[0:-2,0]))
                          / 4
                          - step**2 * step**2 / (2 * (step**2 + step**2))
                          )
                          #* b[1:-1, 0])

            #boundary conditions on wall
            #derivative of pressure =0 at top/bottom walls
            p[-1, :] = p[-2, :]
            p[0, :] = p[1, :]

            p, pp = pp, p


        up[1:-1, 1:-1] = (u[1:-1, 1:-1]
                        - u[1:-1, 1:-1] * dt / step
                        * (u[1:-1, 1:-1] - u[1:-1, 0:-2])
                        - v[1:-1, 1:-1] * dt / step
                        * (u[1:-1, 1:-1] - u[0:-2, 1:-1])
                        - dt / (2 * rho * step)
                        * (p[1:-1, 2:] - p[1:-1, 0:-2])
                        + nu * (dt / step**2
                        * (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, 0:-2])
                        + dt / step**2
                        * (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1]))
                        + P * dt)

        vp[1:-1, 1:-1] = (v[1:-1, 1:-1]
                        - u[1:-1, 1:-1] * dt / step
                        * (v[1:-1, 1:-1] - v[1:-1, 0:-2])
                        - v[1:-1, 1:-1] * dt / step
                        * (v[1:-1, 1:-1] - v[0:-2, 1:-1])
                        - dt / (2 * rho * step)
                        * (p[2:, 1:-1] - p[0:-2, 1:-1])
                        + nu * (dt / step**2
                        * (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, 0:-2])
                        + dt / step**2
                        * (v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[0:-2, 1:-1])))

        """Velocity boundary conditions"""
        #for x = 0
        up[1:-1, 0] = (u[1:-1, 0] - u[1:-1, 0] * dt / step
                    * (u[1:-1, 0] - u[1:-1, -1])
                    - v[1:-1, 0] * dt / step
                    * (u[1:-1, 0] - u[0:-2, 0])
                    - dt / (2 * rho * step)
                    * (p[1:-1, 1] - p[1:-1, -1])
                    + nu * (dt / step**2
                    * (u[1:-1, 1] - 2 * u[1:-1, 0] + u[1:-1, -1])
                    + dt / step**2
                    * (u[2:, 0] - 2 * u[1:-1, 0] + u[0:-2, 0]))
                    + P * dt)

        #for x = L
        up[1:-1, -1] = (u[1:-1, -1] - u[1:-1, -1] * dt / step
                     * (u[1:-1, -1] - u[1:-1, -2])
                     - v[1:-1, -1] * dt / step
                     * (u[1:-1, -1] - u[0:-2, -1])
                     - dt / (2 * rho * step)
                     * (p[1:-1, 0] - p[1:-1, -2])
                     + nu * (dt / step**2
                     * (u[1:-1, 0] - 2 * u[1:-1,-1] + u[1:-1, -2])
                     + dt / step**2
                     * (u[2:, -1] - 2 * u[1:-1, -1] + u[0:-2, -1]))
                     + P * dt)

        #for x = 0
        vp[1:-1, 0] = (v[1:-1, 0] - u[1:-1, 0] * dt / step
                    * (v[1:-1, 0] - v[1:-1, -1])
                    - v[1:-1, 0] * dt / step
                    * (v[1:-1, 0] - v[0:-2, 0])
                    - dt / (2 * rho * step)
                    * (p[2:, 0] - p[0:-2, 0])
                    + nu * (dt / step**2
                    * (v[1:-1, 1] - 2 * v[1:-1, 0] + v[1:-1, -1])
                    + dt / step**2
                    * (v[2:, 0] - 2 * v[1:-1, 0] + v[0:-2, 0])))

        #for x = L
        vp[1:-1, -1] = (v[1:-1, -1] - u[1:-1, -1] * dt / step
                     * (v[1:-1, -1] - v[1:-1, -2])
                     - v[1:-1, -1] * dt / step
                     * (v[1:-1, -1] - v[0:-2, -1])
                     - dt / (2 * rho * step)
                     * (p[2:, -1] - p[0:-2, -1])
                     + nu * (dt / step**2
                     * (v[1:-1, 0] - 2 * v[1:-1, -1] + v[1:-1, -2])
                     + dt / step**2
                     * (v[2:, -1] - 2 * v[1:-1, -1] + v[0:-2, -1])))


        #switching arrays for velocity and pressure
        u, up = up, u
        v, vp = vp, v
        #p, pp = pp, p

        t += dt #advance time step

    return u, v, p

    #fig = plt.figure(figsize = (10,5), dpi=100)
    #plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3]);

    fig1, ax1 = plt.subplots(1, 1, figsize = (10, 5))

    ax1.imshow(p,cmap='hot')
    plt.show()
