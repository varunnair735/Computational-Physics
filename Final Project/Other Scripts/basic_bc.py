def basic_bc(u, v):
    """Sets the basic boundary conditions for the
        channel flow. u and v are 0 on upper and
        lower boundaries"""

    u[0,:] = 0.0 #at y=0
    v[0,:] = 0.0
    u[-1,:] = 0.0 #at y=L/2
    v[-1,:] = 0.0

    return u,v
