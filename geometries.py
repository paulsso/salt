import numpy as np

def CreateArray(properties):
    """ Method for constructing an array of transducers """
    transducer_radius = 4.5*1e-3
    transducer_per_layer = properties["Depth"]
    layers = properties["Layers"]
    socket_radius = properties["Radius"]
    r_c = properties["RadiusCurvature"]
    h = properties["Displacement"]
    s = properties["Orientation"]

    tz0 = s*(h - r_c)

    X, Y = np.mgrid[-socket_radius*1e3:socket_radius*1e3,-socket_radius*1e3:socket_radius*1e3]
    X = X*1e-3
    Y = Y*1e-3
    r = np.array([24, 46, 68, 90, 112, 135, 158, 180])*1e-3
    Vx = []
    Vy = []

    for kk in range(layers):
        n = int(transducer_per_layer[kk])
        beta = np.linspace(2*np.pi/n, 2*np.pi, n)

        for ii in range(n):
            C = np.sqrt((X - 0.5*r[kk]*np.cos(beta[ii]))**2 + (Y - 0.5*r[kk]*np.sin(beta[ii]))**2) <= transducer_radius
            vec1 = X[C]
            vec2 = Y[C]
            Vx = np.append(Vx, vec1)
            Vy = np.append(Vy, vec2)

        if properties["Concave"]:
            Vz = tz0 + s*np.sqrt(r_c**2 - np.power(Vx,2) - np.power(Vy,2))
        elif not properties["Concave"]:
            Vz = s*h*np.ones([len(Vx),1])

    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])
    S = Vx + Vz + Vy

    if np.isnan(S).any():
        Vx = np.delete(Vx,np.argwhere(np.isnan(S)))
        Vy = np.delete(Vy,np.argwhere(np.isnan(S)))
        Vz = np.delete(Vz,np.argwhere(np.isnan(S)))
    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])

    return Vx, Vy, Vz

def CreateMedium(Vx,Vz,Ux,Uz):
    """Function generates points in a plane between components where pressure will be calculated"""
    xMax = np.round(np.max(Ux)+1e-3,3)
    xMin = np.round(-np.max(Ux)-1e-3,3)
    zMax = np.round(np.min(Vz)-1e-3,3)
    zMin = np.round(np.max(Uz)+1e-3,3)

    x_span = np.arange(xMin, xMax+1e-3, 1e-3)
    z_span = np.arange(zMin, zMax+1e-3, 1e-3)

    Mx, Mz = np.meshgrid(x_span, z_span, sparse=False)

    sz = Mx.shape
    My = np.zeros([sz[0],sz[1]])

    szx = Mx.shape
    Mx = Mx.reshape(szx[0]*szx[1],1)

    szy = My.shape
    My = My.reshape(szy[0]*szy[1],1)

    szz = Mz.shape
    Mz = Mz.reshape(szz[0]*szz[1],1)

    return Mx, My, Mz, x_span, z_span

def CreateReflector(properties):
    """ Method for constructing a reflector """
    R = properties["Radius"]*1e3
    r_c = properties["RadiusCurvature"]
    s = properties["Orientation"]
    d = properties["Displacement"]
    z0 = s*d

    R_span = np.arange(-R*1e-3,R*1e-3,1e-3)
    R_length = len(R_span)

    X, Y = np.mgrid[-R:R, -R:R]
    X = X*1e-3
    Y = Y*1e-3
    Z = np.zeros([R_length, R_length])

    rows, cols = np.mgrid[0:R_length, 0:R_length]
    C = np.sqrt((rows-R)**2 + (cols-R)**2)<R

    Vx = X[C]
    Vy = Y[C]
    if properties["Concave"]:
        if s == -1:
            Vz = z0 - np.sqrt(r_c**2 - (Vx)**2 - (Vy)**2) + r_c
        elif s == 1:
            Vz = z0 + np.sqrt(r_c**2 - (Vx)**2 - (Vy)**2) - r_c
        else:
            Vz = 0
    elif not properties["Concave"]:
        if s == -1:
            Vz = z0*np.ones([len(Vx),1])
        elif s == 1:
            Vz = z0*np.ones([len(Vx),1])
        else:
            Vz = 0


    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])
    S = Vx + Vz + Vy

    if np.isnan(S).any():
        Vx = np.delete(Vx,np.argwhere(np.isnan(S)))
        Vy = np.delete(Vy,np.argwhere(np.isnan(S)))
        Vz = np.delete(Vz,np.argwhere(np.isnan(S)))

    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])

    return Vx, Vy, Vz

def CreateTransducer(properties):
    """ Method for constructing a single transducer """
    R = properties["Radius"]*1e3
    r_c = properties["RadiusCurvature"]
    s = properties["Orientation"]
    d = properties["Displacement"]

    z0 = s*d

    R_span = np.arange(-R*1e-3,R*1e-3,1e-3)
    R_length = len(R_span)

    X, Y = np.mgrid[-R:R,-R:R]
    X = X*1e-3
    Y = Y*1e-3
    Z = np.zeros([R_length, R_length])

    rows, cols = np.mgrid[0:R_length, 0:R_length]
    C = np.sqrt((rows-R)**2 + (cols-R)**2)<R

    Vx = X[C]
    Vy = Y[C]

    if properties["Concave"]:
        if s == -1:
            Vz = z0 - np.sqrt(r_c**2 - (Vx)**2 - (Vy)**2) + r_c
        elif s == 1:
            Vz = z0 + np.sqrt(r_c**2 - (Vx)**2 - (Vy)**2) - r_c
        else:
            Vz = 0
    elif not properties["Concave"]:
        if s == -1:
            Vz = z0*np.ones([len(Vx),1])
        elif s == 1:
            Vz = z0*np.ones([len(Vx),1])
        else:
            Vz = 0

    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])
    S = Vx + Vz + Vy

    if np.isnan(S).any():
        Vx = np.delete(Vx,np.argwhere(np.isnan(S)))
        Vy = np.delete(Vy,np.argwhere(np.isnan(S)))
        Vz = np.delete(Vz,np.argwhere(np.isnan(S)))

    Vx = Vx.reshape([len(Vx),1])
    Vy = Vy.reshape([len(Vy),1])
    Vz = Vz.reshape([len(Vz),1])

    return Vx, Vy, Vz
