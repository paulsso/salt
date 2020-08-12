import time, sys
import numpy as np

def CreateTransducer(properties):
    """ Method for constructing a single transducer """
    R = properties["Radius"]*1e3
    r_c = properties["RadiusCurvature"]
    s = properties["Orientation"]
    d = properties["Displacement"]

    z0 = s*d

    R_span = np.arange(-R*1e-3,R*1e-3+1e-3,1e-3)
    R_length = len(R_span)

    X, Y = np.mgrid[-R:R+1, -R:R+1]
    X = X*1e-3
    Y = Y*1e-3
    Z = np.zeros([R_length, R_length])

    rows, cols = np.mgrid[0:R_length, 0:R_length]
    C = np.sqrt((rows-R-1)**2 + (cols-R-1)**2)<=R

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

def CreateMedium(Vx,Vz,Ux,Uz):
    """Function generates points in a plane between components where pressure will be calculated"""
    xMax = np.round(np.max(Vx)+1e-3,3)
    xMin = np.round(-np.max(Vx)-1e-3,3)
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

def DistanceMatrices(Vx, Vy, Vz, Ux, Uy, Uz, Mx, My, Mz):
    """Function computes the distance between measurement plane and components"""
    nT = len(Vx)
    nR = len(Ux)
    nM = len(Mx)

    Ax = np.repeat(Vx,nM,0).reshape(nT,nM).T
    Bx = np.repeat(Mx.reshape(nM,1),nT,axis=1)

    Ay = np.repeat(Vy,nM,0).reshape(nT,nM).T
    By = np.repeat(My.reshape(nM,1),nT,axis=1)

    Az = np.repeat(Vz,nM,0).reshape(nT,nM).T
    Bz = np.repeat(Mz.reshape(nM,1),nT,axis=1)

    r_nm = np.sqrt((Bx-Ax)**2 + (By-Ay)**2 + (Bz-Az)**2)

    Ax = np.repeat(Ux,nM,0).reshape(nR,nM).T
    Bx = np.repeat(Mx.reshape(nM,1),nR,axis=1)

    Ay = np.repeat(Uy,nM,0).reshape(nR,nM).T
    By = np.repeat(My.reshape(nM,1),nR,axis=1)

    Az = np.repeat(Uz,nM,0).reshape(nR,nM).T
    Bz = np.repeat(Mz.reshape(nM,1),nR,axis=1)

    r_im = np.sqrt((Bx-Ax)**2 + (By-Ay)**2 + (Bz-Az)**2)

    Ax = np.repeat(Ux, nT, 1)
    Bx = np.repeat(Vx.T, nR, 0)

    Ay = np.repeat(Uy, nT, 1)
    By = np.repeat(Vy.T, nR, 0)

    Az = np.repeat(Uz, nT, 1)
    Bz = np.repeat(Vz.T, nR, 0)

    r_in = np.sqrt((Bx-Ax)**2 + (By-Ay)**2 + (Bz-Az)**2)
    r_ni = r_in.T

    return r_nm, r_im, r_in, r_ni

def TransferMatrices(mediumProperties, zPosProperties, zNegProperties, r_nm, r_im, r_in, r_ni):
    """ Method for computing transfer matrices """
    Sn = 1e-6
    Si = 1e-6

    c = mediumProperties["SpeedOfSound"]
    f1 = zPosProperties["TransFreq"]
    try:
        f2 = zNegProperties["TransFreq"]
    except:
        f2 = 0

    if f1 != 0:
        wL1 = c/f1
        kk1 = 2*np.pi/wL1
    elif f1 == 0:
        wL1 = 0
        kk1 = 0
    if f2 != 0:
        wL2 = c/f2
        kk2 = 2*np.pi/wL2
    elif f2 == 0:
        wL2 = 0
        kk2 = 0

    if zPosProperties["Type"] == "Array" and zNegProperties["Type"] == "Array":
        T_TR = Sn*np.exp(-1j*kk1*r_in - 1j*kk2*r_ni)/(r_in)
        T_RT = Si*np.exp(-1j*kk1*r_ni - 1j*kk2*r_in)/(r_ni)
        T_RM = Sn*np.exp(-1j*kk1*r_im)/r_im
        T_TM = Si*np.exp(-1j*kk2*r_nm)/r_nm
    elif zPosProperties["Type"] == "Array" and zNegProperties["Type"] =="Reflector":
        T_TR = Sn*np.exp(-1j*kk1*r_in)/(r_in)
        T_RT = Si*np.exp(-1j*kk1*r_ni)/(r_ni)
        T_RM = Sn*np.exp(-1j*kk1*r_im)/(r_im)
        T_TM = Si*np.exp(-1j*kk1*r_nm)/(r_nm)
    elif zPosProperties["Type"] == "Reflector" and zNegProperties["Type"] =="Array":
        T_TR = Sn*np.exp(-1j*kk2*r_in)/(r_in)
        T_RT = Si*np.exp(-1j*kk2*r_ni)/(r_ni)
        T_RM = Sn*np.exp(-1j*kk2*r_im)/(r_im)
        T_TM = Si*np.exp(-1j*kk2*r_nm)/(r_nm)
    elif zPosProperties["Type"] == "Transducer" and zNegProperties["Type"] =="Reflector":
        T_TR = Sn*np.exp(-1j*kk1*r_in)/(r_in)
        T_RT = Si*np.exp(-1j*kk1*r_ni)/(r_ni)
        T_RM = Sn*np.exp(-1j*kk1*r_im)/(r_im)
        T_TM = Si*np.exp(-1j*kk1*r_nm)/(r_nm)

    return T_TR, T_RT, T_RM, T_TM

def ComputePressure(mediumProperties,T_TR,T_RT,T_RM,T_TM,zPosProperties,zNegProperties,nT,nR,nM):
    """ Method computes pressure matrix in accordance with the matrix method """
    t = 0
    f1 = zPosProperties["TransFreq"]

    try:
        f2 = zNegProperties["TransFreq"]
    except:
        f2 = 0

    d1 = 1e-6
    d2 = 1e-6

    rho = mediumProperties["Density"]
    c = mediumProperties["SpeedOfSound"]

    if f1 != 0:
        wL1 = c/f1
        omega1 = 2 * np.pi * f1
        U1 = np.ones([nT, 1])*d1*np.exp(-1j*(omega1*t))
        A1 = (1j/wL1)
        C1 = omega1*rho*c/wL1
    if f1 == 0:
        wL1 = 0
        omega1 = 0
        U1 = np.zeros([nT, 1])
        A1 = 0
        C1 = 0

    if f2 != 0:
        wL2 = c/f2
        omega2 = 2 * np.pi * f2
        U2 = -np.ones([nR, 1])*d2*np.exp(-1j*(omega2*t))
        A2 = (1j/wL2)
        C2 = omega2*rho*c/wL2
    if f2 == 0:
        wL2 = 0
        omega2 = 0
        U2 = np.zeros([nR, 1])
        A2 = 0
        C2 = 0

    if zPosProperties["Type"] == "Transducer" and zNegProperties["Type"] == "Reflector":
        # TODO : COMPUTE PRESSURE WITH MULTIPLE REFLECTIONS
        PT0 = (C1)*T_TM@U1;
        PT1 = (C1)*(A1)*T_RM@T_TR@U1;
        PT2 = (C1)*(A1**2)*T_TM@T_RT@T_TR@U1;
        PT3 = (C1)*(A1**3)*T_RM@T_TR@T_RT@T_TR@U1;
        PT4 = (C1)*(A1**4)*T_TM@T_RT@T_TR@T_RT@T_TR@U1;
        PT5 = (C1)*(A1**5)*T_RM@T_TR@T_RT@T_TR@T_RT@T_TR@U1;
        PT6 = (C1)*(A1**6)*T_TM@T_RT@T_TR@T_RT@T_TR@T_RT@T_TR@U1;
        PT = PT0 + PT1 + PT2 + PT3 + PT4 + PT5 + PT6
        PR = PT*0
    elif (zPosProperties["Type"] == "Array" and zNegProperties["Type"] == "Array") or (zPosProperties["Type"] == "Transducer" and zNegProperties["Type"] == "Array"):
        PT = (C1)*T_TM@U1
        PR = (C2)*T_RM@U2
    elif zPosProperties["Type"] == "Array" and zNegProperties["Type"] == "Reflector":
        PT = (C1)*T_TM@U1
        PR = (C1)*(A1)*T_RM@T_TR@U1
    elif zPosProperties["Type"] == "Reflector" and zNegProperties["Type"] == "Array":
        PR = (C2)*T_RM@U2
        PT = (C2)*(A2)*T_TM@T_RT@U2

    P = PT + PR
    return P

def ComputeRelativePotential(Ptotal, mediumProperties, zPosProperties, zNegProperties):
    """ Method computes the relative acoustic potential """
    c = mediumProperties["SpeedOfSound"]
    rho = mediumProperties["Density"]
    f = zPosProperties["TransFreq"]
    w = c/f
    p = np.real(Ptotal)
    sz = Ptotal.shape
    M = sz[0]*sz[1]
    phi = -Ptotal/(1j*w*rho)

    gradx, grady = np.gradient(phi,5e-4,5e-4)
    gradP = np.sqrt(gradx**2 + grady**2)
    T1 = (np.real(Ptotal*np.conj(Ptotal))/M)/(3*rho*c**2)
    T2 = -0.5*rho*(np.real(np.multiply(gradP,np.conj(gradP)))/M)
    potential = T1+T2

    return potential

def MatrixMethod(mediumProperties,zPosProperties,zNegProperties):

    if zPosProperties["Type"] == "Array":
        print("Creating array...")
        start = time.time()
        Vx, Vy, Vz = CreateArray(zPosProperties)
        end = time.time()
        diff1 = end - start
        print("Create array took %.6f" % diff1, "seconds")

    elif zPosProperties["Type"] == "Transducer":
        print("Creating Transducer...")
        start = time.time()
        Vx, Vy, Vz = CreateTransducer(zPosProperties)
        end = time.time()
        diff1 = end - start
        print("Create transducer took %.6f" % diff1, "seconds")
    else:
        Vx = 0
        Vy = 0
        Vz = 0

    if zNegProperties["Type"] == "Array":
        print("Creating array...")
        start = time.time()
        Ux, Uy, Uz = CreateArray(zNegProperties)
        end = time.time()
        diff2 = end - start
        print("Create array took %.6f" % diff2, "seconds")

    elif zNegProperties["Type"] == "Reflector":
        print("Creating reflector...")
        start = time.time()
        Ux, Uy, Uz = CreateReflector(zNegProperties)
        end = time.time()
        diff2 = end - start
        print("Create reflector took %.6f" % diff2, "seconds")
    else:
        Ux = 0
        Uy = 0
        Uz = 0

    print("Creating medium...")
    start = time.time()
    Mx, My, Mz, x_span, z_span = CreateMedium(Vx, Vz, Ux, Uz)
    end = time.time()
    diff3 = end - start
    print("Creating medium took %.6f" % diff3, "seconds")

    print("Computing distance matrices...")
    start = time.time()
    r_nm, r_im, r_in, r_ni = DistanceMatrices(Vx, Vy, Vz, Ux, Uy, Uz, Mx, My, Mz)
    end = time.time()
    diff4 = end - start
    print("Computing distance matrices took %.6f" % diff4, "seconds")

    print("Computing transfer matrices...")
    start = time.time()
    T_TR, T_RT, T_RM, T_TM = TransferMatrices(mediumProperties, zPosProperties, zNegProperties, r_nm, r_im, r_in, r_ni)
    end = time.time()
    diff5 = end - start
    print("Computing transfer matrices took %.6f" % diff5, "seconds")

    nT = len(Vx)
    nR = len(Ux)
    nM = len(Mx)
    print("Computing pressure matrix...")
    start = time.time()
    pressure = ComputePressure(mediumProperties,T_TR,T_RT,T_RM,T_TM,zPosProperties,zNegProperties,nT,nR,nM)
    end = time.time()
    diff6 = end - start
    print("Computing pressure matrices took %.6f" % diff6, "seconds")

    x = len(x_span)
    z = len(z_span)
    P_ = pressure.reshape([z, x])

    print("Computing relative acoustic potential...")
    start = time.time()
    relative_potential = ComputeRelativePotential(P_, mediumProperties, zPosProperties, zNegProperties)
    end = time.time()
    diff7 = end - start
    print("Computing relative acoustic potential took %.6f" % diff7, "seconds")

    return relative_potential, pressure, x_span, z_span

def CreateGeometry(zPosProperties, zNegProperties):
    if zPosProperties["Type"] == "Array":
        Vx, Vy, Vz = CreateArray(zPosProperties)
    elif zPosProperties["Type"] == "Transducer":
        Vx, Vy, Vz = CreateTransducer(zPosProperties)
    else:
        Vx = 0
        Vy = 0
        Vz = 0

    if zNegProperties["Type"] == "Array":
        Ux, Uy, Uz = CreateArray(zNegProperties)
    elif zNegProperties["Type"] == "Reflector":
        Ux, Uy, Uz = CreateReflector(zNegProperties)
    else:
        Ux = 0
        Uy = 0
        Uz = 0
    return Ux, Uy, Uz, Vx, Vy, Vz
