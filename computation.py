import numpy as np

def ComputePressure(mediumProperties,T_TR,T_RT,T_RM,T_TM,zPosProperties,zNegProperties,nT,nR,nM):
    """ Method computes pressure matrix in accordance with the matrix method """

    t1 = zPosProperties["Phase"]
    f1 = zPosProperties["TransFreq"]

    if "Phase" in zNegProperties:
        t2 = zNegProperties["Phase"]
    else:
        t2 = 0

    if "TransFreq" in zNegProperties:
        f2 = zNegProperties["TransFreq"]
    else:
        f2 = 0

    d1 = 1e-6
    d2 = 1e-6

    rho = mediumProperties["Density"]
    c = mediumProperties["SpeedOfSound"]

    if f1 != 0:
        wL1 = c/f1
        omega1 = 2 * np.pi * f1
        U1 = np.ones([nT, 1])*d1*np.exp(-1j*(omega1*t1))
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
        U2 = -np.ones([nR, 1])*d2*np.exp(-1j*(omega2*t2))
        A2 = (1j/wL2)
        C2 = omega2*rho*c/wL2
    if f2 == 0:
        wL2 = 0
        omega2 = 0
        U2 = np.zeros([nR, 1])
        A2 = 0
        C2 = 0

    if zPosProperties["Type"] == "Transducer" and zNegProperties["Type"] == "Reflector":
        PT0 = (C1)*T_TM@U1;
        PT1 = (C1)*(A1)*T_RM@T_TR@U1;
        PT2 = (C1)*(A1**2)*T_TM@T_RT@T_TR@U1;
        PT3 = (C1)*(A1**3)*T_RM@T_TR@T_RT@T_TR@U1;
        PT4 = (C1)*(A1**4)*T_TM@T_RT@T_TR@T_RT@T_TR@U1;
        PT5 = (C1)*(A1**5)*T_RM@T_TR@T_RT@T_TR@T_RT@T_TR@U1;
        PT6 = (C1)*(A1**6)*T_TM@T_RT@T_TR@T_RT@T_TR@T_RT@T_TR@U1;
        PT = PT0 + PT2 + PT4 + PT6
        PR = PT1 + PT3 + PT5
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
