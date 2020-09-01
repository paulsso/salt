import numpy as np

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

    r_ni = np.sqrt((Bx-Ax)**2 + (By-Ay)**2 + (Bz-Az)**2)
    r_in = r_ni.T

    return r_nm, r_im, r_in, r_ni

def TransferMatrices(mediumProperties, zPosProperties, zNegProperties, r_nm, r_im, r_in, r_ni):
    """ Method for computing transfer matrices """
    Sn = 1e-6
    Si = 1e-6

    c = mediumProperties["SpeedOfSound"]
    f1 = zPosProperties["TransFreq"]

    if "TransFreq" in zNegProperties:
        f2 = zNegProperties["TransFreq"]
    else:
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
        T_TR = Sn*np.exp(-1j*kk1*r_ni - 1j*kk2*r_in)/(r_in)
        T_RT = Si*np.exp(-1j*kk1*r_in - 1j*kk2*r_ni)/(r_ni)
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
