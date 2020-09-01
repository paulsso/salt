import time, sys
import numpy as np
from geometries import CreateTransducer, CreateReflector, CreateArray, CreateMedium
from matrices import DistanceMatrices, TransferMatrices
from computation import ComputePressure, ComputeRelativePotential

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

    c = mediumProperties["SpeedOfSound"]
    rho = mediumProperties["Density"]

    p = np.real(pressure)**2

    acoustic_radiation_pressure = np.real(((p)/(4*rho*c**2)).reshape([z, x]))

    return acoustic_radiation_pressure, relative_potential, pressure, x_span, z_span

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
